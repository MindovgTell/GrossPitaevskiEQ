#include "workflow.hpp"
#include "Core/simple_image.hpp"
#include "Solvers/SplitStep.hpp"
#include "Solvers/BogoliubovdeGennes.hpp"
#include "Solvers/bdg_spectral_analysis.hpp"

namespace {

using namespace workflow;

struct Options {
    gpes::bdg::BdGConfig solver;
    double stationary_tolerance = 1e-4;
    double boundary_tolerance = 1e-8;
    double residual_tolerance = 1e-4;
    bool allow_unconverged = false;
    bool no_mode_fields = false;
};

bool valid_mode(const gpes::bdg::BdGMode& mode, double tolerance) {
    return std::isfinite(mode.omega.real()) && std::isfinite(mode.omega.imag()) &&
        mode.u.allFinite() && mode.v.allFinite() && std::isfinite(mode.relative_residual) &&
        mode.relative_residual <= tolerance;
}

int calculate(const fs::path& input, const fs::path& output, const Options& options) {
    std::ostringstream report;
    const auto& cfg = options.solver;
    report << std::setprecision(17) << std::boolalpha
        << "requested_modes=" << cfg.num_modes << '\n'
        << "requested_ncv=" << cfg.ncv << '\n'
        << "max_iter=" << cfg.max_iter << '\n' << "eigensolver_tolerance=" << cfg.tol << '\n'
        << "eigenpair_residual_tolerance=" << options.residual_tolerance << '\n'
        << "shift_real=" << cfg.shift.real() << '\n' << "shift_imag=" << cfg.shift.imag() << '\n'
        << "linear_tolerance=" << cfg.shift_invert_linear_tol << '\n'
        << "linear_max_iter=" << cfg.shift_invert_linear_max_iter << '\n'
        << "linear_restart=" << cfg.shift_invert_linear_restart << '\n'
        << "dense_max_dimension=" << cfg.shift_invert_dense_max_dim << '\n'
        << "stationary_residual_tolerance=" << options.stationary_tolerance << '\n'
        << "boundary_density_tolerance=" << options.boundary_tolerance << '\n'
        << "allow_unconverged=" << options.allow_unconverged << '\n'
        << "save_mode_fields=" << !options.no_mode_fields << '\n';
    bool attempted = false;
    const auto status = [&](const std::string& state) {
        write_text(output / "bdg_status.txt", "status=" + state + "\nbdg_attempted=" +
            (attempted ? "true\n" : "false\n") + report.str());
    };
    try {
        write_text(output / "RUNNING", "bdg\n");
        write_text(output / "source.txt", fs::absolute(input).string() + '\n');
        status("loading");
        fs::create_directories(output / "log");
        gpes::log::Log::getInstance().setLogDirectory(output / "log");
        const auto params = read_parameters(input);
        model_warning(params);
        fs::copy_file(input / "wavefunction/parameters.txt", output / "parameters.txt");
        const auto grid = params.grid();
        const int dimension = 2 * params.nx * params.ny;
        if (cfg.num_modes >= dimension - 1) throw std::invalid_argument("Too many modes for the saved grid");
        report << "effective_ncv=" << std::min(dimension, std::max(cfg.ncv, 2 * cfg.num_modes + 1)) << '\n';
        Wavefunction psi(grid);
        gpes::io::load_wavefunction_csv(psi, input / "wavefunction/wavefunction.csv", true);
        check_normalization(psi, params);
        gpes::solvers::StationaryDiagnostics diagnostics;
        {
            gpes::solvers::SplitStep<Dim> evaluator(grid, psi, params.physics(), gpes::SimConfig{});
            diagnostics = evaluator.diagnostics(psi);
        }
        const double boundary = boundary_ratio(psi);
        const bool stationary = std::isfinite(diagnostics.relative_residual) &&
            diagnostics.relative_residual <= options.stationary_tolerance &&
            std::isfinite(boundary) && boundary <= options.boundary_tolerance;
        report << "energy=" << diagnostics.energy << '\n'
            << "chemical_potential=" << diagnostics.chemical_potential << '\n'
            << "stationary_relative_residual=" << diagnostics.relative_residual << '\n'
            << "boundary_to_peak_density=" << boundary << '\n'
            << "particle_number=" << particle_number(psi) << '\n'
            << "input_quality_passed=" << stationary << '\n';
        status("validated");
        std::cout << "Input: " << input << "\nBdG output: " << output
                  << "\nStationary residual: " << diagnostics.relative_residual << std::endl;
        if (!stationary && !options.allow_unconverged) {
            status("rejected_input");
            throw std::runtime_error("Saved state fails stationarity or boundary checks; use --allow-unconverged only for diagnostics");
        }
        if (!stationary) std::cerr << "Warning: nonstationary or boundary-contaminated input; results are diagnostic only.\n";
        {
            gpes::bdg::BdGOperator2D op(grid, psi.vec(), params.physics(), diagnostics.chemical_potential);
            report << "goldstone_relative_residual="
                   << gpes::bdg::analysis::goldstone_relative_residual(op, psi.vec()) << '\n';
        }
        attempted = true;
        status("solving");
        gpes::bdg::BdGSolver2D solver(grid, params.physics(), cfg);
        auto result = solver.solve(psi, diagnostics.chemical_potential);
        // Keep the returned order: sorting non-finite candidates would violate
        // the ordering requirements, and mode indices must match validation.csv.
        int valid_count = 0;
        std::ostringstream validation;
        validation << "mode,residual_passed\n";
        std::vector<std::complex<double>> frequencies;
        for (std::size_t m = 0; m < result.modes.size(); ++m) {
            const auto& mode = result.modes[m];
            const bool valid = valid_mode(mode, options.residual_tolerance);
            valid_count += valid ? 1 : 0;
            validation << m << ',' << (valid ? 1 : 0) << '\n';
            if (!std::isfinite(mode.omega.real()) || !std::isfinite(mode.omega.imag()) ||
                !mode.u.allFinite() || !mode.v.allFinite()) continue;
            frequencies.push_back(mode.omega);
            const auto fluctuation = gpes::bdg::analysis::density_fluctuation(psi.vec(), mode.u, mode.v);
            std::ostringstream name;
            name << "mode_" << std::setw(4) << std::setfill('0') << m << ".bmp";
            gpes::image::save_scalar_bmp(fluctuation.cwiseAbs(), params.nx, params.ny,
                output / "images/bdg_modes" / name.str(), true);
        }
        const bool eigensolver_ok = result.info == Spectra::CompInfo::Successful &&
            result.nconv >= cfg.num_modes && valid_count >= cfg.num_modes &&
            valid_count == static_cast<int>(result.modes.size()) && result.linear_failures == 0 &&
            std::isfinite(result.max_linear_residual) && result.max_linear_residual <= cfg.shift_invert_linear_tol;
        report << "converged_modes=" << result.nconv << '\n'
            << "returned_modes=" << result.modes.size() << '\n'
            << "residual_valid_modes=" << valid_count << '\n'
            << "spectra_info=" << static_cast<int>(result.info) << '\n'
            << "arnoldi_iterations=" << result.niter << '\n'
            << "operator_applications=" << result.nops << '\n'
            << "linear_iterations=" << result.linear_iterations << '\n'
            << "linear_failures=" << result.linear_failures << '\n'
            << "max_linear_residual=" << result.max_linear_residual << '\n'
            << "max_eigen_residual=" << result.max_eigen_residual << '\n'
            << "used_dense_direct=" << result.used_dense_direct << '\n'
            << "used_fourier_preconditioner=" << result.used_fourier_preconditioner << '\n'
            << "eigensolver_quality_passed=" << eigensolver_ok << '\n';
        status("saving");
        gpes::bdg::save_bdg_result(result, *grid, output, !options.no_mode_fields);
        write_text(output / "validation.csv", validation.str());
        if (frequencies.size() == result.modes.size()) {
            gpes::image::save_spectrum_svg(frequencies, output / "images/bdg_spectrum.svg");
        }
        const bool complete = stationary && eigensolver_ok;
        status(complete ? "complete" : (stationary ? "incomplete" : "diagnostic_only"));
        fs::remove(output / "RUNNING");
        write_text(output / (complete ? "COMPLETED" : "INCOMPLETE"), "bdg\n");
        std::cout << "BdG: " << result.nconv << '/' << cfg.num_modes << " modes converged; "
                  << valid_count << " pass the direct residual test.\n";
        if (!complete) std::cerr << "BdG is incomplete or diagnostic only. Saved candidates are NOT a validated spectrum.\n";
        return complete ? 0 : 2;
    } catch (const std::exception& error) {
        try { status(attempted ? "failed" : "rejected_input"); }
        catch (const std::exception&) {}
        record_failure(output, error.what());
        std::cerr << "BdG error for " << input << ": " << error.what() << '\n';
        return 1;
    }
}

} // namespace

int main(int argc, char** argv) {
    using namespace workflow;
    Options options;
    auto& cfg = options.solver;
    cfg.num_modes = 24;
    cfg.ncv = 64;
    cfg.max_iter = 1500;
    cfg.tol = 1e-7;
    cfg.use_shift_invert = true;
    cfg.shift_invert_linear_max_iter = 1000;
    cfg.shift_invert_linear_restart = 48;
    cfg.shift_invert_dense_max_dim = 0;
    double shift = 0.05, shift_imag = 0.0;
    std::string input, output;
    bool recursive = false;
    CLI::App cli{"BdG spectra from saved quasi-2D states; no imaginary-time evolution"};
    cli.add_option("--input-dir", input, "Ground-state run directory, or sweep root with --recursive")
        ->required()->check(CLI::ExistingDirectory);
    cli.add_option("--output-dir", output, "New attempt directory; default: INPUT/bdg_attempts/attempt_<unique>");
    cli.add_flag("--recursive", recursive, "Process all saved states below input-dir sequentially");
    cli.add_option("--modes,--bdg-modes", cfg.num_modes, "Eigenpairs near the shift (including negative/phase branches)")
        ->check(CLI::PositiveNumber)->capture_default_str();
    cli.add_option("--ncv,--bdg-ncv", cfg.ncv, "Arnoldi subspace; effective size is at least 2*modes+1")
        ->check(CLI::PositiveNumber)->capture_default_str();
    cli.add_option("--max-iter,--bdg-max-iter", cfg.max_iter, "Maximum Arnoldi iterations")
        ->check(CLI::PositiveNumber)->capture_default_str();
    cli.add_option("--tol,--bdg-tol", cfg.tol, "Outer eigensolver tolerance")->capture_default_str();
    cli.add_option("--residual-tol,--bdg-residual-tol", options.residual_tolerance, "Direct eigenpair residual threshold")
        ->capture_default_str();
    cli.add_option("--shift,--bdg-shift", shift, "Real spectral target in trap-frequency units")->capture_default_str();
    cli.add_option("--shift-imag", shift_imag, "Imaginary part of the spectral target")->capture_default_str();
    cli.add_option("--linear-max-iter,--bdg-linear-max-iter", cfg.shift_invert_linear_max_iter, "Maximum inner GMRES iterations")
        ->check(CLI::PositiveNumber)->capture_default_str();
    cli.add_option("--linear-restart,--bdg-linear-restart", cfg.shift_invert_linear_restart, "GMRES restart length")
        ->check(CLI::PositiveNumber)->capture_default_str();
    cli.add_option("--linear-tol,--bdg-linear-tol", cfg.shift_invert_linear_tol, "Inner GMRES relative tolerance")
        ->capture_default_str();
    cli.add_option("--dense-max-dim,--bdg-dense-max-dim", cfg.shift_invert_dense_max_dim,
        "Dense inverse only when 2*nx*ny <= threshold; 0 disables; tiny-grid testing only")
        ->check(CLI::Range(0, 4096))->capture_default_str();
    cli.add_option("--stationary-residual-tol", options.stationary_tolerance, "Recomputed input stationarity threshold")
        ->capture_default_str();
    cli.add_option("--boundary-density-tol", options.boundary_tolerance, "Maximum input boundary / peak density")
        ->capture_default_str();
    cli.add_flag("--allow-unconverged", options.allow_unconverged, "Permit failed input-quality checks for diagnostics only");
    cli.add_flag("--no-mode-fields", options.no_mode_fields, "Skip full u,v CSV files; retain eigenvalues and images");
    try { cli.parse(argc, argv); }
    catch (const CLI::ParseError& e) { return cli.exit(e); }
    try {
        positive(cfg.tol, "tol");
        positive(cfg.shift_invert_linear_tol, "linear-tol");
        positive(options.residual_tolerance, "residual-tol");
        positive(options.stationary_tolerance, "stationary-residual-tol");
        positive(options.boundary_tolerance, "boundary-density-tol");
        if (!std::isfinite(shift) || !std::isfinite(shift_imag)) throw std::invalid_argument("Shift must be finite");
        if (cfg.num_modes > (std::numeric_limits<int>::max() - 1) / 2) throw std::invalid_argument("Too many requested modes");
        cfg.shift = {shift, shift_imag};
        std::vector<fs::path> states;
        if (recursive) {
            for (const auto& csv : gpes::io::findWavefunctionCsvFiles(input)) {
                if (csv.parent_path().filename() == "wavefunction") states.push_back(csv.parent_path().parent_path());
            }
            std::sort(states.begin(), states.end());
        } else {
            if (!fs::is_regular_file(fs::path(input) / "wavefunction/wavefunction.csv")) {
                throw std::runtime_error("No saved wavefunction in input-dir; use --recursive for a sweep directory");
            }
            states.emplace_back(input);
        }
        if (states.empty()) throw std::runtime_error("No saved wavefunctions found");
        const auto destination = new_output(output, fs::path(input) / "bdg_attempts", "attempt");
        std::cout << "Processing " << states.size() << " state(s) sequentially. Output: " << destination << std::endl;
        int exit_code = 0;
        std::size_t completed = 0;
        for (std::size_t i = 0; i < states.size(); ++i) {
            auto run = destination;
            if (recursive) {
                std::ostringstream name;
                name << "state_" << std::setw(4) << std::setfill('0') << i + 1;
                run /= name.str();
                fs::create_directory(run);
            }
            const int result = calculate(states[i], run, options);
            exit_code = std::max(exit_code, result);
            if (result == 0) ++completed;
        }
        if (recursive) {
            write_text(destination / "batch_status.txt", "total=" + std::to_string(states.size()) +
                "\ncompleted=" + std::to_string(completed) + "\nexit_code=" + std::to_string(exit_code) + '\n');
        }
        return exit_code;
    } catch (const std::exception& error) {
        std::cerr << "BdG error: " << error.what() << '\n';
        return 1;
    }
}
