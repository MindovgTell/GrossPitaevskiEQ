#include "workflow.hpp"
#include "Core/simple_image.hpp"
#include "Solvers/SplitStep.hpp"
#include "Solvers/TDSESolver.hpp"

namespace {

using Peak = workflow::Wavefunction::GaussianComponent;

std::vector<Peak> read_modulation_file(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot read modulation file: " + path);
    std::vector<Peak> peaks;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        std::istringstream row(line);
        row >> std::ws;
        if (row.eof() || row.peek() == '#') continue;
        Peak peak{};
        std::string extra;
        if (!(row >> peak[0] >> peak[1] >> peak[2] >> peak[3] >> peak[4]) ||
            (row >> extra && extra[0] != '#')) {
            throw std::invalid_argument("Invalid modulation row " + std::to_string(line_number) +
                " in " + path + "; expected: x y sigma_x sigma_y relative_amplitude");
        }
        peaks.push_back(peak);
    }
    if (!input.eof() || peaks.empty()) {
        throw std::invalid_argument("Modulation file is unreadable or has no peaks: " + path);
    }
    return peaks;
}

std::string serialize_modulations(const std::vector<Peak>& peaks) {
    std::ostringstream out;
    out << "# x y sigma_x sigma_y relative_amplitude\n" << std::setprecision(17);
    for (const auto& peak : peaks) {
        out << peak[0] << ' ' << peak[1] << ' ' << peak[2] << ' '
            << peak[3] << ' ' << peak[4] << '\n';
    }
    return out.str();
}

} // namespace

int main(int argc, char** argv) {
    using namespace workflow;
    Parameters params;
    gpes::SimConfig sim;
    std::size_t snapshot_interval = 5000;
    double seed_noise = 1e-3;
    unsigned int random_seed = 12345;
    double boundary_tolerance = 1e-8;
    double lambda = 1.0;
    double background_width_x = 2.5, background_width_y = 2.5;
    std::string output, resume, modulation_file;
    CLI::App cli{"One quasi-2D imaginary-time ground-state calculation (no BdG)"};
    auto* resume_option = cli.add_option("--resume", resume, "Load a run's parameters and checkpoint; write to a NEW directory")
        ->check(CLI::ExistingDirectory);
    const auto physical_option = [&](const std::string& name, auto& value, const std::string& help) {
        return cli.add_option(name, value, help)->capture_default_str()->excludes(resume_option);
    };
    physical_option("--a-s", params.a_s, "Scattering length in reference oscillator-length units");
    physical_option("--a-dd", params.a_dd, "Dipolar length in the same units");
    physical_option("--num-particles", params.num_particles, "Particle number N");
    physical_option("--w-x", params.omega_x, "Dimensionless x trap frequency");
    auto* wy_option = physical_option("--w-y", params.omega_y, "Dimensionless y trap frequency; alternative to --lambda");
    physical_option("--w-z", params.omega_z, "Frozen-axis z frequency");
    auto* lambda_option = physical_option("--lambda", lambda, "Aspect ratio omega_x / omega_y; sets omega_y=omega_x/lambda");
    lambda_option->excludes(wy_option);
    physical_option("--rotation-frequency", params.rotation_frequency, "Omega about z in reference-frequency units; |Omega| < min(w-x,w-y)");
    physical_option("--nx", params.nx, "Fourier grid points along x");
    physical_option("--ny", params.ny, "Fourier grid points along y");
    physical_option("--x-half-length", params.x_half_length, "Domain [-Lx,Lx)");
    physical_option("--y-half-length", params.y_half_length, "Domain [-Ly,Ly)");
    physical_option("--background-width-x", background_width_x, "Fresh Gaussian background width along x");
    physical_option("--background-width-y", background_width_y, "Fresh Gaussian background width along y");
    cli.add_option("--modulation-file", modulation_file,
        "Fresh seed peaks: x y sigma_x sigma_y amplitude_relative_to_background_peak")
        ->check(CLI::ExistingFile)->excludes(resume_option);
    cli.add_option("--dt", sim.dt, "Initial imaginary-time step")->capture_default_str();
    cli.add_option("--min-ground-dt", sim.minimum_ground_dt, "Smallest adaptive imaginary-time step")->capture_default_str();
    cli.add_option("--max-ground-steps", sim.max_ground_steps, "Maximum steps in this invocation")
        ->check(CLI::PositiveNumber)->capture_default_str();
    cli.add_option("--convergence-check-interval", sim.convergence_check_interval, "Steps between convergence checks")
        ->check(CLI::PositiveNumber)->capture_default_str();
    cli.add_option("--convergence-checks-required", sim.convergence_checks_required, "Consecutive passing checks")
        ->check(CLI::PositiveNumber)->capture_default_str();
    cli.add_option("--relative-energy-tol", sim.relative_energy_tolerance, "Relative energy stability threshold")
        ->capture_default_str();
    cli.add_option("--absolute-energy-tol", sim.absolute_energy_tolerance, "Absolute energy floor near zero")
        ->capture_default_str();
    cli.add_option("--stationary-residual-tol", sim.stationary_residual_tolerance, "Maximum relative eGPE residual")
        ->capture_default_str();
    cli.add_option("--boundary-density-tol", boundary_tolerance, "Maximum boundary / peak density for success")
        ->capture_default_str();
    cli.add_option("--snapshot-interval", snapshot_interval, "Steps between atomic checkpoints")
        ->check(CLI::PositiveNumber)->capture_default_str();
    cli.add_option("--seed-noise", seed_noise, "Relative complex noise for a fresh Gaussian seed")
        ->capture_default_str()->excludes(resume_option);
    cli.add_option("--random-seed", random_seed, "Random seed for a fresh start")
        ->capture_default_str()->excludes(resume_option);
    cli.add_option("--output-dir,--wavefunction-out", output, "New run directory; default: $SCRATCH/GPES/groundstate_<unique>");
    try { cli.parse(argc, argv); }
    catch (const CLI::ParseError& e) { return cli.exit(e); }

    fs::path run;
    try {
        if (!resume.empty()) params = read_parameters(resume, true);
        if (lambda_option->count()) {
            positive(lambda, "lambda");
            params.omega_y = params.omega_x / lambda;
        }
        validate(params);
        positive(sim.dt, "dt");
        positive(sim.minimum_ground_dt, "min-ground-dt");
        positive(sim.relative_energy_tolerance, "relative-energy-tol");
        positive(sim.absolute_energy_tolerance, "absolute-energy-tol", true);
        positive(sim.stationary_residual_tolerance, "stationary-residual-tol");
        positive(boundary_tolerance, "boundary-density-tol");
        positive(seed_noise, "seed-noise", true);
        positive(background_width_x, "background-width-x");
        positive(background_width_y, "background-width-y");
        if (sim.minimum_ground_dt > sim.dt) throw std::invalid_argument("min-ground-dt cannot exceed dt");
        model_warning(params);
        const auto grid = params.grid();
        Wavefunction psi(grid);
        std::vector<Peak> modulation_peaks;
        if (!resume.empty()) {
            gpes::io::load_wavefunction_csv(psi, fs::path(resume) / "wavefunction/wavefunction.csv", true);
            check_normalization(psi, params);
            const auto saved_modulations = fs::path(resume) / "wavefunction/modulations.tsv";
            if (fs::exists(saved_modulations)) {
                modulation_peaks = read_modulation_file(saved_modulations.string());
            }
        } else {
            psi.set_state_Gauss(0.0, 0.0, background_width_x, background_width_y,
                static_cast<unsigned int>(params.num_particles));
            if (!modulation_file.empty()) {
                modulation_peaks = read_modulation_file(modulation_file);
                psi.add_gaussian_modulations(modulation_peaks);
            }
            if (seed_noise > 0.0) gpes::add_complex_random_noise(psi, seed_noise, random_seed);
        }
        run = new_output(output, default_root(), "groundstate");
        write_parameters(run, params);
        write_text(run / "RUNNING", "groundstate\n");
        if (!resume.empty()) write_text(run / "source.txt", fs::absolute(resume).string() + '\n');
        if (!modulation_peaks.empty()) {
            write_text(run / "wavefunction/modulations.tsv", serialize_modulations(modulation_peaks));
        }
        if (resume.empty()) {
            std::ostringstream seed;
            seed << std::setprecision(17)
                 << "background_width_x=" << background_width_x << '\n'
                 << "background_width_y=" << background_width_y << '\n'
                 << "modulation_count=" << modulation_peaks.size() << '\n'
                 << "seed_noise=" << seed_noise << '\n'
                 << "random_seed=" << random_seed << '\n';
            write_text(run / "wavefunction/seed.txt", seed.str());
        } else {
            const auto saved_seed = fs::path(resume) / "wavefunction/seed.txt";
            if (fs::exists(saved_seed)) fs::copy_file(saved_seed, run / "wavefunction/seed.txt");
        }
        fs::create_directories(run / "log");
        gpes::log::Log::getInstance().setLogDirectory(run / "log");
        std::cout << "Ground-state directory: " << run << "\nGrid: " << params.nx << 'x' << params.ny
                  << "\nOmega: " << params.rotation_frequency << ", lambda=omega_x/omega_y: "
                  << params.omega_x / params.omega_y << std::endl;
        const auto checkpoint = run / "wavefunction/wavefunction.csv";
        gpes::io::save_wavefunction_csv(psi, checkpoint);
        gpes::image::save_scalar_bmp(psi.prob(), params.nx, params.ny, run / "images/initial_density.bmp", true);
        gpes::image::save_phase_bmp(psi.vec(), params.nx, params.ny, run / "images/initial_phase.bmp");

        gpes::solvers::SplitStep<Dim> stepper(grid, psi, params.physics(), sim, params.rotation_frequency);
        gpes::TDSESolver<Dim, decltype(stepper)> solver(grid, stepper, psi, sim);
        solver.run(checkpoint, snapshot_interval);
        const auto& result = solver.ground_state_result();
        const auto& diagnostics = stepper.diagnostic_history().back();
        check_normalization(psi, params);
        const double boundary = boundary_ratio(psi);
        const bool complete = result.converged && std::isfinite(boundary) && boundary <= boundary_tolerance;

        std::ostringstream history;
        history << "step,energy\n" << std::setprecision(17);
        for (std::size_t i = 0; i < stepper.energies().size(); ++i) {
            history << result.diagnostic_steps.at(i) << ',' << stepper.energies()[i] << '\n';
        }
        write_text(run / "wavefunction/energy_history.csv", history.str());
        std::ostringstream diagnostic_history;
        diagnostic_history << "step,rotating_energy,lab_energy,angular_momentum,angular_momentum_per_particle,chemical_potential,relative_residual,scaled_residual\n"
                           << std::setprecision(17);
        for (std::size_t i = 0; i < stepper.diagnostic_history().size(); ++i) {
            const auto& d = stepper.diagnostic_history()[i];
            diagnostic_history << result.diagnostic_steps.at(i) << ',' << d.energy << ',' << d.lab_energy << ','
                << d.angular_momentum << ',' << d.angular_momentum_per_particle << ',' << d.chemical_potential << ','
                << d.relative_residual << ',' << d.scaled_residual << '\n';
        }
        write_text(run / "wavefunction/diagnostics.csv", diagnostic_history.str());
        gpes::image::save_line_svg(stepper.energies(), run / "images/energy_history.svg", "Imaginary-time convergence", "E - Omega <Lz>");
        gpes::image::save_scalar_bmp(psi.prob(), params.nx, params.ny, run / "images/ground_state_density.bmp", true);
        gpes::image::save_phase_bmp(psi.vec(), params.nx, params.ny, run / "images/ground_state_phase.bmp");

        std::ostringstream report;
        report << std::setprecision(17) << std::boolalpha
            << "status=" << (complete ? "complete" : "not_converged") << '\n'
            << "ground_state_converged=" << result.converged << '\n'
            << "boundary_check_passed=" << (boundary <= boundary_tolerance) << '\n'
            << "steps_completed=" << result.steps_completed << '\n'
            << "convergence_checks=" << result.convergence_checks << '\n'
            << "energy=" << result.energy << '\n'
            << "energy_frame=rotating\n"
            << "rotating_energy=" << diagnostics.energy << '\n'
            << "lab_energy=" << diagnostics.lab_energy << '\n'
            << "angular_momentum=" << diagnostics.angular_momentum << '\n'
            << "angular_momentum_per_particle=" << diagnostics.angular_momentum_per_particle << '\n'
            << "rotation_frequency=" << params.rotation_frequency << '\n'
            << "lambda_omega_x_over_omega_y=" << params.omega_x / params.omega_y << '\n'
            << "chemical_potential=" << result.chemical_potential << '\n'
            << "absolute_energy_change=" << result.absolute_energy_change << '\n'
            << "relative_energy_change=" << result.relative_energy_change << '\n'
            << "stationary_relative_residual=" << result.stationary_residual << '\n'
            << "stationary_scaled_residual=" << diagnostics.scaled_residual << '\n'
            << "particle_number=" << particle_number(psi) << '\n'
            << "boundary_to_peak_density=" << boundary << '\n'
            << "initial_dt=" << sim.dt << '\n'
            << "final_dt=" << result.final_time_step << '\n'
            << "minimum_ground_dt=" << sim.minimum_ground_dt << '\n'
            << "relative_energy_tolerance=" << sim.relative_energy_tolerance << '\n'
            << "absolute_energy_tolerance=" << sim.absolute_energy_tolerance << '\n'
            << "stationary_residual_tolerance=" << sim.stationary_residual_tolerance << '\n'
            << "boundary_density_tolerance=" << boundary_tolerance << '\n'
            << "max_ground_steps=" << sim.max_ground_steps << '\n'
            << "convergence_check_interval=" << sim.convergence_check_interval << '\n'
            << "convergence_checks_required=" << sim.convergence_checks_required << '\n'
            << "snapshot_interval=" << snapshot_interval << '\n'
            << "modulation_count=" << modulation_peaks.size() << '\n'
            << "resumed=" << !resume.empty() << '\n';
        if (resume.empty()) {
            report << "background_width_x=" << background_width_x << '\n'
                << "background_width_y=" << background_width_y << '\n'
                << "seed_noise=" << seed_noise << '\n'
                << "random_seed=" << random_seed << '\n';
        } else {
            report << "seed_metadata="
                   << (fs::exists(run / "wavefunction/seed.txt") ? "wavefunction/seed.txt" : "unavailable")
                   << '\n';
        }
        write_text(run / "ground_state.txt", report.str());
        fs::remove(run / "RUNNING");
        write_text(run / (complete ? "COMPLETED" : "NOT_CONVERGED"), "groundstate\n");
        std::cout << "Rotating energy: " << result.energy << ", lab energy: " << diagnostics.lab_energy
                  << ", Lz/N: " << diagnostics.angular_momentum_per_particle << ", stationary residual: " << result.stationary_residual
                  << ", boundary / peak: " << boundary << '\n';
        if (!complete) std::cerr << "Ground-state quality checks failed; checkpoint and diagnostics were saved.\n";
        return complete ? 0 : 2;
    } catch (const std::exception& error) {
        record_failure(run, error.what());
        std::cerr << "Ground-state error: " << error.what() << '\n';
        return 1;
    }
}
