#include "cli/cli.hpp"
#include "gpes.hpp"
#include "Core/simple_image.hpp"
#include "Core/utility.hpp"
#include "Solvers/bdg_spectral_analysis.hpp"
#if GPES_HAS_PLOTTING
#include "Core/plot.hpp"
#endif

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>

namespace {

using Grid2D = gpes::Grid<gpes::Dimension::Two>;
using Wavefunction2D = gpes::WaveFunction<gpes::Dimension::Two>;

void write_json_number(std::ostream& out, double value) {
    if (std::isfinite(value)) out << value;
    else out << "null";
}

double boundary_to_peak_density(const Wavefunction2D& psi) {
    const int nx = static_cast<int>(psi.grid()->size_x());
    const int ny = static_cast<int>(psi.grid()->size_y());
    double peak = 0.0;
    double boundary = 0.0;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            const double density = std::norm(psi.vec()(i * ny + j));
            peak = std::max(peak, density);
            if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1) {
                boundary = std::max(boundary, density);
            }
        }
    }
    return peak > 0.0 ? boundary / peak : std::numeric_limits<double>::infinity();
}

double particle_number(const Wavefunction2D& psi) {
    return psi.vec().squaredNorm() * psi.grid()->step_x() * psi.grid()->step_y();
}

std::string checkpoint_signature(const CliConfig& cfg) {
    std::ostringstream out;
    out << std::setprecision(17)
        << "nx=" << cfg.nx << '\n'
        << "ny=" << cfg.ny << '\n'
        << "x_half_length=" << cfg.x_half_length << '\n'
        << "y_half_length=" << cfg.y_half_length << '\n'
        << "a_s=" << cfg.a_s << '\n'
        << "a_dd=" << cfg.a_dd << '\n'
        << "num_particles=" << cfg.num_particles << '\n'
        << "omega_x=" << cfg.w_x << '\n'
        << "omega_y=" << cfg.w_y << '\n'
        << "omega_z=" << cfg.w_z << '\n';
    return out.str();
}

std::string read_text_file(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot read checkpoint manifest: " + path.string());
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

void write_text_file(const std::filesystem::path& path, const std::string& text) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write checkpoint manifest: " + path.string());
    out << text;
}

void save_run_metadata(
    const std::filesystem::path& path,
    const CliConfig& cfg,
    const gpes::GroundStateResult& ground,
    double norm,
    double boundary_ratio,
    const gpes::bdg::BdGResult* bdg,
    double goldstone_residual,
    bool bdg_valid)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write run metadata: " + path.string());
    out << std::setprecision(17)
        << "{\n"
        << "  \"a_s\": " << cfg.a_s << ",\n"
        << "  \"a_dd\": " << cfg.a_dd << ",\n"
        << "  \"num_particles\": " << cfg.num_particles << ",\n"
        << "  \"omega_x\": " << cfg.w_x << ",\n"
        << "  \"omega_y\": " << cfg.w_y << ",\n"
        << "  \"omega_z\": " << cfg.w_z << ",\n"
        << "  \"omega_y_over_omega_x\": " << cfg.w_y / cfg.w_x << ",\n"
        << "  \"omega_y_over_omega_z\": " << cfg.w_y / cfg.w_z << ",\n"
        << "  \"nx\": " << cfg.nx << ",\n"
        << "  \"ny\": " << cfg.ny << ",\n"
        << "  \"dt\": " << cfg.dt << ",\n"
        << "  \"final_imaginary_time_dt\": " << ground.final_time_step << ",\n"
        << "  \"ground_state_converged\": " << (ground.converged ? "true" : "false") << ",\n"
        << "  \"ground_state_steps\": " << ground.steps_completed << ",\n"
        << "  \"energy\": "; write_json_number(out, ground.energy); out << ",\n"
        << "  \"chemical_potential\": "; write_json_number(out, ground.chemical_potential); out << ",\n"
        << "  \"relative_energy_change\": "; write_json_number(out, ground.relative_energy_change); out << ",\n"
        << "  \"stationary_relative_residual\": "; write_json_number(out, ground.stationary_residual); out << ",\n"
        << "  \"particle_number\": "; write_json_number(out, norm); out << ",\n"
        << "  \"boundary_to_peak_density\": "; write_json_number(out, boundary_ratio); out << ",\n"
        << "  \"bdg_requested\": " << (!cfg.skip_bdg ? "true" : "false") << ",\n"
        << "  \"bdg_attempted\": " << (bdg ? "true" : "false") << ",\n"
        << "  \"bdg_completed\": " << (bdg_valid ? "true" : "false");
    if (bdg) {
        out << ",\n"
            << "  \"bdg_converged_modes\": " << bdg->nconv << ",\n"
            << "  \"bdg_linear_failures\": " << bdg->linear_failures << ",\n"
            << "  \"bdg_max_linear_residual\": "; write_json_number(out, bdg->max_linear_residual); out << ",\n"
            << "  \"bdg_max_eigen_residual\": "; write_json_number(out, bdg->max_eigen_residual); out << ",\n"
            << "  \"bdg_goldstone_residual\": "; write_json_number(out, goldstone_residual);
    }
    out << "\n}\n";
}

void save_mode_images(
    const gpes::bdg::BdGResult& result,
    const Wavefunction2D& psi,
    const std::filesystem::path& directory)
{
    const int nx = static_cast<int>(psi.grid()->size_x());
    const int ny = static_cast<int>(psi.grid()->size_y());
    for (std::size_t m = 0; m < result.modes.size(); ++m) {
        const auto fluctuation = gpes::bdg::analysis::density_fluctuation(
            psi.vec(), result.modes[m].u, result.modes[m].v);
        Eigen::VectorXd magnitude(fluctuation.size());
        for (Eigen::Index i = 0; i < fluctuation.size(); ++i) {
            magnitude(i) = std::abs(fluctuation(i));
        }
        std::ostringstream filename;
        filename << "mode_" << std::setw(4) << std::setfill('0') << m << ".bmp";
        gpes::image::save_scalar_bmp(
            magnitude, nx, ny,
            directory / filename.str(), true);
    }
}

} // namespace

int main(int argc, char **argv) {
    try {
        auto parsed = parse_cli(argc, argv);
        auto phys = parsed.phys;
        auto sim = parsed.sim;
        sim.ground_state = true;
        const auto& cfg = parsed.cli;
        print_parsed_params(parsed);

        const std::filesystem::path run_dir = cfg.wavefunction_out;
        const std::filesystem::path wavefunction_dir = run_dir / "wavefunction";
        const std::filesystem::path bdg_dir = run_dir / "bdg";
        const std::filesystem::path image_dir = run_dir / "images";
        std::filesystem::create_directories(wavefunction_dir);
        std::filesystem::create_directories(bdg_dir);
        std::filesystem::create_directories(image_dir / "bdg_modes");

        std::cout << "Run directory: " << run_dir << "\n";
        if (cfg.w_z <= std::max(cfg.w_x, cfg.w_y)) {
            std::cerr << "Model warning: omega_z=" << cfg.w_z
                      << " is not greater than both in-plane frequencies; "
                         "the frozen-z quasi-2D approximation is not controlled.\n";
        }
        auto grid = std::make_shared<Grid2D>(
            cfg.nx, cfg.ny, -cfg.x_half_length, -cfg.y_half_length,
            cfg.w_x, cfg.w_y, cfg.w_z);

        const std::filesystem::path wavefunction_path = wavefunction_dir / "wavefunction.csv";
        const std::filesystem::path checkpoint_manifest = wavefunction_dir / "parameters.txt";
        const std::string expected_signature = checkpoint_signature(cfg);
        Wavefunction2D psi(grid);
        if (!cfg.fresh_start && std::filesystem::exists(wavefunction_path)) {
            if (!std::filesystem::exists(checkpoint_manifest)) {
                throw std::runtime_error(
                    "Checkpoint has no parameter manifest; use --fresh to replace it safely.");
            }
            if (read_text_file(checkpoint_manifest) != expected_signature) {
                throw std::runtime_error(
                    "Checkpoint parameters differ from this run; use --fresh or a different output directory.");
            }
            gpes::io::load_wavefunction_csv(psi, wavefunction_path);
            const double checkpoint_particles = particle_number(psi);
            const double relative_norm_error = std::abs(checkpoint_particles - phys.num_of_prt)
                / phys.num_of_prt;
            if (relative_norm_error > 1e-8) {
                throw std::runtime_error(
                    "Checkpoint particle number does not match --num-particles; use --fresh or the original parameters.");
            }
            std::cout << "Resuming from checkpoint: " << wavefunction_path << "\n";
        } else {
            write_text_file(checkpoint_manifest, expected_signature);
            psi.set_state_Gauss(0.0, 0.0, 2.5, 2.5,
                static_cast<unsigned int>(std::llround(phys.num_of_prt)));
            if (cfg.seed_noise > 0.0) {
                gpes::add_complex_random_noise(psi, cfg.seed_noise, cfg.random_seed);
            }
        }
        gpes::image::save_scalar_bmp(
            psi.prob(), cfg.nx, cfg.ny, image_dir / "initial_density.bmp", true);

        gpes::solvers::SplitStep<gpes::Dimension::Two> stepper(grid, psi, phys, sim);
        gpes::TDSESolver<gpes::Dimension::Two, decltype(stepper)> solver(grid, stepper, psi, sim);
        solver.run(wavefunction_path, cfg.snapshot_interval);
        const auto ground = solver.ground_state_result();

        gpes::io::save_state_energy_history_csv(
            stepper.energies(), wavefunction_path, "energy_history.csv",
            cfg.convergence_check_interval);
        gpes::image::save_line_svg(
            stepper.energies(), image_dir / "energy_history.svg",
            "Imaginary-time convergence", "energy");
        gpes::image::save_scalar_bmp(
            psi.prob(), cfg.nx, cfg.ny, image_dir / "ground_state_density.bmp", true);

        const double norm = particle_number(psi);
        const double boundary_ratio = boundary_to_peak_density(psi);
        if (boundary_ratio > 1e-8) {
            std::cerr << "Warning: boundary/peak density ratio is " << boundary_ratio
                      << "; enlarge the spatial domain before treating this run as converged.\n";
        }

#if GPES_HAS_PLOTTING
        gpes::plot::save_heatmap(psi, image_dir / "ground_state_density.png", "Ground-state density");
#endif

        save_run_metadata(run_dir / "run_metadata.json", cfg, ground, norm,
            boundary_ratio, nullptr, std::numeric_limits<double>::quiet_NaN(), false);

        if (!ground.converged && !cfg.allow_unconverged) {
            throw std::runtime_error(
                "Ground state failed the convergence gate; BdG was not run. "
                "Inspect run_metadata.json or use --allow-unconverged only for diagnostics.");
        }

        if (!cfg.skip_bdg) {
            gpes::bdg::BdGConfig bdg_config;
            bdg_config.num_modes = cfg.bdg_modes;
            bdg_config.ncv = cfg.bdg_ncv;
            bdg_config.max_iter = cfg.bdg_max_iter;
            bdg_config.tol = cfg.bdg_tolerance;
            bdg_config.use_shift_invert = true;
            bdg_config.shift = {cfg.bdg_shift, 0.0};
            bdg_config.shift_invert_linear_tol = cfg.bdg_linear_tolerance;
            bdg_config.shift_invert_linear_max_iter = cfg.bdg_linear_max_iter;
            bdg_config.shift_invert_linear_restart = cfg.bdg_linear_restart;
            bdg_config.shift_invert_dense_max_dim = cfg.bdg_dense_max_dimension;
            bdg_config.use_fourier_preconditioner = true;
            bdg_config.preconditioner_include_dipolar = true;

            gpes::bdg::BdGSolver2D bdg_solver(grid, phys, bdg_config);
            auto bdg_result = bdg_solver.solve(psi, ground.chemical_potential);
            std::sort(bdg_result.modes.begin(), bdg_result.modes.end(),
                [](const auto& lhs, const auto& rhs) {
                    if (lhs.omega.real() != rhs.omega.real()) return lhs.omega.real() < rhs.omega.real();
                    return lhs.omega.imag() < rhs.omega.imag();
                });

            gpes::bdg::BdGOperator2D validation_operator(
                grid, psi.vec(), phys, ground.chemical_potential);
            const double goldstone_residual =
                gpes::bdg::analysis::goldstone_relative_residual(validation_operator, psi.vec());
            if (bdg_result.nconv < cfg.bdg_modes || bdg_result.linear_failures != 0 ||
                bdg_result.max_eigen_residual > cfg.bdg_residual_tolerance) {
                save_run_metadata(run_dir / "run_metadata.json", cfg, ground, norm,
                    boundary_ratio, &bdg_result, goldstone_residual, false);
                throw std::runtime_error(
                    "BdG failed its convergence gate; unconverged eigenpairs were not saved.");
            }

            gpes::bdg::save_bdg_result(
                bdg_result, grid, bdg_dir, cfg.save_mode_fields);
            std::vector<std::complex<double>> omega;
            omega.reserve(bdg_result.modes.size());
            for (const auto& mode : bdg_result.modes) omega.push_back(mode.omega);
            gpes::image::save_spectrum_svg(omega, image_dir / "bdg_spectrum.svg");
            save_mode_images(bdg_result, psi, image_dir / "bdg_modes");
#if GPES_HAS_PLOTTING
            gpes::plot::save_excitation_spectrum_plot(
                bdg_result, image_dir / "bdg_spectrum.png");
#endif
            save_run_metadata(run_dir / "run_metadata.json", cfg, ground, norm,
                boundary_ratio, &bdg_result, goldstone_residual, true);
        }

        std::ofstream(run_dir / "COMPLETED") << "ok\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return 1;
    }
}
