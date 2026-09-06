#include <cmath>
#include <cstdlib>
#include <limits>
#include <iomanip>
#include <iostream>
#include <filesystem>

#include "CLI11.hpp"

#include "cli.hpp"
#include "Log/Log.hpp"


namespace {

ParsedCli parse_cli_full(int argc, char **argv) {
    CLI::App app{"Gross-Pitaevskii Equation Solver"};

    bool cli_only = false;
    app.add_flag("--cli-only", cli_only,
        "Require all parameters to be provided via CLI (no code defaults)");

    CliConfig cfg;

    double a_s = cfg.a_s;
    double a_dd = cfg.a_dd;
    double num_particles = cfg.num_particles;
    double dt = cfg.dt;
    double duration = cfg.duration;
    double w_x = cfg.w_x;
    double w_y = cfg.w_y;
    double w_z = cfg.w_z;
    int nx = cfg.nx;
    int ny = cfg.ny;
    double x_half_length = cfg.x_half_length;
    double y_half_length = cfg.y_half_length;
    std::size_t max_ground_steps = cfg.max_ground_steps;
    std::size_t convergence_check_interval = cfg.convergence_check_interval;
    std::size_t convergence_checks_required = cfg.convergence_checks_required;
    double relative_energy_tolerance = cfg.relative_energy_tolerance;
    double stationary_residual_tolerance = cfg.stationary_residual_tolerance;
    double minimum_ground_dt = cfg.minimum_ground_dt;
    std::size_t snapshot_interval = cfg.snapshot_interval;
    double seed_noise = cfg.seed_noise;
    unsigned int random_seed = cfg.random_seed;
    bool skip_bdg = false;
    bool allow_unconverged = false;
    bool no_mode_fields = false;
    bool fresh_start = false;
    int bdg_modes = cfg.bdg_modes;
    int bdg_ncv = cfg.bdg_ncv;
    int bdg_max_iter = cfg.bdg_max_iter;
    double bdg_tolerance = cfg.bdg_tolerance;
    double bdg_residual_tolerance = cfg.bdg_residual_tolerance;
    double bdg_shift = cfg.bdg_shift;
    int bdg_linear_max_iter = cfg.bdg_linear_max_iter;
    int bdg_linear_restart = cfg.bdg_linear_restart;
    double bdg_linear_tolerance = cfg.bdg_linear_tolerance;
    int bdg_dense_max_dimension = cfg.bdg_dense_max_dimension;
    std::string log_dir;
    std::string wavefunction_out = cfg.wavefunction_out;

    auto *opt_a_s = app.add_option("--a-s", a_s, "Scattering length")->check(CLI::PositiveNumber);
    auto *opt_a_dd = app.add_option("--a-dd", a_dd, "Dipole-dipole interaction strength")
                         ->check(CLI::PositiveNumber);
    auto *opt_num_particles = app.add_option("--num-particles", num_particles, "Number of particles")
                                  ->check(CLI::PositiveNumber);
    auto *opt_dt = app.add_option("--dt", dt, "Time step")->check(CLI::PositiveNumber);
    auto *opt_duration = app.add_option("--duration", duration, "Simulation duration")
                             ->check(CLI::PositiveNumber);
    auto *opt_w_x = app.add_option("--w-x", w_x, "Trap frequency along x")->check(CLI::PositiveNumber);
    auto *opt_w_y = app.add_option("--w-y", w_y, "Trap frequency along y")->check(CLI::PositiveNumber);
    auto *opt_w_z = app.add_option("--w-z", w_z, "Trap frequency along z")->check(CLI::PositiveNumber);
    app.add_option("--nx", nx, "Grid points along x")->check(CLI::PositiveNumber);
    app.add_option("--ny", ny, "Grid points along y")->check(CLI::PositiveNumber);
    app.add_option("--x-half-length", x_half_length, "Half width of the x domain")->check(CLI::PositiveNumber);
    app.add_option("--y-half-length", y_half_length, "Half width of the y domain")->check(CLI::PositiveNumber);
    app.add_option("--max-ground-steps", max_ground_steps, "Maximum imaginary-time steps")->check(CLI::PositiveNumber);
    app.add_option("--convergence-check-interval", convergence_check_interval, "Steps between convergence checks")->check(CLI::PositiveNumber);
    app.add_option("--convergence-checks-required", convergence_checks_required, "Consecutive passing convergence checks")->check(CLI::PositiveNumber);
    app.add_option("--relative-energy-tol", relative_energy_tolerance, "Relative ground-state energy tolerance")->check(CLI::PositiveNumber);
    app.add_option("--stationary-residual-tol", stationary_residual_tolerance, "Relative stationary eGPE residual tolerance")->check(CLI::PositiveNumber);
    app.add_option("--min-ground-dt", minimum_ground_dt, "Smallest imaginary-time step used during residual refinement")->check(CLI::PositiveNumber);
    app.add_option("--snapshot-interval", snapshot_interval, "Steps between wavefunction checkpoints")->check(CLI::PositiveNumber);
    app.add_option("--seed-noise", seed_noise, "Relative complex noise used to break spatial symmetry")->check(CLI::NonNegativeNumber);
    app.add_option("--random-seed", random_seed, "Random seed for the initial state");
    app.add_flag("--skip-bdg", skip_bdg, "Calculate only the ground state");
    app.add_flag("--allow-unconverged", allow_unconverged, "Run BdG even when the stationary-state quality gate fails");
    app.add_flag("--no-mode-fields", no_mode_fields, "Do not write full u,v mode fields");
    app.add_flag("--fresh", fresh_start, "Ignore an existing checkpoint and start from a Gaussian");
    app.add_option("--bdg-modes", bdg_modes, "Number of BdG eigenpairs")->check(CLI::PositiveNumber);
    app.add_option("--bdg-ncv", bdg_ncv, "Arnoldi subspace size")->check(CLI::PositiveNumber);
    app.add_option("--bdg-max-iter", bdg_max_iter, "Maximum Arnoldi iterations")->check(CLI::PositiveNumber);
    app.add_option("--bdg-tol", bdg_tolerance, "BdG eigensolver tolerance")->check(CLI::PositiveNumber);
    app.add_option("--bdg-residual-tol", bdg_residual_tolerance, "Maximum accepted direct BdG eigenpair residual")->check(CLI::PositiveNumber);
    app.add_option("--bdg-shift", bdg_shift, "Real shift for shift-invert BdG");
    app.add_option("--bdg-linear-max-iter", bdg_linear_max_iter, "Maximum inner GMRES iterations")->check(CLI::PositiveNumber);
    app.add_option("--bdg-linear-restart", bdg_linear_restart, "Inner GMRES restart length")->check(CLI::PositiveNumber);
    app.add_option("--bdg-linear-tol", bdg_linear_tolerance, "Inner GMRES relative tolerance")->check(CLI::PositiveNumber);
    app.add_option("--bdg-dense-max-dim", bdg_dense_max_dimension, "Use dense shift-invert at or below this BdG dimension")->check(CLI::NonNegativeNumber);
    app.add_option("--log-dir", log_dir, "Log directory path");
    auto *opt_wavefunction_out = app.add_option("--wavefunction-out", wavefunction_out,
        "Run output directory");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        std::exit(app.exit(e));
    }

    if (nx < 8 || ny < 8) {
        throw CLI::ValidationError("--nx and --ny must each be at least 8");
    }
    if (minimum_ground_dt > dt) {
        throw CLI::ValidationError("--min-ground-dt cannot exceed --dt");
    }
    if (std::abs(num_particles - std::round(num_particles)) > 1e-10 ||
        num_particles > static_cast<double>(std::numeric_limits<unsigned int>::max())) {
        throw CLI::ValidationError("--num-particles must be an integer representable as unsigned int");
    }

    if (cli_only) {
        if (opt_a_s->count() == 0 || opt_a_dd->count() == 0 || opt_num_particles->count() == 0 ||
            opt_dt->count() == 0 || opt_duration->count() == 0 ||
            opt_w_x->count() == 0 || opt_w_y->count() == 0 || opt_w_z->count() == 0 ||
            opt_wavefunction_out->count() == 0) {
            throw CLI::ValidationError(
                "With --cli-only, all parameters must be provided on the command line.");
        }
    }

    cfg.a_s = a_s;
    cfg.a_dd = a_dd;
    cfg.num_particles = num_particles;
    cfg.dt = dt;
    cfg.duration = duration;
    cfg.w_x = w_x;
    cfg.w_y = w_y;
    cfg.w_z = w_z;
    cfg.nx = nx;
    cfg.ny = ny;
    cfg.x_half_length = x_half_length;
    cfg.y_half_length = y_half_length;
    cfg.max_ground_steps = max_ground_steps;
    cfg.convergence_check_interval = convergence_check_interval;
    cfg.convergence_checks_required = convergence_checks_required;
    cfg.relative_energy_tolerance = relative_energy_tolerance;
    cfg.stationary_residual_tolerance = stationary_residual_tolerance;
    cfg.minimum_ground_dt = minimum_ground_dt;
    cfg.snapshot_interval = snapshot_interval;
    cfg.seed_noise = seed_noise;
    cfg.random_seed = random_seed;
    cfg.skip_bdg = skip_bdg;
    cfg.allow_unconverged = allow_unconverged;
    cfg.save_mode_fields = !no_mode_fields;
    cfg.fresh_start = fresh_start;
    cfg.bdg_modes = bdg_modes;
    cfg.bdg_ncv = bdg_ncv;
    cfg.bdg_max_iter = bdg_max_iter;
    cfg.bdg_tolerance = bdg_tolerance;
    cfg.bdg_residual_tolerance = bdg_residual_tolerance;
    cfg.bdg_shift = bdg_shift;
    cfg.bdg_linear_max_iter = bdg_linear_max_iter;
    cfg.bdg_linear_restart = bdg_linear_restart;
    cfg.bdg_linear_tolerance = bdg_linear_tolerance;
    cfg.bdg_dense_max_dimension = bdg_dense_max_dimension;
    cfg.wavefunction_out = wavefunction_out;
    if (cfg.wavefunction_out.empty()) {
        const char* scratch = std::getenv("SCRATCH");
        cfg.wavefunction_out = scratch && *scratch
            ? std::string(scratch) + "/GPES/manual"
            : "runs/manual";
    }
    cfg.log_dir = log_dir.empty() ? cfg.wavefunction_out + "/log" : log_dir;
    gpes::PhysConfig phys{cfg.a_s, cfg.a_dd, cfg.num_particles};
    gpes::SimConfig sim{cfg.duration, cfg.dt,
        static_cast<std::size_t>(std::ceil(cfg.duration / cfg.dt))};
    sim.max_ground_steps = cfg.max_ground_steps;
    sim.convergence_check_interval = cfg.convergence_check_interval;
    sim.convergence_checks_required = cfg.convergence_checks_required;
    sim.relative_energy_tolerance = cfg.relative_energy_tolerance;
    sim.stationary_residual_tolerance = cfg.stationary_residual_tolerance;
    sim.minimum_ground_dt = cfg.minimum_ground_dt;

    ParsedCli parsed;
    parsed.cli_only = cli_only;
    parsed.cli = cfg;
    parsed.phys = phys;
    parsed.sim = sim;
    return parsed;

}

void print_row(std::ostream &os, const std::string &name, const std::string &value,
    const std::string &desc) {
    os << std::left << std::setw(16) << name << " " << std::setw(16) << value << " " << desc << "\n";
}

} // namespace

ParsedCli parse_cli(int argc, char **argv) {
    ParsedCli parsed = parse_cli_full(argc, argv);
    if (!parsed.cli.log_dir.empty()) {
        std::filesystem::create_directories(parsed.cli.log_dir);
        gpes::log::Log::getInstance().setLogDirectory(parsed.cli.log_dir);
    }
    return parsed;
}

void print_parsed_params(const ParsedCli& parsed) {
    std::cout << "Parsed CLI parameters:\n";
    print_row(std::cout, "name", "value", "description");
    print_row(std::cout, "cli_only", parsed.cli_only ? "true" : "false",
        "Require all parameters from CLI");
    const auto &cfg = parsed.cli;
    print_row(std::cout, "a_s", std::to_string(cfg.a_s), "Scattering length");
    print_row(std::cout, "a_dd", std::to_string(cfg.a_dd), "Dipole-dipole interaction strength");
    print_row(std::cout, "num_particles", std::to_string(cfg.num_particles), "Number of particles");
    print_row(std::cout, "dt", std::to_string(cfg.dt), "Time step");
    print_row(std::cout, "duration", std::to_string(cfg.duration), "Simulation duration");
    print_row(std::cout, "w_x", std::to_string(cfg.w_x), "Trap frequency along x");
    print_row(std::cout, "w_y", std::to_string(cfg.w_y), "Trap frequency along y");
    print_row(std::cout, "w_z", std::to_string(cfg.w_z), "Trap frequency along z");
    print_row(std::cout, "grid", std::to_string(cfg.nx) + "x" + std::to_string(cfg.ny), "Fourier grid");
    print_row(std::cout, "domain", std::to_string(cfg.x_half_length) + " x " + std::to_string(cfg.y_half_length), "Domain half lengths");
    print_row(std::cout, "max_gs_steps", std::to_string(cfg.max_ground_steps), "Maximum imaginary-time steps");
    print_row(std::cout, "bdg_modes", std::to_string(cfg.bdg_modes), "Requested BdG eigenpairs");
    print_row(std::cout, "log_dir", cfg.log_dir, "Log directory path");
    print_row(std::cout, "wavefunction_out", cfg.wavefunction_out, "Run output directory");
}
