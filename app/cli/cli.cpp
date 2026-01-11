// #include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "CLI11.hpp"

#include "cli.hpp"
#include "Log/Log.hpp"


namespace {

struct ParsedCli {
    bool cli_only = false;
    CliConfig cfg;
    gpes::PhysConfig phys{};
    gpes::SimConfig sim{};
};

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
    std::string log_dir;

    auto *opt_a_s = app.add_option("--a-s", a_s, "Scattering length")->check(CLI::PositiveNumber);
    auto *opt_a_dd = app.add_option("--a-dd", a_dd, "Dipole-dipole interaction strength")
                         ->check(CLI::PositiveNumber);
    auto *opt_num_particles = app.add_option("--num-particles", num_particles, "Number of particles")
                                  ->check(CLI::PositiveNumber);
    auto *opt_dt = app.add_option("--dt", dt, "Time step")->check(CLI::PositiveNumber);
    auto *opt_duration = app.add_option("--duration", duration, "Simulation duration")
                             ->check(CLI::PositiveNumber);
    auto *opt_log_dir = app.add_option("--log-dir", log_dir, "Log directory path (required)")
                             ->required();

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        std::exit(app.exit(e));
    }

    if (cli_only) {
        if (opt_a_s->count() == 0 || opt_a_dd->count() == 0 || opt_num_particles->count() == 0 ||
            opt_dt->count() == 0 || opt_duration->count() == 0 || opt_log_dir->count() == 0) {
            throw CLI::ValidationError(
                "With --cli-only, all parameters must be provided on the command line.");
        }
    }

    cfg.a_s = a_s;
    cfg.a_dd = a_dd;
    cfg.num_particles = num_particles;
    cfg.dt = dt;
    cfg.duration = duration;
    cfg.log_dir = log_dir;
    gpes::PhysConfig phys{cfg.a_s, cfg.a_dd, cfg.num_particles};
    gpes::SimConfig sim{cfg.duration, cfg.dt,
        static_cast<std::size_t>(std::ceil(cfg.duration / cfg.dt))};

    ParsedCli parsed;
    parsed.cli_only = cli_only;
    parsed.cfg = cfg;
    parsed.phys = phys;
    parsed.sim = sim;
    return parsed;

}

void print_row(std::ostream &os, const std::string &name, const std::string &value,
    const std::string &desc) {
    os << std::left << std::setw(16) << name << " " << std::setw(16) << value << " " << desc << "\n";
}

} // namespace

std::pair<gpes::PhysConfig, gpes::SimConfig> parse_cli(int argc, char **argv) {
    ParsedCli parsed = parse_cli_full(argc, argv);
    if (!parsed.cfg.log_dir.empty()) {
        gpes::log::Log::getInstance().setLogDirectory(parsed.cfg.log_dir);
    }
    return {parsed.phys, parsed.sim};
}

void print_parsed_params(int argc, char **argv) {
    ParsedCli parsed = parse_cli_full(argc, argv);

    std::cout << "Parsed CLI parameters:\n";
    print_row(std::cout, "name", "value", "description");
    print_row(std::cout, "cli_only", parsed.cli_only ? "true" : "false",
        "Require all parameters from CLI");
    const auto &cfg = parsed.cfg;
    print_row(std::cout, "a_s", std::to_string(cfg.a_s), "Scattering length");
    print_row(std::cout, "a_dd", std::to_string(cfg.a_dd), "Dipole-dipole interaction strength");
    print_row(std::cout, "num_particles", std::to_string(cfg.num_particles), "Number of particles");
    print_row(std::cout, "dt", std::to_string(cfg.dt), "Time step");
    print_row(std::cout, "duration", std::to_string(cfg.duration), "Simulation duration");
    print_row(std::cout, "log_dir", cfg.log_dir, "Log directory path");
}
