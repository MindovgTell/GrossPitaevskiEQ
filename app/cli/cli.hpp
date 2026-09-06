#pragma once

#include <string>
#include <cstddef>
#include <utility>
#include "Core/utility.hpp"

struct CliConfig {
    double a_s = 0.00265;
    double a_dd = 0.008843;
    double num_particles = 80000.0;
    double dt = 0.001;
    double duration = 0.1;
    double w_x = 1.0;
    double w_y = 5.0;
    double w_z = 8.0;
    int nx = 1024;
    int ny = 1024;
    double x_half_length = 25.0;
    double y_half_length = 25.0;
    std::size_t max_ground_steps = 50000;
    std::size_t convergence_check_interval = 100;
    std::size_t convergence_checks_required = 3;
    double relative_energy_tolerance = 1e-10;
    double stationary_residual_tolerance = 1e-4;
    double minimum_ground_dt = 1e-4;
    std::size_t snapshot_interval = 5000;
    double seed_noise = 1e-3;
    unsigned int random_seed = 12345;
    bool skip_bdg = false;
    bool allow_unconverged = false;
    bool save_mode_fields = true;
    bool fresh_start = false;
    int bdg_modes = 24;
    int bdg_ncv = 64;
    int bdg_max_iter = 1500;
    double bdg_tolerance = 1e-7;
    double bdg_residual_tolerance = 1e-4;
    double bdg_shift = 0.05;
    int bdg_linear_max_iter = 1000;
    int bdg_linear_restart = 48;
    double bdg_linear_tolerance = 1e-8;
    int bdg_dense_max_dimension = 0;
    std::string log_dir;
    std::string wavefunction_out;
};

struct ParsedCli {
    bool cli_only = false;
    CliConfig cli;
    gpes::PhysConfig phys{};
    gpes::SimConfig sim{};
};

ParsedCli parse_cli(int argc, char **argv);
void print_parsed_params(const ParsedCli& parsed);
