#pragma once

#include <string>
#include <utility>
#include "Core/utility.hpp"

struct CliConfig {
    double a_s = 0.008;
    double a_dd = 0.008843;
    double num_particles = 80000.0;
    double dt = 0.001;
    double duration = 0.1;
    double w_x = 1.0;
    double w_y = 3.0;
    double w_z = 5.0;
    std::string log_dir;
    std::string wavefunction_out = "wavefunction.csv";
};

struct ParsedCli {
    bool cli_only = false;
    CliConfig cli;
    gpes::PhysConfig phys{};
    gpes::SimConfig sim{};
};

ParsedCli parse_cli(int argc, char **argv);
void print_parsed_params(int argc, char **argv);
