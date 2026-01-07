#pragma once

#include <string>
#include <utility>

#include "Core/definitions.hpp"

struct CliConfig {
    double a_s = 0.00607;
    double a_dd = 0.00882;
    double num_particles = 10000.0;
    double dt = 0.001;
    double duration = 0.1;
    std::string log_file;
};

std::pair<gpes::PhysConfig, gpes::SimConfig> parse_cli(int argc, char **argv);
void print_parsed_params(int argc, char **argv);
