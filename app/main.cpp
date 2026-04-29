#include "cli/cli.hpp"
#include "gpes.hpp"
#include "Core/plot.hpp"
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>

int main(int argc, char **argv) {
    try {
        auto parsed = parse_cli(argc, argv);
        auto phys = parsed.phys;
        auto sim = parsed.sim;
        const auto &cfg = parsed.cli;
        print_parsed_params(argc, argv);

        const gpes::Dimension Dim = gpes::Dimension::Two;
        auto grid_ptr = std::make_shared<gpes::Grid<Dim>>(
            400, 400, -20, -20, cfg.w_x, cfg.w_y, cfg.w_z
        );

        gpes::WaveFunction<Dim> Psi(grid_ptr);
        Psi.set_state_Gauss(0.0, 0.0, 1.5, 1.5, phys.num_of_prt);
        // Psi.set_state_TF(0.0, 0.0,phys.num_of_prt);
        seed_supersolid_state(
            Psi,
            gpes::SeedType::RandomNoise,
            0.05,
            0.0,
            0.0,
            gpes::ModulationAxis::X,
            12345
        );
        {
            const double l_z = 1.0 / std::sqrt(grid_ptr->omega_z());
            gpes::calc_inter_consts<Dim>(phys, l_z);

            std::cout << "Effective simulation parameters:\n";
            std::cout << std::left << std::setw(16) << "name"
                      << std::setw(18) << "value"
                      << "description\n";
            std::cout << std::setw(16) << "l_z"
                      << std::setw(18) << l_z
                      << "harmonic confinement length\n";
            std::cout << std::setw(16) << "eps_dd"
                      << std::setw(18) << (phys.a_dd / phys.a_s)
                      << "dipolar ratio a_dd / a_s\n";
            std::cout << std::setw(16) << "g_scat"
                      << std::setw(18) << phys.g_scat
                      << "contact coupling used in solver\n";
            std::cout << std::setw(16) << "V_dd"
                      << std::setw(18) << phys.V_dd
                      << "dipolar coupling used in solver\n";
            std::cout << std::setw(16) << "g_lhy"
                      << std::setw(18) << phys.g_lhy
                      << "LHY coupling used in solver\n";
        }
        gpes::plot::heatmap(Psi);
        gpes::solvers::SplitStep<Dim> stepper1(grid_ptr, Psi, phys, sim);
        // gpes::solvers::CrankNicolson<Dim> stepper1(grid_ptr, Psi, phys, sim);
        gpes::TDSESolver<Dim, decltype(stepper1)> solver1(grid_ptr, stepper1, Psi, sim);
        const std::filesystem::path wavefunction_out_path(parsed.cli.wavefunction_out);
        if (wavefunction_out_path.has_parent_path()) {
            std::filesystem::create_directories(wavefunction_out_path.parent_path());
        }
        
        solver1.run(wavefunction_out_path, 4000u);
        const auto& energy1 = stepper1.energies();
        gpes::io::save_state_energy_history_csv(energy1, wavefunction_out_path);
        gpes::plot::heatmap(Psi);
        gpes::plot::draw_slices(Psi);

        gpes::plot::draw_energy(energy1);

        return 0;
    } catch (const std::exception& e) {
        std::cout << "Fatal error: " << e.what() << "\n";
        return 1;
    }
}
