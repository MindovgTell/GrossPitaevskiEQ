#include "cli/cli.hpp"
#include "gpes.hpp"
#include "Core/utility.hpp"
#include "Core/plot.hpp"
#include <filesystem>
#include <iostream>


int main(int argc, char **argv) {
    try {
        auto parsed = parse_cli(argc, argv);
        auto phys = parsed.phys;
        auto sim = parsed.sim;
        phys.a_s = 0.0045;
        sim.ground_state = true;
        const auto &cfg = parsed.cli;
        print_parsed_params(argc, argv);

        constexpr gpes::Dimension Dim = gpes::Dimension::Two;
        auto grid_ptr = std::make_shared<gpes::Grid<Dim>>(
            300, 300, -25, -25, cfg.w_x, cfg.w_y, cfg.w_z
        );

        gpes::WaveFunction<Dim> Psi(grid_ptr);
        Psi.set_state_Gauss(0.0, 0.0, 2.5, 2.5, phys.num_of_prt);

        gpes::solvers::SplitStep<Dim> stepper1(grid_ptr, Psi, phys, sim);

        gpes::TDSESolver<Dim, decltype(stepper1)> solver1(grid_ptr, stepper1, Psi, sim);
        const std::filesystem::path wavefunction_out_path(parsed.cli.wavefunction_out);
        if (wavefunction_out_path.has_parent_path()) {
            std::filesystem::create_directories(wavefunction_out_path.parent_path());
        }

        gpes::plot::heatmap(Psi);

        solver1.run(wavefunction_out_path, 4000u);

        gpes::plot::heatmap(Psi);
        gpes::plot::draw_slices(Psi);

        auto temp = Psi.fft_shifted();
        gpes::WaveFunction<Dim> kPsi(grid_ptr, temp);
        gpes::plot::heatmap(kPsi);
        gpes::plot::draw_slices(kPsi);

        // const auto& energy1 = stepper1.energies();

        return 0;
    } catch (const std::exception& e) {
        std::cout << "Fatal error: " << e.what() << "\n";
        return 1;
    }
}
