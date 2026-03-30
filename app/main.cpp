#include "cli/cli.hpp"
#include "gpes.hpp"
#include "Core/plot.hpp"
#include <filesystem>

int main(int argc, char **argv) {
    try {
        auto parsed = parse_cli(argc, argv);
        auto phys = parsed.phys;
        auto sim = parsed.sim;
        const auto &cfg = parsed.cli;
        print_parsed_params(argc, argv);

        const gpes::Dimension Dim = gpes::Dimension::Two;
        auto grid_ptr = std::make_shared<gpes::Grid<Dim>>(
            400, 400, -30, -30, cfg.w_x, cfg.w_y, cfg.w_z
        );

        gpes::WaveFunction<Dim> Psi(grid_ptr);
        Psi.set_state_GaussSum({
            {-6.0, 0.0, 2.8, 2.8, 1.0},
            { 0.0, 0.0, 2.8, 2.8, 1.0},
            { 6.0, 0.0, 2.8, 2.8, 1.0}
        }, phys.num_of_prt);
        seed_supersolid_state(
            Psi,
            gpes::SeedType::Sinusoidal,
            0.08,
            6.0,
            0.0,
            gpes::ModulationAxis::X,
            12345
        );
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

        return 0;
    } catch (const std::exception& e) {
        std::cout << "Fatal error: " << e.what() << "\n";
        return 1;
    }
}
