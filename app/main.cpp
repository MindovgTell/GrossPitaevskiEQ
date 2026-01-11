#include "cli/cli.hpp"
#include "gpes.hpp"
#include "Core/plot.hpp"

int main(int argc, char **argv) {
    try {
        auto [phys, sim] = parse_cli(argc, argv);
        print_parsed_params(argc, argv);

        const gpes::Dimension Dim = gpes::Dimension::One;

        auto grid_ptr = std::make_shared<gpes::Grid<Dim>>(128, -10, 1, 5); // Grid ctor parameters 100 -- grid size; -10 -- start

        gpes::WaveFunction<Dim> Psi(grid_ptr);
        Psi.set_state_Gauss(0, 1, phys.num_of_prt);
        auto Psi2 = Psi;

        // gpes::solvers::CrankNicolson<Dim> stepper(grid_ptr, Psi, phys, sim);

        gpes::solvers::SplitStep<Dim> stepper(grid_ptr, Psi, phys, sim);

        gpes::TDSESolver<Dim, decltype(stepper)> solver(grid_ptr, stepper, Psi, sim);

        solver.run();

        gpes::plot::draw(Psi, Psi2);
        gpes::plot::draw_energy(stepper.energies());

        return 0;
    } catch (const std::exception& e) {
        std::cout << "Fatal error: " << e.what() << "\n";
        return 1;
    }
}
