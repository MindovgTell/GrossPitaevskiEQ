#include "cli/cli.hpp"
#include "gpes.hpp"

int main(int argc, char **argv) {
    auto [phys, sim] = parse_cli(argc, argv);
    print_parsed_params(argc, argv);

    const gpes::Dimension Dim = gpes::Dimension::One;

    auto grid_ptr = std::make_shared<gpes::Grid<Dim>>(100,-10); // Grid ctor parameters 100 -- grid size; -10 -- start 
    grid_ptr->set_harmonic_potential(5);

    gpes::WaveFunction<Dim> Psi(grid_ptr);
    Psi.set_state_Gauss(0, 5, 1000);

    gpes::solvers::CrankNicolson<Dim> stepper(grid_ptr, Psi, phys, sim);

    gpes::TDSESolver<Dim, decltype(stepper)> solver(grid_ptr, stepper, Psi, sim);

    return 0;
}
