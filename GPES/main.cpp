#include <iostream>

#include "grid.hpp"
#include "visualitsation.hpp"
#include "CrankNicolson.hpp"

void simulation1D();
void simulation2D();

int main(){
    // simulation1D();
    simulation2D();

    return 0;
}

void simulation1D(){
    //initialize grid
    GPES::Grid<Dimension::One> grid(100, -20); // grid(Number_of_grid_nodes, starting_position)
    grid.set_harmonic_potential(1);

    //initialize wavefunction
    GPES::WaveFunction<Dimension::One> Psi(grid, 1, 1, 1000); //(Grid, )
    Psi.set_state_Gauss(0,5,0);

    GPES::CrankNicolson<Dimension::One> solver(grid,Psi, 0.05, 10);
    solver.simulation();

    GPES::WaveFunction<Dimension::One> Fin;
    solver.get_final_state(Fin);

    GPES::draw(Psi,Fin);
}


void simulation2D(){
    GPES::Grid<Dimension::Two> grid(100,100, -20, -20);
    grid.set_harmonic_potential(1,100);

    GPES::heatmap(grid);
    draw(grid);
}