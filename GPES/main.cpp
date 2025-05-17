#include <iostream>

#include "grid.hpp"
#include "visualitsation.hpp"

int main(){
    
    GPES::Grid<Dimension::One> grid(100, 0.005);
    grid.set_harmonic_potential(1);
    GPES::draw(grid);
    return 0;
}