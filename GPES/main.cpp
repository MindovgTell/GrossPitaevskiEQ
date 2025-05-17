#include <iostream>

#include "grid.hpp"


int main(){
    
    GPES::Grid<Dimension::One> grid(100, 0.005);
    grid.set_harmonic_potential(1);
    std::cout << grid.get_size_of_grid() << std::endl;
    return 0;
}