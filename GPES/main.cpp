#include <chrono>
#include <regex>

#include "src/grid.hpp"
#include "src/utility.hpp"
#include "src/CrankNicolson.hpp"


void simulation2D(int argc, char const *argv[]);
void simulation1D(int argc, char const *argv[]);


int main(int argc, char const *argv[]){
    auto t0 = std::chrono::steady_clock::now();
    if(argc == 14){
        simulation2D(argc, argv);
    }
    else if(argc == 10){
        simulation1D(argc, argv);
    }
    else {
        std::string executable = argv[0];
        std::cerr << "Error: Wrong number of input parameters" << '\n';
        std::cerr << "Usage:" << executable << " input_file.txt" << '\n';
        return 1;
    }
    auto t1 = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "Simulation have been run for " << elapsed << " s\n";
    return 0;
}


void simulation2D(int argc, char const *argv[]) {
    try {
        GPES::Grid<Dimension::Two> grid(300,300, -10, -10);
        grid.set_harmonic_potential(1,1);
        grid.set_z_freq(5);

        GPES::WaveFunction<Dimension::Two> Phi(grid,0.00607, 0.00882, 10000);
        GPES::CrankNicolson<Dimension::Two> solver(Phi,grid,0.001,0.1, 0.614);
        Phi.print_params();
        solver.print_param_of_eq();
        solver.simulation("directory_for_save_results");
    }
    catch (const std::exception& ex) {           // catches runtime_error and any std::exception
        std::cout << "Error while saving wave-function: " << ex.what() << '\n';
        return;                    // tell the OS something went wrong
    }
    catch (...)                                 // catches non-standard exceptions, just in case
    {
        std::cout << "Unknown error while saving wave-function.\n";
        return;
    }
}

void simulation1D(int argc, char const *argv[]){
    try {
        GPES::Grid<Dimension::One> grid(600, -10);
        grid.set_harmonic_potential(1);
        grid.set_transverse(5);

        GPES::WaveFunction<Dimension::One> Phi(grid,0.00607, 0.00882, 10000);
        GPES::CrankNicolson<Dimension::One> solver(grid, Phi,0.001,0.1);
        Phi.print_params();
        solver.print_param_of_eq();
        solver.simulation("directory_for_save_results");
    }
    catch (const std::exception& ex) {           // catches runtime_error and any std::exception
        std::cout << "Error while saving wave-function: " << ex.what() << '\n';
        return;                    // tell the OS something went wrong
    }
    catch (...)                                 // catches non-standard exceptions, just in case
    {
        std::cout << "Unknown error while saving wave-function.\n";
        return;
    }
}