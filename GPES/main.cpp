#include <iostream>
#include <string>

#include "grid.hpp"
#include "visualitsation.hpp"
#include "CrankNicolson.hpp"

void simulation1D(std::string inputfile);
void simulation2D(std::string inputfile);

int main(int argc, char const *argv[]){

    if(argc != 2){
        std::string executable = argv[0];

        std::cerr << "Error: Wrong number of input parameters" << '\n';
        std::cerr << "Usage:" << executable << " input_file.txt" << '\n';
        return 1;
    }

    // simulation1D(argv[1]);

    simulation2D(argv[1]);

    return 0;
}

void simulation1D(std::string inputfile){

    std::ifstream input_data(inputfile);

    std::string line;
    double Problem,grid_size, deltat, T, start, number_of_mol, a_s, a_dd, x_c, sigma_x;

    std::getline(input_data,line);
    while(std::getline(input_data,line)) {//Skip first line in file

        std::stringstream str_stream(line);

        //str_stream.imbue(std::locale("C"));

        if (str_stream >> Problem >> grid_size >> deltat >> T >> start >> number_of_mol >> a_s >> a_dd >> x_c >> sigma_x) {
                int width = 10;
            std::cout << std::setw(width) << "Problem" << std::setw(width) << "grid_size" << std::setw(width) << "deltat" << std::setw(width) << "T" << std::setw(width) << "start" << std::setw(width)
            << "N_of_mol"  << std::setw(width) << "a_s" << std::setw(width) << "a_dd" << std::setw(width) << "x_c" << std::setw(width) << "sigma_x" << std::endl;

            std::cout << std::setw(width) << Problem << std::setw(width) << grid_size << std::setw(width) << deltat << std::setw(width) << T << std::setw(width) << start << std::setw(width)
            << number_of_mol << std::setw(width) << a_s << std::setw(width) << a_dd << std::setw(width) << x_c << std::setw(width) << sigma_x << std::endl;
        }
        else
            std::cout << "Error reading values from string!" << std::endl;
        
        //initialize grid
        GPES::Grid<Dimension::One> grid(grid_size, start); // grid(Number_of_grid_nodes, starting_position)
        grid.set_harmonic_potential(1);

        //initialize wavefunction
        GPES::WaveFunction<Dimension::One> Psi(grid, a_s, a_dd, number_of_mol);
        GPES::WaveFunction<Dimension::One> Psi2(grid, a_s, a_dd, number_of_mol); //(Grid, )
        Psi.set_state_Gauss(x_c,sigma_x);
        Psi2.set_state_TM(x_c);

        GPES::CrankNicolson<Dimension::One> solver(grid,Psi, deltat, T);
        solver.simulation();

        GPES::WaveFunction<Dimension::One> Fin;
        solver.get_final_state(Fin);

        GPES::draw(Psi,Fin);
        GPES::draw(Psi2,Fin);
    }
}


void simulation2D(std::string inputfile){

    std::ifstream input_data(inputfile);

    std::string line;
    double Problem, grid_size_x, grid_size_y, deltat, T, start_x, start_y, number_of_mol, a_s, a_dd, x_c, y_c, sigma_x, sigma_y;

    std::getline(input_data,line);
    while(std::getline(input_data,line)) {//Skip first line in file

        std::stringstream str_stream(line);

        //str_stream.imbue(std::locale("C"));

        if (str_stream >> Problem >> grid_size_x >> grid_size_y >> deltat >> T >> start_x >> start_y >> number_of_mol >> a_s >> a_dd >> x_c >> y_c >> sigma_x >> sigma_y) {
            int width = 10;
            std::cout << std::setw(width) << "Problem" << std::setw(width) << "grid_size_x" << std::setw(width) << "grid_size_y" << std::setw(width) << "deltat" << std::setw(width) << "T" << std::setw(width)
            << "start_x" << std::setw(width) << "start_y" << std::setw(width) << "N_of_mol"  << std::setw(width) << "a_s" << std::setw(width) << "a_dd" << std::setw(width)
            << "x_c"  << std::setw(width) << "y_c" << std::setw(width) << "sigma_x" << std::setw(width) << "sigma_y" << std::endl;

            std::cout << std::setw(width) << Problem << std::setw(width) << grid_size_x << std::setw(width) << grid_size_y << std::setw(width) << deltat << std::setw(width) << T << std::setw(width)
            << start_x << std::setw(width) << start_y << std::setw(width) << number_of_mol << std::setw(width) << a_s << std::setw(width) << a_dd << std::setw(width) << x_c << std::setw(width)
            << y_c << std::setw(width) << sigma_x  << std::setw(width) << sigma_y << std::endl;
        }
        else
            std::cout << "Error reading values from string!" << std::endl;
        

        GPES::Grid<Dimension::Two> grid(grid_size_x,grid_size_y, start_x, start_y);
        grid.set_harmonic_potential(1,1);

        GPES::WaveFunction<Dimension::Two> Psi(grid,a_s, a_dd,number_of_mol);

        Psi.set_state_Gauss(x_c,y_c,sigma_x,sigma_y);

        GPES::CrankNicolson<Dimension::Two> solver(Psi,grid,deltat,T);

        solver.simulation();

        GPES::WaveFunction<Dimension::Two> Fin;
        solver.get_final_state(Fin);


        GPES::heatmap(Psi);
        GPES::heatmap(Fin);
        GPES::heatmap(grid);

    }
}