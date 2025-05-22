#include <iostream>
#include <string>

#include "src/grid.hpp"
#include "src/visualitsation.hpp"
#include "src/CrankNicolson.hpp"

void simulation1D(std::string inputfile);
void simulation2D(std::string inputfile);

void test();

int main(int argc, char const *argv[]){

    if(argc != 2){
        std::string executable = argv[0];

        std::cerr << "Error: Wrong number of input parameters" << '\n';
        std::cerr << "Usage:" << executable << " input_file.txt" << '\n';
        return 1;
    }

    // simulation1D(argv[1]);

    simulation2D(argv[1]);

    // test();

    return 0;
}

void simulation1D(std::string inputfile){

    std::ifstream input_data(inputfile);

    std::string line;
    double Problem,grid_size, deltat, T, start, number_of_mol, a_s, a_dd, x_c, sigma_x;

    std::getline(input_data,line);
    while(std::getline(input_data,line)) {//Skip first line in file

        std::stringstream str_stream(line);
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
        Psi.set_state_Gauss(x_c, sigma_x);
        Psi2.set_state_TF(x_c);

        GPES::CrankNicolson<Dimension::One> solver(grid,Psi, deltat, T);
        solver.simulation();

        GPES::WaveFunction<Dimension::One> Fin;
        solver.get_final_state(Fin);

        GPES::draw(Psi,Fin);
        GPES::draw(Psi2,Fin);

        //Calculating energy the energy during simulation

        double TM_energy = solver.calc_state_energy(Psi2);

        std::vector<double> E = solver.get_vec_Energy();
        std::vector<double> TM_en(E.size(), TM_energy);

        GPES::draw_energy(E, TM_en);
    }
}


void simulation2D(std::string inputfile){

    std::ifstream input_data(inputfile);

    std::string line;
    double Problem, grid_size_x, grid_size_y, deltat, T, start_x, start_y, number_of_mol, a_s, a_dd, x_c, y_c, sigma_x, sigma_y, omega_x, omega_y;

    std::vector<double> energies_for_diff_grid;
    energies_for_diff_grid.reserve(8);

    std::getline(input_data,line);
    while(std::getline(input_data,line)) {//Skip first line in file

        std::stringstream str_stream(line);

        if (str_stream >> Problem >> grid_size_x >> grid_size_y >> deltat >> T >> start_x >> start_y >> number_of_mol >> a_s >> a_dd >> x_c >> y_c >> sigma_x >> sigma_y >> omega_x >> omega_y) {
            int width = 10;
            std::cout << std::setw(width) << "Problem" << std::setw(width) << "grid_size_x" << std::setw(width) << "grid_size_y" << std::setw(width) << "deltat" << std::setw(width) << "T" << std::setw(width)
            << "start_x" << std::setw(width) << "start_y" << std::setw(width) << "N_of_mol"  << std::setw(width) << "a_s" << std::setw(width) << "a_dd" << std::setw(width)
            << "x_c"  << std::setw(width) << "y_c" << std::setw(width) << "sigma_x" << std::setw(width) << "sigma_y" << std::setw(width) << "omega_x"<< std::setw(width) << "omega_y" << std::endl;

            std::cout << std::setw(width) << Problem << std::setw(width) << grid_size_x << std::setw(width) << grid_size_y << std::setw(width) << deltat << std::setw(width) << T << std::setw(width)
            << start_x << std::setw(width) << start_y << std::setw(width) << number_of_mol << std::setw(width) << a_s << std::setw(width) << a_dd << std::setw(width) << x_c << std::setw(width)
            << y_c << std::setw(width) << sigma_x  << std::setw(width) << sigma_y << std::setw(width) << omega_x << std::setw(width) << omega_y << std::endl;
        }
        else
            std::cout << "Error reading values from string!" << std::endl;
        

        GPES::Grid<Dimension::Two> grid(grid_size_x,grid_size_y, start_x, start_y);
        grid.set_harmonic_potential(omega_x, omega_y);

        GPES::WaveFunction<Dimension::Two> Psi(grid,a_s, a_dd,number_of_mol);
        GPES::WaveFunction<Dimension::Two> Psi2(grid, a_s, a_dd, number_of_mol); //(Grid, )

        Psi.set_state_Gauss(x_c,y_c,sigma_x,sigma_y);
        Psi2.set_state_TF(x_c, y_c);


        GPES::CrankNicolson<Dimension::Two> solver(Psi,grid,deltat,T);

        solver.simulation();

        GPES::WaveFunction<Dimension::Two> Fin;
        solver.get_final_state(Fin);


        GPES::heatmap(Psi);
        GPES::heatmap(Fin);
        GPES::heatmap(Psi2);



        //Calculating energy the energy during simulation

        double TM_energy = solver.calc_state_energy(Psi2);

        std::vector<double> E = solver.get_vec_Energy();
        std::vector<double> TM_en(E.size(), TM_energy);

        energies_for_diff_grid.push_back(E.back());
        GPES::draw_energy(E, TM_en);

    }

    

}


void test(){
    GPES::Grid<Dimension::Two> grid(5,5,-10,-10);
    grid.set_harmonic_potential(100,100);

    GPES::WaveFunction<Dimension::Two> Psi(grid, 1,1, 1000);
    Psi.set_state_Gauss(0,0,10,10);

    GPES::CrankNicolson<Dimension::Two> solver(Psi, grid, 0.5, 1);
    solver.simulation();
    solver.print_Mat_A();
    // std::cout << solver.print_Mat_A() << std::endl;
}