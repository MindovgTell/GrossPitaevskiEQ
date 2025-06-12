#include <chrono>

#include "src/grid.hpp"
#include "src/visualitsation.hpp"
#include "src/CrankNicolson.hpp"

void simulation1D(std::string inputfile);
void simulation2D(std::string inputfile);
void simulation2Dt(std::string inputfile);

void test();

int main(int argc, char const *argv[]){

    if(argc != 2){
        std::string executable = argv[0];

        std::cerr << "Error: Wrong number of input parameters" << '\n';
        std::cerr << "Usage:" << executable << " input_file.txt" << '\n';
        return 1;
    }

    simulation1D(argv[1]);

    // simulation2Dt(argv[1]);

    // test();

    return 0;
}

void simulation1D(std::string inputfile){

    std::ifstream input_data(inputfile);

    std::string line;
    double Problem,grid_size, deltat, T, start, number_of_mol, a_s, e_dd, x_c, sigma_x, omega_t;
    std::vector<double> energies_for_diff_grid;
    energies_for_diff_grid.reserve(8);

    std::getline(input_data,line);

    // Loop for running simulation and getting all the fin states
    while(std::getline(input_data,line)) {//Skip first line in file
        try { 
            std::stringstream str_stream(line);
            if (str_stream >> Problem >> grid_size >> deltat >> T >> start >> number_of_mol >> e_dd >> a_s >> x_c >> sigma_x >> omega_t) {
                    int width = 10;
                std::cout << std::setw(width) << "Problem" << std::setw(width) << "grid_size" << std::setw(width) << "deltat" << std::setw(width) << "T" << std::setw(width) << "start" << std::setw(width)
                << "N_of_mol"  << std::setw(width) << "e_dd" << std::setw(width) << "a_s" << std::setw(width) << "x_c" << std::setw(width) << "sigma_x" << std::setw(width) << "omega_t" << std::endl;

                std::cout << std::setw(width) << Problem << std::setw(width) << grid_size << std::setw(width) << deltat << std::setw(width) << T << std::setw(width) << start << std::setw(width)
                << number_of_mol << std::setw(width) << e_dd << std::setw(width) << a_s << std::setw(width) << x_c << std::setw(width) << sigma_x << std::setw(width) << omega_t << std::endl;
            }
            else {
                std::cout << "Error reading values from string!" << std::endl;
                return;
            }
            double a_dd = e_dd * a_s; 
            //initialize grid
            GPES::Grid<Dimension::One> grid(grid_size, start); // grid(Number_of_grid_nodes, starting_position)
            grid.set_harmonic_potential(1);
            grid.set_transverse(omega_t);
            //initialize wavefunction
            GPES::WaveFunction<Dimension::One> Psi(grid, a_s, a_dd, number_of_mol);
            GPES::WaveFunction<Dimension::One> Psi2(grid, a_s, a_dd, number_of_mol); //(Grid, )
            Psi.set_state_Gauss(x_c, sigma_x);
            Psi2.set_state_TF(x_c);
            //initialize simulator 
            GPES::CrankNicolson<Dimension::One> solver(grid,Psi2, deltat, T);
            //run simulation
            solver.print_param_of_eq();
            solver.simulation();
            //get and save final state to csv
            GPES::WaveFunction<Dimension::One> Fin;
            solver.get_final_state(Fin);
            Fin.savecsv_state("../../res/Fin_states/sim.csv");
            //Calculating energy the energy during simulation
            double TM_energy = solver.calc_state_energy(Psi2);
            std::vector<double> E = solver.get_vec_Energy();
            std::vector<double> TM_en(E.size(), TM_energy);
            energies_for_diff_grid.push_back(E.back());
            std::cout << E.back() << std::endl;

            GPES::draw(Psi, Fin, Psi2);
            GPES::draw_energy(E, TM_en);
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
    // GPES::draw_energy(energies_for_diff_grid);
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


        GPES::CrankNicolson<Dimension::Two> solver(Psi2,grid,deltat,T);

        // std::cout << "g_scattering: " << solver.get_g_s() << std::endl;
        // std::cout << "g_lhy: " << solver.get_g_lhy() << std::endl;
        // std::cout << "C_dd: " << solver.get_C_dd() << std::endl;

        solver.simulation();
        GPES::WaveFunction<Dimension::Two> Fin;
        solver.get_final_state(Fin);
        GPES::heatmap(Psi2, "Initial state");
        GPES::heatmap(Fin, "Final state");
        GPES::heatmap(Psi2, "Thomas-Fermi limit");
        //Calculating energy the energy during simulation
        double TM_energy = solver.calc_state_energy(Psi2);
        std::vector<double> E = solver.get_vec_Energy();
        std::vector<double> TM_en(E.size(), TM_energy);
        energies_for_diff_grid.push_back(E.back());
        GPES::draw_energy(E, TM_en);

    }

    GPES::draw_energy(energies_for_diff_grid);

}
 

void simulation2Dt(std::string inputfile){

    std::ifstream input_data(inputfile);

    std::string line;
    double RunNum, grid_size_x, grid_size_y, deltat, T, start_x, start_y, number_of_mol, e_dd, a_s, x_c, y_c, sigma_x, sigma_y, omega_x, omega_y;

    std::vector<double> energies_for_diff_grid;
    energies_for_diff_grid.reserve(8);

    std::getline(input_data,line);
    while(std::getline(input_data,line)) {//Skip first line in file

        std::stringstream str_stream(line);

        if (str_stream >> RunNum >> grid_size_x >> grid_size_y >> deltat >> T >> start_x >> start_y >> number_of_mol >> e_dd >> a_s >> x_c >> y_c >> sigma_x >> sigma_y >> omega_x >> omega_y) {
            int width = 10;
            std::cout << std::setw(width) << "Problem" << std::setw(width) << "grid_size_x" << std::setw(width) << "grid_size_y" << std::setw(width) << "deltat" << std::setw(width) << "T" << std::setw(width)
            << "start_x" << std::setw(width) << "start_y" << std::setw(width) << "N_of_mol"  << std::setw(width) << "e_dd" << std::setw(width) << "a_s" << std::setw(width)
            << "x_c"  << std::setw(width) << "y_c" << std::setw(width) << "sigma_x" << std::setw(width) << "sigma_y" << std::setw(width) << "omega_x"<< std::setw(width) << "omega_y" << std::endl;

            std::cout << std::setw(width) << RunNum << std::setw(width) << grid_size_x << std::setw(width) << grid_size_y << std::setw(width) << deltat << std::setw(width) << T << std::setw(width)
            << start_x << std::setw(width) << start_y << std::setw(width) << number_of_mol << std::setw(width) << e_dd << std::setw(width) << a_s << std::setw(width) << x_c << std::setw(width)
            << y_c << std::setw(width) << sigma_x  << std::setw(width) << sigma_y << std::setw(width) << omega_x << std::setw(width) << omega_y << std::endl;
        }
        else
            std::cout << "Error reading values from string!" << std::endl;
        
        double a_dd = e_dd * a_s;

        GPES::Grid<Dimension::Two> grid(grid_size_x,grid_size_y, start_x, start_y);
        grid.set_harmonic_potential(omega_x, omega_y);

        GPES::WaveFunction<Dimension::Two> Psi(grid,a_s, a_dd,number_of_mol);
        GPES::WaveFunction<Dimension::Two> Psi2(grid, a_s, a_dd, number_of_mol); //(Grid, )

        Psi.set_state_Gauss(x_c,y_c,sigma_x,sigma_y);
        Psi2.set_state_TF(x_c, y_c);


        GPES::CrankNicolson<Dimension::Two> solver(Psi,grid,deltat,T);

        // // std::cout << "g_scattering: " << solver.get_g_s() << std::endl;
        // // std::cout << "g_lhy: " << solver.get_g_lhy() << std::endl;
        // // std::cout << "C_dd: " << solver.get_C_dd() << std::endl;

        solver.simulation();
        GPES::WaveFunction<Dimension::Two> Fin;
        solver.get_final_state(Fin);
        Eigen::VectorXd x_Fin_slice = Fin.get_x_slice();
        Eigen::VectorXd y_Fin_slice = Fin.get_y_slice();
        GPES::heatmap(Psi, "Initial state");
        GPES::heatmap(Fin, "Final state");
        GPES::draw(x_Fin_slice, "x", "Density");
        GPES::draw(y_Fin_slice, "y", "Density");
        GPES::heatmap(Psi2, "Thomas-Fermi limit");
        //Calculating energy the energy during simulation
        double TM_energy = solver.calc_state_energy(Psi2);
        std::vector<double> E = solver.get_vec_Energy();
        std::vector<double> TM_en(E.size(), TM_energy);
        energies_for_diff_grid.push_back(E.back());
        GPES::draw_energy(E, TM_en);

    }

    GPES::draw_energy(energies_for_diff_grid);

}


void test(){
    try {
        GPES::Grid<Dimension::One> grid(100,-20);
        grid.set_harmonic_potential(1);
        GPES::WaveFunction<Dimension::One> Phi(grid,1,0,1000);
        Phi.set_state_Gauss(0,10);
        
        Phi.savecsv_state("../../res/test/phi.csv");
        Phi.print_params();
        GPES::WaveFunction<Dimension::One> Psi;
        Psi.readcsv("../../res/test/phi.csv");
        Psi.print_params();
        GPES::CrankNicolson<Dimension::One> solver(grid,Psi, 10e-5, 0.01);
        solver.simulation();

        GPES::WaveFunction<Dimension::One> Fin;
        solver.get_final_state(Fin);
        Fin.print_params();
        GPES::draw(Psi, Fin);

        Fin.savecsv_state("../../res/test/psi.csv");
        GPES::WaveFunction<Dimension::One> Psi2;
        Psi2.readcsv("../../res/test/psi.csv");
        Psi2.print_params();
        GPES::draw(Psi2);

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