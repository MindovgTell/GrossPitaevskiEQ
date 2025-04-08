#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <complex>
#include "./CrankNicolson/CrankNicolson.hpp"
#include "./matplotlibcpp.h"

using namespace std::complex_literals;

namespace plt = matplotlibcpp;


//--------//-----------------------------------//--------//
//--------//             Simulators            //--------//
//--------//-----------------------------------//--------//
void simulation1D(std::string inputfile);
void simulation2D(std::string inputfile);

//--------//-----------------------------------//--------//
//--------//Additional drawing functions for 1D//--------//
//--------//-----------------------------------//--------//
void draw1(Eigen::VectorXd& x, Eigen::VectorXd& Psi);

void draw2(Eigen::VectorXd& x, Eigen::VectorXd& Psi, Eigen::VectorXd& Fin);
void draw3(Eigen::VectorXd& x, Eigen::VectorXd& Psi, Eigen::VectorXd& Fin, Eigen::VectorXd& V);

void draw_energy(std::vector<double>& x, std::vector<double>& Psi,  std::vector<double>& TM_en);


//--------//-----------------------------------//--------//
//--------//Additional drawing functions for 2D//--------//
//--------//-----------------------------------//--------//




int main(int argc, char const *argv[]){

    if(argc != 2){
        std::string executable = argv[0];

        std::cerr << "Error: Wrong number of input parameters" << '\n';
        std::cerr << "Usage:" << executable << " input_file.txt" << '\n';
        return 1;
    }

    simulation1D(argv[1]);

    // simulation2D(argv[1]);

    return 0;
}


void simulation1D(std::string inputfile){

    std::ifstream input_data(inputfile);

    std::string line;
    std::getline(input_data,line);//Skip first line in file

    double Problem,h,deltat,T,x_c,sigma_x,p_x,omega,N,a_s, start;

    std::getline(input_data,line);
    std::stringstream str_stream(line);

    //str_stream.imbue(std::locale("C"));

    if (str_stream >> Problem >> h >> deltat >> T >> x_c >> sigma_x >> p_x >> omega >> N >> a_s >> start)
        std::cout << "Values read successfully!" << std::endl;
    else
        std::cout << "Error reading values from string!" << std::endl;
    
    int width = 10;
    std::cout << std::setw(width) << "Problem" << std::setw(width) << "h" << std::setw(width) << "deltat" << std::setw(width) << "T" << std::setw(width) << "x_c" << std::setw(width)
    << "sigma_x" << std::setw(width) << "p_x" << std::setw(width) << "omega" << std::setw(width) << "N" << std::setw(width) << "a_s" << std::setw(width) << "start" << std::endl;

    std::cout << std::setw(width) << Problem << std::setw(width) << h << std::setw(width) << deltat << std::setw(width) << T << std::setw(width) << x_c << std::setw(width)
    << sigma_x << std::setw(width) << p_x << std::setw(width) << omega << std::setw(width) << N << std::setw(width) << a_s << std::setw(width) << start << std::endl;

    CrankNicolson Crank(h, deltat, T, x_c, sigma_x, p_x, omega, N, a_s, start);

    Eigen::VectorXd Psi = Crank.get_m_Psi_prob(); 
    Eigen::VectorXcd Ps = Crank.get_m_Psi();  

    Crank.simulation_1D();

    Eigen::VectorXd Fin = Crank.get_m_Fin_prob(); 
    Eigen::VectorXcd Fn = Crank.get_m_Fin(); 

    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(500, start, -1 * start);

    Eigen::VectorXd V = Crank.get_m_V();
    Eigen::VectorXcd TM = Crank.TM_state();
    Eigen::VectorXd TM_pr = Crank.prob_1D(TM);

    // std::cout << '\n' << '\n' << std::endl;

    // Crank.print_Mat_A();

    std::cout << '\n' << '\n' << std::endl;

    // Crank.print_Mat_B();
    
    std::cout << "Initial state norm: " << Crank.vec_norm(Ps) << std::endl;
    std::cout << "Fin state norm: " << Crank.vec_norm(Fn) << std::endl;
    std::cout << "TM state norm: " << Crank.vec_norm(TM) << std::endl;


    std::cout << '\n' << '\n' << std::endl;

    std::cout << "The energy of the current state: " << Crank.calc_state_energy() << std::endl;
    std::cout << "The energy of the Thomas-Fermi state: " << Crank.calc_state_energy(TM) << std::endl;



    std::cout << '\n' << '\n' << std::endl;

    double TM_energy = Crank.calc_state_energy(TM);

    std::vector<double> E = Crank.get_vec_Energy();
    std::vector<double> en_len(E.size());
    std::vector<double> TM_en(E.size(), TM_energy);
    std::iota(en_len.begin(), en_len.end(), 1); // fills with 1, 2, 3, ..., y.size()


    draw3(x, Psi, Fin, V);
    draw3(x, Psi, Fin, TM_pr);

    draw_energy(en_len, E, TM_en);
}



void simulation2D(std::string inputfile){

    // Functionality for reading all the neccesary for simulation values from source file
    std::ifstream input_data(inputfile);

    std::string line;
    std::getline(input_data,line);//Skip first line in file

    double Problem, h, deltat, T, x_c, sigma_x, y_c, sigma_y, omega_x, omega_y, N, a_s, start;

    std::getline(input_data,line);
    std::stringstream str_stream(line);

    if (str_stream >> Problem >> h >> deltat >> T >> x_c >> sigma_x >> y_c >> sigma_y >> omega_x >> omega_y >> N >> a_s >> start)
        std::cout << "Values read successfully!" << std::endl;
    else
        std::cout << "Error reading values from string!" << std::endl;

    int width = 10;

    //Writing into console all parameters
    std::cout << std::setw(width) << "Problem" << std::setw(width) << "h" << std::setw(width) << "deltat" << std::setw(width) << "T" << std::setw(width) << "x_c" << std::setw(width)
    << "sigma_x" << std::setw(width) << "y_c"<< std::setw(width) << "sigma_y" << std::setw(width) << "omega_x" << std::setw(width) << "omega_y" << std::setw(width) << "N" << std::setw(width) << "a_s" << std::setw(width) << "start" << std::endl;

    std::cout << std::setw(width) << Problem << std::setw(width) << h << std::setw(width) << deltat << std::setw(width) << T << std::setw(width) << x_c << std::setw(width)
    << sigma_x << std::setw(width) << y_c << std::setw(width) << sigma_y << std::setw(width) << omega_x << std::setw(width) << omega_y << std::setw(width) << N << std::setw(width) << a_s << std::setw(width) << start << std::endl;

    //Creating the simulator
    CrankNicolson Crank(h, deltat, T, x_c, sigma_x, y_c, sigma_y, omega_x, omega_y, N, a_s, start);

    

}

//--------//-----------------------------------//--------//
//--------//Additional drawing functions for 1D//--------//
//--------//-----------------------------------//--------//

void draw1(Eigen::VectorXd& x, Eigen::VectorXd& Psi){
    std::vector<double> x_vec(x.data(), x.data() + x.size());
    std::vector<double> y_vec(Psi.data(), Psi.data() + Psi.size());

    plt::figure();
    plt::plot(x_vec,y_vec, std::string("b-"),{{"label", "data trend"}});
    plt::xlabel("time [s]");
    plt::ylabel("observation [m]");
    plt::legend();
    plt::grid();
    plt::show(); 
}

void draw_energy(std::vector<double>& x, std::vector<double>& vec, std::vector<double>& TM_en){
    std::vector<double> x_vec(x.data(), x.data() + x.size());

    plt::figure();
    plt::plot(x_vec, vec, std::string("b-"),{{"label", "data trend"}});
    plt::plot(x_vec, TM_en, std::string("g--"),{{"label", "TM state energy"}});
    plt::xlabel("time [s]");
    plt::ylabel("observation [m]");
    plt::legend();
    plt::grid();
    plt::show(); 
}

void draw2(Eigen::VectorXd& x, Eigen::VectorXd& Psi,  Eigen::VectorXd& Fin){
    std::vector<double> x_vec(x.data(), x.data() + x.size());
    std::vector<double> y_vec(Psi.data(), Psi.data() + Psi.size());
    std::vector<double> z_vec(Fin.data(), Fin.data() + Fin.size());

    plt::figure();
    plt::plot(x_vec,y_vec, std::string("b-"),{{"label", "Initial state"}});
    plt::plot(x_vec,z_vec, std::string("r-"),{{"label", "Final state"}});
    plt::xlabel("");
    plt::ylabel("");
    plt::legend();
    plt::grid();
    plt::show(); 
}

void draw3(Eigen::VectorXd& x, Eigen::VectorXd& Psi,  Eigen::VectorXd& Fin, Eigen::VectorXd& V){
    std::vector<double> x_vec(x.data(), x.data() + x.size());
    std::vector<double> y_vec(Psi.data(), Psi.data() + Psi.size());
    std::vector<double> z_vec(Fin.data(), Fin.data() + Fin.size());
    std::vector<double> t_vec(V.data(), V.data() + V.size());

    plt::figure();
    plt::plot(x_vec,y_vec, std::string("b-"),{{"label", "Initial state"}});
    plt::plot(x_vec,z_vec, std::string("r-"),{{"label", "Final state"}});
    plt::plot(x_vec,t_vec, std::string("g-"),{{"label", "Potential well"}});
    plt::xlabel("");
    plt::ylabel("");
    plt::legend();
    plt::grid();
    plt::show(); 
}



//--------//-----------------------------------//--------//
//--------//Additional drawing functions for 2D//--------//
//--------//-----------------------------------//--------//