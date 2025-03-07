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

void simulation(std::string inputfile);
void draw(Eigen::VectorXd& x, Eigen::VectorXd& Psi);

int main(int argc, char const *argv[]){

    if(argc != 2){
        std::string executable = argv[0];

        std::cerr << "Error: Wrong number of input parameters" << '\n';
        std::cerr << "Usage:" << executable << " input_file.txt" << '\n';
        return 1;
    }

    simulation(argv[1]);

    return 0;
}


void simulation(std::string inputfile){


    std::ifstream input_data(inputfile);

    std::string line;
    std::getline(input_data,line);//Skip first line in file

    double Problem,h,deltat,T,x_c,sigma_x,p_x,omega,N,a_s;

    std::getline(input_data,line);
    std::stringstream str_stream(line);

    //str_stream.imbue(std::locale("C"));

    if (str_stream >> Problem >> h >> deltat >> T >> x_c >> sigma_x >> p_x >> omega >> N >> a_s) {
        std::cout << "Values read successfully!" << std::endl;
    } else {
        std::cout << "Error reading values from string!" << std::endl;
    }



    // str_stream >> Problem >> h >>deltat >>T >> x_c >> sigma_x >> p_x >> y_c >> sigma_y >> p_y >> v_0 >> slits;
    
    int width = 10;
    std::cout << std::setw(width) << "Problem" << std::setw(width) << "h" << std::setw(width) << "deltat" << std::setw(width) << "T" << std::setw(width) << "x_c" << std::setw(width)
    << "sigma_x" << std::setw(width) << "p_x" << std::setw(width) << "omega" << std::setw(width) << "N" << std::setw(width) << "a_s" << std::endl;

    std::cout << std::setw(width) << Problem << std::setw(width) << h << std::setw(width) << deltat << std::setw(width) << T << std::setw(width) << x_c << std::setw(width)
    << sigma_x << std::setw(width) << p_x << std::setw(width) << omega << std::setw(width) << N << std::setw(width) << a_s << std::endl;


    CrankNicolson Crank(h, deltat, T, x_c, sigma_x, p_x, omega, N, a_s);

    //Crank.simulation_1D();
    Crank.print_m_Psi();


    Eigen::VectorXd Psi = Crank.get_m_Psi().real(); 

    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(98, -1,1);

    draw(x, Psi);
}


void draw(Eigen::VectorXd& x, Eigen::VectorXd& Psi){
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