#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <complex>
#include "./CrankNicolson/CrankNicolson.hpp"
//#include "./matplotlibcpp.h"

using namespace std::complex_literals;

//namespace plt = matplotlibcpp;

void simulation(std::string inputfile);

void test(double arg);

int main(int argc, char const *argv[]){

    if(argc != 2){
        std::string executable = argv[0];

        std::cerr << "Error: Wrong number of input parameters" << '\n';
        std::cerr << "Usage:" << executable << " input_file.txt" << '\n';
        return 1;
    }

    //simulation(argv[1]);
    
    test(std::stoi(argv[1]));

    return 0;
}


void simulation(std::string inputfile){
    std::ifstream input_data(inputfile);

    std::string line;
    std::getline(input_data,line);//Skip first line in file

    double Problem,h,deltat,T,x_c,sigma_x,p_x,omega,N , a_s;

    std::getline(input_data,line);
    std::stringstream str_stream(line);
    str_stream >> Problem >> h >>deltat >>T >> x_c >> sigma_x >> p_x >> omega >> N >> a_s;

    int width = 10;
    std::cout << std::setw(width) << "Problem" << std::setw(width) << "h" << std::setw(width) << "deltat" << std::setw(width) << "T" << std::setw(width) << "x_c" << std::setw(width)
    << "sigma_x" << std::setw(width) << "p_x" << std::setw(width) << "omega" << std::setw(width) << "N" << std::setw(width) << "a_s" << std::endl;

    std::cout << std::setw(width) << Problem << std::setw(width) << h << std::setw(width) << deltat << std::setw(width) << T << std::setw(width) << x_c << std::setw(width)
    << sigma_x << std::setw(width) << p_x << std::setw(width) << omega << std::setw(width) << N << std::setw(width) << a_s << std::endl;

    CrankNicolson Crank(h, deltat, T, x_c, sigma_x, p_x, omega, N, a_s);

    //Crank.simulation_1D();
    Crank.print_m_Psi();
}

void test(double arg){
    // CrankNicolson Crank(0.05, 2.5e-04, 0.0075, 0.25,0.01,200, 0.1, 100, 1e-05);
    // Crank.print_m_Psi();

    // std::vector<double> y{1,2.1,2.5,5};
    // std::vector<double> x{2,0.5,3,2.7};

    // plt::figure();
    // plt::plot(x,y, std::string("bo-"),{{"label", "data trend"}});
    // plt::xlabel("time [s]");
    // plt::ylabel("observation [m]");
    // plt::legend();
    // plt::show(); 
}