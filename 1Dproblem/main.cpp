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

void simulation1D(std::string inputfile);
void simulation2D(std::string inputfile);

void draw1(Eigen::VectorXd& x, Eigen::VectorXd& Psi);
void draw2(Eigen::VectorXd& x, Eigen::VectorXd& Psi, Eigen::VectorXd& Fin);
void draw3(Eigen::VectorXd& x, Eigen::VectorXd& Psi, Eigen::VectorXd& Fin, Eigen::VectorXd& V);

int main(int argc, char const *argv[]){

    if(argc != 2){
        std::string executable = argv[0];

        std::cerr << "Error: Wrong number of input parameters" << '\n';
        std::cerr << "Usage:" << executable << " input_file.txt" << '\n';
        return 1;
    }

    simulation1D(argv[1]);

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

    if (str_stream >> Problem >> h >> deltat >> T >> x_c >> sigma_x >> p_x >> omega >> N >> a_s >> start) {
        std::cout << "Values read successfully!" << std::endl;
    } else {
        std::cout << "Error reading values from string!" << std::endl;
    }



    // str_stream >> Problem >> h >>deltat >>T >> x_c >> sigma_x >> p_x >> y_c >> sigma_y >> p_y >> v_0 >> slits;
    
    int width = 10;
    std::cout << std::setw(width) << "Problem" << std::setw(width) << "h" << std::setw(width) << "deltat" << std::setw(width) << "T" << std::setw(width) << "x_c" << std::setw(width)
    << "sigma_x" << std::setw(width) << "p_x" << std::setw(width) << "omega" << std::setw(width) << "N" << std::setw(width) << "a_s" << std::setw(width) << "start" << std::endl;

    std::cout << std::setw(width) << Problem << std::setw(width) << h << std::setw(width) << deltat << std::setw(width) << T << std::setw(width) << x_c << std::setw(width)
    << sigma_x << std::setw(width) << p_x << std::setw(width) << omega << std::setw(width) << N << std::setw(width) << a_s << std::setw(width) << start << std::endl;


    CrankNicolson Crank(h, deltat, T, x_c, sigma_x, p_x, omega, N, a_s, start);

    Eigen::VectorXd Psi = Crank.get_m_Psi_prob();   

    Crank.simulation_1D();

    Eigen::VectorXd Fin = Crank.get_m_Fin_prob(); 

    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(998, start, -1 * start);

    Eigen::VectorXd V = Crank.get_m_V();
    Eigen::VectorXcd TM = Crank.TM_state();
    Eigen::VectorXd TM_pr = Crank.prob_1D(TM);

    // draw2(x, Psi, Fin);
    draw3(x, Psi, Fin, V);
    draw3(x, Psi, Fin, TM_pr);
    // Crank.get_m_V_size();

    std::cout << '\n' << '\n' << std::endl;

    // Crank.print_Mat_A();


    // double R_tf = std::pow(1.5 * 100, 1/3.);
    // double potential = std::pow(( R_tf),2)/ 0.5;

    // std::cout << potential << std::endl;

    // std::cout << '\n' << '\n' << std::endl;

    // std::cout << '\n' << '\n' << std::endl;

    // Crank.print_Mat_B();

}


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



