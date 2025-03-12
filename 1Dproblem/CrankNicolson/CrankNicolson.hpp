#ifndef CRANKNICOLSON_HPP   
#define CRANKNICOLSON_HPP   

#include <iostream>
#include <complex>
#include <fstream>
#include <cmath>
#include <numbers>
#include <Eigen/Dense>
#include <Eigen/Sparse>


class CrankNicolson
{
private:
    Eigen::SparseMatrix<std::complex<double> > m_A;
    Eigen::SparseMatrix<std::complex<double> > m_B;

    Eigen::VectorXcd m_Psi;
    Eigen::VectorXd  m_out;
    Eigen::MatrixXd m_V;
    int m_size, m_T, t_step;
    double m_delta_t, m_h_step, V_0, m_omega, m_N, m_g, m_chem_potential;
    std::complex<double> m_r;

public:

//********************************/1DFunctions/********************************//
    
    //1D constructor
    CrankNicolson(double h, double deltat, double T, double x_c, double sigma_x, double p_x, double omega, double N, double a_s);

    //Thomas-Fermi ansatz
    double thomas_fermi_state(double x);
    //Gauss wave function
    //std::complex<double> CrankNicolson::gauss_wave_packet_1D(double sigma_x, double x, double x_c, double p_x);

    //Initialization of starting 1D state
    void init_start_state_1D(double x_c, double sigma_x, double p_x);

    void init_time_evolution_matrices_1D();
    void init_chem_potential(double omega, double N, double a_s);

    Eigen::VectorXd create_harmonic_potential_1D();
    Eigen::VectorXd create_potential_1D();

    void simulation_1D();
    Eigen::VectorXd  prob_1D(Eigen::VectorXcd &vec);

    void init_Mat_A_1D(std::complex<double> r,Eigen::VectorXcd& d);
    void init_Mat_B_1D(std::complex<double> r,Eigen::VectorXcd& d);
    
//********************************/2DFunctions/********************************//

    //2D constructor
    CrankNicolson(double h, double deltat, double T, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, double v_0, int slits=0);

    // Gauss Wave funciton
    std::complex<double> gauss_wave_packet(double sigma_x, double sigma_y, double x, double y, double x_c, double y_c, double p_x, double p_y = 0);

    //Methods for creating matrixes
    void init_start_state_2D(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y);
    
    int get_m_index(int i,int j, int M);

    void init_time_evolution_matrices_2D();

    void init_Mat_A_2D(std::complex<double> r,Eigen::VectorXcd& d);
    void init_Mat_B_2D(std::complex<double> r,Eigen::VectorXcd& d);

    Eigen::MatrixXd vec_to_mat(const Eigen::VectorXd& vec);

    Eigen::VectorXd prob(Eigen::VectorXcd&  vec);

    void save_matrix_to_csv(std::string filename, Eigen::MatrixXd mat);
    void save_vector_to_csv(std::string filename, Eigen::VectorXd mat);

    void simulation();
    
    // Functions for create potential
    Eigen::MatrixXd create_potential_box();
    Eigen::MatrixXd create_one_slit();
    Eigen::MatrixXd create_double_slit();
    Eigen::MatrixXd create_triple_slit();


//*****************************/Printing and Getter Functions/******************************//
    void print_Mat_A_dim();
    void print_Mat_B_dim();
    void print_m_Psi();
    void print_Mat_A();
    void print_Mat_B();
    Eigen::VectorXcd get_m_Psi();
    Eigen::VectorXd get_m_out();

};


#endif 