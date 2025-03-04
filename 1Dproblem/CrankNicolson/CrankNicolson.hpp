#ifndef CRANKNICOLSON_HPP   
#define CRANKNICOLSON_HPP   

#include <iostream>
#include <complex>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>


class CrankNicolson
{
private:
    Eigen::SparseMatrix<std::complex<double> > m_A;
    Eigen::SparseMatrix<std::complex<double> > m_B;

    Eigen::VectorXcd m_Psi;
    Eigen::MatrixXd m_V;
    int m_size, m_T, t_step;
    double m_delta_t, m_h_step,V_0;
    std::complex<double> m_r;

public:
    //Constructor
    //2D constructor
    CrankNicolson(double h, double deltat, double T, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, double v_0, int slits=0);

    //1D constructor
    CrankNicolson(double h, double deltat, double T, double x_c, double sigma_x, double p_x, double v_0);

    //***************/1DFunctions/***************//
    //Thomas-Fermi ansatz
    std::complex<double> thomas_fermi_state();
    //Initialization of starting 1D state
    void init_start_state_1D(double x_c, double sigma_x, double p_x);

    void init_time_evolution_matrices_1D();

    void init_Mat_A_1D(std::complex<double> r,Eigen::VectorXcd& d);
    void init_Mat_B_1D(std::complex<double> r,Eigen::VectorXcd& d);

    Eigen::VectorXd create_potential_1D();
    
    //***************/2DFunctions/***************//
    // Gauss Wave funciton
    std::complex<double> gauss_wave_packet(double sigma_x, double sigma_y, double x, double y, double x_c, double y_c, double p_x, double p_y = 0);

    //Methods for creating matrixes
    void init_start_state_2D(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y);
    
    int get_m_index(int i,int j, int M);

    void init_time_evolution_matrices();

    void init_Mat_A(std::complex<double> r,Eigen::VectorXcd& d);
    void init_Mat_B(std::complex<double> r,Eigen::VectorXcd& d);

    Eigen::MatrixXd vec_to_mat(const Eigen::VectorXd& vec);

    Eigen::VectorXd prob(Eigen::VectorXcd&  vec);

    void save_to_csv(std::string filename, Eigen::MatrixXd mat);
    void simulation();
    
    // Functions for create potential
    Eigen::MatrixXd create_potential_box();
    Eigen::MatrixXd create_one_slit();
    Eigen::MatrixXd create_double_slit();
    Eigen::MatrixXd create_triple_slit();
};


#endif 