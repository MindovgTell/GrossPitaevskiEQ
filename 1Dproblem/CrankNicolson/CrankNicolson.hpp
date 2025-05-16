#ifndef CRANKNICOLSON_HPP   
#define CRANKNICOLSON_HPP   

#include <iostream>
#include <complex>
#include <fstream>
#include <cmath>
#include <numbers>
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fftw3.h>
#include "DipoleInteraction.hpp"


//---------//GLOBAL COMMENT//----------//
// All functions with same names 
// should be combined into ones without
// _2D or _1D endings by introducing 
// new variable wich will determine 
// properties and functionality of the function.
//---------//--------------//----------//


class CrankNicolson
{
private:
    Eigen::SparseMatrix<std::complex<double> > m_A;
    Eigen::SparseMatrix<std::complex<double> > m_B;

    std::unique_ptr<DipolarInteraction2D> F_ddi;
    Eigen::VectorXd m_V_ddi;

    Eigen::VectorXcd m_Psi;
    Eigen::VectorXcd m_Fin;
    Eigen::MatrixXd m_V;
    int m_size, m_T, t_step;
    double m_g_scattering, m_g_dipole ,m_g_lhy;
    double m_delta_t,m_h_step, V_0, m_omega_x, m_omega_y, m_N, m_chem_potential, _start, step;
    std::complex<double> m_lambda_x, m_lambda_y;

    // vector of the energies

    std::vector<double> vec_Energy;
    std::vector<double> vec_Chem_potential;

public:

//********************************/***********/********************************//
//                                                                             //
//********************************/1DFunctions/********************************//
//                                                                             //
//********************************/***********/********************************//

    //1D constructor
    CrankNicolson(double h, double deltat, double T, double x_c, double sigma_x, double p_x, double omega, double N, double a_s, double start);

    //Thomas-Fermi ansatz
    void init_chem_potential(double omega, double N, double a_s);
    double thomas_fermi_state_1D(double x);

    //Gauss wave function
    std::complex<double> gauss_wave_packet_1D(double sigma_x, double x, double x_c, double p_x);
    double square_func(double x);

    //Initialization of starting 1D state
    void init_start_state_1D(double x_c, double sigma_x, double p_x);
    void init_V_ddi();

    //Function for calculating Dipole-Dipole Interaction
    void calculate_1D_DDI(Eigen::VectorXcd& vec);

    void init_time_evolution_matrices_1D();
    void update_time_evolution_matrices_1D(Eigen::VectorXcd &vec);


    Eigen::VectorXd create_harmonic_potential_1D();
    Eigen::VectorXd create_potential_1D();

    void simulation_1D();
    Eigen::VectorXd  prob_1D(Eigen::VectorXcd &vec);
 
    void init_Mat_A_1D(std::complex<double> r,Eigen::VectorXcd& d);
    void init_Mat_B_1D(std::complex<double> r,Eigen::VectorXcd& d);


    void normalize_1D(Eigen::VectorXcd &vec);

    Eigen::VectorXcd TM_state();
    double vec_norm_1D(Eigen::VectorXcd &vec);
    double vec_norm_1D(Eigen::VectorXd &vec);

    double calc_state_energy();
    double calc_state_energy(Eigen::VectorXcd &vec);

    double calc_state_chem_potential();
    double calc_state_chem_potential(Eigen::VectorXcd &vec);



    //Methods for calculating interaction strength parameters

    void init_g_scatt_1D(double a_s);
    void init_g_dipole_1D(double dipole_moment);
    void init_g_lhy_1D(double dipole_moment, double a_s);
    
//********************************/***********/********************************//
//                                                                             //
//********************************/2DFunctions/********************************//
//                                                                             //
//********************************/***********/********************************//

    //2D constructor
    CrankNicolson(double h, double deltat, double T, double x_c, double sigma_x, double y_c, double sigma_y, double omega_x, double omega_y, double N, double a_s, double start);
    
    void init_g_scatt_2D();
    void init_g_dipole_2D();
    void init_g_lhy_2D();
            
    int get_m_index(int i,int j, int M);

    double thomas_fermi_radius_x();
    double thomas_fermi_radius_y();

    void init_chem_potential(double omega_x, double omega_y, double N, double g);

    // Gauss Wave funciton
    std::complex<double> gauss_wave_packet_2D(double x, double y, double x_c, double y_c, double sigma_x, double sigma_y); //, double p_x, double p_y
    // Thomas-Fermi function
    double thomas_fermi_state_2D(double x, double y);
    Eigen::VectorXcd TM_state_2D();

    //Function for calculating 2D DDI
    Eigen::MatrixXcd calculate_2D_DDI(Eigen::MatrixXcd& vec);

    //Methods for creating matrixes
    void init_start_state_2D(double x_c, double y_c, double sigma_x, double sigma_y); //, double p_x, double p_y


    void init_time_evolution_matrices_2D();
    void update_time_evolution_matrices_2D(Eigen::VectorXcd &vec);

    void init_Mat_A_2D(std::complex<double> r_x, std::complex<double> r_y, Eigen::VectorXcd& d);
    void init_Mat_B_2D(std::complex<double> r_x, std::complex<double> r_y, Eigen::VectorXcd& d);

    Eigen::MatrixXd vec_to_mat(const Eigen::VectorXd& vec);
    Eigen::MatrixXcd vec_to_mat(const Eigen::VectorXcd& vec);

    Eigen::VectorXd prob(Eigen::VectorXcd&  vec);

    void save_matrix_to_csv(std::string filename, Eigen::MatrixXd mat);
    void save_vector_to_csv(std::string filename, Eigen::VectorXd mat);

    void simulation_2D();

    void normalize_2D(Eigen::VectorXcd &vec);
    double vec_norm_2D(Eigen::VectorXcd &vec);
    double vec_norm_2D(Eigen::VectorXd &vec);


    double calc_state_energy_2D();
    double calc_state_energy_2D(Eigen::VectorXcd &vec);
    
    // Functions for create potential
    Eigen::MatrixXd create_potential_box();
    Eigen::MatrixXd create_harmonic_potential_2D();


//********************************/***********/********************************//
//                                                                             //
//*****************************/Getters Functions/*****************************//
//                                                                             //
//********************************/***********/********************************//
    void print_Mat_A_dim();
    void print_Mat_B_dim();
    void m_Psi_len();
    void m_Fin_len();
    void print_Mat_A();
    void print_Mat_B();

    Eigen::VectorXcd get_m_Psi();
    Eigen::VectorXd get_m_Psi_prob();

    Eigen::VectorXd get_m_V_ddi();

    Eigen::VectorXcd get_m_Fin();
    Eigen::VectorXd get_m_Fin_prob();

    Eigen::VectorXd get_m_V();
    void get_m_V_size(){
        std::cout << "Size of V vector: " << this->m_V.size() << std::endl;
    }
    Eigen::VectorXd real(Eigen::VectorXcd& vec);
    Eigen::VectorXd imag(Eigen::VectorXcd& vec);

    std::vector<double> get_vec_Energy(){ return this->vec_Energy;}


};  


#endif 