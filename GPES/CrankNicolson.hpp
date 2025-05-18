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

#include "definitions.hpp"
#include "grid.hpp"
#include "wavefunciton.hpp"
#include "dipoleinteraction.hpp"

namespace GPES{
    
template <Dimension Dim>
class CrankNicolson;

//********************************/***********/********************************//
//                                                                             //
//***************************/One dimensional solver/**************************//
//                                                                             //
//********************************/***********/********************************//


template <>
class CrankNicolson<Dimension::One>{

private:
    Eigen::SparseMatrix<std::complex<double> > _A;
    Eigen::SparseMatrix<std::complex<double> > _B;

    Eigen::VectorXcd _Psi;
    Eigen::VectorXcd _Fin;
    Eigen::VectorXd _V_ext;


    unsigned int _Num, _T, _t_step, _size;
    double _g_scattering, _g_ddi, _g_lhy, _delta_t, _start, _step;

    std::complex<double> _lambda_x;

    std::unique_ptr<DipolarInteraction<Dimension::One>> F_ddi;
    Eigen::VectorXd _U_ddi;
    // vector of the energies

    std::vector<double> vec_Energy;
    std::vector<double> vec_Chem_potential;


public:
    CrankNicolson(Grid<Dimension::One>& grid, WaveFunction<Dimension::One>& Psi, double deltat, double T);

    void init_g_scattering();
    void init_g_ddi();
    void init_g_lhy();

    //Function for calculating Dipole-Dipole Interaction
    void calculate_DDI(Eigen::VectorXcd& vec);

    void init_time_evolution_matrices();
    void update_time_evolution_matrices(Eigen::VectorXcd& vec);


    void simulation();

 
    void init_Mat_A(std::complex<double> r,Eigen::VectorXcd& d);
    void init_Mat_B(std::complex<double> r,Eigen::VectorXcd& d);


    void normalize(Eigen::VectorXcd &vec);


    Eigen::VectorXd  prob(Eigen::VectorXcd &vec);

    double vec_norm(Eigen::VectorXcd &vec);
    double vec_norm(Eigen::VectorXd &vec);

    double calc_state_energy();
    double calc_state_energy(Eigen::VectorXcd &vec);

    // double calc_state_chem_potential();
    // double calc_state_chem_potential(Eigen::VectorXcd &vec);


    //********************//Getters funcitons//*****************//

    void print_Mat_A_dim();
    void print_Mat_B_dim();
    void m_Psi_len();
    void m_Fin_len();
    void print_Mat_A();
    void print_Mat_B();

    Eigen::VectorXcd get_m_Psi();
    Eigen::VectorXd get_m_Psi_prob();

    // Eigen::VectorXcd get_Fin();
    // Eigen::VectorXd get_Fin_prob();
    void get_final_state(WaveFunction<Dimension::One>& fin);

    Eigen::VectorXd get_m_V();
    // void get_m_V_size(){
    //     std::cout << "Size of V vector: " << _V_ext.size() << std::endl;
    // }
    Eigen::VectorXd real(Eigen::VectorXcd& vec);
    Eigen::VectorXd imag(Eigen::VectorXcd& vec);

    std::vector<double> get_vec_Energy(){ return vec_Energy;}

};


#include "CrankNicolson1D.inl"



//********************************/***********/********************************//
//                                                                             //
//**************************/Two dimensional solver/***************************//
//                                                                             //
//********************************/***********/********************************//


template <>
class CrankNicolson<Dimension::Two>{
private:
    Eigen::SparseMatrix<std::complex<double> > _A;
    Eigen::SparseMatrix<std::complex<double> > _B;

    Eigen::VectorXcd _Psi;
    Eigen::VectorXcd _Fin;
    Eigen::MatrixXd _V_ext;
    Eigen::VectorXd _U_ddi;
    std::unique_ptr<DipolarInteraction<Dimension::Two>> F_ddi;

    int _size_x, _size_y, _T, _t_step, _Num;
    double _g_scattering, _g_ddi, _g_lhy;
    double _delta_t, _h_step, _omega_x, _omega_y, _chem_potential, _start_x, _start_y, _step_x, _step_y;
    std::complex<double> _lambda_x, _lambda_y;

    // vector of the energies

    std::vector<double> vec_Energy;
    std::vector<double> vec_Chem_potential;

public: 
    //2D constructor
    CrankNicolson(WaveFunction<Dimension::Two>& Psi, Grid<Dimension::Two>& grid, double T, double delta_t);
    
    int get_index(int i,int j) { return i * _size_x + j;}


    void init_g_scattering();
    void init_g_ddi();
    void init_g_lhy();

    //Function for calculating 2D DDI
    void calculate_DDI(Eigen::VectorXcd &vec);

    void init_time_evolution_matrices();
    void update_time_evolution_matrices(Eigen::VectorXcd &vec);

    void init_Mat_A(std::complex<double> r_x, std::complex<double> r_y, Eigen::VectorXcd& d);
    void init_Mat_B(std::complex<double> r_x, std::complex<double> r_y, Eigen::VectorXcd& d);

    Eigen::MatrixXd vec_to_mat(const Eigen::VectorXd& vec);

    Eigen::VectorXd prob(Eigen::VectorXcd&  vec);

    void save_matrix_to_csv(std::string filename, Eigen::MatrixXd mat);
    void save_vector_to_csv(std::string filename, Eigen::VectorXd mat);

    void simulation();

    void normalize(Eigen::VectorXcd &vec);
    double vec_norm(Eigen::VectorXcd &vec);
    double vec_norm(Eigen::VectorXd &vec);


    double calc_state_energy();
    double calc_state_energy(Eigen::VectorXcd &vec);

    //Gettes Funcitons

    void print_Mat_A_dim();
    void print_Mat_B_dim();
    void m_Psi_len();
    void m_Fin_len();
    void print_Mat_A();
    void print_Mat_B();

    Eigen::VectorXcd get_m_Psi();
    Eigen::VectorXd get_m_Psi_prob();

    Eigen::VectorXcd get_m_Fin();
    Eigen::VectorXd get_m_Fin_prob();

    void get_final_state(WaveFunction<Dimension::Two>& fin);

    Eigen::VectorXd get_m_V();
    void get_V_size(){
        std::cout << "Size of V vector: " << _V_ext.size() << std::endl;
    }
    Eigen::VectorXd real(Eigen::VectorXcd& vec);
    Eigen::VectorXd imag(Eigen::VectorXcd& vec);

    std::vector<double> get_vec_Energy(){ return vec_Energy;}


};


#include "CrankNicolson2D.inl"

};


#endif 