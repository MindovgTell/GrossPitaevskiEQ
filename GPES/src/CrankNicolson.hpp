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
#include <cassert>

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
    double _a_s, _a_dd, _omega, _omega_t, _lam, _g_scattering, _V_dd, _g_lhy, _delta_t, _start, _step, _l_perp;

    std::complex<double> _lambda_x;

    std::unique_ptr<DipolarInteraction<Dimension::One>> F_ddi;
    Eigen::VectorXd _U_ddi;
    // vector of the energies

    std::vector<double> vec_Energy;
    std::vector<double> vec_Chem_potential;


public:
    CrankNicolson(Grid<Dimension::One>& grid, WaveFunction<Dimension::One>& Psi, double deltat, double T);

    void calc_g_scattering(double a_s);
    void calc_V_dd(double a_dd);
    void calc_g_lhy(double a_s, double a_dd);

    //Function for calculating Dipole-Dipole Interaction
    void calculate_DDI(Eigen::VectorXcd& vec);
    std::vector<double> compute_Vdd_k();
    std::vector<double> compute_Phi_dd(const Eigen::VectorXcd& density);

    std::vector<double> calculate_DDI_not_FFT(const Eigen::VectorXcd& vec);
    double V_1DD(double x);
    double R_function(double xi);
    double V_dd_1D(double x);
    std::vector<double> compute_Phi_dd_realspace(const Eigen::VectorXcd& vec);
    // void init_time_evolution_matrices();
    void calc_time_evolution_matrices(const Eigen::VectorXcd& vec);


    void simulation(std::string& outdir);

 
    void init_Mat_A(std::complex<double> r,Eigen::VectorXcd& d);
    void init_Mat_B(std::complex<double> r,Eigen::VectorXcd& d);


    void normalize(Eigen::VectorXcd &vec);


    Eigen::VectorXd  prob(Eigen::VectorXcd &vec);

    double vec_norm(Eigen::VectorXcd &vec);
    double vec_norm(Eigen::VectorXd &vec);

    double calc_state_energy();
    double calc_state_energy(Eigen::VectorXcd &vec);
    double calc_state_energy(WaveFunction<Dimension::One>& vec);

    bool simulation_stop(int i); // Simulation stop function due to the small energy difference between steps


    // double calc_state_chem_potential();
    // double calc_state_chem_potential(Eigen::VectorXcd &vec);


    //********************//Getters funcitons//*****************//

    void print_Mat_A_dim();
    void print_Mat_B_dim();
    void m_Psi_len();
    void m_Fin_len();
    void print_Mat_A();
    void print_Mat_B();
    void print_param_of_eq();

    double get_g_s () {return _g_scattering;} 
    double get_g_lhy () {return _g_lhy;}
    double get_V_dd () {return _V_dd;}  

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

    void save_state(std::string filename, Eigen::VectorXcd& vec);
    void savecsv_wave(std::string file_path, Eigen::VectorXcd& vec);
    void savecsv_vec(std::string file_path, std::vector<double>& vec);
    int num_isnan();
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
    double _a_s, _a_dd, _g_scattering, _V_dd, _g_lhy, _lam_z, _lam_y;
    double _delta_t, _h_step, _omega_x, _omega_y, _omega_z, _chem_potential, _start_x, _start_y, _step_x, _step_y, l_z;
    std::complex<double> _lambda_x, _lambda_y;

    // vector of the energies

    std::vector<double> vec_Energy;

public: 
    //2D constructor
    CrankNicolson(WaveFunction<Dimension::Two>& Psi, Grid<Dimension::Two>& grid, double T, double delta_t);
    
    int get_index(int i,int j) { return i * _size_x + j;}


    void calc_g_scattering(double a_s);
    void calc_V_dd(double a_dd);
    void calc_g_lhy(double a_s, double a_dd);

    //Function for calculating 2D DDI
    void calculate_DDI(Eigen::VectorXcd &vec);

    // void init_time_evolution_matrices();
    void calc_time_evolution_matrices(Eigen::VectorXcd &vec);

    void init_Mat_A(std::complex<double> r_x, std::complex<double> r_y, Eigen::VectorXcd& d);
    void init_Mat_B(std::complex<double> r_x, std::complex<double> r_y, Eigen::VectorXcd& d);

    Eigen::MatrixXd vec_to_mat(const Eigen::VectorXd& vec);

    Eigen::VectorXd prob(Eigen::VectorXcd&  vec);

    void save_matrix_to_csv(std::string filename, Eigen::MatrixXd mat);
    void save_vector_to_csv(std::string filename, Eigen::VectorXd mat);
    void savecsv_wave(std::string file_path, Eigen::VectorXcd& vec);
    void savecsv_vec(std::string file_path, std::vector<double>& vec);

    void simulation(std::string& outdir);
    bool simulation_stop(int i); // Simulation stop function due to the small energy difference between steps

    void normalize(Eigen::VectorXcd &vec);
    double vec_norm(Eigen::VectorXcd &vec);
    double vec_norm(Eigen::VectorXd &vec);


    double calc_state_energy();
    double calc_state_energy(Eigen::VectorXcd &vec);
    double calc_state_energy(WaveFunction<Dimension::Two> &vec);


    Eigen::VectorXd calc_FFT_DDI(const Eigen::VectorXcd& wave);

    //Gettes Funcitons

    double get_g_s () {return _g_scattering;} 
    double get_g_lhy () {return _g_lhy;}
    double get_C_dd () {return _V_dd;}  

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

    Eigen::VectorXd get_U_ddi() {return _U_ddi;}

    void get_final_state(WaveFunction<Dimension::Two>& fin);

    Eigen::VectorXd get_m_V();
    void get_V_size(){
        std::cout << "Size of V vector: " << _V_ext.size() << std::endl;
    }
    Eigen::VectorXd real(Eigen::VectorXcd& vec);
    Eigen::VectorXd imag(Eigen::VectorXcd& vec);

    std::vector<double> get_vec_Energy(){ return vec_Energy;}

    void save_state(std::string filename, Eigen::VectorXcd& vec);
    void save_sim_param(std::string file_path);

    void print_param_of_eq();


};


#include "CrankNicolson2D.inl"

};


#endif 