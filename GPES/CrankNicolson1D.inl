#ifndef CRANKNICOLSON_INL
#define CRANKNICOLSON_INL

#include "CrankNicolson.hpp"


GPES::CrankNicolson<Dimension::One>::CrankNicolson(Grid<Dimension::One>& grid, WaveFunction<Dimension::One>& Psi, double deltat, double T): _delta_t(deltat), _T(T) {
    // Initialize parameters of the grid
    _size = grid.get_size_of_grid(); // number of grid nodes
    _step = grid.get_step_size();
    _start = grid.get_start_position();
    _V_ext = grid.get_potential();

    _t_step = std::round(T/_delta_t) + 1;

    _Psi = Psi.get_wavefunction();
    _Fin = Eigen::VectorXcd::Zero(_size);
    _U_ddi = Eigen::VectorXd::Zero(_size);



    //Physics units
    _Num = Psi.get_Num();
    
    // Initialize interaction strengths

    _g_scattering = Psi.get_g_scat();
    _g_ddi = -0.003; 
    _g_lhy = 0.005; 

    _lambda_x = std::complex<double>(-1.0*_delta_t/(4*std::pow(_step,2)),0);

    double L = std::abs(_start*2);
    double confinment_length = 1;
    F_ddi = std::make_unique<DipolarInteraction<Dimension::One>>(_size, L, confinment_length, _g_ddi);

    init_time_evolution_matrices();
}

void GPES::CrankNicolson<Dimension::One>::init_g_ddi(){
    
}

void GPES::CrankNicolson<Dimension::One>::init_g_lhy(){
    
}


void GPES::CrankNicolson<Dimension::One>::calculate_DDI(Eigen::VectorXcd& vec){
    if(this->F_ddi)
        F_ddi->compute_DDI_term(vec, _U_ddi);
}



//Function for constructing the right-hand and left-hand sides matrices from Crank Nicolson algorithm
void GPES::CrankNicolson<Dimension::One>::init_time_evolution_matrices(){
    Eigen::VectorXcd a(_size);
    Eigen::VectorXcd b(_size);

    calculate_DDI(_Psi);
    for(int i = 0; i < _size; ++i){

        std::complex<double> U_potential = 1.0 * (_delta_t*0.5)*std::complex<double>(_V_ext(i));
        std::complex<double> U_scattering =  1.0 * (_delta_t*0.5 * _g_scattering)*std::norm(_Psi(i));
        std::complex<double> U_dd = 1.0 * (_delta_t*0.5 * _g_ddi) * _U_ddi(i);
        std::complex<double> U_lhy = 1.0 * (_delta_t*0.5) * this->_g_lhy * std::pow(std::norm(_Psi(i)), 2.5);
        
        // // Real time evolution matrices
        // a(i) = (1.0 - 2.0*this->m_lambda_x + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(i)) + 1.0i*(m_delta_t/2)*std::pow(std::abs(m_Psi(i)),2));
        // b(i) = (1.0 + 2.0*this->m_lambda_x + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(i)) + 1.0i*(m_delta_t/2)*std::pow(std::abs(m_Psi(i)),2));
        
        // // Imaginary time evolution matrices
        // a(i) = 1.0 + 2.0*this->m_lambda_x + 1.0 * (m_delta_t*0.5)*std::complex<double>(m_V(i)) + 1.0 * (m_delta_t*0.5 * m_g)*std::norm(m_Psi(i));
        // b(i) = 1.0 - 2.0*this->m_lambda_x + 1.0 * (m_delta_t*0.5)*std::complex<double>(m_V(i)) + 1.0 * (m_delta_t*0.5 * m_g)*std::norm(m_Psi(i));

        // TEST ADDED DDI
        a(i) = 1.0 - 2.0*_lambda_x + U_scattering + U_potential ;//+ U_dd + U_lhy;
        b(i) = 1.0 + 2.0*_lambda_x - U_scattering - U_potential ;//- U_dd - U_lhy;
    }

    this->init_Mat_A(_lambda_x, a);
    this->init_Mat_B(-1.0 * _lambda_x, b); 
}

void GPES::CrankNicolson<Dimension::One>::update_time_evolution_matrices(Eigen::VectorXcd& vec){
    int size = vec.size();
    Eigen::VectorXcd a(size);
    Eigen::VectorXcd b(size);

    calculate_DDI(vec);
    for(int i = 0; i < size; ++i){

        std::complex<double> U_potential = 1.0 * (_delta_t*0.5)*std::complex<double>(_V_ext(i));
        std::complex<double> U_scattering =  1.0 * (_delta_t*0.5 * _g_scattering)*std::norm(vec(i));
        std::complex<double> U_dd = 1.0 * (_delta_t*0.5 * _g_ddi) * _U_ddi(i);
        std::complex<double> U_lhy = 1.0 * (_delta_t*0.5) * this->_g_lhy * std::pow(std::norm(vec(i)), 2.5);
        
        
        // // Real time evolution matrices
        // a(i) = (1.0 - 2.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(i)) + 1.0i*(m_delta_t/2)*std::pow(std::abs(m_Psi(i)),2));
        // b(i) = (1.0 + 2.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(i)) + 1.0i*(m_delta_t/2)*std::pow(std::abs(m_Psi(i)),2));

        // // Imaginary time evolution matrices
        // a(i) = 1.0 + 2.0*this->m_lambda_x + (m_delta_t*0.5)*std::complex<double>(m_V(i)) + (m_delta_t*0.5 * m_g)*std::norm(vec(i));
        // b(i) = 1.0 - 2.0*this->m_lambda_x + (m_delta_t*0.5)*std::complex<double>(m_V(i)) + (m_delta_t*0.5 * m_g)*std::norm(vec(i));
    
        // TEST ADDED DDI
        a(i) = 1.0 - 2.0*_lambda_x + U_scattering + U_potential;// + U_dd + U_lhy;
        b(i) = 1.0 + 2.0*_lambda_x - U_scattering - U_potential;// - U_dd - U_lhy;
    }

    this->init_Mat_A(_lambda_x, a);
    this->init_Mat_B(-1.0 * _lambda_x, b); 
}



//Function for initialization left-hand side matrix according to Crank Nicolson algorithm
void GPES::CrankNicolson<Dimension::One>::init_Mat_A(std::complex<double> r, Eigen::VectorXcd& d){
    int S = d.size();

    typedef Eigen::Triplet<std::complex<double> > T;
    std::vector<T> tripletList;
    tripletList.reserve(std::pow(S,2));
    tripletList.push_back(T(0, 0, d(0)));
    
    for (int i = 1; i < S; ++i) {
        tripletList.push_back(T(i, i,d(i)));
        tripletList.push_back(T(i, i - 1, r));
        tripletList.push_back(T(i - 1, i, r));
    }

    Eigen::SparseMatrix<std::complex<double> > A(S,S);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    _A = A;
}

//Function for initialization right-hand side matrix according to Crank Nicolson algorithm
void GPES::CrankNicolson<Dimension::One>::init_Mat_B(std::complex<double> r, Eigen::VectorXcd& d) {
    int S = d.size();

    typedef Eigen::Triplet<std::complex<double> > T;
    std::vector<T> tripletList;
    tripletList.reserve(std::pow(S,2));

    tripletList.push_back(T(0, 0, d(0)));
    // tripletList.push_back(T(S - 1, S - 1, d(S-1)));

    for (int i = 1; i < S; ++i) {
        tripletList.push_back(T(i, i,d(i)));
        tripletList.push_back(T(i, i - 1, r));
        tripletList.push_back(T(i - 1, i, r));
    }

    Eigen::SparseMatrix<std::complex<double> > B(S,S);
    B.setFromTriplets(tripletList.begin(), tripletList.end());
    _B = B;
}

// Time evolution simulation for 1D Gross-Pitaevskii equation
void GPES::CrankNicolson<Dimension::One>::simulation(){

    Eigen::VectorXcd x(_size);
    Eigen::VectorXcd b(_size);

    x.setZero();
    b.setZero();
    // Prepare the right-hand side for the time-stepping
    b = _Psi;
    
    // Set up the sparse LU solver
    // Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solver;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double>>> solver;
    // Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double>>> solver;

    for (int i = 1; i < _t_step; ++i) {

        normalize(b);
        update_time_evolution_matrices(b);

        solver.compute(_A);
        
        // Update the right-hand side vector b
        b = _B * b;

        // // Solve the system A * x = b
        x = solver.solve(b);
        // Update b for the next iteration
        b = x;

        if(i % 100 == 0)
            std::cout << "Simulation step: #" << i << '\n';

        _Fin = x;
        normalize(_Fin);

        this->vec_Energy.push_back(calc_state_energy(_Fin));
    }
}

void GPES::CrankNicolson<Dimension::One>::normalize(Eigen::VectorXcd &vec){
    int size = vec.size();
    double psum = 0;
    for(int i = 0; i != size; ++i){
        psum += std::norm(vec(i));
    }
    std::complex<double> normalization_factor = std::sqrt(_Num) / std::sqrt(psum * _step); // 

    vec *= std::abs(normalization_factor);
}

double GPES::CrankNicolson<Dimension::One>::vec_norm(Eigen::VectorXcd &vec){
    int size = vec.size();
    double norm = 0;
    for(int i = 0; i < size; ++i){
        norm += std::norm(vec(i));
    }
    return norm * this->_step;
}

double GPES::CrankNicolson<Dimension::One>::vec_norm(Eigen::VectorXd &vec){
    int size = vec.size();
    double norm = 0;
    for(int i = 0; i < size; ++i){
        norm += std::pow(vec(i), 2);
    }
    return norm * this->_step;
}

double GPES::CrankNicolson<Dimension::One>::calc_state_energy(){
    return calc_state_energy(_Fin);
}

double GPES::CrankNicolson<Dimension::One>::calc_state_energy(Eigen::VectorXcd &vec){
    int size = vec.size();
    double energy = 0.0;
    for(int i = 1; i < size-1; ++i){
        std::complex<double> derivative = (vec(i + 1) - vec(i - 1)) * (1. / (2 * _step));
        double kinetic = std::norm(derivative) * 0.5;
        double potential = _V_ext(i) * std::norm(vec(i));
        double interaction = 0.5 * std::norm(vec(i)) * std::norm(vec(i));
        energy += _step * (kinetic + potential + interaction); 
    }
    return energy;
}


// double CrankNicolson::calc_state_chem_potential(){
//     return calc_state_chem_potential(this->m_Fin);
// }

// double CrankNicolson::calc_state_chem_potential(Eigen::VectorXcd &vec){
//     int size = vec.size();
//     double chem_potential = 0.0;
//     for(int i = 1; i < size-1; ++i){
//         std::complex<double> derivative = (vec(i + 1) - vec(i - 1)) / (2 * this->step);
//         double kinetic = std::norm(derivative) * 0.5;
//         double potential = this->m_V(i) * std::norm(vec(i));
//         double interaction = 0.5 * std::norm(vec(i)) * std::norm(vec(i));
//         chem_potential += this->step * (kinetic + potential + interaction); 
//     }
//     return chem_potential;
// }

void GPES::CrankNicolson<Dimension::One>::get_final_state(GPES::WaveFunction<Dimension::One>& fin){ 
    fin.set_Num(_Num);
    fin.set_size_of_grid(_size);
    fin.set_start_position(_start);
    fin.set_step_size(_start);
    fin.set_vec(_Fin);
}



#endif