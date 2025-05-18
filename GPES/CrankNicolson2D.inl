#ifndef CRANKNICOLSON2D_INL
#define CRANKNICOLSON2D_INL

#include "CrankNicolson.hpp"

GPES::CrankNicolson<Dimension::Two>::CrankNicolson(WaveFunction<Dimension::Two>& Psi, Grid<Dimension::Two>& grid, double delta_t,double T):  _delta_t(delta_t), _T(T) {
    _omega_x = grid.get_omega_x();
    _omega_y = grid.get_omega_y();
    _size_x = grid.get_size_of_grid_x();
    _size_y = grid.get_size_of_grid_y();
    _step_x = grid.get_step_size_x();
    _step_y = grid.get_step_size_y();
    _start_x = grid.get_start_position_x();
    _start_y = grid.get_start_position_y();
    
    _t_step = std::round(T/_delta_t) + 1;

    _Psi = Psi.get_wavefunction();
    _Fin = Eigen::VectorXcd::Zero(_size_x * _size_y);
    _U_ddi = Eigen::VectorXd::Zero(_size_x * _size_y);
    _V_ext = grid.get_potential();

    //Physics units
    _Num = Psi.get_Num();

    // Initialize interaction strengths
    init_g_scattering();
    init_g_ddi();
    init_g_lhy();

    _lambda_x = -1.0*_delta_t/(4*std::pow(_step_x,2));
    _lambda_y = -1.0*_delta_t/(4*std::pow(_step_y,2));

    double Lx = std::abs(2*_start_x);
    double Ly = std::abs(2*_start_y);
    F_ddi = std::make_unique<DipolarInteraction<Dimension::Two>>(_size_x,_size_y,Lx,Ly,1,_g_ddi);

    init_time_evolution_matrices();
}

void GPES::CrankNicolson<Dimension::Two>::init_g_scattering(){
    _g_scattering = 1; 
}

void GPES::CrankNicolson<Dimension::Two>::init_g_ddi(){
    _g_ddi = 1; 
}

void GPES::CrankNicolson<Dimension::Two>::init_g_lhy(){
    _g_lhy = 1; 
}

void GPES::CrankNicolson<Dimension::Two>::calculate_DDI(Eigen::VectorXcd& vec)
{
    if(this->F_ddi)
        F_ddi->compute_DDI_term(vec, _U_ddi);

}


//Function for constructing the right-hand and left-hand sides matrices from Crank Nicolson algorithm
void GPES::CrankNicolson<Dimension::Two>::init_time_evolution_matrices(){
    int mat_size = _size_x * _size_y;
    Eigen::VectorXcd a(mat_size);
    Eigen::VectorXcd b(mat_size);

    // Eigen::MatrixXcd mat = vec_to_mat(_Psi);
    // calculate_DDI(_Psi);

    for(int k = 0; k < _size_x; ++k){
        for(int l = 0; l < _size_y; ++l){
            int index = get_index(k,l);

            std::complex<double> U_potential = 1.0*(_delta_t * 0.5) * std::complex<double>(_V_ext(k,l));
            std::complex<double> U_contact = 1.0*(_delta_t * 0.5)*_g_scattering * std::norm(_Psi(index));
            std::complex<double> U_dd = 1.0 * (_delta_t*0.5) * _g_ddi * _U_ddi(index);
            std::complex<double> U_lhy = 1.0 * (_delta_t*0.5) * _g_lhy * std::pow(std::norm(_Psi(index)), 2.5);

            // //Real time evolution matrices
            // a(index) = (1.0 - 4.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));
            // b(index) = (1.0 + 4.0*this->m_lambda - 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));

            //Imaginary time evolution matrices
            a(index) = 1.0 - 2.0*_lambda_x - 2.0*_lambda_y + U_contact; // + U_potential + U_dd + U_lhy;
            b(index) = 1.0 + 2.0*_lambda_x + 2.0*_lambda_y - U_contact; // - U_potential  - U_dd - U_lhy;
        }
    }
    this->init_Mat_A(_lambda_x, _lambda_y,a);
    this->init_Mat_B(-1.0 * _lambda_x,-1.0 * _lambda_y,b); 
}

void GPES::CrankNicolson<Dimension::Two>::update_time_evolution_matrices(Eigen::VectorXcd &vec){
    int mat_size = _size_x * _size_y;
    Eigen::VectorXcd a(mat_size);
    Eigen::VectorXcd b(mat_size);

    // Eigen::MatrixXcd mat = vec_to_mat(vec);

    // calculate_DDI(vec);

    for(int k = 0; k < _size_x; ++k){
        for(int l = 0; l < _size_y; ++l){
            int index = get_index(k,l);

            std::complex<double> U_potential = 1.0 * (_delta_t * 0.5)*std::complex<double>(_V_ext(k,l));
            std::complex<double> U_contact = 1.0 * (_delta_t * 0.5)*std::norm(vec(index));
            std::complex<double> U_dd = 1.0 * (_delta_t*0.5) * _g_ddi * _U_ddi(index);
            std::complex<double> U_lhy = 1.0 * (_delta_t*0.5) * _g_lhy * std::pow(std::norm(vec(index)), 2.5);

            // //Real time evolution matrices
            // a(index) = (1.0 - 4.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));
            // b(index) = (1.0 + 4.0*this->m_lambda - 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));

            //Imaginary time evolution matrices
            a(index) = 1.0 - 2.0*_lambda_x - 2.0*_lambda_y + U_contact; // + U_potential + U_dd + U_lhy;
            b(index) = 1.0 + 2.0*_lambda_x + 2.0*_lambda_y - U_contact; // - U_potential  - U_dd - U_lhy;
        }
    }
    this->init_Mat_A(_lambda_x, _lambda_y,a);
    this->init_Mat_B(-1.0 * _lambda_x,-1.0 * _lambda_y,b); 
}


//Function for initialization left-hand side matrix according to Crank Nicolson algorithm
void GPES::CrankNicolson<Dimension::Two>::init_Mat_A(std::complex<double> r_x, std::complex<double> r_y, Eigen::VectorXcd& d){

    std::complex<double> r = r_x; // First time solution
    
    int S = d.size();
    int s = std::sqrt(S);

    typedef Eigen::Triplet<std::complex<double> > T;
    std::vector<T> tripletList;
    tripletList.reserve(5*s);

    tripletList.push_back(T(0, 0, d(0)));
    tripletList.push_back(T(0, 1, r));
    tripletList.push_back(T(S - 1, S - 2, r));
    tripletList.push_back(T(S - 1, S - 1, d(S-1)));

    for (int i = 1; i < S-1; ++i) {
        if(i + s  < S){
            tripletList.push_back(T(i-1, i+s-1, r));
            tripletList.push_back(T(i+s-1, i-1, r));
        } 
        tripletList.push_back(T(i, i,d(i)));

        if(i%s == 0){
            std::complex<double> z(0.,0.);
            tripletList.push_back(T(i, i - 1, z));
            tripletList.push_back(T(i-1, i, z));
        }
        else {
            tripletList.push_back(T(i, i - 1, r));
            tripletList.push_back(T(i - 1, i, r));
        }
    }

    Eigen::SparseMatrix<std::complex<double> > A(S,S);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    _A = A;
}

//Function for initialization right-hand side matrix according to Crank Nicolson algorithm
void GPES::CrankNicolson<Dimension::Two>::init_Mat_B(std::complex<double> r_x, std::complex<double> r_y, Eigen::VectorXcd& d) {

    std::complex<double> r = r_x; // First time solution


    int S = d.size();
    int s = std::sqrt(S);

    typedef Eigen::Triplet<std::complex<double> > T;
    std::vector<T> tripletList;
    tripletList.reserve(5*s);

    tripletList.push_back(T(0, 0, d(0)));
    tripletList.push_back(T(0, 1, r));
    tripletList.push_back(T(S - 1, S - 2, r));
    tripletList.push_back(T(S - 1, S - 1, d(S-1)));

    for (int i = 1; i < S-1; ++i) {
        if(i + s  < S){
            tripletList.push_back(T(i-1, i+s-1, r));
            tripletList.push_back(T(i+s-1, i-1, r));
        } 
        tripletList.push_back(T(i, i,d(i)));

        if(i%s == 0){
            std::complex<double> z(0.,0.);
            tripletList.push_back(T(i, i - 1, z));
            tripletList.push_back(T(i-1, i, z));
        }
        else {
            tripletList.push_back(T(i, i - 1, r));
            tripletList.push_back(T(i - 1, i, r));
        }
    }

    Eigen::SparseMatrix<std::complex<double> > B(S,S);
    B.setFromTriplets(tripletList.begin(), tripletList.end());
    _B = B;
}


// Time evolution simulation for 2D Gross-Pitaevskii equation
void GPES::CrankNicolson<Dimension::Two>::simulation(){

    int size = _Psi.size();

    Eigen::VectorXcd x(size);
    Eigen::VectorXcd b(size);

    x.setZero();
    b.setZero();
    // Prepare the right-hand side for the time-stepping
    b = _Psi;
    
    // Set up the sparse LU solver
    // Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solver;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double>>> solver;
    // Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double>>> solver;

    for (int i = 1; i < _t_step; ++i) {
        // Add the normalization function for 2D matrices, which would consider different spatial step in x and y directions
        normalize(b);
        update_time_evolution_matrices(b);

        solver.compute(_A);
        
        // Update the right-hand side vector b
        b = _B * b;

        // // Solve the system A * x = b
        x = solver.solve(b);
        // Update b for the next iteration
        b = x;

        if(i % 1 == 0)
            std::cout << "Simulation step: #" << i << '\n';

        _Fin = x;
        normalize(_Fin);

        // this->vec_Energy.push_back(calc_state_energy(this->m_Fin));
    }
}


void GPES::CrankNicolson<Dimension::Two>::normalize(Eigen::VectorXcd &vec){
    int size = vec.size();
    double psum = 0;
    for(int i = 0; i < size; ++i){
        psum += std::norm(vec(i));
    }
    std::complex<double> normalization_factor = std::sqrt(_Num) / std::sqrt(psum * _step_x * _step_y); // 

    vec *= std::abs(normalization_factor);
}

double GPES::CrankNicolson<Dimension::Two>::vec_norm(Eigen::VectorXcd &vec){
    int size = vec.size();
    double norm = 0;
    for(int i = 0; i < size; ++i){
        norm += std::norm(vec(i));
    }
    return norm * _step_x * _step_y;
}

double GPES::CrankNicolson<Dimension::Two>::vec_norm(Eigen::VectorXd &vec){
    int size = vec.size();
    double norm = 0;
    for(int i = 0; i < size; ++i){
        norm += std::pow(vec(i), 2);
    }
    return norm * _step_x * _step_y;
}

// //TODO ENERGY AND CHEMICAL POTENTIAL CALCULATION
// double GPES::CrankNicolson<Dimension::Two>::calc_state_energy(){
//     return calc_state_energy(_Fin);
// }

// double GPES::CrankNicolson<Dimension::Two>::calc_state_energy(Eigen::VectorXcd &vec){
//     int size = vec.size();
//     double energy = 0.0;
//     for(int i = 1; i < size-1; ++i){
//         std::complex<double> derivative = (vec(i + 1) - vec(i - 1)) / (2 * this->step);
//         double kinetic = std::norm(derivative) * 0.5;
//         double potential = this->m_V(i) * std::norm(vec(i));
//         double interaction = 0.5 * std::norm(vec(i)) * std::norm(vec(i));
//         energy += this->step * (kinetic + potential + interaction); 
//     }
//     return energy;
// }


void GPES::CrankNicolson<Dimension::Two>::get_final_state(GPES::WaveFunction<Dimension::Two>& fin){ 
    fin.set_Num(_Num);
    fin.set_size_of_grid_x(_size_x);
    fin.set_size_of_grid_y(_size_y);
    fin.set_start_position_x(_start_x);
    fin.set_start_position_y(_start_y);
    fin.set_start_position_x(_start_x);
    fin.set_start_position_y(_start_y);
    fin.set_step_size_x(_start_x);
    fin.set_step_size_y(_start_y);
    fin.set_vec(_Fin);
}

#endif