#ifndef CRANKNICOLSON2D_INL
#define CRANKNICOLSON2D_INL

#include "CrankNicolson.hpp"

GPES::CrankNicolson<Dimension::Two>::CrankNicolson(WaveFunction<Dimension::Two>& Psi, Grid<Dimension::Two>& grid, double delta_t,double T):  _delta_t(delta_t), _T(T) {
    _omega_x        =   grid.get_omega_x();
    _omega_y        =   grid.get_omega_y();
    _omega_z        =   grid.get_omega_z();
    _size_x         =   grid.get_size_of_grid_x();
    _size_y         =   grid.get_size_of_grid_y();
    _step_x         =   grid.get_step_size_x();
    _step_y         =   grid.get_step_size_y();
    _start_x        =   grid.get_start_position_x();
    _start_y        =   grid.get_start_position_y();
    
    _t_step         =   std::round(T/_delta_t) + 1;

    _Psi            =   Psi.get_wavefunction();
    _Fin            =   Eigen::VectorXcd::Zero(_size_x * _size_y);
    _U_ddi          =   Eigen::VectorXd::Zero(_size_x * _size_y);
    _V_ext          =   grid.get_potential();

    //Physics units
    _Num            =   Psi.get_Num();
    _lam_y          =   _omega_y / _omega_x;
    _lam_z          =   _omega_z / _omega_x;

    // Initialize interaction strengths
    calc_g_scattering(Psi.get_a_s());
    calc_V_dd(Psi.get_a_dd());
    calc_g_lhy(Psi.get_a_s(), Psi.get_a_dd());

    _lambda_x       =   -1.0*_delta_t/(4*std::pow(_step_x,2));
    _lambda_y       =   -1.0*_delta_t/(4*std::pow(_step_y,2));

    double Lx       =   std::abs(2*_start_x);
    double Ly       =   std::abs(2*_start_y);
    F_ddi           =   std::make_unique<DipolarInteraction<Dimension::Two>>(_size_x,_size_y,Lx,Ly,0.1414,_V_dd);

    vec_Energy.reserve(_t_step);

    calc_time_evolution_matrices(_Psi);
}

void  GPES::CrankNicolson<Dimension::Two>::calc_g_scattering(double a_s) {
    _g_scattering = std::sqrt(8. * M_PI) * a_s * std::sqrtf(_lam_z);
}

void GPES::CrankNicolson<Dimension::Two>::calc_V_dd(double a_dd) {
    
}

void GPES::CrankNicolson<Dimension::Two>::calc_g_lhy(double a_s, double a_dd) {
    _g_lhy = (128. / 3) * std::sqrt(M_PI) * std::pow(a_s, 2.5) * (1 + 1.5 * std::pow((a_dd / a_s ), 2.)) ; 
}

void GPES::CrankNicolson<Dimension::Two>::calculate_DDI(Eigen::VectorXcd& vec)
{
    if(this->F_ddi)
        F_ddi->compute_DDI_term(vec, _U_ddi);
}


//Function for constructing the right-hand and left-hand sides matrices from Crank Nicolson algorithm
// void GPES::CrankNicolson<Dimension::Two>::init_time_evolution_matrices(){
//     int mat_size = _size_x * _size_y;
//     Eigen::VectorXcd a(mat_size);
//     Eigen::VectorXcd b(mat_size);
//     a.setZero();
//     b.setZero();
//     // Eigen::MatrixXcd mat = vec_to_mat(_Psi);
//     calculate_DDI(_Psi);
//     for(int k = 0; k < _size_x; ++k){
//         for(int l = 0; l < _size_y; ++l){
//             int index = get_index(k,l);
//             std::complex<double> U_potential = 1.0*(_delta_t * 0.5) * std::complex<double>(_V_ext(k,l));
//             std::complex<double> U_contact = 1.0*(_delta_t * 0.5) * _g_scattering * std::norm(_Psi(index));
//             std::complex<double> U_dd = 1.0 * (_delta_t*0.5) * _U_ddi(index);
//             std::complex<double> U_lhy = 1.0 * (_delta_t*0.5) * _g_lhy * std::pow(std::norm(_Psi(index)), 1.5);
//             // //Real time evolution matrices
//             // a(index) = (1.0 - 4.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));
//             // b(index) = (1.0 + 4.0*this->m_lambda - 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));
//             //Imaginary time evolution matrices
//             a(index) = 1.0 - 2.0*_lambda_x - 2.0*_lambda_y + U_contact + U_potential;// + U_dd + U_lhy;
//             b(index) = 1.0 + 2.0*_lambda_x + 2.0*_lambda_y - U_contact - U_potential;// - U_dd - U_lhy;
//         }
//     }
//     this->init_Mat_A(_lambda_x, _lambda_y,a);
//     this->init_Mat_B(-1.0 * _lambda_x,-1.0 *_lambda_y,b); 
// }

void GPES::CrankNicolson<Dimension::Two>::calc_time_evolution_matrices(Eigen::VectorXcd &vec){
    int mat_size = _size_x * _size_y;
    Eigen::VectorXcd a(mat_size);
    Eigen::VectorXcd b(mat_size);
    
    a.setZero();
    b.setZero();

    // Eigen::MatrixXcd mat = vec_to_mat(vec);

    calculate_DDI(vec);

    for(int k = 0; k < _size_x; ++k){
        for(int l = 0; l < _size_y; ++l){
            int index = get_index(k,l);

            std::complex<double> U_potential = 1.0 * (_delta_t * 0.5)*std::complex<double>(_V_ext(k,l));
            std::complex<double> U_contact = 1.0 * (_delta_t * 0.5)*std::norm(vec(index));
            std::complex<double> U_dd = 1.0 * (_delta_t*0.5) * _U_ddi(index);
            std::complex<double> U_lhy = 1.0 * (_delta_t*0.5) * _g_lhy * std::pow(std::norm(vec(index)), 1.5);

            // //Real time evolution matrices
            // a(index) = (1.0 - 4.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));
            // b(index) = (1.0 + 4.0*this->m_lambda - 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));

            //Imaginary time evolution matrices
            a(index) = 1.0 - 2.0*_lambda_x - 2.0*_lambda_y + U_contact + U_potential + U_dd + U_lhy;
            b(index) = 1.0 + 2.0*_lambda_x + 2.0*_lambda_y - U_contact - U_potential - U_dd - U_lhy;
        }
    }
    this->init_Mat_A(_lambda_x,_lambda_y,a);
    this->init_Mat_B(-1.0 *_lambda_x,-1.0 *_lambda_y,b); 
}

// // // Old version
// // //Function for initialization left-hand side matrix according to Crank Nicolson algorithm
// void GPES::CrankNicolson<Dimension::Two>::init_Mat_A(std::complex<double> r_x, std::complex<double> r_y, Eigen::VectorXcd& d){
//     std::complex<double> r = r_x; // First time solution
//     int S = d.size();
//     int s = std::sqrt(S);
//     typedef Eigen::Triplet<std::complex<double> > T;
//     std::vector<T> tripletList;
//     tripletList.reserve(5*s);
//     tripletList.push_back(T(0, 0, d(0)));
//     tripletList.push_back(T(0, 1, r));
//     tripletList.push_back(T(S - 1, S - 2, r));
//     tripletList.push_back(T(S - 1, S - 1, d(S-1)));
//     for (int i = 1; i < S-1; ++i) {
//         if(i + s  < S){
//             tripletList.push_back(T(i-1, i+s-1, r));
//             tripletList.push_back(T(i+s-1, i-1, r));
//         } 
//         tripletList.push_back(T(i, i,d(i)));
//         if(i%s == 0){
//             std::complex<double> z(0.,0.);
//             tripletList.push_back(T(i, i - 1, z));
//             tripletList.push_back(T(i-1, i, z));
//         }
//         else {
//             tripletList.push_back(T(i, i - 1, r));
//             tripletList.push_back(T(i - 1, i, r));
//         }
//     }
//     Eigen::SparseMatrix<std::complex<double> > A(S,S);
//     A.setFromTriplets(tripletList.begin(), tripletList.end());
//     _A = A;
// }

//Function for initialization left-hand side matrix according to Crank Nicolson algorithm
void GPES::CrankNicolson<Dimension::Two>::init_Mat_A(std::complex<double> r_x, std::complex<double> r_y, Eigen::VectorXcd& d) {
    int S = d.size();
    int s = std::sqrt(S);

    typedef Eigen::Triplet<std::complex<double>> T;
    std::vector<T> tripletList;
    tripletList.reserve(5*S);  // Reserve more space for safety

    for (int i = 0; i < S; ++i) {
        // Main diagonal
        tripletList.push_back(T(i, i, d(i)));

        int x = i % s;
        int y = i / s;

        // x-direction neighbors (left and right)
        if (x > 0) tripletList.push_back(T(i, i-1, r_x));  // left
        if (x < s-1) tripletList.push_back(T(i, i+1, r_x)); // right

        // y-direction neighbors (top and bottom)
        if (y > 0) tripletList.push_back(T(i, i-s, r_y));   // top
        if (y < s-1) tripletList.push_back(T(i, i+s, r_y)); // bottom
    }

    Eigen::SparseMatrix<std::complex<double>> A(S,S);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    _A = A;
}

// // Old version
// //Function for initialization right-hand side matrix according to Crank Nicolson algorithm
// void GPES::CrankNicolson<Dimension::Two>::init_Mat_B(std::complex<double> r_x, std::complex<double> r_y, Eigen::VectorXcd& d) {
//     std::complex<double> r = r_x; // First time solution
//     int S = d.size();
//     int s = std::sqrt(S);
//     typedef Eigen::Triplet<std::complex<double> > T;
//     std::vector<T> tripletList;
//     tripletList.reserve(5*s);
//     tripletList.push_back(T(0, 0, d(0)));
//     tripletList.push_back(T(0, 1, r));
//     tripletList.push_back(T(S - 1, S - 2, r));
//     tripletList.push_back(T(S - 1, S - 1, d(S-1)));
//     for (int i = 1; i < S-1; ++i) {
//         if(i + s  < S){
//             tripletList.push_back(T(i-1, i+s-1, r));
//             tripletList.push_back(T(i+s-1, i-1, r));
//         } 
//         tripletList.push_back(T(i, i,d(i)));
//         if(i%s == 0){
//             std::complex<double> z(0.,0.);
//             tripletList.push_back(T(i, i - 1, z));
//             tripletList.push_back(T(i-1, i, z));
//         }
//         else {
//             tripletList.push_back(T(i, i - 1, r));
//             tripletList.push_back(T(i - 1, i, r));
//         }
//     }
//     Eigen::SparseMatrix<std::complex<double> > B(S,S);
//     B.setFromTriplets(tripletList.begin(), tripletList.end());
//     _B = B;
// }

//Function for initialization right-hand side matrix according to Crank Nicolson algorithm
void GPES::CrankNicolson<Dimension::Two>::init_Mat_B(std::complex<double> r_x, std::complex<double> r_y, Eigen::VectorXcd& d) {
    int S = d.size();
    int s = std::sqrt(S);

    typedef Eigen::Triplet<std::complex<double>> T;
    std::vector<T> tripletList;
    tripletList.reserve(5*S);  // Reserve more space for safety

    for (int i = 0; i < S; ++i) {
        // Main diagonal
        tripletList.push_back(T(i, i, d(i)));

        int x = i % s;
        int y = i / s;

        // x-direction neighbors (left and right)
        if (x > 0) tripletList.push_back(T(i, i-1, r_x));  // left
        if (x < s-1) tripletList.push_back(T(i, i+1, r_x)); // right

        // y-direction neighbors (top and bottom)
        if (y > 0) tripletList.push_back(T(i, i-s, r_y));   // top
        if (y < s-1) tripletList.push_back(T(i, i+s, r_y)); // bottom
    }

    Eigen::SparseMatrix<std::complex<double>> B(S,S);
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
    // Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double>>> solver;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double>>> solver;

    // for (int i = 1; i < _t_step; ++i) {
    //     // Add the normalization function for 2D matrices, which would consider different spatial step in x and y directions
    //     normalize(b);
    //     calc_time_evolution_matrices(b);
    //     solver.compute(_A);    
    //     // Update the right-hand side vector b
    //     b = _B * b;
    //     // // Solve the system A * x = b
    //     x = solver.solve(b);
    //     // Update b for the next iteration
    //     b = x;
    //     if(i % 50 == 0)
    //         std::cout << "Simulation step: #" << i << '\n';
    //     _Fin = x;
    //     normalize(_Fin);
    //     this->vec_Energy.push_back(calc_state_energy(_Fin));
    // }

    int i = 1;

    do {
        normalize(b);
        calc_time_evolution_matrices(b);

        solver.compute(_A);
        
        // Update the right-hand side vector b
        b = _B * b;

        // // Solve the system A * x = b
        x = solver.solve(b);
        // Update b for the next iteration
        b = x;

        if(i % 50 == 0)
            std::cout << "Simulation step: #" << i << '\n';

        _Fin = x;
        normalize(_Fin);
        double current_energy = calc_state_energy(_Fin);
        vec_Energy.push_back(current_energy);
        ++i;
        if(i % 250 == 0){
            // std::string file;
            // std::stringstream ss;
            save_state("../../res/Fin_states/sim.csv", _Fin);
        }
    } while(simulation_stop(i));
}

inline bool GPES::CrankNicolson<Dimension::Two>::simulation_stop(int i)
{
    if(i > 500)
        return false;

    double epsilon = 10e-10, last, before_last;

    if (vec_Energy.size() >= 2) {
        last = *vec_Energy.rbegin();
        before_last = *(vec_Energy.rbegin() + 1);
    }
    else 
        return true;
    
    double diff = std::abs((last - before_last) / last);

    if(diff < epsilon)
        return false;
    return true;
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



double GPES::CrankNicolson<Dimension::Two>::calc_state_energy(){
    return calc_state_energy(_Fin);
}

double GPES::CrankNicolson<Dimension::Two>::calc_state_energy(Eigen::VectorXcd &vec){
    double energy = 0.0;
    for(int i = 1; i < _size_x-1; ++i){
        for(int j = 1; j < _size_y-1; ++j) {

            int index = get_index(i,j);
            int x_index_minus_one = get_index(i - 1, j);
            int x_index_plus_one = get_index(i + 1, j);
            int y_index_minus_one = get_index(i, j - 1);
            int y_index_plus_one = get_index(i, j + 1);

            std::complex<double> derivative_x = (vec(x_index_minus_one) - vec(x_index_plus_one)) * (1. / (2 * _step_x)); 
            std::complex<double> derivative_y = (vec(y_index_minus_one) - vec(y_index_plus_one)) * (1. / (2 * _step_y));

            double kinetic = (std::norm(derivative_x) + std::norm(derivative_y))* 0.5;
            double potential = _V_ext(i,j) * std::norm(vec(index));
            double interaction = 0.5 * _g_scattering * std::norm(vec(index)) * std::norm(vec(index));

            //Dipole-Dipole interaction energy
            double ddi = 0.5 * _U_ddi(index) * std::norm(vec(index)); 

            //LHY correction term energy
            double lhy = _g_lhy * 0.5 * std::pow(std::norm(vec(index)), 3);

            energy += _step_x * _step_y * (kinetic + potential + interaction + ddi + lhy); //  
        }
    }
    return energy;
}

double GPES::CrankNicolson<Dimension::Two>::calc_state_energy(GPES::WaveFunction<Dimension::Two>& vec){
    double energy = 0.0;
    for(int i = 1; i < _size_x-1; ++i){
        for(int j = 1; j < _size_y-1; ++j){
            int index = get_index(i,j);
            int x_index_minus_one = get_index(i - 1, j);
            int x_index_plus_one = get_index(i + 1, j);
            int y_index_minus_one = get_index(i, j - 1);
            int y_index_plus_one = get_index(i, j + 1);

            std::complex<double> derivative_x = (vec(x_index_minus_one) - vec(x_index_plus_one)) * (1. / (2 * _step_x)); 
            std::complex<double> derivative_y = (vec(y_index_minus_one) - vec(y_index_plus_one)) * (1. / (2 * _step_y));

            double kinetic = (std::norm(derivative_x) + std::norm(derivative_y))* 0.5;
            double potential = _V_ext(i,j) * std::norm(vec(index));
            double interaction = 0.5 * _g_scattering * std::norm(vec(index)) * std::norm(vec(index));

            //Dipole-Dipole interaction energy
            double ddi = 0.5 * _U_ddi(index) * std::norm(vec(index)); 

            //LHY correction term energy
            double lhy = _g_lhy * 0.5 * std::pow(std::norm(vec(index)), 2.5);

            energy += _step_x * _step_y * (kinetic + potential + interaction + ddi + lhy); //  
        }   
    }
    return energy;
}


void GPES::CrankNicolson<Dimension::Two>::print_Mat_A_dim(){
    std::cout << "Matrix A number of cols: " << _A.cols() << std::endl;
    std::cout << "Matrix A number of rows: " << _A.rows() << std::endl;
}

void GPES::CrankNicolson<Dimension::Two>::print_Mat_B_dim(){
    std::cout << "Matrix B number of cols: " << _B.cols() << std::endl;
    std::cout << "Matrix B number of rows: " << _B.rows() << std::endl;
}


void GPES::CrankNicolson<Dimension::Two>::print_Mat_A(){
    std::cout << Eigen::MatrixXcd(_A) << std::endl;
}

void GPES::CrankNicolson<Dimension::Two>::print_Mat_B(){
    std::cout << Eigen::MatrixXcd(_B) << std::endl;
}


void GPES::CrankNicolson<Dimension::Two>::save_state(std::string filename, Eigen::VectorXcd& vec){
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ",","\n");
    std::ofstream file(filename);
    if(file.is_open()){
        //file << vec.format(CSVFormat) << '\n';
        file << vec.format(CSVFormat);
        file.close();
    }
    else
        std::cout << "Smth goes wrong with open file for w" << std::endl;
}

#endif