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

    l_z = 1.0 / std::sqrt(_omega_z);

    _a_s = Psi.get_a_s();
    _a_dd = Psi.get_a_dd();

    // Initialize interaction strengths
    calc_g_scattering(_a_s);
    calc_V_dd(_a_dd);
    calc_g_lhy(_a_s, _a_dd);

    _lambda_x       =   -1.0*_delta_t/(4*std::pow(_step_x,2));
    _lambda_y       =   -1.0*_delta_t/(4*std::pow(_step_y,2));

    double l_z      =   std::sqrt(_omega_x / _omega_z);

    double Lx       =   std::abs(2*_start_x);
    double Ly       =   std::abs(2*_start_y);
    F_ddi           =   std::make_unique<DipolarInteraction<Dimension::Two>>(_size_x,_size_y,Lx,Ly,l_z,_V_dd);

    vec_Energy.reserve(_t_step);

    calc_time_evolution_matrices(_Psi);
}

void  GPES::CrankNicolson<Dimension::Two>::calc_g_scattering(double a_s) {
    _g_scattering = std::sqrt(8. * M_PI) * a_s / l_z;
}

void GPES::CrankNicolson<Dimension::Two>::calc_V_dd(double a_dd) {
    _V_dd = std::sqrt(8. * M_PI) * a_dd / l_z;
}

void GPES::CrankNicolson<Dimension::Two>::calc_g_lhy(double a_s, double a_dd) {
    _g_lhy = (128. / ( 3 * std::pow(M_PI, 0.25) ) ) * std::sqrt(0.4) * std::pow(a_s, 2.5) * std::pow(_lam_z, 0.75) * (1 + 1.5 * std::pow((a_dd / a_s ), 2.));
}


Eigen::VectorXd GPES::CrankNicolson<Dimension::Two>::calc_FFT_DDI(const Eigen::VectorXcd& wave){
    const int Nx = _size_x;
    const int Ny = _size_y;
    const int N = Nx*Ny;
    const double pi = M_PI;
    double dx = _step_x;
    double dy = _step_y;
    Eigen::VectorXd phi_dd(N);
    assert(static_cast<int>(wave.size()) == N && "Density size must be Nx * Ny");

    const int Nkx = Nx;
    const int Nky = Ny/2 + 1;

    
    double* in = fftw_alloc_real(N);
    fftw_complex* out = fftw_alloc_complex(Nkx*Nky);

    for(int i = 0; i < N; ++i)
        in[i] = std::pow(wave(i).real(),2) + std::pow(wave(i).imag(),2); 

    fftw_plan forward  = fftw_plan_dft_r2c_2d(Nx, Ny, in, out, FFTW_ESTIMATE);
    fftw_plan backward = fftw_plan_dft_c2r_2d(Nx, Ny, out, in, FFTW_ESTIMATE);

    fftw_execute(forward);

    //Now compute kernel
    for(int i = 0; i < Nx; ++i){
        int ki = (i <= Nx/2)? i : i - Nx;
        double kx = 2.0 * pi * ki / (Nx*dx);
        for(int j = 0; j < Nky; ++j){
            int kj = j;
            double ky = 2.0 * pi * kj / (Ny * dy);
            
            double kappa = std::sqrt(kx*kx + ky*ky);

            int idx = i * Nky + j;
            double V_dd_k = 0.0;

            double k_r = kappa*l_z / std::sqrt(2);

            double K_exp = std::exp(k_r*k_r);
            double K_erfc = std::erfc(k_r);
            if(std::isnan(K_exp) || std::isnan(K_erfc)) std::cout<< "NAN value in erfc or exp" << std::endl;
            double F = 2.0 - 3.0 * std::sqrt(pi) * k_r * K_exp * K_erfc;

            V_dd_k = _V_dd * F;

            // Multiply by kernel
            out[idx][0] *= V_dd_k;
            out[idx][1] *= V_dd_k;
        }
    }

    
    // Inverse FFT
    fftw_execute(backward);

    // Normalize and copy output
    phi_dd.resize(N);
    for (int i = 0; i < N; ++i) {
        phi_dd[i] = in[i]/ N;
    }

    // Cleanup
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
    fftw_free(in);
    fftw_free(out);

    return phi_dd;
}



void GPES::CrankNicolson<Dimension::Two>::calculate_DDI(Eigen::VectorXcd& vec){
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
    // Eigen::VectorXd Phi_dd = calc_FFT_DDI(vec);

    for(int k = 0; k < _size_x; ++k){
        for(int l = 0; l < _size_y; ++l){
            int index = get_index(k,l);

            std::complex<double> U_potential = 1.0 * (_delta_t * 0.5)*std::complex<double>(_V_ext(k,l));
            std::complex<double> U_contact = 1.0 * (_delta_t * 0.5)*std::norm(vec(index));
            std::complex<double> U_dd = 1.0 * (_delta_t*0.5) * _U_ddi(index); //;Phi_dd(index)
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
void GPES::CrankNicolson<Dimension::Two>::simulation(std::string& outdir){

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

        if(i % 200 == 0){
            std::cout << "Simulation step: #" << i << '\n';
            auto s = std::to_string(i);
            std::string outfile = outdir + "/fin.csv";
            savecsv_wave(outfile,x);
        }
        _Fin = x;
        normalize(_Fin);
        double current_energy = calc_state_energy(_Fin);
        vec_Energy.push_back(current_energy);
        ++i;
    } while(simulation_stop(i));
    std::string outfile_en = outdir + "/energy_history.csv";
    savecsv_vec(outfile_en, vec_Energy);
}

inline bool GPES::CrankNicolson<Dimension::Two>::simulation_stop(int i)
{
    if(i > 4000)
        return false;

    double epsilon = 10e-7, last, before_last;

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
    fin.set_step_size_x(_step_x);
    fin.set_step_size_y(_step_y);
    fin.set_omega_x(_omega_x);
    fin.set_omega_y(_omega_y);
    fin.set_a_dd(_a_dd);
    fin.set_a_s(_a_s);
    fin.set_vec(_Fin);
}


double GPES::CrankNicolson<Dimension::Two>::calc_state_energy(){
    return calc_state_energy(_Fin);
}

double GPES::CrankNicolson<Dimension::Two>::calc_state_energy(Eigen::VectorXcd &vec){
    double energy = 0.0;
    // Eigen::VectorXd Phi_dd = calc_FFT_DDI(vec);
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
            double ddi = 0.5 * std::norm(vec(index))* _U_ddi(index);// * Phi_dd(index); //

            //LHY correction term energy
            double lhy = _g_lhy * 0.4 * std::pow(std::norm(vec(index)), 2.5);

            energy += _step_x * _step_y * (kinetic + potential + interaction + ddi + lhy); //  
        }
    }
    return energy;
}

double GPES::CrankNicolson<Dimension::Two>::calc_state_energy(GPES::WaveFunction<Dimension::Two>& vec){
    double energy = 0.0;
    // Eigen::VectorXd Phi_dd = calc_FFT_DDI(vec.get_wavefunction());
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
            double ddi = 0.5 * std::norm(vec(index))* _U_ddi(index) ;// *  Phi_dd(index); //

            //LHY correction term energy
            double lhy = _g_lhy * 0.4 * std::pow(std::norm(vec(index)), 2.5);

            energy += _step_x * _step_y * (kinetic + potential + interaction + ddi + lhy); //  
        }   
    }
    return energy;
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

void GPES::CrankNicolson<Dimension::Two>::savecsv_wave(std::string file_path, Eigen::VectorXcd& v){
    // 1) Ensure parent directory exists
    std::string dir = parent_dir(file_path);
    if (!dir.empty()) {
        // mkdir -p dir
        std::string cmd = "mkdir -p '" + dir + "'";
        if (std::system(cmd.c_str()) != 0) {
            throw std::runtime_error("Failed to create directory: " + dir);
        }
    }

    std::ofstream file(file_path, std::ios::out /*| std::ios::trunc is implicit*/);
    if (!file) {
        std::cerr << "Error: cannot open " << file_path << " for writing\n";
        return;
    }

    ParamList params{
        {"a_s",                _a_s                        },
        {"a_dd",               _a_dd                       },
        {"omega_x",            _omega_x                    },
        {"omega_y",            _omega_y                    },
        {"Num_of_particle",    static_cast<double>(_Num)   },
        {"size_of_grid_x",     static_cast<double>(_size_x)},
        {"size_of_grid_y",     static_cast<double>(_size_y)},
        {"step_size_x",        _step_x                     },
        {"step_size_y",        _step_y                     },
        {"grid_start_point_x", _start_x                    },
        {"grid_start_point_y", _start_y                    }
    };

    /* ---------- 1) metadata names ---------- */
    bool first = true;
    for (const auto& kv : params) {
        if (!first) file << ',';
        file << kv.first;
        first = false;
    }
    file << '\n';

    /* ---------- 2) metadata values ---------- */
    first = true;
    file << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const auto& kv : params) {
        if (!first) file << ',';
        file << kv.second;
        first = false;
    }
    file << "\n\n";                           // blank line before the data block

    /* ---------- 3) the actual wave-function ---------- */
    file << "index,real,imag\n";
    for (Eigen::Index i = 0; i < v.size(); ++i) {
        file << i << ',' << v[i].real() << ',' << v[i].imag() << '\n';
    }

    std::cout << "State have been saved" << std::endl;
}

void GPES::CrankNicolson<Dimension::Two>::savecsv_vec(std::string file_path, std::vector<double>& v){
    // 1) Ensure parent directory exists
    std::string dir = parent_dir(file_path);
    if (!dir.empty()) {
        // mkdir -p dir
        std::string cmd = "mkdir -p '" + dir + "'";
        if (std::system(cmd.c_str()) != 0) {
            throw std::runtime_error("Failed to create directory: " + dir);
        }
    }

    std::ofstream file(file_path, std::ios::out /*| std::ios::trunc is implicit*/);
    if (!file) {
        std::cerr << "Error: cannot open " << file_path << " for writing\n";
        return;
    }

    /* ---------- 3) the actual wave-function ---------- */
    file << "index,value\n";
    for (Eigen::Index i = 0; i < v.size(); ++i) {
        file << i << ',' << v[i] << '\n';
    }

    std::cout << "Vector have been saved" << std::endl;
}


void GPES::CrankNicolson<Dimension::Two>::print_param_of_eq(){
    int width = 15;
    std::cout << std::setw(width) << "a_s"<< std::setw(width) << "a_dd" << std::setw(width) << "g_scattering" << std::setw(width) << "V_dd" << std::setw(width) << "g_lhy" << std::setw(width) << "lambda" << std::endl;
    std::cout << std::setw(width) << _a_s << std::setw(width) << _a_dd << std::setw(width) << _g_scattering << std::setw(width) << _V_dd << std::setw(width) << _g_lhy << std::setw(width) << l_z << std::endl;
}


#endif