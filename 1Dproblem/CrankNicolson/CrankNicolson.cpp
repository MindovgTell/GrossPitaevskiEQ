#include "CrankNicolson.hpp"


// #define _IMAGINARY_UNIT std::complex i(0,1);

using namespace std::complex_literals;

//********************************/***********/********************************//
//                                                                             //
//********************************/1DFunctions/********************************//
//                                                                             //
//********************************/***********/********************************//

CrankNicolson::CrankNicolson(double h, double deltat, double T, double x_c, double sigma_x, double p_x, double omega, double N, double a_s, double start){
    int M = 1/h;
    this->m_h_step = h;
    this->m_delta_t = deltat;
    this->m_T = T;
    this->V_0 = 1e+10;
    this->_start = start;
    this->step = (-2*_start)/M;
    this->t_step = std::round(T/deltat) + 1;
    this->m_lambda_x = -1.0*m_delta_t/(4*std::pow(this->step,2));
    this->m_lambda_y = std::complex(0,0);
    this->m_size = M;
    this->m_omega_x = omega;
    this->m_omega_y = 0.0;

    this->m_Fin = Eigen::VectorXcd(m_size).setZero();
    //Physical parameters
    this->m_N = N;
    //this->m_g = 2 * M_PI * a_s;
    this->m_g_scattering = 1;
    this->m_g_lhy = 0.0005;
    this->m_g_dipole = 0.0005;

    this->init_chem_potential(omega, N, a_s);

    this->m_V = this->create_harmonic_potential_1D();
    // this->m_V = this->create_potential_1D();

    this->init_start_state_1D(x_c,sigma_x,p_x);

    this->init_time_evolution_matrices_1D();
}

void CrankNicolson::init_chem_potential(double omega, double N, double a_s){
    double R_tf = std::cbrt(1.5 * this->m_N);
    double potential = R_tf * R_tf * 0.5;
    this->m_chem_potential = potential;
}

double CrankNicolson::thomas_fermi_state_1D(double x){
    double R_tf = std::cbrt(1.5 * this->m_N);
    double out = this->m_chem_potential * (1 - std::pow(x/R_tf,2.));
    if(out > 0)
        return std::sqrt(out);
    else 
        return 0;
}

double CrankNicolson::square_func(double x){
    double R_tf = std::cbrt(1.5 * this->m_N);
    double high = this->m_N / R_tf;

    double out = (this->m_chem_potential) * (1 - std::pow(x/R_tf,2.));

    if(out > 0)
        return high;
    else 
        return 0;
    
}

std::complex<double> CrankNicolson::gauss_wave_packet_1D(double sigma_x,  double x,  double x_c, double p_x){
    std::complex<double> i(0, 1); // Define the imaginary unit
    double exponent = -(pow(x-x_c,2) / (2 * pow(sigma_x,2))); //this->m_size/2 //-x_c
    std::complex<double> phase = i * p_x * (x - x_c); 
    return std::exp(exponent); //+ phase
}

void CrankNicolson::init_start_state_1D(double x_c, double sigma_x, double p_x){
    int size = this->m_size;
    Eigen::VectorXcd U(size);
    std::complex<double> psum = 0;
    double x = this->_start;
    for(int i = 0; i < size; ++i){
        // std::complex<double> c = thomas_fermi_state_1D(x - x_c); // Add parameters for Thomas-Fermi function
        std::complex<double> c = gauss_wave_packet_1D(sigma_x, x, x_c, p_x);
        // std::complex<double> c = square_func(x);

        U(i) = c;
        psum += std::norm(c);

        x += this->step;
    }

    std::complex<double> normalization_factor = std::sqrt(m_N) / std::sqrt(psum * this->step);
    this->m_Psi = U * std::abs(normalization_factor);
}

Eigen::VectorXcd CrankNicolson::TM_state(){
    int size = this->m_size;
    Eigen::VectorXcd U(size);

    double x = this->_start;
    std::complex<double> psum = 0;

    for(int i = 0; i < size; ++i){
        x +=  this->step;
        std::complex<double> c = thomas_fermi_state_1D(x); // Add parameters for Thomas-Fermi function
        U(i) = c;
        psum += std::norm(c);
    }

    std::complex<double> normalization_factor = std::sqrt(m_N) / std::sqrt(psum * this->step);

    U *= std::abs(normalization_factor);
    return U;
}

//Function for calculating dipole-dipole interaction

std::complex<double> CrankNicolson::calculate_1D_DDI(int grid_position, Eigen::VectorXcd& vec){
    std::complex<double> answ(0,0);
    std::complex<double> V_1d;

    double l_transv =  7.;
    double dipole_mom = 1.;
    double coeff = 1.; //Its supposed to be coefficient which containe both dipole_moment etc.

    int size = this->m_size;
    double x_prime = this->_start;

    double x = this->_start + this->step * grid_position;


    for(int i = 0; i < size; ++i){
        double pos = x - x_prime;
        double alfa = std::abs(pos) / (std::sqrt(2) * l_transv);

        //Firstly calculating V_1d
        V_1d = -0.001 * (2 * l_transv * std::abs(pos) - std::sqrt(2 * M_PI) * (std::pow(l_transv,2) + std::pow(pos, 2)) * std::erfc(alfa) * std::exp(alfa*alfa));

        answ += V_1d * std::norm(vec(grid_position)) * this->step;

        x_prime += this->step;
    }

    // std::cout << answ << std::endl;

    return answ;
}



//Function for constructing the right-hand and left-hand sides matrices from Crank Nicolson algorithm
void CrankNicolson::init_time_evolution_matrices_1D(){
    Eigen::VectorXcd a(this->m_size);
    Eigen::VectorXcd b(this->m_size);

    // Eigen::VectorXcd V_dd = ;

    for(int i = 0; i < this->m_size; ++i){

        std::complex<double> U_potential =  1.0 * (m_delta_t*0.5)*std::complex<double>(m_V(i)) + 1.0 * (m_delta_t*0.5 * m_g_scattering)*std::norm(m_Psi(i));
        std::complex<double> U_dd = 1.0 * (m_delta_t*0.5 * m_g_dipole) * calculate_1D_DDI(i, m_Psi);
        std::complex<double> U_lhy = 1.0 * (m_delta_t*0.5) * this->m_g_lhy * std::pow(std::norm(m_Psi(i)), 2.5);

        // // Real time evolution matrices
        // a(i) = (1.0 - 2.0*this->m_lambda_x + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(i)) + 1.0i*(m_delta_t/2)*std::pow(std::abs(m_Psi(i)),2));
        // b(i) = (1.0 + 2.0*this->m_lambda_x + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(i)) + 1.0i*(m_delta_t/2)*std::pow(std::abs(m_Psi(i)),2));
        
        // // Imaginary time evolution matrices
        // a(i) = 1.0 + 2.0*this->m_lambda_x + 1.0 * (m_delta_t*0.5)*std::complex<double>(m_V(i)) + 1.0 * (m_delta_t*0.5 * m_g)*std::norm(m_Psi(i));
        // b(i) = 1.0 - 2.0*this->m_lambda_x + 1.0 * (m_delta_t*0.5)*std::complex<double>(m_V(i)) + 1.0 * (m_delta_t*0.5 * m_g)*std::norm(m_Psi(i));

        // TEST ADDED DDI
        a(i) = 1.0 - 2.0*this->m_lambda_x + U_potential + U_dd + U_lhy;
        b(i) = 1.0 + 2.0*this->m_lambda_x - U_potential - U_dd - U_lhy;
    }

    this->init_Mat_A_1D(m_lambda_x,a);
    this->init_Mat_B_1D(-1.0 * m_lambda_x,b); 
}

void CrankNicolson::update_time_evolution_matrices_1D(Eigen::VectorXcd &vec){

    int size = this->m_size; 
    Eigen::VectorXcd a(size);
    Eigen::VectorXcd b(size);
    for(int i = 0; i < size; ++i){

        std::complex<double> U_potential = 1.0 * (m_delta_t*0.5)*std::complex<double>(m_V(i)) + 1.0 * (m_delta_t*0.5 * m_g_scattering)*std::norm(vec(i));
        std::complex<double> U_dd = 1.0 * (m_delta_t*0.5 * m_g_dipole) * calculate_1D_DDI(i, vec);
        std::complex<double> U_lhy = 1.0 * (m_delta_t*0.5) * this->m_g_lhy * std::pow(std::norm(vec(i)), 2.5);
        
        
        // // Real time evolution matrices
        // a(i) = (1.0 - 2.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(i)) + 1.0i*(m_delta_t/2)*std::pow(std::abs(m_Psi(i)),2));
        // b(i) = (1.0 + 2.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(i)) + 1.0i*(m_delta_t/2)*std::pow(std::abs(m_Psi(i)),2));

        // // Imaginary time evolution matrices
        // a(i) = 1.0 + 2.0*this->m_lambda_x + (m_delta_t*0.5)*std::complex<double>(m_V(i)) + (m_delta_t*0.5 * m_g)*std::norm(vec(i));
        // b(i) = 1.0 - 2.0*this->m_lambda_x + (m_delta_t*0.5)*std::complex<double>(m_V(i)) + (m_delta_t*0.5 * m_g)*std::norm(vec(i));
    
        // TEST ADDED DDI
        a(i) = 1.0 + 2.0*this->m_lambda_x + U_potential + U_dd + U_lhy;
        b(i) = 1.0 - 2.0*this->m_lambda_x + U_potential + U_dd + U_lhy;

        // // //Test
        // a(i) = 1.0 - 2.0*this->m_lambda_x + 1.0 * (m_delta_t*0.5)*std::complex<double>(m_V(i)) + 1.0 * (m_delta_t*0.5 * m_g)*std::norm(vec(i));
        // b(i) = 1.0 + 2.0*this->m_lambda_x - 1.0 * (m_delta_t*0.5)*std::complex<double>(m_V(i)) - 1.0 * (m_delta_t*0.5 * m_g)*std::norm(vec(i));
        // //
        // a(i) = 1.0 - 2.0*this->m_lambda_x + A;
        // b(i) = 1.0 + 2.0*this->m_lambda_x - A;
    }

    this->init_Mat_A_1D(-1.0 *  m_lambda_x,a);
    this->init_Mat_B_1D(m_lambda_x,b); 
}

//Function for initialization left-hand side matrix according to Crank Nicolson algorithm
void CrankNicolson::init_Mat_A_1D(std::complex<double> r, Eigen::VectorXcd& d){
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
    this->m_A = A;
}

//Function for initialization right-hand side matrix according to Crank Nicolson algorithm
void CrankNicolson::init_Mat_B_1D(std::complex<double> r, Eigen::VectorXcd& d) {
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
    this->m_B = B;
}

Eigen::VectorXd CrankNicolson::create_harmonic_potential_1D(){
    Eigen::VectorXd V(this->m_size);
    V.setZero();
    double x = this->_start;
    for(int i = 0; i < this->m_size; ++i){
        V(i) = 0.5 * std::pow(x,2.);
        x += this->step;
    }
    return V;
}

Eigen::VectorXd CrankNicolson::create_potential_1D(){
    Eigen::VectorXd V(this->m_size);
    V.setZero();
    V(0) = V_0;
    V(this->m_size - 1) = V_0;
    return V;
}

// Time evolution simulation for 1D Gross-Pitaevskii equation
void CrankNicolson::simulation_1D(){
    int size = this->m_Psi.size();

    Eigen::VectorXcd x(size);
    Eigen::VectorXcd b(size);

    x.setZero();
    b.setZero();
    // Prepare the right-hand side for the time-stepping
    b = (this->m_Psi);
    
    // Set up the sparse LU solver
    // Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solver;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double>>> solver;
    // Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double>>> solver;

    for (int i = 1; i < this->t_step; ++i) {

        normalize_1D(b);
        update_time_evolution_matrices_1D(b);

        solver.compute(m_A);
        
        // Update the right-hand side vector b
        b = (this->m_B) * b;

        // // Solve the system A * x = b
        x = solver.solve(b);
        // Update b for the next iteration
        b = x;

        if(i % 100 == 0)
            std::cout << "Simulation step: #" << i << '\n';

        this->m_Fin = x;
        normalize_1D(this->m_Fin);

        this->vec_Energy.push_back(calc_state_energy(this->m_Fin));
    }
}


// Returning the probability density vector
Eigen::VectorXd CrankNicolson::prob_1D(Eigen::VectorXcd &vec)
{
    int size = vec.size();
    Eigen::VectorXd pr(size);
    for(int i = 0; i != size; ++i){
        pr(i) = std::norm(vec(i));
    }
    return pr;
}

void CrankNicolson::normalize_1D(Eigen::VectorXcd &vec){
    int size = vec.size();

    std::complex<double> psum = 0;
    for(int i = 0; i != size; ++i){
        psum += std::norm(vec(i));
    }
    std::complex<double> normalization_factor = std::sqrt(m_N) / std::sqrt(psum * this->step); // 

    vec *= std::abs(normalization_factor);
}

//Function for saving Eigen Matrices as csv tabels 
void CrankNicolson::save_vector_to_csv(std::string filename, Eigen::VectorXd v){

    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ",","\n");
    int size = v.size();
    Eigen::VectorXd vector(size);

    std::ofstream file(filename);
    if(file.is_open()){
        //file << vec.format(CSVFormat) << '\n';
        file << vector.format(CSVFormat);
        file.close();
    }
}

double CrankNicolson::vec_norm_1D(Eigen::VectorXcd &vec){
    int size = vec.size();
    double norm = 0;
    for(int i = 0; i < size; ++i){
        norm += std::norm(vec(i));
    }
    return norm * this->step;
}

double CrankNicolson::vec_norm_1D(Eigen::VectorXd &vec){
    int size = vec.size();
    double norm = 0;
    for(int i = 0; i < size; ++i){
        norm += std::pow(vec(i), 2);
    }
    return norm * this->step;
}

double CrankNicolson::calc_state_energy(){
    return calc_state_energy(this->m_Fin);
}

double CrankNicolson::calc_state_energy(Eigen::VectorXcd &vec){
    int size = vec.size();
    double energy = 0.0;
    for(int i = 1; i < size-1; ++i){
        std::complex<double> derivative = (vec(i + 1) - vec(i - 1)) / (2 * this->step);
        double kinetic = std::norm(derivative) * 0.5;
        double potential = this->m_V(i) * std::norm(vec(i));
        double interaction = 0.5 * m_g_scattering * std::norm(vec(i)) * std::norm(vec(i));
        double dipole ;
        energy += this->step * (kinetic + potential + interaction); 
    }
    return energy;
}


double CrankNicolson::calc_state_chem_potential(){
    return calc_state_chem_potential(this->m_Fin);
}

double CrankNicolson::calc_state_chem_potential(Eigen::VectorXcd &vec){
    int size = vec.size();
    double chem_potential = 0.0;
    for(int i = 1; i < size-1; ++i){
        std::complex<double> derivative = (vec(i + 1) - vec(i - 1)) / (2 * this->step);
        double kinetic = std::norm(derivative) * 0.5;
        double potential = this->m_V(i) * std::norm(vec(i));
        double interaction = 0.5 * std::norm(vec(i)) * std::norm(vec(i));
        chem_potential += this->step * (kinetic + potential + interaction); 
    }
    return chem_potential;
}


//********************************/***********/********************************//
//                                                                             //
//********************************/2DFunctions/********************************//
//                                                                             //
//********************************/***********/********************************//


//Constructor 2D
CrankNicolson::CrankNicolson(double h, double deltat, double T, double x_c, double sigma_x, double y_c, double sigma_y, double omega_x, double omega_y, double N, double a_s, double start){ 
    int M = 1/h;
    this->m_h_step = h;
    this->m_delta_t = deltat;
    this->m_T = T;
    this->V_0 = 1e+10;

    this->_start = start;
    double L = std::abs(_start*2);

    this->step = (-2*_start)/M;

    this->t_step = std::round(T/deltat) + 1;

    this->m_lambda_x = -1.0*m_delta_t/(4*std::pow(this->step,2));
    this->m_lambda_y = -1.0*m_delta_t/(4*std::pow(this->step,2));

    this->m_size = M;
    this->m_omega_x = omega_x;
    this->m_omega_y = omega_y;

    this->m_Fin = Eigen::VectorXcd((this->m_size - 2)* (this->m_size - 2)).setZero();

    //Physical parameters
    this->m_N = N;
    this->m_g_scattering = 0.3;
    this->m_g_lhy = 0.005;
    this->m_g_dipole = 0.005;

    this->init_chem_potential(omega_x, omega_y, N, a_s);

    //Choosing the potential for system 
    this->m_V = this->create_harmonic_potential_2D();
    // this->m_V = this->create_potential_box();

    this->init_start_state_2D(x_c,y_c,sigma_x,sigma_y);
    this->init_time_evolution_matrices_2D();

    F_ddi = std::make_unique<DipolarInteraction2D>(M-2,M-2,L,L,1,m_g_dipole);
}

//Indexes in vector representation of wave function
int CrankNicolson::get_m_index(int i, int j, int M){
    return (i-1)*(M-2) + j-1;
}

double CrankNicolson::thomas_fermi_radius_x(){   
    double numerator = 4.0 * this->m_g_scattering * this->m_omega_y * this->m_N;
    double denumerator = std::pow(this->m_omega_x, 3.) * M_PI;
    return std::pow(( numerator / denumerator), 0.25);
}

double CrankNicolson::thomas_fermi_radius_y(){
    double numerator = 4.0 * this->m_g_scattering * this->m_omega_x * this->m_N;
    double denumerator = std::pow(this->m_omega_y, 3.) * M_PI;
    return std::pow(( numerator / denumerator), 0.25);
}

// Function for calculating the chemical potential for Thomas-Fermi limit with associated parameters
void CrankNicolson::init_chem_potential(double omega_x, double omega_y, double N, double g){
    double numerator = g * omega_x * omega_y * N;
    double potential = std::sqrt(numerator / M_PI);
    this->m_chem_potential = potential;
}

//The value of the wave funciton in the specific point on grid
std::complex<double> CrankNicolson::gauss_wave_packet_2D(double x, double y, double x_c, double y_c, double sigma_x, double sigma_y){ //, double p_x, double p_y
    //std::complex<double> i(0, 1); // Define the imaginary unit
    double exponent = -1.0 * (std::pow(x - x_c ,2.) / (2 * std::pow(sigma_x,2.))) - (std::pow(y - y_c, 2.) / (2*std::pow(sigma_y, 2.)));

    //std::complex<double> phase = i * (p_x * (x - x_c) + p_y * (y - y_c));
    return std::exp(exponent); // + phase
}

double CrankNicolson::thomas_fermi_state_2D(double x, double y){
    double out, R_x, R_y;
    R_x = thomas_fermi_radius_x();
    R_y = thomas_fermi_radius_y();

    out = this->m_chem_potential * this->m_g_scattering * (1 - std::pow(x/R_x, 2.) - std::pow(y/R_y, 2.)); //  /,  

    if (out > 0)
        return std::sqrt(out);
    else 
        return 0;
}

Eigen::VectorXcd CrankNicolson::TM_state_2D(){
    int size = std::pow(this->m_size-2,2);
    Eigen::VectorXcd U(size);
    std::complex<double> psum = 0;

    for(int i = 1; i != m_size-1; ++i){
        double x = this->_start + i * this->step;
        for(int j = 1; j != m_size-1; ++j){
            double y = this->_start + j * this->step;
            //Initial state function
            std::complex<double> c = thomas_fermi_state_2D(x,y);

            int index = get_m_index(i,j,this->m_size);
            U(index) = c;
            psum += std::norm(c);
        }
    }
    //std::sqrt(this->m_N)
    std::complex<double> normalization_factor = std::sqrt(this->m_N) / std::sqrt(psum * this->step * this->step);
    U *= std::abs(normalization_factor);
    return U;
}

// Function for initializing wave function
void CrankNicolson::init_start_state_2D(double x_c, double y_c, double sigma_x, double sigma_y){
    int size = std::pow(this->m_size-2,2);
    Eigen::VectorXcd U(size);
    std::complex<double> psum = 0;

    for(int i = 1; i != m_size-1; ++i){
        double x = this->_start + i * this->step;
        for(int j = 1; j != m_size-1; ++j){
            double y = this->_start + j * this->step;
            //Initial state function

            std::complex<double> c = gauss_wave_packet_2D(x, y, x_c, y_c, sigma_x, sigma_y); // ,p_x, p_y
            // std::complex<double> c = thomas_fermi_state_2D(x,y);

            int index = get_m_index(i,j,this->m_size);
            U(index) = c;
            psum += std::norm(c);
        }
    }

    std::complex<double> normalization_factor = std::sqrt(this->m_N) / std::sqrt(psum * this->step * this->step);
    U *= std::abs(normalization_factor);
    this->m_Psi = U;
}


Eigen::MatrixXcd CrankNicolson::calculate_2D_DDI(Eigen::MatrixXcd &mat)
{
    int rows = mat.rows();
    int cols = mat.cols();
    Eigen::MatrixXd Phi_DDI(rows,cols);
    Phi_DDI.setZero();

    if(this->F_ddi){
        F_ddi->compute_DDI_term(mat, Phi_DDI);
    }

    return Phi_DDI;
}

//Function for constructing the right-hand and left-hand sides matrices from Crank Nicolson algorithm
void CrankNicolson::init_time_evolution_matrices_2D(){
    int mat_size = pow(this->m_size-2,2);
    Eigen::VectorXcd a(mat_size);
    Eigen::VectorXcd b(mat_size);

    Eigen::MatrixXcd mat = vec_to_mat(this->m_Psi);
    Eigen::MatrixXcd V_DDI = calculate_2D_DDI(mat);

    for(int k = 1; k < m_size-1; ++k){
        for(int l = 1; l < m_size-1; ++l){
            int index = get_m_index(k,l,m_size);

            std::complex<double> U_potential = 1.0*(m_delta_t * 0.5)*std::complex<double>(m_V(k,l));
            std::complex<double> U_contact = 1.0*(m_delta_t * 0.5)*std::norm(this->m_Psi(index));
            std::complex<double> U_dd = 1.0 * (m_delta_t*0.5 * m_g_dipole) * V_DDI(k-1,l-1);
            std::complex<double> U_lhy = 1.0 * (m_delta_t*0.5) * this->m_g_lhy * std::pow(std::norm(m_Psi(index)), 2.5);



            // //Real time evolution matrices
            // a(index) = (1.0 - 4.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));
            // b(index) = (1.0 + 4.0*this->m_lambda - 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));

            //Imaginary time evolution matrices
            a(index) = 1.0 - 2.0*this->m_lambda_x - 2.0*this->m_lambda_y + U_potential + U_contact + U_dd + U_lhy;
            b(index) = 1.0 + 2.0*this->m_lambda_x + 2.0*this->m_lambda_y - U_potential - U_contact - U_dd - U_lhy;
        }
    }
    this->init_Mat_A_2D(m_lambda_x, m_lambda_y,a);
    this->init_Mat_B_2D(-1.0 * m_lambda_x,-1.0 * m_lambda_y,b); 
}

void CrankNicolson::update_time_evolution_matrices_2D(Eigen::VectorXcd &vec){
    int mat_size = pow(this->m_size-2,2);
    Eigen::VectorXcd a(mat_size);
    Eigen::VectorXcd b(mat_size);

    Eigen::MatrixXcd mat = vec_to_mat(vec);

    Eigen::MatrixXcd V_DDI = calculate_2D_DDI(mat);

    for(int k = 1; k < m_size-1; ++k){
        for(int l = 1; l < m_size-1; ++l){
            int index = get_m_index(k,l,m_size);

            std::complex<double> U_potential = 1.0*(m_delta_t * 0.5)*std::complex<double>(m_V(k,l));
            std::complex<double> U_contact = 1.0*(m_delta_t * 0.5)*std::norm(vec(index));
            std::complex<double> U_dd = 1.0 * (m_delta_t*0.5 * m_g_dipole) * V_DDI(k-1,l-1);
            std::complex<double> U_lhy = 1.0 * (m_delta_t*0.5) * this->m_g_lhy * std::pow(std::norm(m_Psi(index)), 2.5);

            // //Real time evolution matrices
            // a(index) = (1.0 - 4.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));
            // b(index) = (1.0 + 4.0*this->m_lambda - 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));

            //Imaginary time evolution matrices
            a(index) = 1.0 - 2.0*this->m_lambda_x - 2.0*this->m_lambda_y + U_potential + U_contact + U_dd + U_lhy;
            b(index) = 1.0 + 2.0*this->m_lambda_x + 2.0*this->m_lambda_y - U_potential - U_contact - U_dd - U_lhy;
        }
    }
    this->init_Mat_A_2D(m_lambda_x, m_lambda_y,a);
    this->init_Mat_B_2D(-1.0 * m_lambda_x,-1.0 * m_lambda_y,b); 
}

//Function which return Eigen::MatrixXd based on Eigen::VectorXd
Eigen::MatrixXd CrankNicolson::vec_to_mat(const Eigen::VectorXd& vec) {   
    int size = std::sqrt(vec.size());
    Eigen::MatrixXd mat(size,size);
    Eigen::Map<const Eigen::MatrixXd> mat_map(vec.data(), size, size);
    mat = mat_map;
    
    return mat;
}


//Function which return Eigen::MatrixXcd based on Eigen::VectorXcd
Eigen::MatrixXcd CrankNicolson::vec_to_mat(const Eigen::VectorXcd& vec) {   
    int size = std::sqrt(vec.size());
    Eigen::MatrixXcd mat(size,size);
    Eigen::Map<const Eigen::MatrixXcd> mat_map(vec.data(), size, size);
    mat = mat_map;
    
    return mat;
}

//Function for calculate the density of probability in each point of space
Eigen::VectorXd CrankNicolson::prob(Eigen::VectorXcd &vec)
{
    int size = std::pow(this->m_size-2,2);
    Eigen::VectorXd pr(size);
    for(int i = 0; i < size; ++i){
        pr(i) = std::pow(vec(i).real(),2) + std::pow(vec(i).imag(),2);
    }
    return pr;
}

//sFunction for saving Eigen Matrices as csv tabels 
void CrankNicolson::save_matrix_to_csv(std::string filename, Eigen::MatrixXd mat){

    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ",","\n");
    int size = std::sqrt(mat.size());
    Eigen::MatrixXd vec(1,size);

    std::ofstream file(filename);
    if(file.is_open()){
        //file << vec.format(CSVFormat) << '\n';
        file << mat.format(CSVFormat);
        file.close();
    }
}

//Function for creating potential box without the wall
Eigen::MatrixXd CrankNicolson::create_potential_box(){
    int size = this->m_size;
    Eigen::MatrixXd V(size, size);
    V.setZero();
    Eigen::VectorXd S(size);
    S.fill(V_0);
    V.col(0) = S;
    V.col(size-1) = S;
    V.row(0) = S;
    V.row(size-1) = S;

    return V;
}

Eigen::MatrixXd CrankNicolson::create_harmonic_potential_2D(){
    int size = this->m_size;
    Eigen::MatrixXd V(size, size);
    V.setZero();
    for(int i = 1; i != size-1; ++i){
        double x = this->_start + i * this->step;
        for(int j = 1; j != size-1; ++j){
            double y = this->_start + j * this->step;
            V(i,j) = 0.5 * (std::pow(this->m_omega_x * x,2.) + std::pow(this->m_omega_y * y,2.));
        }
    }
    return V;
}

//Function for initialization left-hand side matrix according to Crank Nicolson algorithm
void CrankNicolson::init_Mat_A_2D(std::complex<double> r_x, std::complex<double> r_y, Eigen::VectorXcd& d){

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
    this->m_A = A;
}

//Function for initialization right-hand side matrix according to Crank Nicolson algorithm
void CrankNicolson::init_Mat_B_2D(std::complex<double> r_x, std::complex<double> r_y, Eigen::VectorXcd& d) {

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
    this->m_B = B;
}


// Time evolution simulation for 1D Gross-Pitaevskii equation
void CrankNicolson::simulation_2D(){

    int size = this->m_Psi.size();

    Eigen::VectorXcd x(size);
    Eigen::VectorXcd b(size);

    x.setZero();
    b.setZero();
    // Prepare the right-hand side for the time-stepping
    b = (this->m_Psi);
    
    // Set up the sparse LU solver
    // Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solver;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double>>> solver;
    // Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double>>> solver;

    for (int i = 1; i < this->t_step; ++i) {
        // Add the normalization function for 2D matrices, which would consider different spatial step in x and y directions
        normalize_2D(b);
        update_time_evolution_matrices_2D(b);

        solver.compute(m_A);
        
        // Update the right-hand side vector b
        b = (this->m_B) * b;

        // // Solve the system A * x = b
        x = solver.solve(b);
        // Update b for the next iteration
        b = x;

        if(i % 1 == 0)
            std::cout << "Simulation step: #" << i << '\n';

        this->m_Fin = x;
        normalize_2D(this->m_Fin);

        // this->vec_Energy.push_back(calc_state_energy(this->m_Fin));
    }
}

void CrankNicolson::normalize_2D(Eigen::VectorXcd &vec){
    int size = vec.size();

    double psum = 0;
    for(int i = 0; i < size; ++i){
        psum += std::norm(vec(i));
    }
    
    std::complex<double> normalization_factor = std::sqrt(this->m_N) / std::sqrt(psum * this->step * this->step); // 

    vec *= std::abs(normalization_factor);
}

double CrankNicolson::vec_norm_2D(Eigen::VectorXcd &vec){
    int size = vec.size();
    double norm = 0;
    for(int i = 0; i < size; ++i){
        norm += std::norm(vec(i));
    }
    return norm * this->step * this->step;
}

double CrankNicolson::vec_norm_2D(Eigen::VectorXd &vec){
    int size = vec.size();
    double norm = 0;
    for(int i = 0; i < size; ++i){
        norm += std::pow(vec(i), 2);
    }
    return norm * this->step * this->step;
}

//TODO ENERGY AND CHEMICAL POTENTIAL CALCULATION
double CrankNicolson::calc_state_energy_2D(){
    return calc_state_energy(this->m_Fin);
}

double CrankNicolson::calc_state_energy_2D(Eigen::VectorXcd &vec){
    int size = vec.size();
    double energy = 0.0;
    for(int i = 1; i < size-1; ++i){
        std::complex<double> derivative = (vec(i + 1) - vec(i - 1)) / (2 * this->step);
        double kinetic = std::norm(derivative) * 0.5;
        double potential = this->m_V(i) * std::norm(vec(i));
        double interaction = 0.5 * std::norm(vec(i)) * std::norm(vec(i));
        energy += this->step * (kinetic + potential + interaction); 
    }
    return energy;
}


// // Function for solving systems of equations for each time step dependently on start conditions

// Old version
// void CrankNicolson::simulation(){
//     int size = this->m_Psi.size();
//     Eigen::VectorXcd x(size);
//     Eigen::VectorXcd b(size);
//     x.setZero();
//     b.setZero();
//     // // Save initial data before the loop
//     // Eigen::VectorXd p_init = prob(m_Psi);
//     // save_matrix_to_csv("./Matrice/matrix0.csv", vec_to_mat(p_init));
//     // Prepare the right-hand side for the time-stepping
//     b = (this->m_Psi);
//     // Set up the sparse LU solver
//     Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> lg;
//     lg.compute(m_A);
//     for (int i = 1; i < this->t_step; ++i) {
//         // Update the right-hand side vector b
//         b = (this->m_B) * b;
//         // Solve the system A * x = b
//         x = lg.solve(b);
//         // Update b for the next iteration
//         b = x;
//         // // Calculate probability and save to CSV
//         // Eigen::VectorXd p = prob(x);
//         // Eigen::MatrixXd solution = vec_to_mat(p);
//         // save_matrix_to_csv("./Matrice/matrix" + std::to_string(i) + ".csv", solution);   
//         if(i % 4 == 0)
//             std::cout << "Simulation step: #" << i << '\n';
//     }
// }


//********************************/***********/********************************//
//                                                                             //
//*****************************/Getters Functions/*****************************//
//                                                                             //
//********************************/***********/********************************//


void CrankNicolson::m_Psi_len(){
    std::cout << "Psi length is: " << this->m_Psi.size() << std::endl;
}   
void CrankNicolson::m_Fin_len(){
    std::cout << "Fin length is: " << this->m_Fin.size() << std::endl;
}   

Eigen::VectorXcd CrankNicolson::get_m_Psi(){
    Eigen::VectorXcd output = this->m_Psi;
    return output;
}

Eigen::VectorXd CrankNicolson::get_m_Psi_prob(){
    int size = this->m_Psi.size();
    Eigen::VectorXd output(size);
    for(int i = 0; i < size; ++i){
        output(i) = std::pow(m_Psi(i).real(),2) + std::pow(m_Psi(i).imag(),2);
    }
    return output;
}

Eigen::VectorXcd CrankNicolson::get_m_Fin(){
    Eigen::VectorXcd output = this->m_Fin;
    return output;
}

Eigen::VectorXd CrankNicolson::get_m_Fin_prob(){
    int size = this->m_Fin.size();
    Eigen::VectorXd output(size);
    for(int i = 0; i < size; ++i){
        output(i) = std::pow(m_Fin(i).real(),2) + std::pow(m_Fin(i).imag(),2);
    }
    return output;
}

Eigen::VectorXd CrankNicolson::get_m_V(){
    Eigen::VectorXd output = this->m_V;
    return output;
}

void CrankNicolson::print_Mat_A_dim(){
    std::cout << "Matrix A number of cols: " << this->m_A.cols() << std::endl;
    std::cout << "Matrix A number of rows: " << this->m_A.rows() << std::endl;
}

void CrankNicolson::print_Mat_B_dim(){
    std::cout << "Matrix B number of cols: " << this->m_B.cols() << std::endl;
    std::cout << "Matrix B number of rows: " << this->m_B.rows() << std::endl;
}

void CrankNicolson::print_Mat_A(){
    std::cout << Eigen::MatrixXcd(this->m_A) << std::endl;
}

void CrankNicolson::print_Mat_B(){
    std::cout << Eigen::MatrixXcd(this->m_B) << std::endl;
}

Eigen::VectorXd CrankNicolson::real(Eigen::VectorXcd& vec){
    int size = vec.size();
    Eigen::VectorXd answ(size);
    for(int i = 0; i < size; ++i){
        answ(i) = vec(i).real();
    }
    return answ;
}

Eigen::VectorXd CrankNicolson::imag(Eigen::VectorXcd& vec){
    int size = vec.size();
    Eigen::VectorXd answ(size);
    for(int i = 0; i < size; ++i){
        answ(i) = vec(i).imag();
    }
    return answ;
}