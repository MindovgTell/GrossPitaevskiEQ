#include "CrankNicolson.hpp"

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
    this->m_lambda = -1.0*m_delta_t/(4*std::pow(this->step,2));
    this->m_size = M;
    this->m_omega_x = omega;
    this->m_omega_y = 0.0;

    this->m_Fin = Eigen::VectorXcd(m_size).setZero();
    //Physical parameters
    this->m_N = N;
    //this->m_g = 2 * M_PI * a_s;
    this->m_g = 1;

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

double CrankNicolson::thomas_fermi_state(double x){
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
        x += this->step;

        // std::complex<double> c = thomas_fermi_state(x - x_c); // Add parameters for Thomas-Fermi function
        // std::complex<double> c = gauss_wave_packet_1D(sigma_x, x, x_c, p_x);
        std::complex<double> c = square_func(x);

        U(i) = c;
        psum += std::norm(c);
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
        std::complex<double> c = thomas_fermi_state(x); // Add parameters for Thomas-Fermi function
        U(i) = c;
        psum += std::norm(c);
    }

    std::complex<double> normalization_factor = std::sqrt(m_N) / std::sqrt(psum * this->step);

    U = U * std::abs(normalization_factor);
    return U;
}

//Function for constructing the right-hand and left-hand sides matrices from Crank Nicolson algorithm
void CrankNicolson::init_time_evolution_matrices_1D(){
    Eigen::VectorXcd a(this->m_size);
    Eigen::VectorXcd b(this->m_size);

    for(int i = 0; i < this->m_size; ++i){
            // // Real time evolution matrices
            // a(i) = (1.0 - 2.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(i)) + 1.0i*(m_delta_t/2)*std::pow(std::abs(m_Psi(i)),2));
            // b(i) = (1.0 + 2.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(i)) + 1.0i*(m_delta_t/2)*std::pow(std::abs(m_Psi(i)),2));
            // Imaginary time evolution matrices
            a(i) = 1.0 + 2.0*this->m_lambda + 1.0 * (m_delta_t*0.5)*std::complex<double>(m_V(i)) + 1.0 * (m_delta_t*0.5 * m_g)*std::norm(m_Psi(i));
            b(i) = 1.0 - 2.0*this->m_lambda + 1.0 * (m_delta_t*0.5)*std::complex<double>(m_V(i)) + 1.0 * (m_delta_t*0.5 * m_g)*std::norm(m_Psi(i));

            // //Test
            // a(i) = 1.0 - 2.0*this->m_lambda + 1.0 * (m_delta_t*0.5)*std::complex<double>(m_V(i)) + 1.0 * (m_delta_t*0.5 * m_g)*std::norm(m_Psi(i));
            // b(i) = 1.0 + 2.0*this->m_lambda - 1.0 * (m_delta_t*0.5)*std::complex<double>(m_V(i)) - 1.0 * (m_delta_t*0.5 * m_g)*std::norm(m_Psi(i));
        }

    this->init_Mat_A_1D( -1.0 * m_lambda,a);
    this->init_Mat_B_1D(m_lambda,b); 
}

void CrankNicolson::update_time_evolution_matrices_1D(Eigen::VectorXcd &vec){

    int size = this->m_size; 
    Eigen::VectorXcd a(size);
    Eigen::VectorXcd b(size);
    for(int i = 0; i < size; ++i){
        // // Real time evolution matrices
        // a(i) = (1.0 - 2.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(i)) + 1.0i*(m_delta_t/2)*std::pow(std::abs(m_Psi(i)),2));
        // b(i) = (1.0 + 2.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(i)) + 1.0i*(m_delta_t/2)*std::pow(std::abs(m_Psi(i)),2));

        // Imaginary time evolution matrices
        a(i) = 1.0 + 2.0*this->m_lambda + (m_delta_t*0.5)*std::complex<double>(m_V(i)) + (m_delta_t*0.5 * m_g)*std::norm(vec(i));
        b(i) = 1.0 - 2.0*this->m_lambda + (m_delta_t*0.5)*std::complex<double>(m_V(i)) + (m_delta_t*0.5 * m_g)*std::norm(vec(i));
    
        // //Test
        // a(i) = 1.0 - 2.0*this->m_lambda + 1.0 * (m_delta_t*0.5)*std::complex<double>(m_V(i)) + 1.0 * (m_delta_t*0.5 * m_g)*std::norm(m_Psi(i));
        // b(i) = 1.0 + 2.0*this->m_lambda - 1.0 * (m_delta_t*0.5)*std::complex<double>(m_V(i)) - 1.0 * (m_delta_t*0.5 * m_g)*std::norm(m_Psi(i));

    }

    this->init_Mat_A_1D( -1.0 * m_lambda,a);
    this->init_Mat_B_1D(m_lambda,b); 
}

//Function for initialization left-hand side matrix according to Crank Nicolson algorithm
void CrankNicolson::init_Mat_A_1D(std::complex<double> r, Eigen::VectorXcd& d){
    int S = d.size();

    typedef Eigen::Triplet<std::complex<double> > T;
    std::vector<T> tripletList;
    tripletList.reserve(std::pow(S,2));

    tripletList.push_back(T(0, 0, d(0)));
    // tripletList.push_back(T(0, 1, r));
    // tripletList.push_back(T(S - 1, S - 2, r));
    //tripletList.push_back(T(S - 1, S - 1, d(S-1)));

    //std::cout << "the a vector is equal: " << '\n' << d << std::endl;

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
        x += this->step;
        V(i) = 0.5 * std::pow(x,2.);
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

        normalize(b);
        update_time_evolution_matrices_1D(b);

        solver.compute(m_A);
        
        // Update the right-hand side vector b
        b = (this->m_B) * b;

        // // Solve the system A * x = b
        x = solver.solve(b);
        // Update b for the next iteration
        b = x;
        // Calculate probability and save to CSV
        // Eigen::VectorXd solution = prob(x);
        // save_vector_to_csv("./Matrice/matrix" + std::to_string(i) + ".csv", solution);   
        if(i % 1000 == 0)
            std::cout << "Simulation step: #" << i << '\n';

        this->m_Fin = x;
        normalize(this->m_Fin);

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

void CrankNicolson::normalize(Eigen::VectorXcd &vec){
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

double CrankNicolson::vec_norm(Eigen::VectorXcd &vec){
    int size = vec.size();
    double norm = 0;
    for(int i = 0; i < size; ++i){
        norm += std::norm(vec(i));
    }
    return norm * this->step;
}

double CrankNicolson::vec_norm(Eigen::VectorXd &vec){
    int size = vec.size();
    double norm = 0;
    for(int i = 0; i < size; ++i){
        norm += pow(vec(i), 2);
    }
    return norm*this->step;
}

double CrankNicolson::calc_state_energy(){
    int size = this->m_Fin.size();
    double energy = 0.0;
    for(int i = 1; i < size-1; ++i){
        std::complex<double> derivative = (this->m_Fin(i + 1) - this->m_Fin(i - 1)) / (2 * this->step);
        double kinetic = std::norm(derivative) * 0.5;
        double potential = this->m_V(i) * std::norm(this->m_Fin(i));
        double interaction = 0.5 * std::norm(this->m_Psi(i)) * std::norm(this->m_Fin(i));
        energy += this->step * (kinetic + potential + interaction); 
    }

    return energy;
}

double CrankNicolson::calc_state_energy(Eigen::VectorXcd &vec){
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



//********************************/***********/********************************//
//                                                                             //
//********************************/2DFunctions/********************************//
//                                                                             //
//********************************/***********/********************************//


//Constructor 2D
CrankNicolson::CrankNicolson(double h, double deltat, double T, double x_c, double y_c, double sigma_x, double sigma_y, double omega_x, double omega_y, double N, double a_s, double start){ 
    int M = 1/h;
    this->m_h_step = h;
    this->m_delta_t = deltat;
    this->m_T = T;
    this->V_0 = 1e+10;
    this->_start = start;
    this->step = (-2*_start)/M;
    this->t_step = std::round(T/deltat) + 1;


    this->m_lambda = -1.0*m_delta_t/(4*std::pow(this->step,2));
    this->m_size = M;
    this->m_omega_x = omega_x;
    this->m_omega_y = omega_y;

    this->m_Fin = Eigen::VectorXcd(m_size).setZero();

    //Physical parameters
    this->m_N = N;
    this->m_g = 1;

    this->init_chem_potential(omega_x, omega_y, N, a_s);

    //Choosing the potential for system 
    // this->m_V = this->create_harmonic_potential_2D();
    // this->m_V = this->create_potential_box();

    // this->init_start_state_2D(x_c,y_c,sigma_x,sigma_y,0,0);
    // this->init_time_evolution_matrices_2D();
}

//Indexes in vector representation of wave function
int CrankNicolson::get_m_index(int i, int j, int M){
    return (i-1)*(M-2) + j-1;
}

double CrankNicolson::thomas_fermi_radius_x(){   
    double numerator = 4.0 * this->m_g * this->m_omega_y * this->m_N;
    double denumerator = std::pow(this->m_omega_x, 3.) * M_PI;
    return std::pow(( numerator / denumerator), 0.25);
}

double CrankNicolson::thomas_fermi_radius_y(){
    double numerator = 4.0 * this->m_g * this->m_omega_x * this->m_N;
    double denumerator = std::pow(this->m_omega_y, 3.) * M_PI;
    return std::pow(( numerator / denumerator), 0.25);
}

void CrankNicolson::init_chem_potential(double omega_x, double omega_y, double N, double a_s){
    
}

//The value of the wave funciton in the specific point on grid
std::complex<double> CrankNicolson::gauss_wave_packet_2D(double x, double y, double x_c, double y_c, double sigma_x, double sigma_y){ //, double p_x, double p_y
    std::complex<double> i(0, 1); // Define the imaginary unit
    double exponent = -(pow(x - x_c ,2) / (2 * pow(sigma_x,2))) - (pow(y - y_c, 2) / (2*pow(sigma_y, 2)));
    //std::complex<double> phase = i * (p_x * (x - x_c) + p_y * (y - y_c));
    return std::exp(exponent); // + phase
}



double CrankNicolson::thomas_fermi_state_2D(double x, double y){
    double out;

}


//Function for initializing wave function
void CrankNicolson::init_start_state_2D(double x_c, double y_c, double sigma_x, double sigma_y, double p_x=0, double p_y=0){
    int size = std::pow(this->m_size-2,2);
    Eigen::VectorXcd U(size);
    std::complex<double> psum = 0;

    for(int i = 1; i != m_size-1; ++i){
        double x = this->_start + i * this->step;
        for(int j = 1; j != m_size-1; ++j){
            double y = this->_start + this->step;
            //Initial state function
            std::complex<double> c = gauss_wave_packet_2D(x, y, x_c, y_c, sigma_x, sigma_y); // ,p_x, p_y
            // std::complex<double> c = thomas_fermi_state(x,y);

            int index = get_m_index(i,j,this->m_size);
            U(index) = c;
            psum += std::norm(c);
        }
    }

    std::complex<double> normalization_factor = std::sqrt(this->m_N) / std::sqrt(psum * std::pow(this->step, 2.));
    U *= std::abs(normalization_factor);
    this->m_Psi = U;
}

//Function for constructing the right-hand and left-hand sides matrices from Crank Nicolson algorithm
void CrankNicolson::init_time_evolution_matrices_2D(){
    int mat_size = pow(this->m_size-2,2);
    Eigen::VectorXcd a(mat_size);
    Eigen::VectorXcd b(mat_size);

    for(int k = 1; k < this->m_size-1; ++k){
        for(int l = 1; l < this->m_size-1; ++l){
            int index = get_m_index(k,l,m_size);
            a(index) = (1.0 - 4.0*this->m_lambda+ 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));
            b(index) = (1.0 + 4.0*this->m_lambda- 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));
        }
    }
    this->init_Mat_A_2D(m_lambda,a);
    this->init_Mat_B_2D(-m_lambda,b); 
}

//Function which return Eigen::MatrixXd based on Eigen::VectorXd
Eigen::MatrixXd CrankNicolson::vec_to_mat(const Eigen::VectorXd& vec) {   
    int size = std::sqrt(vec.size());
    Eigen::MatrixXd mat(size,size);
    Eigen::Map<const Eigen::MatrixXd> mat_map(vec.data(), size, size);
    mat = mat_map;
    
    return mat;
}

//Function for calculate the density of probability in each point of space
Eigen::VectorXd CrankNicolson::prob(Eigen::VectorXcd &vec)
{
    int size = std::pow(this->m_size-2,2);
    Eigen::VectorXd pr(size);
    for(int i = 0; i != size; ++i){
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

// Function for solving systems of equations for each time step dependently on start conditions
void CrankNicolson::simulation(){
    int size = this->m_Psi.size();

    Eigen::VectorXcd x(size);
    Eigen::VectorXcd b(size);

    x.setZero();
    b.setZero();

    // Save initial data before the loop
    Eigen::VectorXd p_init = prob(m_Psi);
    save_matrix_to_csv("./Matrice/matrix0.csv", vec_to_mat(p_init));

    // Prepare the right-hand side for the time-stepping
    b = (this->m_Psi);

    // Set up the sparse LU solver
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> lg;
    lg.compute(m_A);
    

    for (int i = 1; i < this->t_step; ++i) {
        // Update the right-hand side vector b
        b = (this->m_B) * b;

        // Solve the system A * x = b
        x = lg.solve(b);

        // Update b for the next iteration
        b = x;

        // Calculate probability and save to CSV
        Eigen::VectorXd p = prob(x);
        Eigen::MatrixXd solution = vec_to_mat(p);
        save_matrix_to_csv("./Matrice/matrix" + std::to_string(i) + ".csv", solution);   
        
        if(i % 4 == 0)
            std::cout << "Simulation step: #" << i << '\n';
    }
}

//Function for creating potential box without the wall
Eigen::MatrixXd CrankNicolson::create_potential_box(){
    Eigen::MatrixXd V(this->m_size,this->m_size);
    V.setZero();
    Eigen::VectorXd S(this->m_size);
    S.fill(V_0);
    V.col(0) = S;
    V.col(this->m_size-1) = S;
    V.row(0) = S;
    V.row(this->m_size-1) = S;

    return V;
}

//Function for initialization left-hand side matrix according to Crank Nicolson algorithm
void CrankNicolson::init_Mat_A_2D(std::complex<double> r, Eigen::VectorXcd& d){
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
void CrankNicolson::init_Mat_B_2D(std::complex<double> r, Eigen::VectorXcd& d) {
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







//********************************/***********/********************************//
//                                                                             //
//*****************************/Getters Functions/*****************************//
//                                                                             //
//********************************/***********/********************************//


void CrankNicolson::print_m_Psi(){
    std::cout << "Psi length is: " << this->m_Psi.size() << std::endl;
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