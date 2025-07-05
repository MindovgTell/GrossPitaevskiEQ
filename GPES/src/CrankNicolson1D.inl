#ifndef CRANKNICOLSON_INL
#define CRANKNICOLSON_INL

#include "CrankNicolson.hpp"


GPES::CrankNicolson<Dimension::One>::CrankNicolson(Grid<Dimension::One>& grid, WaveFunction<Dimension::One>& Psi, double deltat, double T): _delta_t(deltat), _T(T) {
    // Initialize parameters of the grid
    _size       =   grid.get_size_of_grid(); // number of grid nodes
    _step       =   grid.get_step_size();
    _start      =   grid.get_start_position();
    _V_ext      =   grid.get_potential();    
    _t_step     =   std::round(1/_delta_t) + 1;

    _Psi        =   Psi.get_wavefunction();
    _Fin        =   Eigen::VectorXcd::Zero(_size);
    _U_ddi      =   Eigen::VectorXd::Zero(_size);
    _lambda_x   =   std::complex<double>(-1.0*_delta_t/(4.0*std::pow(_step,2)),0);    
     
    //Physics units
    _a_dd       =   Psi.get_a_dd();
    _a_s        =   Psi.get_a_s();
    _Num        =   Psi.get_Num();
    _omega      =   grid.get_omega();
    _omega_t    =   grid.get_omega_t();
    _lam        =   _omega_t / _omega;


    _l_perp = 1.0 / std::sqrt(_omega_t);//1.0/std::sqrt(1);

    // Initialize interaction strengths

    calc_V_dd(_a_dd); 
    calc_g_scattering(_a_s);
    calc_g_lhy(_a_s, _a_dd); 


    // double L = std::abs(_start*2);
    // double confinment_length = 1;
    // F_ddi = std::make_unique<DipolarInteraction<Dimension::One>>(_size, L, confinment_length, _V_dd);

    vec_Energy.reserve(_t_step);
    calc_time_evolution_matrices(_Psi);
}

void GPES::CrankNicolson<Dimension::One>::calc_g_scattering(double a_s) {
    double C = 1.4603; // riemann -zeta(1/2)
    _g_scattering =  2 * a_s / ((_l_perp*_l_perp) * (1 - C*(a_s/_l_perp))) + 8./ 3 * _V_dd;
}

void GPES::CrankNicolson<Dimension::One>::calc_V_dd(double a_dd){
    // double cosTheta = 0; // Theta = 90 deg
    // _V_dd = 0.375 * a_dd / ( std::pow(_l_perp,3));
    // double cosTheta = 1; // Theta = 0 deg
    _V_dd = -0.75 * a_dd / ( std::pow(_l_perp,3.)); 
    // _V_dd = 1.5 * a_dd / std::pow(_l_perp,3); 
}
 
void GPES::CrankNicolson<Dimension::One>::calc_g_lhy(double a_s, double a_dd){
    _g_lhy = (256. / (15 * M_PI) ) * std::pow(a_s, 2.5) / std::pow(_l_perp, 3.) * (1 + 1.5 * std::pow((a_dd / a_s ), 2.));
}

void GPES::CrankNicolson<Dimension::One>::calculate_DDI(Eigen::VectorXcd& vec){
    if(this->F_ddi)
        F_ddi->compute_DDI_term(vec, _U_ddi);
}

double GPES::CrankNicolson<Dimension::One>::V_1DD(double x){
    double V;

    double xi = std::abs(x) / _l_perp;
    double xi2 = xi*xi;

    double alfa = std::abs(xi) / (std::sqrt(2.0));
    double sqalfa = 0.5 * xi2;

    if (xi2 / 2.0 < 700.0){
        double erfc_val = std::erfc(alfa);
        double exp_val = std::exp(sqalfa);
        V = _V_dd * (2.0 * std::abs(xi) - (std::sqrt(2.0 * M_PI) * (1.0 + xi2) * exp_val * erfc_val));
    } 
    else 
        V = _V_dd * 4 / std::pow(xi,3);

    if(std::isnan(V)) std::cout << x << '\t' << xi<< '\t'  << xi2 << '\t' << alfa << '\t' << sqalfa<< '\t' << std::endl;

    return V;
}

std::vector<double> GPES::CrankNicolson<Dimension::One>::calculate_DDI_not_FFT(const Eigen::VectorXcd &vec) {

    int size = vec.size(); //vec
    std::vector<double> U(size, 0.0);

    for (int i = 0; i < size; ++i) {
        double x = _start + _step * i; 
        double Sum = 0.0;
        for(int j = 0; j < size; ++j){
            double x_prime = _start + _step * j;
            
            double dist = std::abs(x-x_prime);

            double V_1d = V_1DD(dist);

            // if(std::isnan(std::norm(vec(j)))) std::cout << vec(i) << '\t';

            // std::cout << "Distance: " << dist << "\t Potential: "<< V_1d << std::endl;

            Sum += V_1d * std::norm(vec(j)) * _step;
        }
        // if(std::isnan(Sum)) std::cout << vec(i) << '\t';
        U[i] = Sum;
    }
    // std::cout << std::endl;

    return U;
}

void GPES::CrankNicolson<Dimension::One>::calc_time_evolution_matrices(const Eigen::VectorXcd& vec){
    int size = vec.size();
    Eigen::VectorXcd a(size);
    Eigen::VectorXcd b(size);

    // calculate_DDI(vec);
    std::vector<double> U = calculate_DDI_not_FFT(vec);
    // std::vector<double> U = compute_Phi_dd_realspace(vec);
    // std::vector<double> U = compute_Phi_dd(vec);
    for(int i = 0; i < size; ++i){

        std::complex<double> U_potential = 1.0 * (_delta_t*0.5)*std::complex<double>(_V_ext(i));
        std::complex<double> U_scattering =  1.0 * (_delta_t*0.5 * _g_scattering)*std::norm(vec(i));
        std::complex<double> U_dd = 1.0 * (_delta_t*0.5) * U[i];  //;
        std::complex<double> U_lhy = 1.0 * (_delta_t*0.5) * _g_lhy * std::pow(std::norm(vec(i)), 1.5);
        
        
        // // Real time evolution matrices
        // a(i) = (1.0 - 2.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(i)) + 1.0i*(m_delta_t/2)*std::pow(std::abs(m_Psi(i)),2));
        // b(i) = (1.0 + 2.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(i)) + 1.0i*(m_delta_t/2)*std::pow(std::abs(m_Psi(i)),2));
        // // Imaginary time evolution matrices
        // a(i) = 1.0 + 2.0*this->m_lambda_x + (m_delta_t*0.5)*std::complex<double>(m_V(i)) + (m_delta_t*0.5 * m_g)*std::norm(vec(i));
        // b(i) = 1.0 - 2.0*this->m_lambda_x + (m_delta_t*0.5)*std::complex<double>(m_V(i)) + (m_delta_t*0.5 * m_g)*std::norm(vec(i));
    
        // TEST ADDED DDI
        a(i) = 1.0 - 2.0*_lambda_x + U_scattering + U_potential + U_dd + U_lhy;
        b(i) = 1.0 + 2.0*_lambda_x - U_scattering - U_potential - U_dd - U_lhy;
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
    tripletList.reserve(S*3);

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
void GPES::CrankNicolson<Dimension::One>::simulation(std::string& outdir){

    Eigen::VectorXcd x(_size);
    Eigen::VectorXcd b(_size);

    x.setZero();
    b.setZero();
    // Prepare the right-hand side for the time-stepping
    b = _Psi;
    
    // Set up the sparse LU solver
    // Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solver;
    // Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double>>> solver;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double>>> solver;

    int i = 0;

    // do {
    //     normalize(b);
    //     calc_time_evolution_matrices(b);
    //     solver.compute(_A);
    //     // Update the right-hand side vector b
    //     b = _B * b;
    //     // // Solve the system A * x = b
    //     x = solver.solve(b);
    //     // Update b for the next iteration
    //     b = x;
    //     if(i % 100 == 0)
    //         std::cout << "Simulation step: #" << i << '\n';
    //     _Fin = x;
    //     normalize(_Fin);
    //     double current_energy = calc_state_energy(_Fin);
    //     vec_Energy.push_back(current_energy);
    //     ++i;
    // } while(simulation_stop(i));

    do {
        normalize(b);
        calc_time_evolution_matrices(b);

        // for(auto e: b)
        //     std::cout << e << std::endl;

        solver.compute(_A);
        
        // Update the right-hand side vector b
        b = _B * b;

        // // Solve the system A * x = b
        x = solver.solve(b);
        // Update b for the next iteration
        b = x;


        if(i % 1000 == 0){
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


    std::string outfile_en = outdir + "/energy.csv";
    savecsv_vec(outfile_en, vec_Energy);
}

inline bool GPES::CrankNicolson<Dimension::One>::simulation_stop(int i)
{
    if(i > 300000)
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
    return norm * _step;
}

double GPES::CrankNicolson<Dimension::One>::vec_norm(Eigen::VectorXd &vec){
    int size = vec.size();
    double norm = 0;
    for(int i = 0; i < size; ++i){
        norm += std::pow(vec(i), 2);
    }
    return norm * _step;
}

double GPES::CrankNicolson<Dimension::One>::calc_state_energy(){
    return calc_state_energy(_Fin);
}

double GPES::CrankNicolson<Dimension::One>::calc_state_energy(Eigen::VectorXcd &vec){
    int size = vec.size();
    double energy = 0.0;

    // std::vector<double> U = compute_Phi_dd_realspace(vec);
    // std::vector<double> U = compute_Phi_dd(vec);
    std::vector<double> U = calculate_DDI_not_FFT(vec);

    for(int i = 1; i < size-1; ++i){
        std::complex<double> derivative = (vec(i + 1) - vec(i - 1)) * (1. / (2 * _step));
        double kinetic = std::norm(derivative) * 0.5;
        double potential = _V_ext(i) * std::norm(vec(i));
        double interaction = 0.5 * _g_scattering * std::norm(vec(i)) * std::norm(vec(i));

        //Dipole-Dipole interaction energy
        // double ddi = 0.5 * _U_ddi(i) * std::norm(vec(i)); 
        double ddi = 0.5 * U[i] * std::norm(vec(i)); 

        //LHY correction term energy
        double lhy = _g_lhy * 0.4 * std::pow(std::norm(vec(i)), 2.5); //

        energy += _step * (kinetic + potential + interaction + ddi + lhy); // 
    }
    return energy;
}

double GPES::CrankNicolson<Dimension::One>::calc_state_energy(GPES::WaveFunction<Dimension::One>& vec){
    int size = vec.size();
    double energy = 0.0;

    // std::vector<double> U = compute_Phi_dd_realspace(vec.get_wavefunction());
    // std::vector<double> U = compute_Phi_dd(vec.get_wavefunction());
    std::vector<double> U = calculate_DDI_not_FFT(vec.get_wavefunction());

    for(int i = 1; i < size-1; ++i){
        std::complex<double> derivative = (vec(i + 1) - vec(i - 1)) * (1. / (2 * _step));
        double kinetic = std::norm(derivative) * 0.5;
        double potential = _V_ext(i) * std::norm(vec(i));
        double interaction = 0.5 * _g_scattering * std::norm(vec(i)) * std::norm(vec(i));

        //Dipole-Dipole interaction energy
        // double ddi = 0.5 * _U_ddi(i) * std::norm(vec(i)); 
        double ddi = 0.5 * U[i] * std::norm(vec(i)); 

        //LHY correction term energy
        double lhy = _g_lhy * 0.4 * std::pow(std::norm(vec(i)), 2.5);

        energy += _step * (kinetic + potential + interaction + ddi + lhy); //
    }
    return energy;
}

void GPES::CrankNicolson<Dimension::One>::get_final_state(GPES::WaveFunction<Dimension::One>& fin){ 
    fin.set_Num(_Num);
    fin.set_size_of_grid(_size);
    fin.set_start_position(_start);
    fin.set_step_size(_step);
    fin.set_a_s(_a_s);
    fin.set_a_dd(_a_dd);
    fin.set_omega(_omega);
    fin.set_vec(_Fin);
}

void GPES::CrankNicolson<Dimension::One>::print_param_of_eq(){
    int width = 15;
    std::cout << std::setw(width) << "a_s"<< std::setw(width) << "a_dd" << std::setw(width) << "g_scattering" << std::setw(width) << "V_dd" << std::setw(width) << "g_lhy" << std::setw(width) << "lambda" << std::endl;
    std::cout << std::setw(width) << _a_s << std::setw(width) << _a_dd << std::setw(width) << _g_scattering << std::setw(width) << _V_dd << std::setw(width) << _g_lhy << std::setw(width) << _lam << std::endl;
}

void GPES::CrankNicolson<Dimension::One>::save_state(std::string filename, Eigen::VectorXcd& vec){
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

void GPES::CrankNicolson<Dimension::One>::savecsv_wave(std::string file_path, Eigen::VectorXcd& v){
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
        {"omega",              _omega                      },
        {"Num_of_particle",    static_cast<double>(_Num)   },
        {"size_of_grid",       static_cast<double>(_size)  },
        {"step_size",          _step                       },
        {"grid_start_point",   _start                      }
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

void GPES::CrankNicolson<Dimension::One>::savecsv_vec(std::string file_path, std::vector<double>& v){
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

int GPES::CrankNicolson<Dimension::One>::num_isnan(){
    int size = _Psi.size();
    int count = 0;

    for (int i = 0; i < size; ++i) {
        if (!std::isfinite(_Psi(i).real()) || !std::isfinite(_Psi(i).imag())) {
            ++count;
        }
    }
    return count;
}

#endif