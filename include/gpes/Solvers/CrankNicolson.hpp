#ifndef CRANKNICOLSON_HPP   
#define CRANKNICOLSON_HPP  


#include <cmath>
#include <stdexcept>
#include <limits>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include "Core/utility.hpp"

#include "wavefunction/wavefunction.hpp"
#include "grid/grid.hpp"
#include "dipoleinteraction/dipoleinteraction.hpp"
#include "Core/traits.hpp"

namespace gpes::solvers{


    
template <Dimension D, 
          typename SolverType = Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>>
class CrankNicolson;

//********************************/***********/********************************//
//********************************/***********/********************************//
//***************************/One dimensional solver/**************************//
//********************************/***********/********************************//
//********************************/***********/********************************//


template <typename SolverType>
class CrankNicolson<Dimension::One, SolverType> {
public:
    static constexpr Dimension Dim = Dimension::One;
    using SparseMat = Eigen::SparseMatrix<std::complex<double>>;
    using WaveFuncType  = WaveFunction<Dim>;
    using GridType = Grid<Dim>;
    using Tag = tags::CrankNicolson;
    using ShrdPtrGrid = std::shared_ptr<const Grid<Dim>>;
private:
    ShrdPtrGrid grid_; // shared_ptr to the grid object

    SparseMat A_; // Left-hand side matrix from CN algorithm 
    SparseMat B_; // Right-hand side matrix from CN algorithm 


    PhysConfig ph_config_;
    SimConfig sim_config_;

    double lambda_x_;

    std::unique_ptr<DipolarInteraction<Dim>> F_ddi_;

    // vector of the energies
    std::vector<double> vec_energy_;
    std::vector<double> vec_chem_potential_;

public:
    CrankNicolson(
        ShrdPtrGrid grid,
        WaveFuncType& vec, 
        PhysConfig phcnfg, 
        SimConfig simcnfg
    ) : 
        grid_(std::move(grid)),
        ph_config_(phcnfg),
        sim_config_(simcnfg) 
    {
        if (!grid_) {
                throw std::runtime_error("SplitStep requires a valid grid");
            }
        if (vec.size() != static_cast<Eigen::Index>(grid_->size())) {
            throw std::runtime_error("SplitStep wavefunction and grid size mismatch");
        }

        detail::validate_phys_config(ph_config_);
        // TODO: I'm not sure about this form, need to varify 

        if (!std::isfinite(grid_->omega_t()) || grid_->omega_t() <= 0.0) {
            throw std::invalid_argument("CrankNicolson<Dimension::One>: omega_t must be positive and finite");
        }
        double l_perp = 1 / std::sqrt(grid_->omega_t());

        calc_inter_consts<Dim>(ph_config_,l_perp);
        detail::validate_inter_consts(ph_config_, "CrankNicolson<Dimension::One>");

        lambda_x_ = -1.0*sim_config_.dt/(4.0*std::pow(grid_->step(),2.));
        vec_energy_.reserve(sim_config_.num_of_steps);
        calc_time_evolution_matrices(vec.vec());
    }

    double last_energy() const {
        if (vec_energy_.empty()) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return vec_energy_.back();
    }

    double V_1DD(double x) {
        double V;

        double xi = std::abs(x) / ph_config_.l_; // TODO: check correctness of the code
        double xi2 = xi*xi;

        double alfa = std::abs(xi) / (std::sqrt(2.0));
        double sqalfa = 0.5 * xi2;

        if (xi2 / 2.0 < 700.0){
            double erfc_val = std::erfc(alfa);
            double exp_val = std::exp(sqalfa);
            V = ph_config_.V_dd * (2.0 * std::abs(xi) - (std::sqrt(2.0 * M_PI) * (1.0 + xi2) * exp_val * erfc_val));
        } 
        else 
            V = ph_config_.V_dd * 4 / std::pow(xi,3.);

        // TODO: in case of error throw exception
        if(std::isnan(V)) std::cout << x << '\t' << xi<< '\t'  << xi2 << '\t' << alfa << '\t' << sqalfa<< '\t' << std::endl;

        return V;
    }

    // Function for calculating Dipole-Dipole Interaction using FFT
    std::vector<double> calculate_DDI_not_FFT(const Eigen::VectorXcd& vec) {
        int size = vec.size(); //vec
        std::vector<double> U(size, 0.0);
        // TODO add checking correct numerical values
        for (int i = 0; i < size; ++i) {
            double x = grid_->start() + grid_->step() * i; 
            double Sum = 0.0;
            for(int j = 0; j < size; ++j){
                double x_prime = grid_->start() + grid_->step() * j;
                
                double dist = std::abs(x-x_prime);

                double V_1d = V_1DD(dist);

                Sum += V_1d * std::norm(vec(j)) * grid_->step();
            }
            
            U[i] = Sum;
        }

        return U;
    }

    void init_Mat_A(std::complex<double> r, Eigen::VectorXcd& d) {
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
        A_ = A;
    }

    void init_Mat_B(std::complex<double> r, Eigen::VectorXcd& d) {
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
        B_ = B;
    }

    void calc_time_evolution_matrices(Eigen::VectorXcd& vec) {
        int size = vec.size();
        Eigen::VectorXcd a(size);
        Eigen::VectorXcd b(size);

        // Possible to calculate via FFT, but not necessarily 
        std::vector<double> U = calculate_DDI_not_FFT(vec);

        double _delta_t = sim_config_.dt;

        for(int i = 0; i < size; ++i){

            std::complex<double> U_potential = (_delta_t*0.5)*std::complex<double>(grid_->potential()(i));
            std::complex<double> U_scattering =  (_delta_t*0.5 * ph_config_.g_scat)*std::norm(vec(i));
            std::complex<double> U_dd = (_delta_t*0.5) * U[i];
            std::complex<double> U_lhy = (_delta_t*0.5) * ph_config_.g_lhy * std::pow(std::norm(vec(i)), 1.5);
            
            // // Real time evolution matrices
            // a(i) = (1.0 - 2.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(i)) + 1.0i*(m_delta_t/2)*std::pow(std::abs(m_Psi(i)),2));
            // b(i) = (1.0 + 2.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(i)) + 1.0i*(m_delta_t/2)*std::pow(std::abs(m_Psi(i)),2));
            // // Imaginary time evolution matrices
            // a(i) = 1.0 + 2.0*this->m_lambda_x + (m_delta_t*0.5)*std::complex<double>(m_V(i)) + (m_delta_t*0.5 * m_g)*std::norm(vec(i));
            // b(i) = 1.0 - 2.0*this->m_lambda_x + (m_delta_t*0.5)*std::complex<double>(m_V(i)) + (m_delta_t*0.5 * m_g)*std::norm(vec(i));
        
            a(i) = 1.0 - 2.0*lambda_x_ + U_scattering + U_potential + U_dd + U_lhy;
            b(i) = 1.0 + 2.0*lambda_x_ - U_scattering - U_potential - U_dd - U_lhy;
        }

        init_Mat_A(lambda_x_, a);
        init_Mat_B(-1.0 * lambda_x_, b); 
    }


    void step(WaveFuncType& vec) {

        Eigen::VectorXcd x(vec.size());
        Eigen::VectorXcd b(vec.size());

        x.setZero();
        b.setZero();

        // Prepare the right-hand side for the time-stepping
        b = vec.vec();
        // assert(decltype(b) == decltype(std::declval<WF&>().vec()));
        
        // Set up the sparse LU solver
        // Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solver;
        // Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double>>> solver;
        Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double>>> solver;
        // SolverType solver;

        normalize(b);
        calc_time_evolution_matrices(b);

        solver.compute(A_);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("CrankNicolson<Dimension::One>: solver factorization failed");
        }

        b = B_ * b;

        // // Solve the system A * x = b
        x = solver.solve(b);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("CrankNicolson<Dimension::One>: solver solve failed");
        }

        normalize(x);
        vec.vec() = x;
        double current_energy = calc_state_energy(x);
        vec_energy_.push_back(current_energy);

    }

    void normalize(Eigen::VectorXcd &vec) {
        int size = vec.size();
        double psum = 0;
        for(int i = 0; i != size; ++i){
            psum += std::norm(vec(i));
        }
        std::complex<double> normalization_factor = std::sqrt(ph_config_.num_of_prt) / std::sqrt(psum * grid_->step()); // 

        vec *= std::abs(normalization_factor);
    }

    double vec_norm(Eigen::VectorXcd &vec) {
        int size = vec.size();
        double norm = 0;
        for(int i = 0; i < size; ++i){
            norm += std::norm(vec(i));
        }
        return norm * grid_->step();
    }

    double vec_norm(Eigen::VectorXd &vec) {
        int size = vec.size();
        double norm = 0;
        for(int i = 0; i < size; ++i){
            norm += std::pow(vec(i), 2.);
        }
        return norm * grid_->step();
    }

    // double calc_state_energy();
    double calc_state_energy(Eigen::VectorXcd &vec) {
        int size = vec.size();
        double energy = 0.0;
        std::vector<double> U = calculate_DDI_not_FFT(vec);

        for(int i = 1; i < size-1; ++i){
            std::complex<double> derivative = (vec(i + 1) - vec(i - 1)) * (1. / (2 * grid_->step()));
            double kinetic = std::norm(derivative) * 0.5;
            double potential = grid_->potential()(i) * std::norm(vec(i));
            double interaction = 0.5 * ph_config_.g_scat * std::norm(vec(i)) * std::norm(vec(i));
            double ddi = 0.5 * U[i] * std::norm(vec(i)); 
            double lhy = ph_config_.g_lhy * 0.4 * std::pow(std::norm(vec(i)), 2.5);

            energy += grid_->step() * (kinetic + potential + interaction + ddi + lhy); 
        }
        return energy;
    }

    double calc_state_energy(WaveFuncType& vec) {
        calc_state_energy(vec.vec());
    }

    const std::vector<double>& energies(){ return vec_energy_;}
};


//********************************/***********/********************************//
//********************************/***********/********************************//
//**************************/Two dimensional solver/***************************//
//********************************/***********/********************************//
//********************************/***********/********************************//


template <typename SolverType>
class CrankNicolson<Dimension::Two, SolverType> {
public:
    static constexpr Dimension Dim = Dimension::Two;
    using SparseMat = Eigen::SparseMatrix<std::complex<double>>;
    using WaveFuncType  = WaveFunction<Dim>;
    using GridType = Grid<Dim>;
    using Tag = tags::CrankNicolson;
    using ShrdPtrGrid = std::shared_ptr<const Grid<Dim>>;
private:
    SparseMat A_;
    SparseMat B_;

    ShrdPtrGrid grid_;

    PhysConfig& ph_config_;
    SimConfig& sim_config_;

    double lambda_x_, lambda_y_;

    std::unique_ptr<DipolarInteraction<Dim>> F_ddi_;

    std::vector<double> vec_energy_;

public: 
    CrankNicolson(ShrdPtrGrid grid, const WaveFuncType& vec, PhysConfig& phcnfg, SimConfig& simcnfg) : 
        grid_(std::move(grid)),
        ph_config_(phcnfg),
        sim_config_(simcnfg) 
    {
        detail::validate_phys_config(ph_config_);
        double l_z = 1.0 / std::sqrt(grid->omega_z());
        calc_inter_consts<Dim>(ph_config_, l_z);
        detail::validate_inter_consts(ph_config_, "CrankNicolson<Dimension::Two>");
        vec_energy_.reserve(sim_config_.num_of_steps);

        lambda_x_       =   -1.0*sim_config_.dt/(4*std::pow(grid_->step_x(),2.));
        lambda_y_       =   -1.0*sim_config_.dt/(4*std::pow(grid_->step_y(),2.));

        F_ddi_ =  std::make_unique<DipolarInteraction<Dim>>(grid_, l_z, ph_config_.V_dd);

        calc_time_evolution_matrices(vec.vec());
    }

    double last_energy() const {
        if (vec_energy_.empty()) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return vec_energy_.back();
    }
    
    int get_index(int i,int j) {
        return i * grid_->size_x() + j;
    }
    // Old version previosly used, now prefere DipolarInteraction class
    // //Function for calculating 2D DDI
    // Eigen::VectorXd calc_FFT_DDI(const Eigen::VectorXcd& wave) {
    //     const int Nx = _Psi.grid_size_x();
    //     const int Ny = _Psi.grid_size_y();
    //     const int N = Nx*Ny;
    //     const double pi = M_PI;
    //     double dx = _Psi.step_x();
    //     double dy = _Psi.step_y();
    //     Eigen::VectorXd phi_dd(N);
    //     assert(static_cast<int>(wave.size()) == N && "Density size must be Nx * Ny");
    //     const int Nkx = Nx;
    //     const int Nky = Ny/2 + 1;
    //     double* in = fftw_alloc_real(N);
    //     fftw_complex* out = fftw_alloc_complex(Nkx*Nky);
    //     for(int i = 0; i < N; ++i)
    //         in[i] = std::pow(wave(i).real(),2) + std::pow(wave(i).imag(),2); 
    //     fftw_plan forward  = fftw_plan_dft_r2c_2d(Nx, Ny, in, out, FFTW_ESTIMATE);
    //     fftw_plan backward = fftw_plan_dft_c2r_2d(Nx, Ny, out, in, FFTW_ESTIMATE);
    //     fftw_execute(forward);
    //     //Now compute kernel
    //     for(int i = 0; i < Nx; ++i){
    //         int ki = (i <= Nx/2)? i : i - Nx;
    //         double kx = 2.0 * pi * ki / (Nx*dx);
    //         for(int j = 0; j < Nky; ++j){
    //             int kj = j;
    //             double ky = 2.0 * pi * kj / (Ny * dy);
    //             double kappa = std::sqrt(kx*kx + ky*ky);
    //             int idx = i * Nky + j;
    //             double V_dd_k = 0.0;
    //             double k_r = kappa*l_z / std::sqrt(2);
    //             double K_exp = std::exp(k_r*k_r);
    //             double K_erfc = std::erfc(k_r);
    //             if(std::isnan(K_exp) || std::isnan(K_erfc)) std::cout<< "NAN value in erfc or exp" << std::endl;
    //             double F = 2.0 - 3.0 * std::sqrt(pi) * k_r * K_exp * K_erfc;
    //             V_dd_k = _V_dd * F;
    //             // Multiply by kernel
    //             out[idx][0] *= V_dd_k;
    //             out[idx][1] *= V_dd_k;
    //         }
    //     }
    //     // Inverse FFT
    //     fftw_execute(backward);
    //     // Normalize and copy output
    //     phi_dd.resize(N);
    //     for (int i = 0; i < N; ++i) {
    //         phi_dd[i] = in[i]/ N;
    //     }
    //     // Cleanup
    //     fftw_destroy_plan(forward);
    //     fftw_destroy_plan(backward);
    //     fftw_free(in);
    //     fftw_free(out);
    //     return phi_dd;
    // }
    
    void calculate_DDI(const Eigen::VectorXcd &vec, Eigen::VectorXd& U_ddi) {
        if(F_ddi_)
            F_ddi_->compute_DDI_term(vec, U_ddi);
    }

    // void init_time_evolution_matrices();
    void calc_time_evolution_matrices(const Eigen::VectorXcd &vec) {
        double _delta_t = sim_config_.dt;
        int _size_x = grid_->size_x();
        int _size_y = grid_->size_y();

        // TODO: check correct sizes, throw exception and log fatal error
        // if(vec.size() != _size_x * _size_y)
        //     throw();

        int mat_size = _size_x * _size_y;
        Eigen::VectorXcd a(mat_size);
        Eigen::VectorXcd b(mat_size);
        
        a.setZero();
        b.setZero();

        // Eigen::MatrixXcd mat = vec_to_mat(vec);

        Eigen::VectorXd U_ddi(mat_size);
        calculate_DDI(vec, U_ddi);
        // Eigen::VectorXd Phi_dd = calc_FFT_DDI(vec);

        for(int k = 0; k < _size_x; ++k){
            for(int l = 0; l < _size_y; ++l){
                int index = get_index(k,l);

                std::complex<double> U_potential = (_delta_t * 0.5) * grid_->potential()(k,l); 
                std::complex<double> U_contact =  ph_config_.g_scat * (_delta_t * 0.5)*std::norm(vec(index));
                std::complex<double> U_dd =  (_delta_t*0.5) * U_ddi(index);
                std::complex<double> U_lhy =  (_delta_t*0.5) * ph_config_.g_lhy * std::pow(std::norm(vec(index)), 1.5);

                // //Real time evolution matrices
                // a(index) = (1.0 - 4.0*this->m_lambda + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));
                // b(index) = (1.0 + 4.0*this->m_lambda - 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));

                //Imaginary time evolution matrices
                a(index) = 1.0 - 2.0*lambda_x_ - 2.0*lambda_y_ + U_contact + U_potential + U_dd + U_lhy; //
                b(index) = 1.0 + 2.0*lambda_x_ + 2.0*lambda_y_ - U_contact - U_potential - U_dd - U_lhy; // 
            }
        }
        init_Mat_A(lambda_x_,lambda_y_,a);
        init_Mat_B(-1.0 *lambda_x_,-1.0 * lambda_y_,b); 
    }

    void init_Mat_A(std::complex<double> r_x, std::complex<double> r_y, Eigen::VectorXcd& d) {
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
        A_ = A;
    }
    void init_Mat_B(std::complex<double> r_x, std::complex<double> r_y, Eigen::VectorXcd& d) {
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
        B_ = B;
    }

    // Eigen::MatrixXd vec_to_mat(const Eigen::VectorXd& vec);

    // TODO make a template method for checking is vec is correct WF with correct size
    void step(WaveFuncType& vec) {
        int size = vec.size();

        //TODO check the correct sizes
        // assert(size == grid_->size_x() * grid_->size_y());

        Eigen::VectorXcd x(size);
        Eigen::VectorXcd b(size);

        x.setZero();
        // Prepare the right-hand side for the time-stepping
        b = vec.vec();
        
        // Set up the sparse LU solver
        // Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solver;
        // Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double>>> solver;
        // Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double>>> solver;
        SolverType solver;

        normalize(b);
        calc_time_evolution_matrices(b);

        solver.compute(A_);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("CrankNicolson<Dimension::Two>: solver factorization failed");
        }
        
        // Update the right-hand side vector b
        b = B_ * b;

        // // Solve the system A * x = b
        x = solver.solve(b);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("CrankNicolson<Dimension::Two>: solver solve failed");
        }

        normalize(x);
        vec.vec() = x;
        double current_energy = calc_state_energy(x);
        vec_energy_.push_back(current_energy);
    }

    // Eigen::VectorXd prob(Eigen::VectorXcd&  vec);

    void normalize(Eigen::VectorXcd &vec) {
        int size = vec.size();
        double psum = 0;
        for(int i = 0; i < size; ++i){
            psum += std::norm(vec(i));
        }
        std::complex<double> normalization_factor = std::sqrt(ph_config_.num_of_prt) / std::sqrt(psum * grid_->step_x() * grid_->step_y()); // 

        vec *= std::abs(normalization_factor);
    }

    void normalize(WaveFuncType& vec) {
        int size = vec.size();
        double psum = 0;
        for(int i = 0; i < size; ++i){
            psum += std::norm(vec(i));
        }
        std::complex<double> normalization_factor = std::sqrt(ph_config_.num_of_prt) / std::sqrt(psum * grid_->step_x() * grid_->step_y()); // 
        vec *= std::abs(normalization_factor);
    }

    double vec_norm(Eigen::VectorXcd &vec) {
        int size = vec.size();
        double norm = 0;
        for(int i = 0; i < size; ++i){
            norm += std::norm(vec(i));
        }
        return norm * grid_->step_x() * grid_->step_y();
    }

    double calc_state_energy(Eigen::VectorXcd &vec) {
        double energy = 0.0;
        double _step_x = grid_->step_x();
        double _step_y = grid_->step_y();

        // Eigen::VectorXd Phi_dd = calc_FFT_DDI(vec);
        // TODO: optimize work with DDI calculations 
        Eigen::VectorXd U_ddi(vec.size());
        calculate_DDI(vec, U_ddi);

        for(int i = 1; i < grid_->size_x()-1; ++i){
            for(int j = 1; j < grid_->size_y()-1; ++j) {

                int index = get_index(i,j);
                int x_index_minus_one = get_index(i - 1, j);
                int x_index_plus_one = get_index(i + 1, j);
                int y_index_minus_one = get_index(i, j - 1);
                int y_index_plus_one = get_index(i, j + 1);

                std::complex<double> derivative_x = (vec(x_index_minus_one) - vec(x_index_plus_one)) * (1. / (2 * _step_x)); 
                std::complex<double> derivative_y = (vec(y_index_minus_one) - vec(y_index_plus_one)) * (1. / (2 * _step_y));

                double kinetic = (std::norm(derivative_x) + std::norm(derivative_y))* 0.5;
                double potential = grid_->potential()(i,j) * std::norm(vec(index));
                double interaction = 0.5 * ph_config_.g_scat * std::norm(vec(index)) * std::norm(vec(index));

                //Dipole-Dipole interaction energy
                double ddi = 0.5 * std::norm(vec(index)) * U_ddi(index); // Check how it defined 

                //LHY correction term energy
                double lhy = ph_config_.g_lhy * 0.4 * std::pow(std::norm(vec(index)), 2.5);

                energy += _step_x * _step_y * (kinetic + potential + interaction + ddi + lhy); // 
            }
        }
        return energy;
    }

    double calc_state_energy(WaveFuncType &vec) {
        return calc_state_energy(vec.vec());
    }

};

};


#endif 
