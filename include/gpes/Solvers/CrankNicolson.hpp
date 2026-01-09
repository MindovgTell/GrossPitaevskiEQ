#ifndef CRANKNICOLSON_HPP   
#define CRANKNICOLSON_HPP  


#include <vector>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>


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
    using WaveFuncType  = gpes::WaveFunction<Dim>;
    using GridType = gpes::Grid<Dim>;
    using Tag = gpes::tags::CrankNicolson;
    using ShrdPtrGrid = std::shared_ptr<const Grid<Dim>>;
private:
    ShrdPtrGrid grid_; // shared_ptr to the grid object

    SparseMat A_; // Left-hand side matrix from CN algorithm 
    SparseMat B_; // Right-hand side matrix from CN algorithm 


    PhysConfig ph_config_;
    SimConfig sim_config_;

    double _lambda_x;

    std::unique_ptr<DipolarInteraction<Dim>> F_ddi_;
    Eigen::VectorXd U_ddi_;
    // vector of the energies

    std::vector<double> vec_energy;
    std::vector<double> vec_chem_potential;
public:
    CrankNicolson(ShrdPtrGrid grid, const WaveFuncType& vec, PhysConfig phcnfg, SimConfig simcnfg) : grid_(std::move(grid)),
                                                                                        ph_cnfg_(phcnfg),
                                                                                        sim_cnfg_(simcnfg) 
    {
        calc_V_dd(); 
        calc_g_scattering();
        calc_g_lhy(); 

        // vec_Energy.reserve(sim_config_.num_of_steps);
        calc_time_evolution_matrices(vec);
    }

    // void calc_g_scattering() {
    //     double a_s = 
    //     double C = 1.4603; // riemann -zeta(1/2)
    //     double g_scattering =  2 * a_s / ((_l_perp*_l_perp) * (1 - C*(a_s/_l_perp))) + 8./ 3 * _V_dd;
    // }
    // void calc_V_dd(double a_dd) {
    //     // double cosTheta = 0; // Theta = 90 deg
    //     // _V_dd = 0.375 * a_dd / ( std::pow(_l_perp,3));
    //     // double cosTheta = 1; // Theta = 0 deg
    //     double _V_dd = -0.75 * a_dd / ( std::pow(_l_perp,3.)); 
    //     // _V_dd = 1.5 * a_dd / std::pow(_l_perp,3); 
    // }

    // void calc_g_lhy(double a_s, double a_dd) {
    //     double _g_lhy = (256. / (15 * M_PI) ) * std::pow(a_s, 2.5) / std::pow(_l_perp, 3.) * (1 + 1.5 * std::pow((a_dd / a_s ), 2.));
    // }

    //Function for calculating Dipole-Dipole Interaction
    void calculate_DDI(Eigen::VectorXcd& vec) {
        if(F_ddi)
            F_ddi->compute_DDI_term(vec, _U_ddi);
    }


    double V_1DD(double x) {
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

    std::vector<double> calculate_DDI_not_FFT(const Eigen::VectorXcd& vec) {
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

    void calc_time_evolution_matrices(const Eigen::VectorXcd& vec) {
        int size = vec.size();
        Eigen::VectorXcd a(size);
        Eigen::VectorXcd b(size);

        std::vector<double> U = calculate_DDI_not_FFT(vec);

        double _delta_t = sim_config_.dt;

        for(int i = 0; i < size; ++i){

            std::complex<double> U_potential = 1.0 * (_delta_t*0.5)*std::complex<double>(_V_ext(i));
            std::complex<double> U_scattering =  1.0 * (_delta_t*0.5 * _g_scattering)*std::norm(vec(i));
            std::complex<double> U_dd = 1.0 * (_delta_t*0.5) * U[i];
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

        init_Mat_A(_lambda_x, a);
        init_Mat_B(-1.0 * _lambda_x, b); 
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
        // Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double>>> solver;
        SolverType solver;

        normalize(b);
        calc_time_evolution_matrices(b);

        solver.compute(A_);

        b = B_ * b;

        // // Solve the system A * x = b
        x = solver.solve(b);
        // Update b for the next iteration
        b = x;

        normalize(x);
        double current_energy = calc_state_energy(x);
        vec_Energy.push_back(current_energy);

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
        return norm * _step;
    }

    double vec_norm(Eigen::VectorXd &vec) {
        int size = vec.size();
        double norm = 0;
        for(int i = 0; i < size; ++i){
            norm += std::pow(vec(i), 2);
        }
        return norm * _step;
    }

    // double calc_state_energy();
    double calc_state_energy(Eigen::VectorXcd &vec) {
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

    double calc_state_energy(WaveFuncType& vec) {
        calc_state_energy(vec.vec());
    }


    // int num_isnan() {
    //     int size = _Psi.size();
    //     int count = 0;
    //     for (int i = 0; i < size; ++i) {
    //         if (!std::isfinite(_Psi(i).real()) || !std::isfinite(_Psi(i).imag())) {
    //             ++count;
    //         }
    //     }
    //     return count;
    // }
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
    using WaveFuncType  = gpes::WaveFunction<Dim>;
    using GridType = gpes::Grid<Dim>;
    using Tag = gpes::tags::CrankNicolson;
    using ShrdPtrGrid = std::shared_ptr<const Grid<Dim>>;
private:
    SparseMat A_;
    SparseMat B_;

    ShrdPtrGrid grid_;

    PhysConfig& ph_config_;
    SimConfig& sim_config_;

    double lambda_x_, lambda_y_;

    // Eigen::MatrixXd _V_ext;
    // Eigen::VectorXd _U_ddi;

    std::unique_ptr<DipolarInteraction<Dim>> F_ddi_;

    std::vector<double> vec_energy_;

public: 
    CrankNicolson(ShrdPtrGrid grid, const WaveFuncType& vec, PhysConfig& phcnfg, SimConfig& simcnfg) : grid_(std::move(grid)),\
                                                                                               ph_config_(phcnfg),
                                                                                               sim_config_(simcnfg) 
    {
        // _U_ddi          =   Eigen::VectorXd::Zero(_size_x * _size_y);
        // _V_ext          =   _Psi._grid.potential();


        // l_z = 1.0 / std::sqrt(omega_z);
        // F_ddi           =   std::make_unique<DipolarInteraction<Dimension::Two>>(_size_x,_size_y,Lx,Ly,l_z,_V_dd);

        vec_Energy.reserve(sim_config_.num_of_steps);

        calc_time_evolution_matrices(vec.vec());
    }
    
    int get_index(int i,int j, ) {
        return i * grid_->size_x() + j;
    }

    // inline double calc_g_scattering(double a_s) {
    //     return std::sqrt(8. * M_PI) * a_s / l_z;
    // }
    // void calc_V_dd(double a_dd) {
    //     _V_dd = std::sqrt(8. * M_PI) * a_dd / l_z;
    // }
    // void calc_g_lhy(double a_s, double a_dd) {
    //     _g_lhy = (128. / ( 3 * std::pow(M_PI, 0.25) ) ) * std::sqrt(0.4) * std::pow(a_s, 2.5) / std::pow(l_z,1.5) * (1 + 1.5 * std::pow((a_dd / a_s ), 2.));
    // }

    //Function for calculating 2D DDI
    Eigen::VectorXd calc_FFT_DDI(const Eigen::VectorXcd& wave) {
        const int Nx = _Psi.grid_size_x();
        const int Ny = _Psi.grid_size_y();
        const int N = Nx*Ny;
        const double pi = M_PI;
        double dx = _Psi.step_x();
        double dy = _Psi.step_y();
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
    
    void calculate_DDI(const Eigen::VectorXcd &vec) {
        if(this->F_ddi)
            F_ddi->compute_DDI_term(vec, _U_ddi);
    }

    // void init_time_evolution_matrices();
    void calc_time_evolution_matrices(const Eigen::VectorXcd &vec) {

        double a_s = ph_config_.a_s;
        double g_scat = calc_g_scattering(a_s);
        double _delta_t = sim_config_.dt;

        int _size_x = grid_->size_x();
        int _size_y = grid_->size_y();

        int mat_size = _size_x * _size_y;
        Eigen::VectorXcd a(mat_size);
        Eigen::VectorXcd b(mat_size);
        
        a.setZero();
        b.setZero();

        // Eigen::MatrixXcd mat = vec_to_mat(vec);

        calculate_DDI(vec);
        // Eigen::VectorXd Phi_dd = calc_FFT_DDI(vec);

        double g_scat = calc_g_scattering(a_s);

        for(int k = 0; k < _size_x; ++k){
            for(int l = 0; l < _size_y; ++l){
                int index = get_index(k,l);

                std::complex<double> U_potential = (_delta_t * 0.5) * _V_ext(k,l); // std::complex<double>(
                std::complex<double> U_contact =  g_scat * (_delta_t * 0.5)*std::norm(vec(index));
                std::complex<double> U_dd =  (_delta_t*0.5) *_U_ddi(index);  //;Phi_dd(index);
                std::complex<double> U_lhy =  (_delta_t*0.5) * _g_lhy * std::pow(std::norm(vec(index)), 1.5);

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

        solver.compute(_A);
        
        // Update the right-hand side vector b
        b = B_ * b;

        // // Solve the system A * x = b
        x = solver.solve(b);
        // Update b for the next iteration
        b = x;

        normalize(x);
        double current_energy = calc_state_energy(x);
        vec_Energy.push_back(current_energy);
    }

    Eigen::VectorXd prob(Eigen::VectorXcd&  vec);

    void normalize(Eigen::VectorXcd &vec) {
         double _step_x = _Psi.step_x();
        double _step_y = _Psi.step_y();
        int _Num = _Psi.get_Num();

        int size = vec.size();
        double psum = 0;
        for(int i = 0; i < size; ++i){
            psum += std::norm(vec(i));
        }
        std::complex<double> normalization_factor = std::sqrt(_Num) / std::sqrt(psum * _step_x * _step_y); // 

        vec *= std::abs(normalization_factor);
    }

    void normalize(WaveFuncType& vec) {
        double _step_x = _Psi.step_x();
        double _step_y = _Psi.step_y();
        int _Num = _Psi.get_Num();

        int size = vec.size();
        double psum = 0;
        for(int i = 0; i < size; ++i){
            psum += std::norm(vec(i));
        }
        std::complex<double> normalization_factor = std::sqrt(_Num) / std::sqrt(psum * _step_x * _step_y); // 

        vec *= std::abs(normalization_factor);
    }

    double vec_norm(Eigen::VectorXcd &vec) {
        double _step_x = _Psi.step_x();
        double _step_y = _Psi.step_y();
        int size = vec.size();
        double norm = 0;
        for(int i = 0; i < size; ++i){
            norm += std::norm(vec(i));
        }
        return norm * _step_x * _step_y;
    }

    // double vec_norm(WaveFuncType &vec) {
    //     double _step_x = _Psi.step_x();
    //     double _step_y = _Psi.step_y();
    //     int size = vec.size();
    //     double norm = 0;
    //     for(int i = 0; i < size; ++i){
    //         norm += std::pow(vec(i), 2);
    //     }
    //     return norm * _step_x * _step_y;
    // }


    // double calc_state_energy();
    double calc_state_energy(Eigen::VectorXcd &vec) {
        double energy = 0.0;

        double a_s = _Psi.get_a_s();
        double g_scat = calc_g_scattering(a_s);
        int _size_x = _Psi.grid_size_x();
        int _size_y = _Psi.grid_size_y();

        double _step_x = _Psi.step_x();
        double _step_y = _Psi.step_y();

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
                double interaction = 0.5 * g_scat * std::norm(vec(index)) * std::norm(vec(index));

                //Dipole-Dipole interaction energy
                double ddi = 0.5 * std::norm(vec(index))* _U_ddi(index);// * Phi_dd(index); //

                //LHY correction term energy
                double lhy = _g_lhy * 0.4 * std::pow(std::norm(vec(index)), 2.5);

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