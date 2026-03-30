#pragma once

#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include <fftw3.h>
#include <Eigen/Dense>
#include <Eigen/SparseLU>

#include "Core/definitions.hpp"
#include "Core/traits.hpp"
#include "Log/Log.hpp"
#include "Solvers/experimental.hpp"
#include "dipoleinteraction/dipoleinteraction.hpp"
#include "grid/grid.hpp"
#include "wavefunction/wavefunction.hpp"

#ifndef M_PI
constexpr double M_PI = 3.14159265358979323846;
#endif

namespace gpes::solvers {

    template <Dimension Dim,
                //typename SolverType = Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>>
                typename SolverType =Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double>>>>
    class SplitStep;

    //********************************/***********/********************************//
    //                                                                             //
    //**************************/One dimensional solver/***************************//
    //                                                                             //
    //********************************/***********/********************************//

    template<typename SolverType>
    class SplitStep<Dimension::One, SolverType> {
    public:
        static constexpr Dimension Dim = Dimension::One;
        using WaveFuncType = WaveFunction<Dim>;
        using GridType = Grid<Dim>;
        using Tag = tags::SSFM;
        using ShrdPtrGrid = std::shared_ptr<const Grid<Dim>>;
        inline static const gpes::log::LogCategory kLogCategory{"SplitStep1D"};
    private:

        Eigen::VectorXcd U_k_;
        Eigen::VectorXcd U_v_;
        Eigen::VectorXd U_ddi_;

        // FFT setup (in-place on vec memory)
        fftw_plan plan_fwd_ = nullptr;
        fftw_plan plan_bwd_ = nullptr;
        fftw_complex* plan_data_ptr_ = nullptr;
        int plan_size_ = 0;
        const std::complex<double>* base_ptr_ = nullptr;
        int base_size_ = 0;

        ShrdPtrGrid grid_;
        PhysConfig ph_config_;
        SimConfig sim_config_;

        std::vector<double> vec_energy_;

    public:
        
        SplitStep(
            ShrdPtrGrid grid,
            WaveFuncType& vec,
            PhysConfig phcnfg,
            SimConfig simcnfg
        ) : grid_(std::move(grid)),
            ph_config_(phcnfg),
            sim_config_(simcnfg) 
        {
            if (!grid_) {
                GPES_LOG(kLogCategory, Error, "SplitStep requires a valid grid");
                throw std::runtime_error("SplitStep requires a valid grid");
            }
            if (vec.size() != static_cast<Eigen::Index>(grid_->size())) {
                GPES_LOG(kLogCategory, Error, "SplitStep wavefunction and grid size mismatch");
                throw std::runtime_error("SplitStep wavefunction and grid size mismatch");
            }

            try {
                detail::validate_phys_config(ph_config_);
            } catch (const std::exception& ex) {
                GPES_LOG(kLogCategory, Error, "{}", ex.what());
                throw;
            }

            if (!std::isfinite(grid_->step()) || grid_->step() == 0.0) {
                GPES_LOG(kLogCategory, Error, "SplitStep<Dimension::One>: step must be non-zero and finite");
                throw std::invalid_argument("SplitStep<Dimension::One>: step must be non-zero and finite");
            }
            if (!std::isfinite(sim_config_.dt) || sim_config_.dt == 0.0) {
                GPES_LOG(kLogCategory, Error, "SplitStep<Dimension::One>: dt must be non-zero and finite");
                throw std::invalid_argument("SplitStep<Dimension::One>: dt must be non-zero and finite");
            }
            if (!std::isfinite(grid_->omega_t()) || grid_->omega_t() <= 0.0) {
                GPES_LOG(kLogCategory, Error, "CrankNicolson<Dimension::One>: omega_t must be positive and finite");
                throw std::invalid_argument("CrankNicolson<Dimension::One>: omega_t must be positive and finite");
            }
            double l_perp = 1 / std::sqrt(grid_->omega_t());

            calc_inter_consts<Dim>(ph_config_,l_perp);
            try {
                detail::validate_inter_consts(ph_config_, "CrankNicolson<Dimension::One>");
            } catch (const std::exception& ex) {
                GPES_LOG(kLogCategory, Error, "{}", ex.what());
                throw;
            }

            vec_energy_.reserve(sim_config_.num_of_steps);
        
            base_ptr_ = vec.vec().data();
            base_size_ = static_cast<int>(vec.size());
            init_kinetic_op();
        }

        ~SplitStep() {
            destroy_fft();
        }

        void init_fft(fftw_complex* data_ptr, int N) {
            destroy_fft();
            plan_fwd_ = fftw_plan_dft_1d(N, data_ptr, data_ptr, FFTW_FORWARD, FFTW_ESTIMATE);
            plan_bwd_ = fftw_plan_dft_1d(N, data_ptr, data_ptr, FFTW_BACKWARD, FFTW_ESTIMATE);
            plan_data_ptr_ = data_ptr;
            plan_size_ = N;
        }

        void destroy_fft() {
            if (plan_fwd_) fftw_destroy_plan(plan_fwd_);
            if (plan_bwd_) fftw_destroy_plan(plan_bwd_);
            plan_fwd_ = nullptr;
            plan_bwd_ = nullptr;
            plan_data_ptr_ = nullptr;
            plan_size_ = 0;
        }
        
        void step(WaveFuncType& vec) {
            const int N = static_cast<int>(vec.size()); // Size of the grid
            auto* data = reinterpret_cast<fftw_complex*>(vec.vec().data());

            if (vec.vec().data() != base_ptr_ || N != base_size_) {
                GPES_LOG(kLogCategory, Error, "SplitStep<Dimension::One>: wavefunction resized or reallocated during solver lifetime");
                throw std::runtime_error("SplitStep<Dimension::One>: wavefunction resized or reallocated during solver lifetime");
            }
            if (!plan_fwd_ || !plan_bwd_ || plan_data_ptr_ != data || plan_size_ != N) {
                init_fft(data, N);
            }

            fftw_execute_dft(plan_fwd_, data, data);
            for (int i = 0; i < N; ++i) {
                vec.vec()(i) *= U_k_(i);
            }
            fftw_execute_dft(plan_bwd_, data, data);
            vec.vec() /= static_cast<double>(N);

            calc_potential_op(vec.vec());
            vec.vec().array() *= U_v_.array();

            fftw_execute_dft(plan_fwd_, data, data);
            for (int i = 0; i < N; ++i) {
                vec.vec()(i) *= U_k_(i);
            }
            fftw_execute_dft(plan_bwd_, data, data);
            vec.vec() /= static_cast<double>(N);

            if (sim_config_.ground_state) {
                normalize(vec.vec());
            }

            const int size = static_cast<int>(vec.size());
            std::vector<double> U = calculate_DDI_not_FFT(vec.vec());
            Eigen::VectorXd U_ddi(size);
            for (int i = 0; i < size; ++i) {
                U_ddi(i) = U[i];
            }
            const auto energy_terms = experimental::energy_1d(*grid_, vec.vec(), ph_config_, &U_ddi);
            double current_energy = energy_terms.total();
            vec_energy_.push_back(current_energy);
        }

        const std::vector<double>& energies() const {
            return vec_energy_;
        }

        double last_energy() const {
            if (vec_energy_.empty()) {
                return std::numeric_limits<double>::quiet_NaN();
            }
            return vec_energy_.back();
        }

    private:
        // Repeated functions for calculating DDI for CN and SSFM algorithm in future should be in DipolarInteraction class
        std::complex<double> derivative_x_at(const Eigen::VectorXcd& vec, int i) const {
            const int size = static_cast<int>(vec.size());
            const double dx = grid_->step();

            if (size <= 1) {
                return std::complex<double>(0.0, 0.0);
            }
            if (size == 2) {
                return (vec(1) - vec(0)) / dx;
            }

            if (i <= 0) {
                return (-3.0 * vec(0) + 4.0 * vec(1) - vec(2)) / (2.0 * dx);
            }
            if (i >= size - 1) {
                return (3.0 * vec(size - 1) - 4.0 * vec(size - 2) + vec(size - 3)) / (2.0 * dx);
            }

            return (vec(i + 1) - vec(i - 1)) / (2.0 * dx);
        }

        double integration_weight_1d(int i, int size) const {
            if (size <= 1) {
                return 1.0;
            }
            return (i == 0 || i == size - 1) ? 0.5 : 1.0;
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


            
        void init_kinetic_op() {
            int N = static_cast<int>(grid_->size());
            double dx = grid_->step();
            double L = static_cast<double>(N) * dx;
            const double dk = 2.0 * M_PI / L;

            U_k_.resize(N);

            for (int m = 0; m < N; ++m) {
                int q = (m <= N / 2) ? m : m - N;
                double k = dk * q;
                if (sim_config_.ground_state) {
                    U_k_(m) = std::exp(std::complex<double>(-0.25 * k * k * sim_config_.dt, 0.0));
                } else {
                    U_k_(m) = std::exp(std::complex<double>(0.0, -0.25 * k * k * sim_config_.dt));
                }
            }
        }


        void calc_potential_op(const Eigen::VectorXcd& vec){
            int size = vec.size();
            U_v_.resize(size);
            std::vector<double> U = calculate_DDI_not_FFT(vec);


            //TODO check the correctness of this definition
            for (int i = 0; i < size; ++i) {
                double U_trap = grid_->potential()(i);
                double U_scat = ph_config_.g_scat * std::norm(vec(i));
                double U_dd = U[i];
                double U_lhy = ph_config_.g_lhy * std::pow(std::norm(vec(i)), 1.5);

                double total = U_trap + U_scat + U_dd + U_lhy;
                if (sim_config_.ground_state) {
                    U_v_(i) = std::exp(std::complex<double>(-sim_config_.dt * total, 0.0));
                } else {
                    U_v_(i) = std::exp(std::complex<double>(0.0, -sim_config_.dt * total));
                }
            }
        }

        double calc_state_energy(const Eigen::VectorXcd& vec) {
            int size = vec.size();
            double energy = 0.0;

            std::vector<double> U = calculate_DDI_not_FFT(vec);

            for (int i = 0; i < size; ++i) {
                std::complex<double> derivative = derivative_x_at(vec, i);
                double kinetic = std::norm(derivative) * 0.5;
                double potential = grid_->potential()(i) * std::norm(vec(i));
                double interaction = 0.5 * ph_config_.g_scat * std::norm(vec(i)) * std::norm(vec(i));
                double ddi = 0.5 * U[i] * std::norm(vec(i));
                double lhy = ph_config_.g_lhy * 0.4 * std::pow(std::norm(vec(i)), 2.5);

                const double weight = integration_weight_1d(i, size);
                energy += grid_->step() * weight * (kinetic + potential + interaction + ddi + lhy);
            }
            return energy;
        }

        void normalize(Eigen::VectorXcd& vec) {
            int size = vec.size();
            double psum = 0.0;
            for (int i = 0; i < size; ++i) {
                psum += std::norm(vec(i));
            }
            std::complex<double> norm_factor =
                std::sqrt(ph_config_.num_of_prt) / std::sqrt(psum * grid_->step());
            vec *= std::abs(norm_factor);
        }
    };






    //********************************/***********/********************************//
    //                                                                             //
    //**************************/Two dimensional solver/***************************//
    //                                                                             //
    //********************************/***********/********************************//

    template<typename SolverType>
    class SplitStep<Dimension::Two, SolverType> {
    public:
        static constexpr Dimension Dim = Dimension::Two;
        using WaveFuncType = gpes::WaveFunction<Dim>;
        using GridType = gpes::Grid<Dim>;
        using Tag = gpes::tags::SSFM;
        using ShrdPtrGrid = std::shared_ptr<const Grid<Dim>>;
        inline static const gpes::log::LogCategory kLogCategory{"SplitStep2D"};
    private:
    
        Eigen::VectorXcd U_k_;
        Eigen::VectorXcd U_v_;
        Eigen::VectorXd U_ddi_;

        // FFT setup (in-place on vec memory)
        fftw_plan plan_fwd_ = nullptr;
        fftw_plan plan_bwd_ = nullptr;
        fftw_complex* plan_data_ptr_ = nullptr;
        int plan_size_ = 0;
        const std::complex<double>* base_ptr_ = nullptr;
        int base_size_ = 0;

        ShrdPtrGrid grid_;
        PhysConfig ph_config_;
        SimConfig sim_config_;

        std::unique_ptr<DipolarInteraction<Dim>> F_ddi_;
        std::vector<double> vec_energy_;

    public:

        SplitStep(
            ShrdPtrGrid grid,
            WaveFuncType& vec,
            PhysConfig phcnfg,
            SimConfig simcnfg
        ) : grid_(std::move(grid)),
            ph_config_(phcnfg),
            sim_config_(simcnfg)
        {
            if (!grid_) {
                GPES_LOG(kLogCategory, Error, "SplitStep requires a valid grid");
                throw std::runtime_error("SplitStep requires a valid grid");
            }
            if (vec.size() != static_cast<Eigen::Index>(grid_->size_x() * grid_->size_y())) {
                GPES_LOG(kLogCategory, Error, "SplitStep wavefunction and grid size mismatch");
                throw std::runtime_error("SplitStep wavefunction and grid size mismatch");
            }

            try {
                detail::validate_phys_config(ph_config_);
            } catch (const std::exception& ex) {
                GPES_LOG(kLogCategory, Error, "{}", ex.what());
                throw;
            }

            if (grid_->size_x() <= 0 || grid_->size_y() <= 0) {
                GPES_LOG(kLogCategory, Error, "SplitStep<Dimension::Two>: grid sizes must be positive");
                throw std::invalid_argument("SplitStep<Dimension::Two>: grid sizes must be positive");
            }
            if (!std::isfinite(grid_->step_x()) || grid_->step_x() == 0.0 ||
                !std::isfinite(grid_->step_y()) || grid_->step_y() == 0.0) {
                GPES_LOG(kLogCategory, Error, "SplitStep<Dimension::Two>: step_x and step_y must be non-zero and finite");
                throw std::invalid_argument("SplitStep<Dimension::Two>: step_x and step_y must be non-zero and finite");
            }
            if (!std::isfinite(sim_config_.dt) || sim_config_.dt == 0.0) {
                GPES_LOG(kLogCategory, Error, "SplitStep<Dimension::Two>: dt must be non-zero and finite");
                throw std::invalid_argument("SplitStep<Dimension::Two>: dt must be non-zero and finite");
            }
            if (!std::isfinite(grid_->omega_z()) || grid_->omega_z() <= 0.0) {
                GPES_LOG(kLogCategory, Error, "SplitStep<Dimension::Two>: omega_z must be positive and finite");
                throw std::invalid_argument("SplitStep<Dimension::Two>: omega_z must be positive and finite");
            }
            double l_z = 1.0 / std::sqrt(grid_->omega_z());

            calc_inter_consts<Dim>(ph_config_, l_z);
            try {
                detail::validate_inter_consts(ph_config_, "SplitStep<Dimension::Two>");
            } catch (const std::exception& ex) {
                GPES_LOG(kLogCategory, Error, "{}", ex.what());
                throw;
            }

            F_ddi_ = std::make_unique<DipolarInteraction<Dim>>(grid_, l_z, ph_config_.V_dd);

            vec_energy_.reserve(sim_config_.num_of_steps);

            base_ptr_ = vec.vec().data();
            base_size_ = static_cast<int>(vec.size());
            init_kinetic_op();
        }

        ~SplitStep() {
            destroy_fft();
        }

        void init_fft(fftw_complex* data_ptr, int size_x, int size_y) {
            destroy_fft();
            plan_fwd_ = fftw_plan_dft_2d(size_x, size_y, data_ptr, data_ptr, FFTW_FORWARD, FFTW_ESTIMATE);
            plan_bwd_ = fftw_plan_dft_2d(size_x, size_y, data_ptr, data_ptr, FFTW_BACKWARD, FFTW_ESTIMATE);
            plan_data_ptr_ = data_ptr;
            plan_size_ = size_x * size_y;
        }

        void destroy_fft() {
            if (plan_fwd_) fftw_destroy_plan(plan_fwd_);
            if (plan_bwd_) fftw_destroy_plan(plan_bwd_);
            plan_fwd_ = nullptr;
            plan_bwd_ = nullptr;
            plan_data_ptr_ = nullptr;
            plan_size_ = 0;
        }

        void step(WaveFuncType& vec) {
            const int size_x = static_cast<int>(grid_->size_x());
            const int size_y = static_cast<int>(grid_->size_y());
            const int size = size_x * size_y;

            auto* data = reinterpret_cast<fftw_complex*>(vec.vec().data());
            if (vec.vec().data() != base_ptr_ || size != base_size_) {
                GPES_LOG(kLogCategory, Error, "SplitStep<Dimension::Two>: wavefunction resized or reallocated during solver lifetime");
                throw std::runtime_error("SplitStep<Dimension::Two>: wavefunction resized or reallocated during solver lifetime");
            }
            if (!plan_fwd_ || !plan_bwd_ || plan_data_ptr_ != data || plan_size_ != size) {
                init_fft(data, size_x, size_y);
            }

            fftw_execute_dft(plan_fwd_, data, data);
            for (int i = 0; i < size; ++i) {
                vec.vec()(i) *= U_k_(i);
            }
            fftw_execute_dft(plan_bwd_, data, data);
            vec.vec() /= static_cast<double>(size);

            calc_potential_op(vec.vec());
            vec.vec().array() *= U_v_.array();

            fftw_execute_dft(plan_fwd_, data, data);
            for (int i = 0; i < size; ++i) {
                vec.vec()(i) *= U_k_(i);
            }
            fftw_execute_dft(plan_bwd_, data, data);
            vec.vec() /= static_cast<double>(size);

            if (sim_config_.ground_state) {
                normalize(vec.vec());
            }

            Eigen::VectorXd U_ddi(vec.size());
            if (F_ddi_) {
                F_ddi_->compute_DDI_term(vec.vec(), U_ddi);
            } else {
                U_ddi.setZero();
            }
            const auto energy_terms = experimental::energy_2d(*grid_, vec.vec(), ph_config_, &U_ddi);
            double current_energy = energy_terms.total();
            vec_energy_.push_back(current_energy);
        }

        const std::vector<double>& energies() const {
            return vec_energy_;
        }

        double last_energy() const {
            if (vec_energy_.empty()) {
                return std::numeric_limits<double>::quiet_NaN();
            }
            return vec_energy_.back();
        }

    private:
        int get_index(int i, int j) const {
            return i * static_cast<int>(grid_->size_y()) + j;
        }

        double integration_weight_1d(int idx, int size) const {
            if (size <= 1) {
                return 1.0;
            }
            return (idx == 0 || idx == size - 1) ? 0.5 : 1.0;
        }

        void init_kinetic_op() {
            const int size_x = static_cast<int>(grid_->size_x());
            const int size_y = static_cast<int>(grid_->size_y());
            const int size = size_x * size_y;
            const double dx = grid_->step_x();
            const double dy = grid_->step_y();
            const double Lx = static_cast<double>(size_x) * dx;
            const double Ly = static_cast<double>(size_y) * dy;
            const double dkx = 2.0 * M_PI / Lx;
            const double dky = 2.0 * M_PI / Ly;

            U_k_.resize(size);

            for (int i = 0; i < size_x; ++i) {
                int qx = (i <= size_x / 2) ? i : i - size_x;
                double kx = dkx * qx;
                for (int j = 0; j < size_y; ++j) {
                    int qy = (j <= size_y / 2) ? j : j - size_y;
                    double ky = dky * qy;
                    const int index = get_index(i, j);
                    double k2 = (kx * kx) + (ky * ky);
                    if (sim_config_.ground_state) {
                        U_k_(index) = std::exp(std::complex<double>(-0.25 * k2 * sim_config_.dt, 0.0));
                    }
                    else { 
                        U_k_(index) = std::exp(std::complex<double>(0.0, -0.25 * k2 * sim_config_.dt));
                    }
                }
            }
        }

        void calc_potential_op(const Eigen::VectorXcd& vec) {
            const int size_x = static_cast<int>(grid_->size_x());
            const int size_y = static_cast<int>(grid_->size_y());
            const int size = size_x * size_y;

            U_v_.resize(size);
            U_ddi_.resize(size);

            if (F_ddi_) {
                F_ddi_->compute_DDI_term(vec, U_ddi_);
            } else {
                U_ddi_.setZero();
            }

            for (int i = 0; i < size_x; ++i) {
                for (int j = 0; j < size_y; ++j) {
                    const int index = get_index(i, j);
                    double density = std::norm(vec(index));
                    double U_trap = grid_->potential()(i, j);
                    double U_scat = ph_config_.g_scat * density;
                    double U_dd = U_ddi_(index);
                    double U_lhy = ph_config_.g_lhy * std::pow(density, 1.5);

                    double total = U_trap + U_scat + U_dd + U_lhy;

                    if (sim_config_.ground_state) {
                        U_v_(index) = std::exp(std::complex<double>(-sim_config_.dt * total, 0.0));
                    }
                    else {
                        U_v_(index) = std::exp(std::complex<double>(0.0, -sim_config_.dt * total));
                    }
                }
            }
        }

        // Calculating energy of the state

        std::complex<double> derivative_x_at(const Eigen::VectorXcd& vec, int i, int j) const {
            const int size_x = static_cast<int>(grid_->size_x());
            const double dx = grid_->step_x();

            if (size_x <= 1) {
                return std::complex<double>(0.0, 0.0);
            }
            if (size_x == 2) {
                return (vec(get_index(1, j)) - vec(get_index(0, j))) / dx;
            }
            if (i <= 0) {
                return (-3.0 * vec(get_index(0, j)) + 4.0 * vec(get_index(1, j)) - vec(get_index(2, j))) / (2.0 * dx);
            }
            if (i >= size_x - 1) {
                return (3.0 * vec(get_index(size_x - 1, j)) - 4.0 * vec(get_index(size_x - 2, j)) + vec(get_index(size_x - 3, j))) / (2.0 * dx);
            }
            return (vec(get_index(i + 1, j)) - vec(get_index(i - 1, j))) / (2.0 * dx);
        }

        std::complex<double> derivative_y_at(const Eigen::VectorXcd& vec, int i, int j) const {
            const int size_y = static_cast<int>(grid_->size_y());
            const double dy = grid_->step_y();

            if (size_y <= 1) {
                return std::complex<double>(0.0, 0.0);
            }
            if (size_y == 2) {
                return (vec(get_index(i, 1)) - vec(get_index(i, 0))) / dy;
            }
            if (j <= 0) {
                return (-3.0 * vec(get_index(i, 0)) + 4.0 * vec(get_index(i, 1)) - vec(get_index(i, 2))) / (2.0 * dy);
            }
            if (j >= size_y - 1) {
                return (3.0 * vec(get_index(i, size_y - 1)) - 4.0 * vec(get_index(i, size_y - 2)) + vec(get_index(i, size_y - 3))) / (2.0 * dy);
            }
            return (vec(get_index(i, j + 1)) - vec(get_index(i, j - 1))) / (2.0 * dy);
        }

        double calc_state_energy(const Eigen::VectorXcd& vec) {
            const int size_x = static_cast<int>(grid_->size_x());
            const int size_y = static_cast<int>(grid_->size_y());
            const double dx = grid_->step_x();
            const double dy = grid_->step_y();
            double energy = 0.0;

            Eigen::VectorXd U_ddi(vec.size());
            if (F_ddi_) {
                F_ddi_->compute_DDI_term(vec, U_ddi);
            } else {
                U_ddi.setZero();
            }

            for (int i = 0; i < size_x; ++i) {
                const double wx = integration_weight_1d(i, size_x);
                for (int j = 0; j < size_y; ++j) {
                    const double wy = integration_weight_1d(j, size_y);
                    const int index = get_index(i, j);

                    std::complex<double> derivative_x = derivative_x_at(vec, i, j);
                    std::complex<double> derivative_y = derivative_y_at(vec, i, j);

                    double density = std::norm(vec(index));
                    double kinetic = (std::norm(derivative_x) + std::norm(derivative_y)) * 0.5;
                    double potential = grid_->potential()(i, j) * density;
                    double interaction = 0.5 * ph_config_.g_scat * density * density;
                    double ddi = 0.5 * U_ddi(index) * density;
                    double lhy = ph_config_.g_lhy * 0.4 * std::pow(density, 2.5);

                    energy += dx * dy * wx * wy * (kinetic + potential + interaction + ddi + lhy);
                }
            }
            return energy;
        }


        void normalize(Eigen::VectorXcd &vec) {
            int size = vec.size();
            double psum = 0;
            for(int i = 0; i < size; ++i){
                psum += std::norm(vec(i));
            }
            std::complex<double> normalization_factor = std::sqrt(ph_config_.num_of_prt) / std::sqrt(psum * grid_->step_x() * grid_->step_y()); // 
            vec *= std::abs(normalization_factor);
        }

    };

} // namespace gpes::solvers
