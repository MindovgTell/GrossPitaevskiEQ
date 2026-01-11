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
#include "dipoleinteraction/dipoleinteraction.hpp"
#include "grid/grid.hpp"
#include "wavefunction/wavefunction.hpp"

#ifndef M_PI
constexpr double M_PI = 3.14159265358979323846;
#endif

namespace gpes::solvers {

    template <Dimension Dim,
                typename SolverType = Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>>
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
    private:

        Eigen::VectorXcd U_k_;
        Eigen::VectorXcd U_v_;
        Eigen::VectorXd U_ddi_;

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

            vec_energy_.reserve(sim_config_.num_of_steps);

            init_kinetic_op();
        }
        
        void step(WaveFuncType& vec) {
            int N = static_cast<int>(vec.size());// Size of the grid

            Eigen::VectorXcd space_wf = vec.vec();
            Eigen::VectorXcd momentum_wf(N);

            fftw_plan plan_fwd = fftw_plan_dft_1d(
                N,
                reinterpret_cast<fftw_complex*>(space_wf.data()),
                reinterpret_cast<fftw_complex*>(momentum_wf.data()),
                FFTW_FORWARD,
                FFTW_ESTIMATE
            );
            fftw_plan plan_bwd = fftw_plan_dft_1d(
                N,
                reinterpret_cast<fftw_complex*>(momentum_wf.data()),
                reinterpret_cast<fftw_complex*>(space_wf.data()),
                FFTW_BACKWARD,
                FFTW_ESTIMATE
            );

            fftw_execute(plan_fwd);
            for (int i = 0; i < N; ++i) {
                momentum_wf(i) *= U_k_(i);
            }
            fftw_execute(plan_bwd);
            space_wf /= static_cast<double>(N);

            calc_potential_op(space_wf);
            for (int i = 0; i < N; ++i) {
                space_wf(i) *= U_v_(i);
            }

            fftw_execute(plan_fwd);
            for (int i = 0; i < N; ++i) {
                momentum_wf(i) *= U_k_(i);
            }
            fftw_execute(plan_bwd);
            space_wf /= static_cast<double>(N);

            vec.vec() = space_wf;

            fftw_destroy_plan(plan_fwd);
            fftw_destroy_plan(plan_bwd);

            double current_energy = calc_state_energy(vec.vec());
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
                U_k_(m) = std::exp(std::complex<double>(0.0, -0.25 * k * k * sim_config_.dt));
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
                double U_dd = 0.5 * U[i] * std::norm(vec(i));
                double U_lhy = ph_config_.g_lhy * std::pow(std::norm(vec(i)), 1.5);

                double total = U_trap + U_scat + U_dd + U_lhy;
                U_v_(i) = std::exp(std::complex<double>(0.0, -sim_config_.dt * total));
            }
        }

        double calc_state_energy(const Eigen::VectorXcd& vec) {
            int size = vec.size();
            double energy = 0.0;

            std::vector<double> U = calculate_DDI_not_FFT(vec);

            for (int i = 1; i < size - 1; ++i) {
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
    private:
    
        Eigen::VectorXcd U_k_;
        Eigen::VectorXcd U_v_;
        Eigen::VectorXd U_ddi_;

        ShrdPtrGrid grid_;
        PhysConfig ph_config_;
        SimConfig sim_config_;

        double g_scattering_ = 0.0;
        double g_lhy_ = 0.0;
        double l_perp_ = 1.0;
        double V_dd_ = 0.0;

        std::unique_ptr<DipolarInteraction<Dim>> F_ddi_;

    public:

        void step();

    private:
        void init_kinetic_op();
        void calc_potential_op(const Eigen::VectorXcd& vec);

    };

} // namespace gpes::solvers
