#ifndef SPLITSTEP_HPP
#define SPLITSTEP_HPP

#include <cmath>
#include <complex>
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
        
        SplitStep(
            ShrdPtrGrid grid,
            const WaveFuncType& vec,
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

            init_physics_params();
            init_kinetic_op();

            if (V_dd_ != 0.0) {
                int N = static_cast<int>(grid_->size());
                double L = static_cast<double>(N) * grid_->step();
                F_ddi_ = std::make_unique<DipolarInteraction<Dim>>(N, L, l_perp_, V_dd_);
                U_ddi_.resize(N);
            }
        }
        
        void step(WaveFuncType& vec) {
            int N = static_cast<int>(vec.size());

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
        }

    private:
        void init_physics_params() {
            double a_s = ph_config_.a_s;
            double a_dd = ph_config_.a_dd;

            double omega_t = grid_->get_omega_t();
            if (omega_t > 0.0) {
                l_perp_ = 1.0 / std::sqrt(omega_t);
            }

            if (a_dd != 0.0) {
                V_dd_ = -0.75 * a_dd / std::pow(l_perp_, 3.0);
            }

            g_scattering_ = (8.0 / 3.0) * V_dd_;

            if (a_s != 0.0) {
                double C = 1.4603;
                double denom = l_perp_ * l_perp_ * (1.0 - C * (a_s / l_perp_));
                if (denom != 0.0) {
                    g_scattering_ += 2.0 * a_s / denom;
                }

                g_lhy_ = (256.0 / (15.0 * M_PI)) * std::pow(a_s, 2.5) / std::pow(l_perp_, 3.0)
                    * (1.0 + 1.5 * std::pow(a_dd / a_s, 2.0));
            }
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


        void SplitStep<Dimension::One>::calc_potential_op(const Eigen::VectorXcd& vec){
            int N = static_cast<int>(vec.size());

            if (F_ddi_) {
                F_ddi_->compute_DDI_term(vec, U_ddi_);
            }

            U_v_.resize(N);

            for (int i = 0; i < N; ++i) {
                double U_trap = grid_->potential()(i);
                double U_scattering = g_scattering_ * std::norm(vec(i));
                double U_dd = (U_ddi_.size() == N) ? U_ddi_(i) : 0.0;
                double U_lhy = g_lhy_ * std::pow(std::norm(vec(i)), 1.5);

                double total = U_trap + U_scattering + U_dd + U_lhy;
                U_v_(i) = std::exp(std::complex<double>(0.0, -sim_config_.dt * total));
            }
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
#endif
