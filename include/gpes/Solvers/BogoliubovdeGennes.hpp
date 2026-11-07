#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>
#include <fftw3.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/Util/SelectionRule.h>

#include "Core/definitions.hpp"
#include "Core/traits.hpp"
#include "Core/utility.hpp"
#include "Log/Log.hpp"
#include "dipoleinteraction/dipoleinteraction.hpp"
#include "grid/grid.hpp"
#include "wavefunction/wavefunction.hpp"

namespace gpes::bdg {

    struct BdGMode {
        Eigen::VectorXcd u;
        Eigen::VectorXcd v;
        std::complex<double> omega;
        double relative_residual = std::numeric_limits<double>::quiet_NaN();
    };

    struct BdGConfig {
        int num_modes = 20;
        int ncv = 80;
        int max_iter = 2000;
        double tol = 1e-10;

        // Usually useful choices:
        // Spectra::SortRule::LargestMagn
        // Spectra::SortRule::LargestReal
        // Spectra::SortRule::SmallestMagn may catch Goldstone/low modes poorly
        // without shift-invert, depending on Spectra version/problem conditioning.
        Spectra::SortRule sort_rule = Spectra::SortRule::SmallestMagn;

        bool project_phase_mode = false;
        bool verbose = false;

        bool use_shift_invert = false;
        std::complex<double> shift = {0.05, 0.0};
        int shift_invert_linear_max_iter = 300;
        int shift_invert_linear_restart = 40;
        double shift_invert_linear_tol = 1e-8;
        int shift_invert_dense_max_dim = 2048;
        bool use_fourier_preconditioner = true;
        bool preconditioner_include_dipolar = true;
        double preconditioner_regularization = 1e-8;
    };

    struct BdGResult {
        std::vector<BdGMode> modes;
        int nconv = 0;
        int niter = 0;
        int nops = 0;
        Spectra::CompInfo info = Spectra::CompInfo::NotComputed;
        bool used_shift_invert = false;
        std::complex<double> shift = {0.0, 0.0};
        int linear_iterations = 0;
        int linear_failures = 0;
        bool used_fourier_preconditioner = false;
        bool used_dense_direct = false;
        double max_linear_residual = 0.0;
        double max_eigen_residual = 0.0;
    };

    template<Dimension dim>
    class BdGOperator;

    // Class for constructing an operator which will be later used for calculations of excitation spectra

    template<>
    class BdGOperator<Dimension::Two> {
    public:
        using Scalar = std::complex<double>;
        using ShrdPtrGrid = std::shared_ptr<const Grid<Dimension::Two>>;
        static constexpr Dimension Dim = Dimension::Two;
        inline static const gpes::log::LogCategory kLogCategory{"BdGOperator2D"};
    private:
        ShrdPtrGrid grid_;
        PhysConfig  config_;

        Eigen::VectorXcd psi0_;
        Eigen::VectorXd  density_;
        Eigen::VectorXd  phi_ddi_;
        Eigen::VectorXd  diag_A_;
        Eigen::VectorXcd diag_B_;
        Eigen::VectorXd  k2_;

        double mu_ = 0.0;
        double dx_ = 0.0;
        double dy_ = 0.0;
        int nx_ = 0;
        int ny_ = 0;
        int n_ = 0;

        std::unique_ptr<DipolarInteraction<Dimension::Two>> ddi_;

        mutable Eigen::VectorXcd fft_buffer_;
        mutable fftw_plan plan_fwd_ = nullptr;
        mutable fftw_plan plan_bwd_ = nullptr;

    public:

        BdGOperator(
            ShrdPtrGrid grid,
            const Eigen::VectorXcd& psi0,
            PhysConfig phys,
            double chemical_potential
        ) :
            grid_(std::move(grid)),
            config_(phys),
            psi0_(psi0),
            mu_(chemical_potential)
        {
            if (!grid_) {
                GPES_LOG(kLogCategory, Error, "BdGOperator2D requires a valid grid");
                throw std::runtime_error("BdGOperator2D requires a valid grid");
            }

            nx_ = static_cast<int>(grid_->size_x());
            ny_ = static_cast<int>(grid_->size_y());
            n_ = nx_ * ny_;

            dx_ = grid_->step_x();
            dy_ = grid_->step_y();

            if (nx_ <= 0 || ny_ <= 0) {
                throw std::invalid_argument("BdGOperator2D: grid sizes must be positive");
            }
            if (!std::isfinite(dx_) || dx_ == 0.0 ||
                !std::isfinite(dy_) || dy_ == 0.0) {
                throw std::invalid_argument("BdGOperator2D: grid steps must be finite and non-zero");
            }
            if (psi0_.size() != n_) {
                throw std::invalid_argument("BdGOperator2D: psi0 and grid size mismatch");
            }
            if (!std::isfinite(grid_->omega_z()) || grid_->omega_z() <= 0.0) {
                throw std::invalid_argument("BdGOperator2D: omega_z must be positive and finite");
            }
            if (!std::isfinite(mu_)) {
                throw std::invalid_argument("BdGOperator2D: chemical potential must be finite");
            }

            const double lz = 1.0 / std::sqrt(grid_->omega_z());

            calc_inter_consts<Dimension::Two>(config_, lz);
            gpes::solvers::detail::validate_inter_consts(config_, "BdGOperator2D");

            ddi_ = std::make_unique<DipolarInteraction<Dimension::Two>>(
                grid_,
                lz,
                config_.V_dd
            );

            init_spectral_kinetic();
            precompute_stationary_terms();
        }

        ~BdGOperator() {
            destroy_fft();
        }

        BdGOperator(const BdGOperator&) = delete;
        BdGOperator& operator=(const BdGOperator&) = delete;

        BdGOperator(BdGOperator&&) = delete;
        BdGOperator& operator=(BdGOperator&&) = delete;

        int rows() const {
            return 2 * n_;
        }

        int cols() const {
            return 2 * n_;
        }

        void perform_op(const Scalar* x_in, Scalar* y_out) const {
            Eigen::Map<const Eigen::VectorXcd> x(x_in, 2 * n_);
            Eigen::Map<Eigen::VectorXcd> y(y_out, 2 * n_);

            const Eigen::VectorXcd u = x.segment(0, n_);
            const Eigen::VectorXcd v = x.segment(n_, n_);

            Eigen::VectorXcd Au(n_);
            Eigen::VectorXcd Av(n_);
            Eigen::VectorXcd ddi_source(n_);
            Eigen::VectorXcd ddi_response(n_);

            apply_A(u, Au);
            apply_A_conjugate(v, Av);

            ddi_source = psi0_.conjugate().array() * u.array()
                    + psi0_.array() * v.array();

            ddi_->compute_DDI_from_complex_source(ddi_source, ddi_response);

            Eigen::VectorXcd yu(n_);
            Eigen::VectorXcd yv(n_);

            yu = Au;
            yv = -Av;

            for (int i = 0; i < n_; ++i) {
                yu(i) += diag_B_(i) * v(i) + psi0_(i) * ddi_response(i);
                yv(i) -= std::conj(diag_B_(i)) * u(i) + std::conj(psi0_(i)) * ddi_response(i);
            }

            y.segment(0, n_) = yu;
            y.segment(n_, n_) = yv;
        }

        int size_x() const { return nx_; }
        int size_y() const { return ny_; }
        double step_x() const { return dx_; }
        double step_y() const { return dy_; }
        const Eigen::VectorXd& kinetic_k2() const { return k2_; }
        double dipolar_strength() const { return config_.V_dd; }
        double confinement_length() const {
            return 1.0 / std::sqrt(grid_->omega_z());
        }

        // A density-weighted homogeneous approximation is used by the
        // Fourier preconditioner. It captures the bulk interaction scale while
        // retaining the exact spectral kinetic-energy dependence.
        double preconditioner_density() const {
            double weight = 0.0;
            double weighted_density = 0.0;
            for (int i = 0; i < n_; ++i) {
                weight += density_(i);
                weighted_density += density_(i) * density_(i);
            }
            return weight > 0.0 ? weighted_density / weight : 0.0;
        }

        double preconditioner_diag_a() const {
            double weight = 0.0;
            double weighted_value = 0.0;
            for (int i = 0; i < n_; ++i) {
                weight += density_(i);
                weighted_value += density_(i) * diag_A_(i);
            }
            return weight > 0.0 ? weighted_value / weight : 0.0;
        }

        Scalar preconditioner_diag_b() const {
            double weight = 0.0;
            Scalar weighted_value = Scalar(0.0, 0.0);
            for (int i = 0; i < n_; ++i) {
                weight += density_(i);
                weighted_value += density_(i) * diag_B_(i);
            }
            return weight > 0.0 ? weighted_value / weight : Scalar(0.0, 0.0);
        }

    private:
        int idx(int i, int j) const {
            return i * ny_ + j;
        }

        void precompute_stationary_terms() {
            density_.resize(n_);
            phi_ddi_.resize(n_);
            diag_A_.resize(n_);
            diag_B_.resize(n_);

            for (int i = 0; i < n_; ++i) {
                density_(i) = std::norm(psi0_(i));
            }

            ddi_->compute_DDI_term(psi0_, phi_ddi_);

            for (int i = 0; i < nx_; ++i) {
                for (int j = 0; j < ny_; ++j) {
                    const int k = idx(i, j);
                    const double n = density_(k);
                    const double sqrt_n = std::sqrt(std::max(n, 0.0));

                    const double V = grid_->potential()(i, j);

                    const double lhy_A = 2.5 * config_.g_lhy * std::pow(n, 1.5);
                    const std::complex<double> lhy_B =
                        1.5 * config_.g_lhy * sqrt_n * psi0_(k) * psi0_(k);

                    diag_A_(k) =
                        V
                        - mu_
                        + 2.0 * config_.g_scat * n
                        + phi_ddi_(k)
                        + lhy_A;

                    diag_B_(k) =
                        config_.g_scat * psi0_(k) * psi0_(k)
                        + lhy_B;
                }
            }
        }

        void init_spectral_kinetic() {
            k2_.resize(n_);

            const double Lx = static_cast<double>(nx_) * dx_;
            const double Ly = static_cast<double>(ny_) * dy_;
            const double dkx = 2.0 * M_PI / Lx;
            const double dky = 2.0 * M_PI / Ly;

            for (int i = 0; i < nx_; ++i) {
                const int qx = (i <= nx_ / 2) ? i : i - nx_;
                const double kx = dkx * static_cast<double>(qx);
                for (int j = 0; j < ny_; ++j) {
                    const int qy = (j <= ny_ / 2) ? j : j - ny_;
                    const double ky = dky * static_cast<double>(qy);
                    k2_(idx(i, j)) = kx * kx + ky * ky;
                }
            }

            fft_buffer_ = Eigen::VectorXcd::Zero(n_);
            auto* data = reinterpret_cast<fftw_complex*>(fft_buffer_.data());

            plan_fwd_ = fftw_plan_dft_2d(
                nx_,
                ny_,
                data,
                data,
                FFTW_FORWARD,
                FFTW_ESTIMATE
            );

            plan_bwd_ = fftw_plan_dft_2d(
                nx_,
                ny_,
                data,
                data,
                FFTW_BACKWARD,
                FFTW_ESTIMATE
            );

            if (!plan_fwd_ || !plan_bwd_) {
                destroy_fft();
                throw std::runtime_error("BdGOperator2D: FFTW plan creation failed");
            }
        }

        void destroy_fft() {
            if (plan_fwd_) {
                fftw_destroy_plan(plan_fwd_);
                plan_fwd_ = nullptr;
            }
            if (plan_bwd_) {
                fftw_destroy_plan(plan_bwd_);
                plan_bwd_ = nullptr;
            }
        }

        void apply_spectral_kinetic(const Eigen::VectorXcd& in, Eigen::VectorXcd& out) const {
            if (in.size() != n_) {
                throw std::invalid_argument("BdGOperator2D::apply_spectral_kinetic: input size mismatch");
            }
            if (!plan_fwd_ || !plan_bwd_) {
                throw std::runtime_error("BdGOperator2D::apply_spectral_kinetic: FFTW plans are not initialized");
            }

            fft_buffer_ = in;
            auto* data = reinterpret_cast<fftw_complex*>(fft_buffer_.data());

            fftw_execute_dft(plan_fwd_, data, data);
            fft_buffer_.array() *= (0.5 * k2_.array()).cast<std::complex<double>>();
            fftw_execute_dft(plan_bwd_, data, data);
            fft_buffer_ /= static_cast<double>(n_);

            out = fft_buffer_;
        }

        void apply_A(const Eigen::VectorXcd& in, Eigen::VectorXcd& out) const {
            apply_spectral_kinetic(in, out);

            out.array() += diag_A_.array().cast<std::complex<double>>() * in.array();
        }

        void apply_A_conjugate(const Eigen::VectorXcd& in, Eigen::VectorXcd& out) const {
            // For real trap/density/DDI/LHY diag_A_ is real.
            // If you later allow rotating frames, gauge fields, or complex stationary
            // background operators, replace this with the actual conjugate block.
            apply_A(in, out);
        }


    };

    using BdGOperator2D = BdGOperator<Dimension::Two>;

    class BdGFourierPreconditioner2D {
    public:
        using Scalar = std::complex<double>;

    private:
        const BdGOperator2D& op_;
        Scalar shift_;
        bool include_dipolar_ = true;
        double regularization_ = 1e-8;
        int nx_ = 0;
        int ny_ = 0;
        int n_ = 0;
        double dx_ = 0.0;
        double dy_ = 0.0;

        Eigen::VectorXcd inv11_;
        Eigen::VectorXcd inv12_;
        Eigen::VectorXcd inv21_;
        Eigen::VectorXcd inv22_;

        mutable Eigen::VectorXcd work_u_;
        mutable Eigen::VectorXcd work_v_;
        mutable fftw_plan plan_fwd_ = nullptr;
        mutable fftw_plan plan_bwd_ = nullptr;

    public:
        BdGFourierPreconditioner2D(
            const BdGOperator2D& op,
            Scalar shift,
            bool include_dipolar,
            double regularization
        ) :
            op_(op),
            shift_(shift),
            include_dipolar_(include_dipolar),
            regularization_(regularization),
            nx_(op.size_x()),
            ny_(op.size_y()),
            n_(nx_ * ny_),
            dx_(op.step_x()),
            dy_(op.step_y())
        {
            if (regularization_ < 0.0 || !std::isfinite(regularization_)) {
                throw std::invalid_argument(
                    "BdGFourierPreconditioner2D: invalid regularization"
                );
            }

            build_inverse_blocks();

            work_u_ = Eigen::VectorXcd::Zero(n_);
            work_v_ = Eigen::VectorXcd::Zero(n_);
            auto* data_u = reinterpret_cast<fftw_complex*>(work_u_.data());

            plan_fwd_ = fftw_plan_dft_2d(
                nx_, ny_, data_u, data_u, FFTW_FORWARD, FFTW_ESTIMATE
            );
            plan_bwd_ = fftw_plan_dft_2d(
                nx_, ny_, data_u, data_u, FFTW_BACKWARD, FFTW_ESTIMATE
            );

            if (!plan_fwd_ || !plan_bwd_) {
                destroy_fft();
                throw std::runtime_error(
                    "BdGFourierPreconditioner2D: FFTW plan creation failed"
                );
            }
        }

        ~BdGFourierPreconditioner2D() {
            destroy_fft();
        }

        BdGFourierPreconditioner2D(const BdGFourierPreconditioner2D&) = delete;
        BdGFourierPreconditioner2D& operator=(const BdGFourierPreconditioner2D&) = delete;

        void apply_inverse(
            const Eigen::VectorXcd& input,
            Eigen::VectorXcd& output
        ) const {
            if (input.size() != 2 * n_) {
                throw std::invalid_argument(
                    "BdGFourierPreconditioner2D::apply_inverse: input size mismatch"
                );
            }

            work_u_ = input.head(n_);
            work_v_ = input.tail(n_);

            fftw_execute_dft(
                plan_fwd_,
                reinterpret_cast<fftw_complex*>(work_u_.data()),
                reinterpret_cast<fftw_complex*>(work_u_.data())
            );
            fftw_execute_dft(
                plan_fwd_,
                reinterpret_cast<fftw_complex*>(work_v_.data()),
                reinterpret_cast<fftw_complex*>(work_v_.data())
            );

            for (int p = 0; p < n_; ++p) {
                const Scalar u_hat = work_u_(p);
                const Scalar v_hat = work_v_(p);
                work_u_(p) = inv11_(p) * u_hat + inv12_(p) * v_hat;
                work_v_(p) = inv21_(p) * u_hat + inv22_(p) * v_hat;
            }

            fftw_execute_dft(
                plan_bwd_,
                reinterpret_cast<fftw_complex*>(work_u_.data()),
                reinterpret_cast<fftw_complex*>(work_u_.data())
            );
            fftw_execute_dft(
                plan_bwd_,
                reinterpret_cast<fftw_complex*>(work_v_.data()),
                reinterpret_cast<fftw_complex*>(work_v_.data())
            );

            const double normalization = static_cast<double>(n_);
            work_u_ /= normalization;
            work_v_ /= normalization;

            output.resize(2 * n_);
            output.head(n_) = work_u_;
            output.tail(n_) = work_v_;
        }

    private:
        static double dipolar_kernel(double k, double lz, double strength) {
            if (strength == 0.0) {
                return 0.0;
            }

            const double q = k * lz / std::sqrt(2.0);
            const double form_factor =
                2.0 - 3.0 * std::sqrt(M_PI) * q * scaled_erfc_nonnegative(q);
            return strength * form_factor;
        }

        void build_inverse_blocks() {
            inv11_.resize(n_);
            inv12_.resize(n_);
            inv21_.resize(n_);
            inv22_.resize(n_);

            const double density_ref = op_.preconditioner_density();
            const double a_ref = op_.preconditioner_diag_a();
            const Scalar b_ref = op_.preconditioner_diag_b();
            const double lz = op_.confinement_length();
            const double dd_strength = op_.dipolar_strength();
            const Eigen::VectorXd& k2 = op_.kinetic_k2();

            for (int i = 0; i < nx_; ++i) {
                const int qx = (i <= nx_ / 2) ? i : i - nx_;
                const double kx = 2.0 * M_PI * static_cast<double>(qx) /
                                  (static_cast<double>(nx_) * dx_);

                for (int j = 0; j < ny_; ++j) {
                    const int qy = (j <= ny_ / 2) ? j : j - ny_;
                    const double ky = 2.0 * M_PI * static_cast<double>(qy) /
                                      (static_cast<double>(ny_) * dy_);
                    const int p = i * ny_ + j;

                    const double ddi = include_dipolar_
                        ? density_ref * dipolar_kernel(
                              std::hypot(kx, ky), lz, dd_strength)
                        : 0.0;

                    const Scalar a = Scalar(0.5 * k2(p) + a_ref + ddi, 0.0);
                    const Scalar b = b_ref + Scalar(ddi, 0.0);

                    const Scalar c11 = a - shift_;
                    const Scalar c12 = b;
                    const Scalar c21 = -std::conj(b);
                    const Scalar c22 = -std::conj(a) - shift_;
                    const Scalar determinant = c11 * c22 - c12 * c21;

                    Scalar stable_determinant = determinant;
                    if (std::abs(stable_determinant) < regularization_) {
                        stable_determinant += Scalar(regularization_, 0.0);
                    }
                    if (std::abs(stable_determinant) == 0.0) {
                        throw std::runtime_error(
                            "BdGFourierPreconditioner2D: singular Fourier block"
                        );
                    }

                    inv11_(p) = c22 / stable_determinant;
                    inv12_(p) = -c12 / stable_determinant;
                    inv21_(p) = -c21 / stable_determinant;
                    inv22_(p) = c11 / stable_determinant;
                }
            }
        }

        void destroy_fft() const {
            if (plan_fwd_) {
                fftw_destroy_plan(plan_fwd_);
                plan_fwd_ = nullptr;
            }
            if (plan_bwd_) {
                fftw_destroy_plan(plan_bwd_);
                plan_bwd_ = nullptr;
            }
        }
    };

    class BdGShiftInvertOperator2D {
    public:
        using Scalar = std::complex<double>;

    private:
        const BdGOperator2D& op_;
        Scalar shift_;
        int max_iter_ = 300;
        int restart_ = 40;
        double tol_ = 1e-8;
        int dense_max_dim_ = 2048;
        bool use_dense_direct_ = false;
        Eigen::FullPivLU<Eigen::MatrixXcd> dense_lu_;
        std::unique_ptr<BdGFourierPreconditioner2D> preconditioner_;

        mutable int total_linear_iterations_ = 0;
        mutable int linear_failures_ = 0;
        mutable double max_linear_residual_ = 0.0;

    public:
        BdGShiftInvertOperator2D(
            const BdGOperator2D& op,
            Scalar shift,
            int max_iter,
            int restart,
            double tol,
            int dense_max_dim,
            bool use_fourier_preconditioner,
            bool include_dipolar,
            double preconditioner_regularization
        ) :
            op_(op),
            shift_(shift),
            max_iter_(max_iter),
            restart_(restart),
            tol_(tol),
            dense_max_dim_(dense_max_dim)
        {
            if (op_.rows() != op_.cols()) {
                throw std::invalid_argument("BdGShiftInvertOperator2D: BdG operator must be square");
            }
            if (!std::isfinite(shift_.real()) || !std::isfinite(shift_.imag())) {
                throw std::invalid_argument("BdGShiftInvertOperator2D: shift must be finite");
            }
            if (max_iter_ <= 0 || restart_ <= 0 || tol_ <= 0.0) {
                throw std::invalid_argument("BdGShiftInvertOperator2D: invalid GMRES controls");
            }
            restart_ = std::min(restart_, op_.rows());

            if (dense_max_dim_ > 0 && op_.rows() <= dense_max_dim_) {
                dense_lu_.compute(assemble_shifted_matrix());
                if (!dense_lu_.isInvertible()) {
                    throw std::runtime_error("BdGShiftInvertOperator2D: dense shifted matrix is singular");
                }
                use_dense_direct_ = true;
            } else if (use_fourier_preconditioner) {
                preconditioner_ = std::make_unique<BdGFourierPreconditioner2D>(
                    op_, shift_, include_dipolar, preconditioner_regularization
                );
            }
        }

        int rows() const {
            return op_.rows();
        }

        int cols() const {
            return op_.cols();
        }

        int total_linear_iterations() const {
            return total_linear_iterations_;
        }

        int linear_failures() const {
            return linear_failures_;
        }

        double max_linear_residual() const {
            return max_linear_residual_;
        }

        bool uses_fourier_preconditioner() const {
            return static_cast<bool>(preconditioner_);
        }

        bool uses_dense_direct() const {
            return use_dense_direct_;
        }

        void perform_op(const Scalar* x_in, Scalar* y_out) const {
            Eigen::Map<const Eigen::VectorXcd> b(x_in, rows());
            Eigen::Map<Eigen::VectorXcd> y(y_out, rows());

            if (use_dense_direct_) {
                y = dense_lu_.solve(b);
                ++total_linear_iterations_;
                return;
            }

            int iterations = 0;
            bool converged = false;
            y = solve_shifted_system(b, iterations, converged);
            total_linear_iterations_ += iterations;

            Eigen::VectorXcd residual;
            apply_shifted(y, residual);
            const double b_norm = std::max(b.norm(), 1e-300);
            const double relative_residual = (residual - b).norm() / b_norm;
            max_linear_residual_ = std::max(
                max_linear_residual_, relative_residual
            );

            if (!converged || relative_residual > tol_) {
                ++linear_failures_;
            }
        }

    private:
        void apply_preconditioner(
            const Eigen::VectorXcd& input,
            Eigen::VectorXcd& output
        ) const {
            if (preconditioner_) {
                preconditioner_->apply_inverse(input, output);
            } else {
                output = input;
            }
        }

        void apply_preconditioned_shifted(
            const Eigen::VectorXcd& input,
            Eigen::VectorXcd& output
        ) const {
            Eigen::VectorXcd physical_input;
            apply_preconditioner(input, physical_input);
            apply_shifted(physical_input, output);
        }

        Eigen::MatrixXcd assemble_shifted_matrix() const {
            const int n = rows();
            Eigen::MatrixXcd matrix(n, n);
            Eigen::VectorXcd basis = Eigen::VectorXcd::Zero(n);
            Eigen::VectorXcd column(n);

            for (int j = 0; j < n; ++j) {
                basis(j) = 1.0;
                apply_shifted(basis, column);
                matrix.col(j) = column;
                basis(j) = 0.0;
            }

            return matrix;
        }

        void apply_shifted(
            const Eigen::VectorXcd& x,
            Eigen::VectorXcd& y
        ) const {
            if (x.size() != rows()) {
                throw std::invalid_argument("BdGShiftInvertOperator2D::apply_shifted: input size mismatch");
            }

            y.resize(rows());
            op_.perform_op(x.data(), y.data());
            y.array() -= shift_ * x.array();
        }

        Eigen::VectorXcd solve_shifted_system(
            const Eigen::VectorXcd& b,
            int& iterations,
            bool& converged
        ) const {
            const int n = rows();
            if (b.size() != n) {
                throw std::invalid_argument("BdGShiftInvertOperator2D::solve_shifted_system: RHS size mismatch");
            }

            iterations = 0;
            converged = false;

            const double b_norm = b.norm();
            Eigen::VectorXcd x = Eigen::VectorXcd::Zero(n);
            if (b_norm == 0.0) {
                converged = true;
                return x;
            }

            Eigen::VectorXcd ax(n);
            Eigen::VectorXcd r = b;
            double rel_residual = r.norm() / b_norm;
            if (rel_residual <= tol_) {
                converged = true;
                return x;
            }

            while (iterations < max_iter_) {
                const int cycle_dim = std::min(restart_, max_iter_ - iterations);

                std::vector<Eigen::VectorXcd> basis;
                basis.reserve(static_cast<std::size_t>(cycle_dim + 1));

                const double beta = r.norm();
                if (beta == 0.0) {
                    converged = true;
                    Eigen::VectorXcd physical_x;
                    apply_preconditioner(x, physical_x);
                    return physical_x;
                }
                basis.push_back(r / beta);

                Eigen::MatrixXcd hessenberg =
                    Eigen::MatrixXcd::Zero(cycle_dim + 1, cycle_dim);
                Eigen::VectorXcd rhs =
                    Eigen::VectorXcd::Zero(cycle_dim + 1);
                rhs(0) = beta;

                Eigen::VectorXcd best_x = x;
                double best_rel_residual = rel_residual;

                for (int j = 0; j < cycle_dim; ++j) {
                    Eigen::VectorXcd w(n);
                    apply_preconditioned_shifted(
                        basis[static_cast<std::size_t>(j)], w
                    );

                    for (int i = 0; i <= j; ++i) {
                        hessenberg(i, j) = basis[static_cast<std::size_t>(i)].dot(w);
                        w -= hessenberg(i, j) * basis[static_cast<std::size_t>(i)];
                    }

                    const double next_norm = w.norm();
                    hessenberg(j + 1, j) = next_norm;
                    if (next_norm > 1e-14) {
                        basis.push_back(w / next_norm);
                    } else {
                        basis.push_back(Eigen::VectorXcd::Zero(n));
                    }

                    const Eigen::MatrixXcd h_active =
                        hessenberg.block(0, 0, j + 2, j + 1);
                    const Eigen::VectorXcd rhs_active =
                        rhs.segment(0, j + 2);
                    const Eigen::VectorXcd y_small =
                        h_active.colPivHouseholderQr().solve(rhs_active);

                    Eigen::VectorXcd candidate = x;
                    for (int i = 0; i <= j; ++i) {
                        candidate += y_small(i) * basis[static_cast<std::size_t>(i)];
                    }

                    Eigen::VectorXcd candidate_physical;
                    apply_preconditioner(candidate, candidate_physical);
                    apply_shifted(candidate_physical, ax);
                    const double candidate_rel_residual = (b - ax).norm() / b_norm;
                    ++iterations;

                    if (candidate_rel_residual < best_rel_residual) {
                        best_rel_residual = candidate_rel_residual;
                        best_x = std::move(candidate);
                    }

                    if (candidate_rel_residual <= tol_) {
                        converged = true;
                        return candidate_physical;
                    }

                    if (iterations >= max_iter_) {
                        break;
                    }
                }

                x = best_x;
                Eigen::VectorXcd physical_x;
                apply_preconditioner(x, physical_x);
                apply_shifted(physical_x, ax);
                r = b - ax;
                rel_residual = r.norm() / b_norm;
                if (rel_residual <= tol_) {
                    converged = true;
                    return physical_x;
                }
            }

            Eigen::VectorXcd physical_x;
            apply_preconditioner(x, physical_x);
            return physical_x;
        }
    };

    template<Dimension Dim, typename SolverType = void>
    class BdGSolver;


    template<typename SolverType>
    class BdGSolver<Dimension::Two, SolverType>
    {
    public:
        static constexpr Dimension Dim = Dimension::Two;

        using WaveFuncType = gpes::WaveFunction<Dim>;
        using GridType = gpes::Grid<Dim>;
        using ShrdPtrGrid = std::shared_ptr<const GridType>;
        using Tag = gpes::tags::BdG;

        inline static const gpes::log::LogCategory kLogCategory{"BdGSolver2D"};

    private:
        ShrdPtrGrid grid_;
        PhysConfig phys_;
        BdGConfig config_;

        double mu_ = std::numeric_limits<double>::quiet_NaN();

    public:
        BdGSolver(
            ShrdPtrGrid grid,
            PhysConfig phys,
            BdGConfig config = {}
        ) :
            grid_(std::move(grid)),
            phys_(phys),
            config_(config)
        {
            if (!grid_) {
                GPES_LOG(kLogCategory, Error, "BdGSolver requires a valid grid");
                throw std::runtime_error("BdGSolver requires a valid grid");
            }

            if (grid_->size_x() <= 0 || grid_->size_y() <= 0) {
                throw std::invalid_argument("BdGSolver<Dimension::Two>: grid sizes must be positive");
            }

            if (!std::isfinite(grid_->step_x()) || grid_->step_x() == 0.0 ||
                !std::isfinite(grid_->step_y()) || grid_->step_y() == 0.0) {
                throw std::invalid_argument("BdGSolver<Dimension::Two>: invalid grid spacing");
            }

            if (config_.num_modes <= 0) {
                throw std::invalid_argument("BdGSolver<Dimension::Two>: num_modes must be positive");
            }
            if (config_.ncv <= config_.num_modes) {
                config_.ncv = 2 * config_.num_modes + 1;
            }
        }

        BdGResult solve(const WaveFuncType& psi0) {
            if (!psi0.grid()) {
                throw std::invalid_argument("BdGSolver<Dimension::Two>::solve: psi0 has no grid");
            }

            const int n = static_cast<int>(grid_->size_x() * grid_->size_y());

            if (psi0.size() != n) {
                throw std::invalid_argument("BdGSolver<Dimension::Two>::solve: psi0 and grid size mismatch");
            }

            mu_ = estimate_chemical_potential(psi0.vec());
            GPES_LOG(
                kLogCategory,
                Log,
                "BdG chemical potential estimated: mu={}",
                mu_
            );

            return solve(psi0, mu_);
        }

        BdGResult solve(const WaveFuncType& psi0, double chemical_potential) {
            const int n = static_cast<int>(grid_->size_x() * grid_->size_y());
            const int bdg_dim = 2 * n;

            if (psi0.size() != n) {
                throw std::invalid_argument("BdGSolver<Dimension::Two>::solve: psi0 and grid size mismatch");
            }

            if (!std::isfinite(chemical_potential)) {
                throw std::invalid_argument("BdGSolver<Dimension::Two>::solve: chemical potential must be finite");
            }

            const int nev = config_.num_modes;
            int ncv = config_.ncv;

            if (nev >= bdg_dim - 1) {
                throw std::invalid_argument("BdGSolver<Dimension::Two>::solve: too many requested modes");
            }

            ncv = std::max(ncv, 2 * nev + 1);
            ncv = std::min(ncv, bdg_dim);
            if (ncv <= nev) {
                throw std::invalid_argument("BdGSolver<Dimension::Two>::solve: ncv must be larger than num_modes");
            }

            GPES_LOG(
                kLogCategory,
                Log,
                "BdG run start: grid={}x{}, dimension={}, modes={}, ncv={}, mu={}, method={}",
                grid_->size_x(),
                grid_->size_y(),
                bdg_dim,
                nev,
                ncv,
                chemical_potential,
                config_.use_shift_invert ? "shift-invert" : "standard Arnoldi"
            );

            BdGOperator2D op(grid_, psi0.vec(), phys_, chemical_potential);

            BdGResult result;
            Eigen::VectorXcd evals;
            Eigen::MatrixXcd evecs;

            if (config_.use_shift_invert) {
                BdGShiftInvertOperator2D shift_op(
                    op,
                    config_.shift,
                    config_.shift_invert_linear_max_iter,
                    config_.shift_invert_linear_restart,
                    config_.shift_invert_linear_tol,
                    config_.shift_invert_dense_max_dim,
                    config_.use_fourier_preconditioner,
                    config_.preconditioner_include_dipolar,
                    config_.preconditioner_regularization
                );

                GPES_LOG(
                    kLogCategory,
                    Log,
                    "BdG shift-invert setup: shift=({}, {}), linear_solver={}, "
                    "fourier_preconditioner={}, linear_tol={}",
                    config_.shift.real(),
                    config_.shift.imag(),
                    shift_op.uses_dense_direct() ? "dense direct" : "GMRES",
                    shift_op.uses_fourier_preconditioner() ? "enabled" : "disabled",
                    config_.shift_invert_linear_tol
                );

                Spectra::GenEigsSolver<BdGShiftInvertOperator2D> eigs(shift_op, nev, ncv);

                eigs.init();
                const int nconv = eigs.compute(
                    Spectra::SortRule::LargestMagn,
                    config_.max_iter,
                    config_.tol
                );

                result.nconv = nconv;
                result.niter = eigs.num_iterations();
                result.nops = eigs.num_operations();
                result.info = eigs.info();
                result.used_shift_invert = true;
                result.shift = config_.shift;
                result.linear_iterations = shift_op.total_linear_iterations();
                result.linear_failures = shift_op.linear_failures();
                result.max_linear_residual = shift_op.max_linear_residual();
                result.used_fourier_preconditioner =
                    shift_op.uses_fourier_preconditioner();
                result.used_dense_direct = shift_op.uses_dense_direct();

                if (eigs.info() != Spectra::CompInfo::Successful) {
                    GPES_LOG(
                        kLogCategory,
                        Warning,
                        "BdG shift-invert Spectra solve did not fully converge: nconv = {}",
                        nconv
                    );
                }
                if (result.linear_failures > 0) {
                    GPES_LOG(
                        kLogCategory,
                        Warning,
                        "BdG shift-invert GMRES had {} unconverged linear solves",
                        result.linear_failures
                    );
                }

                const Eigen::VectorXcd inverse_evals = eigs.eigenvalues();
                evals.resize(inverse_evals.size());
                for (int m = 0; m < inverse_evals.size(); ++m) {
                    if (std::abs(inverse_evals(m)) == 0.0) {
                        evals(m) = std::numeric_limits<double>::infinity();
                    } else {
                        evals(m) = config_.shift + std::complex<double>(1.0, 0.0) / inverse_evals(m);
                    }
                }
                evecs = eigs.eigenvectors();
            } else {
                Spectra::GenEigsSolver<BdGOperator2D> eigs(op, nev, ncv);

                eigs.init();
                const int nconv = eigs.compute(
                    config_.sort_rule,
                    config_.max_iter,
                    config_.tol
                );

                result.nconv = nconv;
                result.niter = eigs.num_iterations();
                result.nops = eigs.num_operations();
                result.info = eigs.info();

                if (eigs.info() != Spectra::CompInfo::Successful) {
                    GPES_LOG(
                        kLogCategory,
                        Warning,
                        "BdG Spectra solve did not fully converge: nconv = {}",
                        nconv
                    );
                }

                evals = eigs.eigenvalues();
                evecs = eigs.eigenvectors();
            }

            result.modes.reserve(static_cast<std::size_t>(evals.size()));

            for (int m = 0; m < evals.size(); ++m) {
                BdGMode mode;
                mode.omega = evals(m);
                mode.u = evecs.col(m).segment(0, n);
                mode.v = evecs.col(m).segment(n, n);

                normalize_bdg_mode(mode.u, mode.v);

                Eigen::VectorXcd full_mode(2 * n);
                full_mode.head(n) = mode.u;
                full_mode.tail(n) = mode.v;
                Eigen::VectorXcd applied_mode(2 * n);
                op.perform_op(full_mode.data(), applied_mode.data());
                const double denominator = std::max({
                    applied_mode.norm(),
                    std::abs(mode.omega) * full_mode.norm(),
                    1e-300
                });
                mode.relative_residual =
                    (applied_mode - mode.omega * full_mode).norm() / denominator;
                result.max_eigen_residual = std::max(
                    result.max_eigen_residual, mode.relative_residual);

                result.modes.push_back(std::move(mode));
            }

            mu_ = chemical_potential;
            if (result.used_shift_invert) {
                GPES_LOG(
                    kLogCategory,
                    Log,
                    "BdG run finish: converged={}/{}, iterations={}, operator_applications={}, "
                    "linear_iterations={}, linear_failures={}, max_linear_residual={}, "
                    "max_eigen_residual={}",
                    result.nconv,
                    nev,
                    result.niter,
                    result.nops,
                    result.linear_iterations,
                    result.linear_failures,
                    result.max_linear_residual,
                    result.max_eigen_residual
                );
            } else {
                GPES_LOG(
                    kLogCategory,
                    Log,
                    "BdG run finish: converged={}/{}, iterations={}, operator_applications={}, "
                    "max_eigen_residual={}",
                    result.nconv,
                    nev,
                    result.niter,
                    result.nops,
                    result.max_eigen_residual
                );
            }
            return result;
        }

        double chemical_potential() const {
            return mu_;
        }

    private:
        double estimate_chemical_potential(const Eigen::VectorXcd& psi) const {
            const int nx = static_cast<int>(grid_->size_x());
            const int ny = static_cast<int>(grid_->size_y());
            const int n = nx * ny;

            if (psi.size() != n) {
                throw std::invalid_argument("estimate_chemical_potential: psi size mismatch");
            }

            const double dx = grid_->step_x();
            const double dy = grid_->step_y();
            const double measure = dx * dy;

            const double lz = 1.0 / std::sqrt(grid_->omega_z());

            PhysConfig phys = phys_;
            calc_inter_consts<Dim>(phys, lz);
            gpes::solvers::detail::validate_inter_consts(phys, "BdGSolver<Dimension::Two>");

            DipolarInteraction<Dim> ddi(grid_, lz, phys.V_dd);

            Eigen::VectorXd phi_ddi;
            ddi.compute_DDI_term(psi, phi_ddi);

            Eigen::VectorXcd kinetic;
            apply_spectral_kinetic(psi, kinetic);

            double numerator = 0.0;
            double norm = 0.0;

            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    const int k = i * ny + j;

                    const double density = std::norm(psi(k));
                    const double V = grid_->potential()(i, j);

                    const std::complex<double> Hpsi =
                        kinetic(k)
                        + (V
                        + phys.g_scat * density
                        + phi_ddi(k)
                        + phys.g_lhy * std::pow(density, 1.5)) * psi(k);

                    numerator += measure * std::real(std::conj(psi(k)) * Hpsi);
                    norm += measure * density;
                }
            }

            if (norm <= 0.0) {
                throw std::runtime_error("estimate_chemical_potential: zero norm state");
            }

            return numerator / norm;
        }

        void apply_spectral_kinetic(
            const Eigen::VectorXcd& in,
            Eigen::VectorXcd& out
        ) const {
            const int nx = static_cast<int>(grid_->size_x());
            const int ny = static_cast<int>(grid_->size_y());
            const int n = nx * ny;

            if (in.size() != n) {
                throw std::invalid_argument("BdGSolver<Dimension::Two>::apply_spectral_kinetic: input size mismatch");
            }

            const double dx = grid_->step_x();
            const double dy = grid_->step_y();
            const double Lx = static_cast<double>(nx) * dx;
            const double Ly = static_cast<double>(ny) * dy;
            const double dkx = 2.0 * M_PI / Lx;
            const double dky = 2.0 * M_PI / Ly;

            Eigen::VectorXd k2(n);

            for (int i = 0; i < nx; ++i) {
                const int qx = (i <= nx / 2) ? i : i - nx;
                const double kx = dkx * static_cast<double>(qx);
                for (int j = 0; j < ny; ++j) {
                    const int qy = (j <= ny / 2) ? j : j - ny;
                    const double ky = dky * static_cast<double>(qy);
                    k2(i * ny + j) = kx * kx + ky * ky;
                }
            }

            Eigen::VectorXcd work = in;
            auto* data = reinterpret_cast<fftw_complex*>(work.data());

            fftw_plan plan_fwd = fftw_plan_dft_2d(
                nx,
                ny,
                data,
                data,
                FFTW_FORWARD,
                FFTW_ESTIMATE
            );

            fftw_plan plan_bwd = fftw_plan_dft_2d(
                nx,
                ny,
                data,
                data,
                FFTW_BACKWARD,
                FFTW_ESTIMATE
            );

            if (!plan_fwd || !plan_bwd) {
                if (plan_fwd) {
                    fftw_destroy_plan(plan_fwd);
                }
                if (plan_bwd) {
                    fftw_destroy_plan(plan_bwd);
                }
                throw std::runtime_error("BdGSolver<Dimension::Two>::apply_spectral_kinetic: FFTW plan creation failed");
            }

            fftw_execute_dft(plan_fwd, data, data);
            work.array() *= (0.5 * k2.array()).cast<std::complex<double>>();
            fftw_execute_dft(plan_bwd, data, data);
            work /= static_cast<double>(n);

            fftw_destroy_plan(plan_fwd);
            fftw_destroy_plan(plan_bwd);

            out = work;
        }

        void normalize_bdg_mode(Eigen::VectorXcd& u, Eigen::VectorXcd& v) const {
            const double dx = grid_->step_x();
            const double dy = grid_->step_y();

            double qnorm = 0.0;

            for (int i = 0; i < u.size(); ++i) {
                qnorm += std::norm(u(i)) - std::norm(v(i));
            }

            qnorm *= dx * dy;

            if (std::abs(qnorm) < 1e-14) {
                return;
            }

            const double scale = 1.0 / std::sqrt(std::abs(qnorm));
            u *= scale;
            v *= scale;
        }



    };

    using BdGSolver2D = BdGSolver<Dimension::Two>;

    inline double bdg_mode_norm(
        const BdGMode& mode,
        double dx,
        double dy
    ) {
        if (mode.u.size() != mode.v.size()) {
            throw std::invalid_argument("bdg_mode_norm: u and v sizes differ");
        }

        double qnorm = 0.0;
        for (Eigen::Index i = 0; i < mode.u.size(); ++i) {
            qnorm += std::norm(mode.u(i)) - std::norm(mode.v(i));
        }

        return qnorm * dx * dy;
    }

    inline std::string mode_filename(std::size_t mode_index) {
        std::ostringstream name;
        name << "mode_" << std::setw(4) << std::setfill('0') << mode_index << ".csv";
        return name.str();
    }

    inline std::filesystem::path save_bdg_result(
        const BdGResult& result,
        const Grid<Dimension::Two>& grid,
        const std::filesystem::path& destination_folder,
        bool save_mode_fields = true
    ) {
        const int nx = static_cast<int>(grid.size_x());
        const int ny = static_cast<int>(grid.size_y());
        const int n = nx * ny;
        const double dx = grid.step_x();
        const double dy = grid.step_y();

        std::filesystem::create_directories(destination_folder);
        const std::filesystem::path modes_folder = destination_folder / "excited_states";
        std::filesystem::create_directories(modes_folder);

        const std::filesystem::path spectrum_path = destination_folder / "spectrum.csv";
        std::ofstream spectrum(spectrum_path);
        if (!spectrum) {
            throw std::runtime_error("save_bdg_result: cannot open file: " + spectrum_path.string());
        }

        spectrum << "mode,omega_real,omega_imag,bdg_norm,relative_residual,state_file\n";
        spectrum << std::setprecision(17);

        for (std::size_t m = 0; m < result.modes.size(); ++m) {
            const BdGMode& mode = result.modes[m];
            if (mode.u.size() != n || mode.v.size() != n) {
                throw std::invalid_argument("save_bdg_result: mode size does not match grid size");
            }

            const std::string state_filename = mode_filename(m);
            const std::filesystem::path state_path = modes_folder / state_filename;
            const double qnorm = bdg_mode_norm(mode, dx, dy);

            spectrum << m << ','
                     << mode.omega.real() << ','
                     << mode.omega.imag() << ','
                     << qnorm << ','
                     << mode.relative_residual << ','
                     << (save_mode_fields ? "excited_states/" + state_filename : "") << '\n';

            if (!save_mode_fields) {
                continue;
            }

            std::ofstream state(state_path);
            if (!state) {
                throw std::runtime_error("save_bdg_result: cannot open file: " + state_path.string());
            }

            state << "ix,iy,x,y,u_real,u_imag,u_abs2,v_real,v_imag,v_abs2,local_bdg_weight\n";
            state << std::setprecision(17);

            for (int i = 0; i < nx; ++i) {
                const double x = grid.start_pos_x() + static_cast<double>(i) * dx;
                for (int j = 0; j < ny; ++j) {
                    const double y = grid.start_pos_y() + static_cast<double>(j) * dy;
                    const int k = i * ny + j;
                    const double local_bdg_weight =
                        std::norm(mode.u(k)) - std::norm(mode.v(k));

                    state << i << ','
                          << j << ','
                          << x << ','
                          << y << ','
                          << mode.u(k).real() << ','
                          << mode.u(k).imag() << ','
                          << std::norm(mode.u(k)) << ','
                          << mode.v(k).real() << ','
                          << mode.v(k).imag() << ','
                          << std::norm(mode.v(k)) << ','
                          << local_bdg_weight << '\n';
                }
            }
            state.close();
            if (!state) throw std::runtime_error("save_bdg_result: write failed: " + state_path.string());
        }

        spectrum.close();
        if (!spectrum) throw std::runtime_error("save_bdg_result: write failed: " + spectrum_path.string());

        return spectrum_path;
    }

    inline std::filesystem::path save_bdg_result(
        const BdGResult& result,
        const std::shared_ptr<const Grid<Dimension::Two>>& grid,
        const std::filesystem::path& destination_folder,
        bool save_mode_fields = true
    ) {
        if (!grid) {
            throw std::invalid_argument("save_bdg_result: grid is null");
        }

        return save_bdg_result(result, *grid, destination_folder, save_mode_fields);
    }

} // namespace gpes::bdg
