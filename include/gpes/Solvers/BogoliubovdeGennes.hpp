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
    };

    struct BdGResult {
        std::vector<BdGMode> modes;
        int nconv = 0;
        int niter = 0;
        int nops = 0;
        Spectra::CompInfo info = Spectra::CompInfo::NotComputed;
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

            BdGOperator2D op(grid_, psi0.vec(), phys_, chemical_potential);

            Spectra::GenEigsSolver<BdGOperator2D> eigs(op, nev, ncv);

            eigs.init();
            const int nconv = eigs.compute(
                config_.sort_rule,
                config_.max_iter,
                config_.tol
            );

            BdGResult result;
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

            const Eigen::VectorXcd evals = eigs.eigenvalues();
            const Eigen::MatrixXcd evecs = eigs.eigenvectors();

            result.modes.reserve(static_cast<std::size_t>(evals.size()));

            for (int m = 0; m < evals.size(); ++m) {
                BdGMode mode;
                mode.omega = evals(m);
                mode.u = evecs.col(m).segment(0, n);
                mode.v = evecs.col(m).segment(n, n);

                normalize_bdg_mode(mode.u, mode.v);

                result.modes.push_back(std::move(mode));
            }

            mu_ = chemical_potential;
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
        const std::filesystem::path& destination_folder
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

        spectrum << "mode,omega_real,omega_imag,bdg_norm,state_file\n";
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
                     << ("excited_states/" + state_filename) << '\n';

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
        }

        return spectrum_path;
    }

    inline std::filesystem::path save_bdg_result(
        const BdGResult& result,
        const std::shared_ptr<const Grid<Dimension::Two>>& grid,
        const std::filesystem::path& destination_folder
    ) {
        if (!grid) {
            throw std::invalid_argument("save_bdg_result: grid is null");
        }

        return save_bdg_result(result, *grid, destination_folder);
    }

} // namespace gpes::bdg
