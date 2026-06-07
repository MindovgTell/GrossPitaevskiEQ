#pragma once

#include <fftw3.h>
#include <Eigen/Dense>
#include <memory>

#include "Core/definitions.hpp"
#include "grid/grid.hpp"

namespace gpes {

    template <Dimension Dim>
    class DipolarInteraction;

    #ifndef M_PI
    constexpr double M_PI       =   3.14159265358979323846;
    #endif
    constexpr double SQRT_PI    =   1.77245385090551602729;
    constexpr double SQRT2      =   1.41421356237309504880;


//********************************/***********/********************************//
//********************************/***********/********************************//
//*************************/One dimensional FFT sim/***************************//
//********************************/***********/********************************//
//********************************/***********/********************************//


    template <>
    class DipolarInteraction<Dimension::One> {
    public:
        inline static constexpr Dimension Dim = Dimension::One;
        using ShrdPtrGrid = std::shared_ptr<const Grid<Dim>>;

    private:
        int N;
        double L;
        double dx;
        double lz;
        double Cdd;

        fftw_plan plan_forward, plan_backward;
        double* density_real;
        fftw_complex* density_fourier;
        fftw_complex* potential_fourier;
        double* potential_real;

        Eigen::VectorXd U_tilde; // kernel in Fourier space

    public:
        explicit DipolarInteraction(ShrdPtrGrid& grid);

        DipolarInteraction(int size_of_grid, double length, double confinement_length, double dipolar_interaction_strength)
            : N(size_of_grid), L(length), lz(confinement_length), Cdd(dipolar_interaction_strength) {

            dx = L / N;

            // FFTW memory allocation
            density_real = (double*)fftw_malloc(sizeof(double) * N);
            density_fourier = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N/2 + 1));
            potential_fourier = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N/2 + 1));
            potential_real = (double*)fftw_malloc(sizeof(double) * N);

            // FFTW plans
            plan_forward = fftw_plan_dft_r2c_1d(N, density_real, density_fourier, FFTW_MEASURE);
            plan_backward = fftw_plan_dft_c2r_1d(N, potential_fourier, potential_real, FFTW_MEASURE);

            precompute_U_tilde();
        }

        ~DipolarInteraction() {
            fftw_destroy_plan(plan_forward);
            fftw_destroy_plan(plan_backward);
            fftw_free(density_real);
            fftw_free(density_fourier);
            fftw_free(potential_fourier);
            fftw_free(potential_real);
        }

    private:
        void precompute_U_tilde() {
            U_tilde.resize(N/2 + 1);
            for (int k = 0; k <= N/2; ++k) {
                double kx = 2.0 * M_PI * k / L;
                double kappa = std::abs(kx) * lz / std::sqrt(2.0);

                if (kx == 0.0) {
                    U_tilde[k] = 2.0 * Cdd / 3.0;
                } else {
                    U_tilde[k] = (Cdd / 3.0) * (2.0 - 3.0 * M_PI * kappa * std::exp(kappa * kappa) * std::erfc(kappa));
                }
            }
        }


    public:
        void compute_DDI_term(const Eigen::VectorXcd& psi, Eigen::VectorXd& Phi_DDI) {
            // Fill density array
            for (int i = 0; i < N; ++i) {
                density_real[i] = std::norm(psi(i));
            }

            // FFT
            fftw_execute(plan_forward);

            // Multiply by kernel in Fourier space
            for (int k = 0; k <= N/2; ++k) {
                potential_fourier[k][0] = U_tilde[k] * density_fourier[k][0]; // Real
                potential_fourier[k][1] = U_tilde[k] * density_fourier[k][1]; // Imag
            }

            // Inverse FFT
            fftw_execute(plan_backward);

            // Normalize and copy result
            Phi_DDI.resize(N);
            for (int i = 0; i < N; ++i) {
                Phi_DDI[i] = potential_real[i] / N;
            }
        }


    };


//********************************/***********/********************************//
//********************************/***********/********************************//
//*************************/Two dimensional FFT sim/***************************//
//********************************/***********/********************************//
//********************************/***********/********************************//

    template<>
    class DipolarInteraction<Dimension::Two> {
    public:
        inline static constexpr Dimension Dim = Dimension::Two;
        using ShrdPtrGrid = std::shared_ptr<const Grid<Dim>>;

    private:
        int Nx = 0;
        int Ny = 0;
        int Nx_pad = 0;
        int Ny_pad = 0;
        int Nk_y = 0;

        double Lx = 0.0;
        double Ly = 0.0;
        double Lx_pad = 0.0;
        double Ly_pad = 0.0;
        double dx = 0.0;
        double dy = 0.0;
        double lz = 0.0;
        double V_dd = 0.0;

        fftw_plan plan_forward = nullptr;
        fftw_plan plan_backward = nullptr;

        double* source_real = nullptr;
        fftw_complex* source_fourier = nullptr;
        fftw_complex* potential_fourier = nullptr;
        double* potential_real = nullptr;

        Eigen::MatrixXd U_tilde;

    public:
        DipolarInteraction(
            ShrdPtrGrid grid,
            double confinement_length,
            double interaction_strength
        ) :
            Nx(static_cast<int>(grid->size_x())),
            Ny(static_cast<int>(grid->size_y())),
            Nx_pad(2 * Nx),
            Ny_pad(2 * Ny),
            Nk_y(Ny_pad / 2 + 1),
            dx(grid->step_x()),
            dy(grid->step_y()),
            lz(confinement_length),
            V_dd(interaction_strength)
        {
            if (!grid) {
                throw std::runtime_error("DipolarInteraction<Dimension::Two>: grid is null");
            }

            if (Nx <= 0 || Ny <= 0) {
                throw std::invalid_argument("DipolarInteraction<Dimension::Two>: grid sizes must be positive");
            }

            if (!std::isfinite(dx) || dx == 0.0 ||
                !std::isfinite(dy) || dy == 0.0) {
                throw std::invalid_argument("DipolarInteraction<Dimension::Two>: grid steps must be finite and non-zero");
            }

            if (!std::isfinite(lz) || lz <= 0.0) {
                throw std::invalid_argument("DipolarInteraction<Dimension::Two>: confinement length must be positive and finite");
            }

            if (!std::isfinite(V_dd)) {
                throw std::invalid_argument("DipolarInteraction<Dimension::Two>: V_dd must be finite");
            }

            Lx = static_cast<double>(Nx) * dx;
            Ly = static_cast<double>(Ny) * dy;

            Lx_pad = static_cast<double>(Nx_pad) * dx;
            Ly_pad = static_cast<double>(Ny_pad) * dy;

            allocate_fftw();
            precompute_U_tilde();
        }

        ~DipolarInteraction() {
            destroy_fftw();
        }

        DipolarInteraction(const DipolarInteraction&) = delete;
        DipolarInteraction& operator=(const DipolarInteraction&) = delete;

        DipolarInteraction(DipolarInteraction&&) = delete;
        DipolarInteraction& operator=(DipolarInteraction&&) = delete;

        // -------------------------------------------------------------------------
        // Main GPE API:
        //
        //     Phi_DDI = U_dd * |psi|^2
        //
        // -------------------------------------------------------------------------
        void compute_DDI_term(
            const Eigen::VectorXcd& psi,
            Eigen::VectorXd& Phi_DDI
        ) {
            validate_complex_field(psi, "compute_DDI_term");

            Eigen::VectorXd density(Nx * Ny);

            for (int i = 0; i < Nx * Ny; ++i) {
                density(i) = std::norm(psi(i));
            }

            compute_DDI_from_real_source(density, Phi_DDI);
        }

        // -------------------------------------------------------------------------
        // Real-source convolution:
        //
        //     Phi_DDI = U_dd * source
        //
        // -------------------------------------------------------------------------
        void compute_DDI_from_real_source(
            const Eigen::VectorXd& source,
            Eigen::VectorXd& Phi_DDI
        ) {
            validate_real_field(source, "compute_DDI_from_real_source");

            clear_source();

            for (int i = 0; i < Nx; ++i) {
                for (int j = 0; j < Ny; ++j) {
                    const int index = get_index(i, j);
                    const int pindex = get_padded_index(i, j);

                    source_real[pindex] = source(index);
                }
            }

            fftw_execute(plan_forward);

            multiply_kernel();

            fftw_execute(plan_backward);

            copy_real_output(Phi_DDI);
        }

        // -------------------------------------------------------------------------
        // Complex-source convolution:
        //
        //     Phi_DDI = U_dd * source
        //
        // This is needed for BdG because the linearized density perturbation
        //
        //     delta_n = conj(psi0) * u + psi0 * v
        //
        // is complex in general.
        //
        // This implementation reuses the real r2c/c2r FFTW plans twice:
        // one convolution for Re(source), one for Im(source).
        // -------------------------------------------------------------------------
        void compute_DDI_from_complex_source(
            const Eigen::VectorXcd& source,
            Eigen::VectorXcd& Phi_DDI
        ) {
            validate_complex_field(source, "compute_DDI_from_complex_source");

            Eigen::VectorXd source_re(Nx * Ny);
            Eigen::VectorXd source_im(Nx * Ny);

            for (int i = 0; i < Nx * Ny; ++i) {
                source_re(i) = source(i).real();
                source_im(i) = source(i).imag();
            }

            Eigen::VectorXd phi_re;
            Eigen::VectorXd phi_im;

            compute_DDI_from_real_source(source_re, phi_re);
            compute_DDI_from_real_source(source_im, phi_im);

            Phi_DDI.resize(Nx * Ny);

            for (int i = 0; i < Nx * Ny; ++i) {
                Phi_DDI(i) = std::complex<double>(phi_re(i), phi_im(i));
            }
        }

        // -------------------------------------------------------------------------
        // BdG convenience API:
        //
        //     delta_Phi_DDI = U_dd * delta_n
        //
        // where
        //
        //     delta_n = conj(psi0) * u + psi0 * v
        //
        // -------------------------------------------------------------------------
        void compute_BdG_DDI_response(
            const Eigen::VectorXcd& psi0,
            const Eigen::VectorXcd& u,
            const Eigen::VectorXcd& v,
            Eigen::VectorXcd& delta_Phi_DDI
        ) {
            validate_complex_field(psi0, "compute_BdG_DDI_response: psi0");
            validate_complex_field(u, "compute_BdG_DDI_response: u");
            validate_complex_field(v, "compute_BdG_DDI_response: v");

            Eigen::VectorXcd delta_density(Nx * Ny);

            for (int i = 0; i < Nx * Ny; ++i) {
                delta_density(i) = std::conj(psi0(i)) * u(i) + psi0(i) * v(i);
            }

            compute_DDI_from_complex_source(delta_density, delta_Phi_DDI);
        }

        // Optional helper if you want to explicitly build the BdG perturbation.
        Eigen::VectorXcd compute_BdG_density_perturbation(
            const Eigen::VectorXcd& psi0,
            const Eigen::VectorXcd& u,
            const Eigen::VectorXcd& v
        ) const {
            validate_complex_field(psi0, "compute_BdG_density_perturbation: psi0");
            validate_complex_field(u, "compute_BdG_density_perturbation: u");
            validate_complex_field(v, "compute_BdG_density_perturbation: v");

            Eigen::VectorXcd delta_density(Nx * Ny);

            for (int i = 0; i < Nx * Ny; ++i) {
                delta_density(i) = std::conj(psi0(i)) * u(i) + psi0(i) * v(i);
            }

            return delta_density;
        }

    private:
        int get_index(int i, int j) const {
            return i * Ny + j;
        }

        int get_padded_index(int i, int j) const {
            return i * Ny_pad + j;
        }

        int get_fourier_index(int i, int j) const {
            return i * Nk_y + j;
        }

        void allocate_fftw() {
            const int real_size = Nx_pad * Ny_pad;
            const int fourier_size = Nx_pad * Nk_y;

            source_real = static_cast<double*>(
                fftw_malloc(sizeof(double) * real_size)
            );

            source_fourier = static_cast<fftw_complex*>(
                fftw_malloc(sizeof(fftw_complex) * fourier_size)
            );

            potential_fourier = static_cast<fftw_complex*>(
                fftw_malloc(sizeof(fftw_complex) * fourier_size)
            );

            potential_real = static_cast<double*>(
                fftw_malloc(sizeof(double) * real_size)
            );

            if (!source_real || !source_fourier || !potential_fourier || !potential_real) {
                destroy_fftw();
                throw std::runtime_error("DipolarInteraction<Dimension::Two>: FFTW allocation failed");
            }

            plan_forward = fftw_plan_dft_r2c_2d(
                Nx_pad,
                Ny_pad,
                source_real,
                source_fourier,
                FFTW_MEASURE
            );

            plan_backward = fftw_plan_dft_c2r_2d(
                Nx_pad,
                Ny_pad,
                potential_fourier,
                potential_real,
                FFTW_MEASURE
            );

            if (!plan_forward || !plan_backward) {
                destroy_fftw();
                throw std::runtime_error("DipolarInteraction<Dimension::Two>: FFTW plan creation failed");
            }
        }

        void destroy_fftw() {
            if (plan_forward) {
                fftw_destroy_plan(plan_forward);
                plan_forward = nullptr;
            }

            if (plan_backward) {
                fftw_destroy_plan(plan_backward);
                plan_backward = nullptr;
            }

            if (source_real) {
                fftw_free(source_real);
                source_real = nullptr;
            }

            if (source_fourier) {
                fftw_free(source_fourier);
                source_fourier = nullptr;
            }

            if (potential_fourier) {
                fftw_free(potential_fourier);
                potential_fourier = nullptr;
            }

            if (potential_real) {
                fftw_free(potential_real);
                potential_real = nullptr;
            }
        }

        void clear_source() {
            std::fill(
                source_real,
                source_real + Nx_pad * Ny_pad,
                0.0
            );
        }

        void multiply_kernel() {
        for (int i = 0; i < Nx_pad; ++i) {
            for (int j = 0; j < Nk_y; ++j) {
                const int k = get_fourier_index(i, j);
                const double U = U_tilde(i, j);

                potential_fourier[k][0] = U * source_fourier[k][0];
                potential_fourier[k][1] = U * source_fourier[k][1];
            }
        }
    }

        void copy_real_output(Eigen::VectorXd& Phi_DDI) const {
            Phi_DDI.resize(Nx * Ny);

            const double norm = static_cast<double>(Nx_pad * Ny_pad);

            for (int i = 0; i < Nx; ++i) {
                for (int j = 0; j < Ny; ++j) {
                    const int index = get_index(i, j);
                    const int pindex = get_padded_index(i, j);

                    Phi_DDI(index) = potential_real[pindex] / norm;
                }
            }
        }

        void validate_real_field(
            const Eigen::VectorXd& field,
            const char* where
        ) const {
            if (field.size() != Nx * Ny) {
                throw std::invalid_argument(
                    std::string("DipolarInteraction<Dimension::Two>::") +
                    where +
                    ": field size mismatch"
                );
            }
        }

        void validate_complex_field(
            const Eigen::VectorXcd& field,
            const char* where
        ) const {
            if (field.size() != Nx * Ny) {
                throw std::invalid_argument(
                    std::string("DipolarInteraction<Dimension::Two>::") +
                    where +
                    ": field size mismatch"
                );
            }
        }

        double F_perp(double q) const {
            const double q2 = q * q;

            if (q2 < 700.0) {
                const double erfc_val = std::erfc(q);
                const double exp_val = std::exp(q2);

                const double answ =
                    2.0 - 3.0 * SQRT_PI * q * exp_val * erfc_val;

                if (!std::isfinite(answ)) {
                    throw std::runtime_error("DipolarInteraction<Dimension::Two>::F_perp produced non-finite value");
                }

                return answ;
            }

            return -1.0;
        }

        double F_parallel(double q_x, double q_y) const {
            const double q = std::sqrt(q_x * q_x + q_y * q_y);

            if (q == 0.0) {
                return -1.0;
            }

            const double q2 = q * q;

            if (q2 < 700.0) {
                const double erfc_val = std::erfc(q);
                const double exp_val = std::exp(q2);

                const double answ =
                    -1.0 + 3.0 * SQRT_PI * (q_x * q_x / q) * exp_val * erfc_val;

                if (!std::isfinite(answ)) {
                    throw std::runtime_error("DipolarInteraction<Dimension::Two>::F_parallel produced non-finite value");
                }

                return answ;
            }

            return 2.0;
        }

        void precompute_U_tilde() {
            U_tilde.resize(Nx_pad, Nk_y);

            for (int i = 0; i < Nx_pad; ++i) {
                const double kx =
                    (i < Nx_pad / 2)
                        ? (2.0 * M_PI * static_cast<double>(i) / Lx_pad)
                        : (2.0 * M_PI * static_cast<double>(i - Nx_pad) / Lx_pad);

                for (int j = 0; j < Nk_y; ++j) {
                    const double ky =
                        2.0 * M_PI * static_cast<double>(j) / Ly_pad;

                    const double k_perp = std::sqrt(kx * kx + ky * ky);

                    if (k_perp == 0.0) {
                        U_tilde(i, j) = 2.0 * V_dd;
                    } else {
                        const double kappa = k_perp * lz / SQRT2;
                        U_tilde(i, j) = V_dd * F_perp(kappa);
                    }

                    // If you later want in-plane polarization, replace the block
                    // above with something like:
                    //
                    // const double q_x = kx * lz / SQRT2;
                    // const double q_y = ky * lz / SQRT2;
                    // U_tilde(i, j) = V_dd * F_parallel(q_x, q_y);
                }
            }
        }

    };
}
