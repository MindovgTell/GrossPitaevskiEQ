#pragma once

// #include "ReadSimParams.hpp"
// #include "visualitsation.hpp"
// #include "SaveSimResults.hpp"

#include <cmath>
#include <cstddef>
#include <complex>
#include <stdexcept>
#include <random>

#include "definitions.hpp"
#include "wavefunction/wavefunction.hpp"
#include "wavefunction_io.hpp"

namespace gpes {

    class NonCopyable {
    protected:
        NonCopyable() = default;
        ~NonCopyable() = default;

        NonCopyable(const NonCopyable&) = delete;
        NonCopyable& operator=(const NonCopyable&) = delete;

        NonCopyable(NonCopyable&&) = delete;
        NonCopyable& operator=(NonCopyable&&) = delete;
    }; 


    struct SimConfig {
        double duration;
        double dt;
        size_t num_of_steps;
        bool ground_state = true;
    };

    struct PhysConfig {
        double a_s;
        double a_dd;
        double num_of_prt;
        double g_scat = 0.0;
        double V_dd = 0.0;
        double g_lhy = 0.0;
        double l_ = 0.0; // variable which correspond to wavelength in reduced dimensions
    };

    template <Dimension D>
    void calc_inter_consts(PhysConfig&, double);

    enum class ModulationAxis {
        X,
        Y,
    };

    enum class SeedType {
        Sinusoidal,
        RandomNoise
    };

    namespace detail {
        inline void validate_modulation_inputs(
            double amplitude,
            double wavelength,
            double phase,
            const char* context
        ) {
            if (!std::isfinite(amplitude)) {
                throw std::invalid_argument(std::string(context) + ": amplitude must be finite");
            }
            if (!std::isfinite(wavelength) || wavelength == 0.0) {
                throw std::invalid_argument(std::string(context) + ": wavelength must be non-zero and finite");
            }
            if (!std::isfinite(phase)) {
                throw std::invalid_argument(std::string(context) + ": phase must be finite");
            }
        }

        inline double current_norm(const WaveFunction<Dimension::One>& psi) {
            if (!psi.grid()) {
                throw std::invalid_argument("current_norm(1D): wavefunction has no grid");
            }

            const double dx = psi.grid()->step();
            double norm = 0.0;
            for (Eigen::Index i = 0; i < psi.size(); ++i) {
                norm += std::norm(psi.vec()(i));
            }
            return norm * dx;
        }

        inline double current_norm(const WaveFunction<Dimension::Two>& psi) {
            if (!psi.grid()) {
                throw std::invalid_argument("current_norm(2D): wavefunction has no grid");
            }

            const double dA = psi.grid()->step_x() * psi.grid()->step_y();
            double norm = 0.0;
            for (Eigen::Index i = 0; i < psi.size(); ++i) {
                norm += std::norm(psi.vec()(i));
            }
            return norm * dA;
        }


        inline void validate_noise_inputs(
            double amplitude,
            unsigned int /*seed*/,
            const char* function_name
        ) {
            if (!std::isfinite(amplitude) || amplitude < 0.0) {
                throw std::invalid_argument(std::string(function_name) + ": amplitude must be non-negative and finite");
            }
        }





        // ----------------------------------------------------------------
        //  64-point Gauss-Legendre quadrature on [-1, 1].
        //  Only 32 positive nodes stored; negative nodes follow by symmetry.
        // ----------------------------------------------------------------
        namespace gl64 {

        inline constexpr std::array<double, 32> kNodes = {{
            0.02435069280380807, 0.07299312178779904, 0.12146281929612057,
            0.16964441042399283, 0.21742364374000720, 0.26468716220876750,
            0.31132287199021144, 0.35722015833737524, 0.40227015796399163,
            0.44636601725346409, 0.48940314570705264, 0.53127946401989455,
            0.57189564620263404, 0.61115565518529115, 0.64896547125465734,
            0.68523631305423258, 0.71988185017161082, 0.75281990726053190,
            0.78397296107852894, 0.81326531512279756, 0.84062929625258036,
            0.86599939815409282, 0.88931544599511410, 0.91052213707850214,
            0.92956917213193958, 0.94641137485840268, 0.96100879965205372,
            0.97332682778991096, 0.98333625388462596, 0.99101337147674432,
            0.99634011677195528, 0.99930504173577214
        }};

        inline constexpr std::array<double, 32> kWeights = {{
            0.04869095700913972, 0.04857546744129617, 0.04834476223480296,
            0.04799399306119873, 0.04754053634182810, 0.04696818281621017,
            0.04628427260225052, 0.04549051395354865, 0.04458975198667646,
            0.04358438723094860, 0.04247673042249869, 0.04126929292348930,
            0.03996564116617411, 0.03856971190749052, 0.03708459440549989,
            0.03551493579369814, 0.03386448403584631, 0.03213766176091653,
            0.03033858995124842, 0.02847186956396420, 0.02654155115373869,
            0.02455182230400890, 0.02250676695539950, 0.02041081169392176,
            0.01826846790897158, 0.01608436367056735, 0.01386254088856210,
            0.01160768608166252, 0.00932455454895537, 0.00701669807153073,
            0.00468856751676779, 0.00234406987817898
        }};

        // Integrate f on [a, b].  Precondition: f is C-infinity on [a, b].
        template <typename Func>
        [[nodiscard]] double Integrate(Func&& f, double a, double b) noexcept(noexcept(f(0.0))) {
            const double half = 0.5 * (b - a);
            const double mid  = 0.5 * (b + a);
            double result = 0.0;
            for (int i = 0; i < 32; ++i)
                result += kWeights[i] * (f(mid + half * kNodes[i])
                                    + f(mid - half * kNodes[i]));
            return result * half;
        }

        } // namespace gl64
                
        // Case A: eps in (0, 1] — smooth integrand, direct GL.
        [[nodiscard]] inline double Q5Smooth(double eps) noexcept {
            return gl64::Integrate(
                [eps](double u) noexcept -> double {
                    const double v = 1.0 - eps + 3.0 * eps * u * u;
                    return (v > 0.0) ? std::pow(v, 2.5) : 0.0;
                }, 0.0, 1.0);
        }

        // Case B: eps in (1, +inf) — branch point at u*; substitution regularises.
        [[nodiscard]] inline double Q5Singular(double eps) noexcept {
            const double u_star = std::sqrt((eps - 1.0) / (3.0 * eps));
            const double one_mu = 1.0 - u_star;

            // Prefactor: 2*(3*eps)^(5/2)*(1-u*)^(7/2)
            //   (5/2) from f^(5/2) factorisation, +1 from Jacobian => (7/2) on (1-u*)
            const double pre = 2.0 * std::pow(3.0 * eps, 2.5) * std::pow(one_mu, 3.5);

            // Regularised integrand after u = u_star + (1-u_star) s^2.
            // The Jacobian contributes one additional power of s.
            return pre * gl64::Integrate(
                [u_star, one_mu](double s) noexcept -> double {
                    const double bracket = 2.0 * u_star + one_mu * s * s;
                    return s * s * s * s * s * s * std::pow(bracket, 2.5);
                }, 0.0, 1.0);
        }

        // Public entry point inside gpes::detail.
        // @param eps_dd  a_dd / a_ss — must be positive
        // @throws std::domain_error on out-of-range input
        [[nodiscard]] inline double ComputeQ5(double eps_dd) {
            if (eps_dd <= 0.0)
                throw std::domain_error(
                    "gpes::detail::ComputeQ5: eps_dd must be positive, got "
                    + std::to_string(eps_dd));

            return (eps_dd <= 1.0) ? Q5Smooth(eps_dd) : Q5Singular(eps_dd);
        }


    } // namespace detail



    inline void add_sinusoidal_modulation(
        WaveFunction<Dimension::One>& psi,
        double amplitude,
        double wavelength,
        double phase = 0.0
    ) {
        detail::validate_modulation_inputs(
            amplitude,
            wavelength,
            phase,
            "add_sinusoidal_modulation(1D)"
        );

        if (!psi.grid()) {
            throw std::invalid_argument("add_sinusoidal_modulation(1D): wavefunction has no grid");
        }

        const double initial_norm = detail::current_norm(psi);
        if (initial_norm <= 0.0) {
            throw std::invalid_argument("add_sinusoidal_modulation(1D): wavefunction norm must be positive");
        }

        const double k = 2.0 * M_PI / wavelength;
        const double x0 = psi.grid()->start();
        const double dx = psi.grid()->step();

        for (Eigen::Index i = 0; i < psi.size(); ++i) {
            const double x = x0 + static_cast<double>(i) * dx;
            const double envelope = 1.0 + amplitude * std::sin(k * x + phase);
            psi.vec()(i) *= envelope;
        }

        const double modulated_norm = detail::current_norm(psi);
        if (modulated_norm <= 0.0) {
            throw std::runtime_error("add_sinusoidal_modulation(1D): modulation collapsed the wavefunction norm");
        }

        psi.vec() *= std::sqrt(initial_norm / modulated_norm);
    }

    inline void add_sinusoidal_modulation(
        WaveFunction<Dimension::Two>& psi,
        double amplitude,
        double wavelength,
        double phase = 0.0,
        ModulationAxis axis = ModulationAxis::X
    ) {
        detail::validate_modulation_inputs(
            amplitude,
            wavelength,
            phase,
            "add_sinusoidal_modulation(2D)"
        );

        if (!psi.grid()) {
            throw std::invalid_argument("add_sinusoidal_modulation(2D): wavefunction has no grid");
        }

        const double initial_norm = detail::current_norm(psi);
        if (initial_norm <= 0.0) {
            throw std::invalid_argument("add_sinusoidal_modulation(2D): wavefunction norm must be positive");
        }

        const auto& grid = *psi.grid();
        const double k = 2.0 * M_PI / wavelength;
        const int size_x = static_cast<int>(grid.size_x());
        const int size_y = static_cast<int>(grid.size_y());

        for (int i = 0; i < size_x; ++i) {
            const double x = grid.start_pos_x() + static_cast<double>(i) * grid.step_x();
            for (int j = 0; j < size_y; ++j) {
                const double y = grid.start_pos_y() + static_cast<double>(j) * grid.step_y();
                const double coordinate = axis == ModulationAxis::X ? x : y;
                const double envelope = 1.0 + amplitude * std::sin(k * coordinate + phase);
                psi.vec()(i * size_y + j) *= envelope;
            }
        }

        const double modulated_norm = detail::current_norm(psi);
        if (modulated_norm <= 0.0) {
            throw std::runtime_error("add_sinusoidal_modulation(2D): modulation collapsed the wavefunction norm");
        }

        psi *= std::sqrt(initial_norm / modulated_norm);
    }

    inline void add_complex_random_noise(
        WaveFunction<Dimension::Two>& psi,
        double amplitude,
        unsigned int seed = std::random_device{}()
    ) {
        detail::validate_noise_inputs(
            amplitude,
            seed,
            "add_complex_random_noise(2D)"
        );

        if (!psi.grid()) {
            throw std::invalid_argument("add_complex_random_noise(2D): wavefunction has no grid");
        }

        const double initial_norm = detail::current_norm(psi);
        if (initial_norm <= 0.0) {
            throw std::invalid_argument("add_complex_random_noise(2D): wavefunction norm must be positive");
        }

        const auto& grid = *psi.grid();
        const int size_x = static_cast<int>(grid.size_x());
        const int size_y = static_cast<int>(grid.size_y());

        std::mt19937 rng(seed);
        std::uniform_real_distribution<double> dist(-amplitude, amplitude);

        for (int i = 0; i < size_x; ++i) {
            for (int j = 0; j < size_y; ++j) {
                const double nr = dist(rng);
                const double ni = dist(rng);
                psi.vec()(i * size_y + j) *= std::complex<double>(1.0 + nr, ni);
            }
        }

        const double noisy_norm = detail::current_norm(psi);
        if (noisy_norm <= 0.0) {
            throw std::runtime_error("add_complex_random_noise(2D): noise collapsed the wavefunction norm");
        }

        psi *= std::sqrt(initial_norm / noisy_norm);
    }

    inline void seed_supersolid_state(
        WaveFunction<Dimension::Two>& psi,
        SeedType seed_type,
        double amplitude,
        double wavelength = 0.0,
        double phase = 0.0,
        ModulationAxis axis = ModulationAxis::X,
        unsigned int noise_seed = std::random_device{}()
    ) {
        switch (seed_type) {
            case SeedType::Sinusoidal:
                add_sinusoidal_modulation(psi, amplitude, wavelength, phase, axis);
                break;

            case SeedType::RandomNoise:
                add_complex_random_noise(psi, amplitude, noise_seed);
                break;

            default:
                throw std::invalid_argument("seed_supersolid_state(2D): unknown seed type");
        }
    }

    // Function calculate the value of interactions constants for scattering, 
    // dipole-dipole interaction, LHY correction in 1D case.
    template <>
    inline void calc_inter_consts<Dimension::One>(PhysConfig& cnf, double l_perp){

        cnf.l_ = l_perp;
        
        // There are also other possibilities to compute V_dd
        // dependently on the simulation conditions 
        // for now assume it has the form as bellow
        // double cosTheta = 0; // Theta = 90 deg
        // _V_dd = 0.375 * a_dd / ( std::pow(_l_perp,3));
        // double cosTheta = 1; // Theta = 0 deg
        // _V_dd = 1.5 * a_dd / std::pow(_l_perp,3); 
        double V_dd = -0.75 * cnf.a_dd / ( std::pow(l_perp,3.)); 

        double C = 1.4603; // riemann -zeta(1/2)
        double g_scat =  2 * cnf.a_s / ((l_perp*l_perp) * (1 - C*(cnf.a_s/l_perp))) + 8./ 3 * V_dd;
        double g_lhy = (256. / (15 * M_PI) ) * std::pow(cnf.a_s, 2.5) / std::pow(l_perp, 3.) * (1 + 1.5 * std::pow((cnf.a_dd / cnf.a_s ), 2.));

        cnf.g_scat = g_scat;
        cnf.g_lhy = g_lhy;
        cnf.V_dd = V_dd;
    }


    inline double compute_Q5(double eps_dd) {
        if (eps_dd <= 1.0)
            return detail::Q5Smooth(eps_dd);   // path A — no singularity
        else
            return detail::Q5Singular(eps_dd); // path B — has branch point at u*
    }

    // Function calculate the value of interactions constants for scattering, 
    // dipole-dipole interaction, LHY correction in 2D case.    
    template <>
    inline void calc_inter_consts<Dimension::Two>(PhysConfig& cnf, double l_z) {
        if (cnf.a_s  <= 0.0) throw std::invalid_argument(
            "calc_inter_consts<2D>: a_s must be positive, got "      + std::to_string(cnf.a_s));
        if (cnf.a_dd <  0.0) throw std::invalid_argument(
            "calc_inter_consts<2D>: a_dd must be non-negative, got " + std::to_string(cnf.a_dd));
        if (l_z      <= 0.0) throw std::invalid_argument(
            "calc_inter_consts<2D>: l_z must be positive, got "      + std::to_string(l_z));

        cnf.l_ = l_z;

        const double eps_dd = cnf.a_dd / cnf.a_s;
        const double Q5     = gpes::detail::ComputeQ5(eps_dd);

        cnf.g_scat = std::sqrt(8.0 * M_PI) * cnf.a_s / l_z;
        cnf.V_dd   = std::sqrt(8.0 * M_PI) * cnf.a_dd / l_z;
        cnf.g_lhy  = (128.0 / (3.0 * std::pow(M_PI, 0.25))) //1.25
                    * std::sqrt(0.4) // 2.5
                    * std::pow(cnf.a_s, 2.5)
                    / std::pow(l_z,     1.5)
                    * Q5;

        //(1 + 1.5 * std::pow((cnf.a_dd / cnf.a_s ), 2.));            

        // GPES_LOG(kLogCategory, Info,
        //     "calc_inter_consts<2D>: a_s={:.8f} a_dd={:.8f} l_z={:.8f} "
        //     "eps_dd={:.6f} Q5={:.8f} (approx err={:.1f}%) "
        //     "g_scat={:.8f} V_dd={:.8f} g_lhy={:.10f}",
        //     cnf.a_s, cnf.a_dd, l_z, eps_dd, Q5,
        //     100.0 * ((1.0 + 1.5 * eps_dd * eps_dd) / Q5 - 1.0),
        //     cnf.g_scat, cnf.V_dd, cnf.g_lhy);
    }


    namespace solvers::detail {
        inline void validate_phys_config(const PhysConfig& cfg) {
            if (!std::isfinite(cfg.a_s) || cfg.a_s <= 0.0) {
                throw std::invalid_argument("PhysConfig.a_s must be positive and finite");
            }
            if (!std::isfinite(cfg.a_dd) || cfg.a_dd < 0.0) {
                throw std::invalid_argument("PhysConfig.a_dd must be non-negative and finite");
            }
            if (!std::isfinite(cfg.num_of_prt) || cfg.num_of_prt <= 0.0) {
                throw std::invalid_argument("PhysConfig.num_of_prt must be positive and finite");
            }
        }

        inline void validate_inter_consts(const PhysConfig& cfg, const char* context) {
            if (!std::isfinite(cfg.l_) || cfg.l_ == 0.0) {
                throw std::invalid_argument(std::string(context) + ": l_ must be non-zero and finite");
            }
            if (!std::isfinite(cfg.g_scat) || cfg.g_scat == 0.0) {
                throw std::invalid_argument(std::string(context) + ": g_scat must be non-zero and finite");
            }
            if (!std::isfinite(cfg.g_lhy) || cfg.g_lhy == 0.0) {
                throw std::invalid_argument(std::string(context) + ": g_lhy must be non-zero and finite");
            }
            if (!std::isfinite(cfg.V_dd) || cfg.V_dd == 0.0) {
                throw std::invalid_argument(std::string(context) + ": V_dd must be non-zero and finite");
            }
        }

    } // namespace detail


    namespace operators {
        template <typename GridType>
        class KineticOperator {

        };
    }
}
