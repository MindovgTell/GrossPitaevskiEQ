#pragma once

// #include "ReadSimParams.hpp"
// #include "visualitsation.hpp"
// #include "SaveSimResults.hpp"

#include <cmath>
#include <cstddef>

#include "definitions.hpp"

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

    // Function calculate the value of interactions constants for scattering, 
    // dipole-dipole interaction, LHY correction in 2D case.
    template <>
    inline void calc_inter_consts<Dimension::Two>(PhysConfig& cnf, double l_z){

        cnf.l_ = l_z;

        double g_scat = std::sqrt(8. * M_PI) * cnf.a_s / l_z;
        double V_dd = std::sqrt(8. * M_PI) * cnf.a_dd / l_z;
        double g_lhy = (128. / ( 3 * std::pow(M_PI, 0.25) ) ) * std::sqrt(0.4) * std::pow(cnf.a_s, 2.5) / std::pow(l_z,1.5) * (1 + 1.5 * std::pow((cnf.a_dd / cnf.a_s ), 2.));        

        cnf.g_scat = g_scat;
        cnf.g_lhy = g_lhy;
        cnf.V_dd = V_dd;
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
}


