#pragma once

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

#include "grid/grid.hpp"
#include "wavefunction/wavefunction.hpp"
#include "Core/traits.hpp"
#include "Log/Log.hpp"

namespace gpes {

    template <Dimension D, typename Algorythm,
                typename GridT = Grid<D>,
                typename WFT = WaveFunction<D>,
                typename = std::enable_if_t<
                    trait::is_grid_v<GridT> && trait::is_wavefunction_v<WFT> && 
                    trait::is_tstep_algo_v<Algorythm, WFT>
            >>    
    class TDSESolver {
    public:
        inline static constexpr Dimension Dim = D;
        using GridType = GridT;
        using WFType = WFT;
        using ShrdPtrGrid = std::shared_ptr<const Grid<Dim>>;
        inline static const gpes::log::LogCategory kLogCategory{"TDSESolver"};
        inline static constexpr double kEnergyTol = 1e-3;

        TDSESolver(ShrdPtrGrid grid, Algorythm& algo, WFType& psi, SimConfig simcnfg)
            : Psi_(psi),
              gridptr_(std::move(grid)),
              algo_(algo),
              sim_config_(simcnfg) {
            if (!gridptr_) {
                throw std::runtime_error("TDSESolver requires a valid grid");
            }
        }

        void run() {
            sim_config_.num_of_steps = sim_config_.duration / sim_config_.dt;
            // Logging info about given solver algorithm
            if constexpr (requires { typename Algorythm::Tag; }) {
                using AlgoTag = typename Algorythm::Tag;
                GPES_LOG(kLogCategory, Log,
                    "TDSE config: dim={}, algo={}, implicit={}, requires_fft={}",
                    (Dim == Dimension::One ? "1D" : "2D"),
                    algo_traits<AlgoTag>::name(),
                    algo_traits<AlgoTag>::implicit ? "true" : "false",
                    algo_traits<AlgoTag>::requires_fft ? "true" : "false");
            } else {
                GPES_LOG(kLogCategory, Log, "TDSE config: dim={}, algo=unknown",
                    (Dim == Dimension::One ? "1D" : "2D"));
            }


            if (sim_config_.ground_state) {
                GPES_LOG(kLogCategory, Log, "TDSE ground state enabled: energy_tol={}", kEnergyTol);
            }
            // Logging info before simulation start
            GPES_LOG(kLogCategory, Log, "TDSE run start: steps={}, dt={}, duration={}",
                sim_config_.num_of_steps, sim_config_.dt, sim_config_.duration);


            double prev_energy = std::numeric_limits<double>::quiet_NaN();
            // main simulation loop
            for (size_t step = 0; step < sim_config_.num_of_steps; ++step) {
                algo_.step(Psi_);

                if (step > 0 && (step % 1000u) == 0u) {
                    GPES_LOG(kLogCategory, Log, "TDSE progress: step={}", step);
                }
                if (sim_config_.ground_state) {
                    const double current_energy = algo_.last_energy();
                    if (std::isfinite(prev_energy) &&
                        std::abs(current_energy - prev_energy) < kEnergyTol) {
                        GPES_LOG(kLogCategory, Log,
                            "TDSE ground state reached: step={}, energy={}, delta_energy={}",
                            step, current_energy, std::abs(current_energy - prev_energy));
                        break;
                    }
                    prev_energy = current_energy;
                }
            }

            GPES_LOG(kLogCategory, Log, "TDSE run finish: steps_completed={}", sim_config_.num_of_steps);
        }

        const WFType& wavefunction() const { return Psi_; }
        WFType& wavefunction() { return Psi_; }

    private:
        WFType& Psi_;
        ShrdPtrGrid gridptr_;
        Algorythm& algo_;
        SimConfig sim_config_;
    };
};
