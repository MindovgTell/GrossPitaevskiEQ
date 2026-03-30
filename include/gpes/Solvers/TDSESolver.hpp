#pragma once

#include <cmath>
#include <filesystem>
#include <limits>
#include <memory>
#include <stdexcept>

#include "Core/wavefunction_io.hpp"
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
                    trait::is_tstep_algo_v<std::decay_t<Algorythm>, WFT>
            >>    
    class TDSESolver {
    public:
        inline static constexpr Dimension Dim = D;
        using GridType = GridT;
        using WFType = WFT;
        using ShrdPtrGrid = std::shared_ptr<const Grid<Dim>>;
        inline static const gpes::log::LogCategory kLogCategory{"TDSESolver"};
        inline static constexpr double kEnergyTol = 1e-8;
        inline static constexpr std::size_t kMaxGroundSteps = 40000;

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
            run_impl(nullptr, 0);
        }

        void run(const std::filesystem::path& snapshot_path, std::size_t snapshot_interval = 4000u) {
            if (snapshot_interval == 0u) {
                throw std::invalid_argument("TDSESolver: snapshot_interval must be greater than zero");
            }
            run_impl(&snapshot_path, snapshot_interval);
        }

    private:
        void run_impl(const std::filesystem::path* snapshot_path, std::size_t snapshot_interval) {
            if (!std::isfinite(sim_config_.dt) || sim_config_.dt <= 0.0) {
                throw std::invalid_argument("TDSESolver: dt must be positive and finite");
            }
            if (!std::isfinite(sim_config_.duration) || sim_config_.duration <= 0.0) {
                throw std::invalid_argument("TDSESolver: duration must be positive and finite");
            }
            sim_config_.num_of_steps = sim_config_.duration / sim_config_.dt;
            prepare_snapshot_path(snapshot_path);

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
                GPES_LOG(kLogCategory, Log, "TDSE ground state enabled: energy_tol={}, max_steps={}", kEnergyTol, kMaxGroundSteps);
                double prev_energy = std::numeric_limits<double>::quiet_NaN();
                std::size_t steps_completed = 0;
                for (std::size_t step = 0; step < kMaxGroundSteps; ++step) {
                    algo_.step(Psi_);
                    steps_completed = step + 1;
                    maybe_save_snapshot(snapshot_path, snapshot_interval, steps_completed);

                    if (step > 0 && (step % 1000u) == 0u) {
                        GPES_LOG(kLogCategory, Log, "TDSE progress: step={}", step);
                    }

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
                save_snapshot(snapshot_path);
                GPES_LOG(kLogCategory, Log, "TDSE run finish: steps_completed={}", steps_completed);
            } else {
                // Logging info before simulation start
                GPES_LOG(kLogCategory, Log, "TDSE run start: steps={}, dt={}, duration={}",
                    sim_config_.num_of_steps, sim_config_.dt, sim_config_.duration);

                for (std::size_t step = 0; step < sim_config_.num_of_steps; ++step) {
                    algo_.step(Psi_);
                    maybe_save_snapshot(snapshot_path, snapshot_interval, step + 1);
                    if (step > 0 && (step % 1000u) == 0u) {
                        GPES_LOG(kLogCategory, Log, "TDSE progress: step={}", step);
                    }
                }
                save_snapshot(snapshot_path);
                GPES_LOG(kLogCategory, Log, "TDSE run finish: steps_completed={}", sim_config_.num_of_steps);
            }
        }

        void prepare_snapshot_path(const std::filesystem::path* snapshot_path) const {
            if (snapshot_path == nullptr) {
                return;
            }
            if (snapshot_path->has_parent_path()) {
                std::filesystem::create_directories(snapshot_path->parent_path());
            }
        }

        void maybe_save_snapshot(
            const std::filesystem::path* snapshot_path,
            std::size_t snapshot_interval,
            std::size_t completed_steps
        ) const {
            if (snapshot_path == nullptr || snapshot_interval == 0u) {
                return;
            }
            if ((completed_steps % snapshot_interval) != 0u) {
                return;
            }
            save_snapshot(snapshot_path);
        }

        void save_snapshot(const std::filesystem::path* snapshot_path) const {
            if (snapshot_path == nullptr) {
                return;
            }
            io::save_wavefunction_csv(Psi_, *snapshot_path);
        }

    public:

        const WFType& wavefunction() const { return Psi_; }
        WFType& wavefunction() { return Psi_; }

    private:
        WFType& Psi_;
        ShrdPtrGrid gridptr_;
        Algorythm& algo_;
        SimConfig sim_config_;
    };
};
