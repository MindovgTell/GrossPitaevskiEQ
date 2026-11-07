#pragma once

#include <cmath>
#include <algorithm>
#include <filesystem>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include "Core/utility.hpp"
#include "Core/wavefunction_io.hpp"
#include "grid/grid.hpp"
#include "wavefunction/wavefunction.hpp"
#include "Core/traits.hpp"
#include "Log/Log.hpp"

namespace gpes {

    struct GroundStateResult {
        bool converged = false;
        std::size_t steps_completed = 0;
        std::size_t convergence_checks = 0;
        double energy = std::numeric_limits<double>::quiet_NaN();
        double relative_energy_change = std::numeric_limits<double>::infinity();
        double chemical_potential = std::numeric_limits<double>::quiet_NaN();
        double stationary_residual = std::numeric_limits<double>::infinity();
        double final_time_step = std::numeric_limits<double>::quiet_NaN();
        double absolute_energy_change = std::numeric_limits<double>::infinity();
        std::vector<std::size_t> diagnostic_steps;
    };

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
                validate_ground_state_config();
                GPES_LOG(kLogCategory, Log,
                    "Ground-state search: relative_energy_tol={}, absolute_energy_tol={}, residual_tol={}, check_interval={}, max_steps={}",
                    sim_config_.relative_energy_tolerance,
                    sim_config_.absolute_energy_tolerance,
                    sim_config_.stationary_residual_tolerance,
                    sim_config_.convergence_check_interval,
                    sim_config_.max_ground_steps);
                double prev_energy = std::numeric_limits<double>::quiet_NaN();
                std::size_t consecutive_converged = 0;
                std::size_t energy_plateau_checks = 0;
                ground_result_ = {};
                ground_result_.final_time_step = sim_config_.dt;
                for (std::size_t step = 0; step < sim_config_.max_ground_steps; ++step) {
                    algo_.step(Psi_);
                    ground_result_.steps_completed = step + 1;
                    maybe_save_snapshot(snapshot_path, snapshot_interval, ground_result_.steps_completed);

                    if ((ground_result_.steps_completed % sim_config_.convergence_check_interval) != 0u) {
                        continue;
                    }

                    update_ground_state_diagnostics();
                    ++ground_result_.convergence_checks;
                    update_energy_change(prev_energy);
                    const bool energy_plateau = ground_energy_converged(
                        ground_result_.energy, prev_energy,
                        sim_config_.relative_energy_tolerance, sim_config_.absolute_energy_tolerance);
                    const bool check_converged =
                        energy_plateau &&
                        ground_result_.stationary_residual <= sim_config_.stationary_residual_tolerance;
                    consecutive_converged = check_converged ? consecutive_converged + 1 : 0;
                    energy_plateau_checks = energy_plateau ? energy_plateau_checks + 1 : 0;

                    GPES_LOG(kLogCategory, Log,
                        "Ground-state check: step={}, energy={}, relative_delta={}, residual={}, consecutive={}",
                        ground_result_.steps_completed,
                        ground_result_.energy,
                        ground_result_.relative_energy_change,
                        ground_result_.stationary_residual,
                        consecutive_converged);

                    if (consecutive_converged >= sim_config_.convergence_checks_required) {
                        ground_result_.converged = true;
                        GPES_LOG(kLogCategory, Log,
                            "Ground state converged at step={}", ground_result_.steps_completed);
                        break;
                    }
                    if (energy_plateau_checks >= 3 &&
                        ground_result_.stationary_residual > sim_config_.stationary_residual_tolerance) {
                        if constexpr (requires { algo_.set_time_step(1.0); algo_.time_step(); }) {
                            const double old_dt = algo_.time_step();
                            const double new_dt = std::max(sim_config_.minimum_ground_dt, 0.5 * old_dt);
                            if (new_dt < old_dt * (1.0 - 1e-12)) {
                                algo_.set_time_step(new_dt);
                                ground_result_.final_time_step = new_dt;
                                prev_energy = std::numeric_limits<double>::quiet_NaN();
                                energy_plateau_checks = 0;
                                consecutive_converged = 0;
                                GPES_LOG(kLogCategory, Log,
                                    "Refining imaginary-time step: {} -> {}", old_dt, new_dt);
                                continue;
                            }
                        }
                    }
                    prev_energy = ground_result_.energy;
                }
                // The final checkpoint and diagnostics must describe the same state,
                // including runs ending between regular convergence checks.
                if (ground_result_.steps_completed % sim_config_.convergence_check_interval != 0) {
                    update_ground_state_diagnostics();
                    update_energy_change(prev_energy);
                }
                if (!ground_result_.converged) {
                    GPES_LOG(kLogCategory, Warning,
                        "Ground state did not converge within {} steps; final residual={}",
                        ground_result_.steps_completed, ground_result_.stationary_residual);
                }
                save_snapshot(snapshot_path);
                GPES_LOG(kLogCategory, Log, "TDSE run finish: steps_completed={}", ground_result_.steps_completed);
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

        void validate_ground_state_config() const {
            if (sim_config_.max_ground_steps == 0 ||
                sim_config_.convergence_check_interval == 0 ||
                sim_config_.convergence_checks_required == 0 ||
                !(sim_config_.relative_energy_tolerance > 0.0) ||
                !std::isfinite(sim_config_.relative_energy_tolerance) ||
                !(sim_config_.absolute_energy_tolerance >= 0.0) ||
                !std::isfinite(sim_config_.absolute_energy_tolerance) ||
                !(sim_config_.stationary_residual_tolerance > 0.0) ||
                !std::isfinite(sim_config_.stationary_residual_tolerance) ||
                !(sim_config_.minimum_ground_dt > 0.0) ||
                sim_config_.minimum_ground_dt > sim_config_.dt) {
                throw std::invalid_argument("TDSESolver: invalid ground-state convergence controls");
            }
        }

        void update_energy_change(double previous) {
            ground_result_.absolute_energy_change = std::numeric_limits<double>::infinity();
            ground_result_.relative_energy_change = std::numeric_limits<double>::infinity();
            if (std::isfinite(previous) && std::isfinite(ground_result_.energy)) {
                const double delta = std::abs(ground_result_.energy - previous);
                const double scale = std::max(std::abs(ground_result_.energy), std::abs(previous));
                ground_result_.absolute_energy_change = delta;
                ground_result_.relative_energy_change = scale > 0.0 ? delta / scale : 0.0;
            }
        }

        void update_ground_state_diagnostics() {
            if constexpr (requires { algo_.diagnostics(Psi_); }) {
                const auto diagnostics = algo_.diagnostics(Psi_);
                ground_result_.energy = diagnostics.energy;
                ground_result_.chemical_potential = diagnostics.chemical_potential;
                ground_result_.stationary_residual = diagnostics.relative_residual;
            } else {
                ground_result_.energy = algo_.last_energy();
                ground_result_.stationary_residual = 0.0;
            }
            ground_result_.diagnostic_steps.push_back(ground_result_.steps_completed);
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
        const GroundStateResult& ground_state_result() const { return ground_result_; }

    private:
        WFType& Psi_;
        ShrdPtrGrid gridptr_;
        Algorythm& algo_;
        SimConfig sim_config_;
        GroundStateResult ground_result_;
    };
};
