#include "Solvers/TDSESolver.hpp"

#include <iostream>

namespace {

constexpr auto Dim = gpes::Dimension::Two;
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

struct ScriptedStepper {
    static constexpr auto Dim = gpes::Dimension::Two;
    struct Diagnostics { double energy, chemical_potential, relative_residual; };
    int step_count = 0;
    double dt = 0.001;
    bool refining = false;
    void step(gpes::WaveFunction<Dim>&) { ++step_count; }
    double last_energy() const { return refining ? 1.0 : step_count; }
    Diagnostics diagnostics(const gpes::WaveFunction<Dim>&) {
        return {last_energy(), 2.0 * last_energy(), refining && step_count < 5 ? 1.0 : 0.0};
    }
    void set_time_step(double value) { dt = value; }
    double time_step() const { return dt; }
};

} // namespace

int main() {
    try {
        require(gpes::ground_energy_converged(1000.000005, 1000.0, 1e-8, 1e-12), "Relative large-energy tolerance failed");
        require(gpes::ground_energy_converged(-1000.000005, -1000.0, 1e-8, 1e-12), "Negative-energy tolerance failed");
        require(!gpes::ground_energy_converged(1.01e-6, 1e-6, 1e-8, 1e-12), "Energy scale was incorrectly clamped to one");
        require(gpes::ground_energy_converged(5e-13, 0.0, 1e-8, 1e-12), "Absolute tolerance near zero failed");
        require(!gpes::ground_energy_converged(0.0, std::numeric_limits<double>::quiet_NaN(), 1e-8, 1e-12), "First diagnostic must not pass");
        require(!gpes::ground_energy_converged(std::numeric_limits<double>::infinity(), 1.0, 1e-8, 1e-12), "Nonfinite energy passed");
        auto grid = std::make_shared<gpes::Grid<Dim>>(8, 8, -4.0, -4.0, 1.0, 1.0, 8.0);
        gpes::WaveFunction<Dim> psi(grid);
        psi.set_state_Gauss(0.0, 0.0, 1.0, 1.0, 1);
        gpes::SimConfig sim;
        sim.max_ground_steps = 5;
        sim.convergence_check_interval = 3;
        ScriptedStepper stepper;
        gpes::TDSESolver<Dim, ScriptedStepper> solver(grid, stepper, psi, sim);
        solver.run();
        const auto& result = solver.ground_state_result();
        require(!result.converged && result.energy == 5.0 && result.chemical_potential == 10.0,
                "Final diagnostics do not describe the final checkpoint");
        require(result.diagnostic_steps == std::vector<std::size_t>{3, 5}, "Diagnostic step labels are wrong");

        sim.max_ground_steps = 6;
        sim.convergence_check_interval = 1;
        sim.convergence_checks_required = 2;
        ScriptedStepper refining;
        refining.refining = true;
        gpes::TDSESolver<Dim, ScriptedStepper> refinement_solver(grid, refining, psi, sim);
        refinement_solver.run();
        require(refining.dt == 0.0005, "Timestep was not refined");
        require(!refinement_solver.ground_state_result().converged, "Stale energy change passed after timestep refinement");

        sim.absolute_energy_tolerance = std::numeric_limits<double>::infinity();
        bool rejected = false;
        try {
            gpes::TDSESolver<Dim, ScriptedStepper> invalid(grid, refining, psi, sim);
            invalid.run();
        } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "Infinite convergence tolerance was accepted");
        std::cout << "Energy convergence and checkpoint diagnostics tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
