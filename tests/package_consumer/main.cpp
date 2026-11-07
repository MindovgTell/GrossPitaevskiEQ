#include "gpes.hpp"
#ifdef GPES_TEST_PLOT
#include "Core/plot.hpp"
#endif

#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>

int main() {
    // Exercise installed headers, compiled logging, FFTW, and Spectra linkage.
    // This checks the package interface, not physical convergence on this grid.
    constexpr auto dim = gpes::Dimension::Two;
    auto grid = std::make_shared<gpes::Grid<dim>>(8, 8, -4.0, -4.0, 1.0, 1.0, 8.0);
    gpes::WaveFunction<dim> psi(grid);
    psi.set_state_Gauss(0.0, 0.0, 1.0, 1.0, 1);
    const gpes::PhysConfig phys{0.01, 0.005, 1.0};
    const gpes::SimConfig sim{};
    gpes::solvers::SplitStep<dim> stepper(grid, psi, phys, sim);
    stepper.step(psi);
    const auto diagnostics = stepper.diagnostics(psi);
    const double norm = psi.vec().squaredNorm() * grid->step_x() * grid->step_y();
    if (!std::isfinite(diagnostics.energy) || std::abs(norm - 1.0) > 1e-10) {
        std::cerr << "Installed split-step calculation failed\n";
        return 1;
    }

    gpes::bdg::BdGConfig cfg;
    cfg.num_modes = 2;
    cfg.ncv = 10;
    cfg.use_shift_invert = true;
    cfg.shift_invert_dense_max_dim = 128;
    gpes::bdg::BdGSolver2D solver(grid, phys, cfg);
    const auto result = solver.solve(psi, diagnostics.chemical_potential);
    if (result.nconv < cfg.num_modes || result.modes.empty() ||
        !std::isfinite(result.modes.front().omega.real())) {
        std::cerr << "Installed BdG calculation failed\n";
        return 1;
    }
#ifdef GPES_TEST_PLOT
    const std::filesystem::path image_path = "package-energy.png";
    std::filesystem::remove(image_path);
    gpes::plot::save_energy_plot({diagnostics.energy, diagnostics.energy}, image_path);
    if (!std::filesystem::exists(image_path) || std::filesystem::file_size(image_path) == 0) {
        std::cerr << "Installed plotting component failed\n";
        return 1;
    }
#endif
    std::cout << "Installed GPES: energy=" << diagnostics.energy
              << ", norm=" << norm << ", BdG modes=" << result.nconv << '\n';
}
