#include "cli/cli.hpp"
#include "gpes.hpp"
#include "Core/utility.hpp"
#if GPES_HAS_PLOTTING
#include "Core/plot.hpp"
#endif
#include "Log/Log.hpp"
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iomanip>
#include <iostream>


int main(int argc, char **argv) {
    try {
        auto parsed = parse_cli(argc, argv);
        auto phys = parsed.phys;
        auto sim = parsed.sim;
        sim.ground_state = true;
        const auto &cfg = parsed.cli;
        print_parsed_params(argc, argv);

        const std::filesystem::path run_dir = cfg.wavefunction_out;
        const std::filesystem::path wavefunction_dir = run_dir / "wavefunction";
        const std::filesystem::path bdg_dir = run_dir / "bdg";

        std::filesystem::create_directories(wavefunction_dir);
        std::filesystem::create_directories(bdg_dir);
#if GPES_HAS_PLOTTING
        const std::filesystem::path plot_dir = run_dir / "plot";
        const std::filesystem::path excited_plot_dir = plot_dir / "excited_states";
        std::filesystem::create_directories(plot_dir);
        std::filesystem::create_directories(excited_plot_dir);
#endif

        std::cout << "Run directory: " << run_dir << "\n";

        constexpr gpes::Dimension Dim = gpes::Dimension::Two;
        auto grid_ptr = std::make_shared<gpes::Grid<Dim>>(
            128, 128, -25, -25, cfg.w_x, cfg.w_y, cfg.w_z
        );

        gpes::WaveFunction<Dim> Psi(grid_ptr);
        Psi.set_state_Gauss(0.0, 0.0, 2.5, 2.5, phys.num_of_prt);

        gpes::solvers::SplitStep<Dim> stepper1(grid_ptr, Psi, phys, sim);

        gpes::TDSESolver<Dim, decltype(stepper1)> solver1(grid_ptr, stepper1, Psi, sim);
        const std::filesystem::path wavefunction_out_path = wavefunction_dir / "wavefunction.csv";

#if GPES_HAS_PLOTTING
        gpes::plot::save_heatmap(Psi, plot_dir / "initial_wavefunction_heatmap.png", "Initial wavefunction");
#endif

        solver1.run(wavefunction_out_path, 4000u);

        const auto& energy_history = stepper1.energies();
        const std::filesystem::path energy_history_path =
            gpes::io::save_state_energy_history_csv(
                energy_history,
                wavefunction_out_path,
                "energy_history.csv"
            );
#if GPES_HAS_PLOTTING
        gpes::plot::save_energy_plot(energy_history, plot_dir / "energy_history.png");

        gpes::plot::save_heatmap(Psi, plot_dir / "wavefunction_heatmap.png", "Wavefunction");
        gpes::plot::save_slices(Psi, plot_dir / "wavefunction_slices.png", "Wavefunction");

        auto temp = Psi.fft_shifted();
        gpes::WaveFunction<Dim> kPsi(grid_ptr, temp);
        gpes::plot::save_heatmap(kPsi, plot_dir / "momentum_wavefunction_heatmap.png", "Momentum wavefunction");
        gpes::plot::save_slices(kPsi, plot_dir / "momentum_wavefunction_slices.png", "Momentum wavefunction");
#endif

        // 128x128 shift-invert test. Increase num_modes only after this
        // configuration converges reliably on the target machine.
        gpes::bdg::BdGConfig bdg_config;
        bdg_config.num_modes = 20;
        bdg_config.ncv = 60;
        bdg_config.max_iter = 1500;
        bdg_config.tol = 1e-6;
        bdg_config.use_shift_invert = true;
        bdg_config.shift = {0.05, 0.0};
        bdg_config.shift_invert_linear_tol = 1e-9;
        bdg_config.shift_invert_linear_max_iter = 1500;
        bdg_config.shift_invert_linear_restart = 60;
        bdg_config.shift_invert_dense_max_dim = 0;
        bdg_config.use_fourier_preconditioner = true;
        bdg_config.preconditioner_include_dipolar = true;
        bdg_config.preconditioner_regularization = 1e-8;

        gpes::bdg::BdGSolver2D bdg_solver(grid_ptr, phys, bdg_config);
        const auto bdg_result = bdg_solver.solve(Psi);

        const std::filesystem::path bdg_spectrum_path =
            gpes::bdg::save_bdg_result(bdg_result, grid_ptr, bdg_dir);

        std::cout << "BdG chemical potential: "
                  << std::setprecision(16)
                  << bdg_solver.chemical_potential()
                  << "\n";

        std::cout << "BdG solver diagnostics: nconv="
                  << bdg_result.nconv
                  << ", niter=" << bdg_result.niter
                  << ", nops=" << bdg_result.nops
                  << "\n";
        std::cout << "BdG shift-invert diagnostics: shift="
                  << bdg_result.shift
                  << ", linear_iterations=" << bdg_result.linear_iterations
                  << ", linear_failures=" << bdg_result.linear_failures
                  << ", max_linear_residual=" << bdg_result.max_linear_residual
                  << ", max_eigen_residual=" << bdg_result.max_eigen_residual
                  << ", fourier_preconditioner="
                  << (bdg_result.used_fourier_preconditioner ? "yes" : "no")
                  << ", dense_direct="
                  << (bdg_result.used_dense_direct ? "yes" : "no")
                  << "\n";

        std::cout << "BdG spectrum saved to " << bdg_spectrum_path << "\n";

#if GPES_HAS_PLOTTING
        gpes::plot::save_excitation_spectrum_plot(
            bdg_result,
            plot_dir / "excitation_spectrum.png"
        );
#endif


        return 0;
    } catch (const std::exception& e) {
        std::cout << "Fatal error: " << e.what() << "\n";
        return 1;
    }
}
