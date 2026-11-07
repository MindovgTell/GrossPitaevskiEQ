#include "wavefunction/wavefunction.hpp"

#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

int main() {
    try {
        using Grid = gpes::Grid<gpes::Dimension::Two>;
        using State = gpes::WaveFunction<gpes::Dimension::Two>;
        auto grid = std::make_shared<Grid>(20, 24, -5.0, -6.0, 1.0, 1.25, 5.0);
        State psi(grid);
        psi.set_state_Gauss(0.0, 0.0, 1.2, 1.5, 10);
        const auto original = psi.vec();
        const double norm = original.squaredNorm();
        const double scale = original.cwiseAbs().maxCoeff();
        const std::vector<State::GaussianComponent> peaks{
            {1.0, -1.0, 0.5, 0.8, 0.75},
            {-1.5, 1.25, 0.6, 0.4, 0.3}
        };
        psi.add_gaussian_modulations(peaks);
        Eigen::VectorXcd expected = original;
        for (int i = 0; i < 20; ++i) {
            const double x = grid->start_pos_x() + i * grid->step_x();
            for (int j = 0; j < 24; ++j) {
                const double y = grid->start_pos_y() + j * grid->step_y();
                for (const auto& peak : peaks) {
                    expected(i * 24 + j) += scale * peak[4] * std::exp(
                        -0.5 * (std::pow((x - peak[0]) / peak[2], 2)
                              + std::pow((y - peak[1]) / peak[3], 2)));
                }
            }
        }
        expected *= std::sqrt(norm / expected.squaredNorm());
        if ((psi.vec() - expected).norm() > 1e-13 * expected.norm() ||
            std::abs(psi.vec().squaredNorm() - norm) > 1e-12 * norm) {
            throw std::runtime_error("Gaussian peak sum or normalization is incorrect");
        }
        const auto before_failure = psi.vec();
        bool rejected = false;
        try { psi.add_gaussian_modulations({{0.0, 0.0, 0.0, 1.0, 1.0}}); }
        catch (const std::invalid_argument&) { rejected = true; }
        if (!rejected || (psi.vec() - before_failure).norm() != 0.0) {
            throw std::runtime_error("Invalid peak changed the wavefunction");
        }
        psi.add_gaussian_modulations({});
        if ((psi.vec() - before_failure).norm() != 0.0) {
            throw std::runtime_error("Empty modulation list should preserve the state");
        }
        std::cout << "Gaussian modulation and norm-preservation checks passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
