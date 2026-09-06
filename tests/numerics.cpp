#include "Core/utility.hpp"
#include "Solvers/BogoliubovdeGennes.hpp"
#include "dipoleinteraction/dipoleinteraction.hpp"
#include "grid/grid.hpp"

#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <stdexcept>

#include <Eigen/Core>
#include <fftw3.h>

namespace {

using Complex = std::complex<double>;
constexpr gpes::Dimension Dim = gpes::Dimension::Two;

Eigen::VectorXcd spectral_kinetic(
    const Eigen::VectorXcd& psi,
    const gpes::Grid<Dim>& grid)
{
    const int nx = static_cast<int>(grid.size_x());
    const int ny = static_cast<int>(grid.size_y());
    const int n = nx * ny;
    Eigen::VectorXcd out = psi;
    auto* data = reinterpret_cast<fftw_complex*>(out.data());
    fftw_plan forward = fftw_plan_dft_2d(nx, ny, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan backward = fftw_plan_dft_2d(nx, ny, data, data, FFTW_BACKWARD, FFTW_ESTIMATE);
    if (!forward || !backward) throw std::runtime_error("test FFTW plan creation failed");
    fftw_execute_dft(forward, data, data);
    for (int i = 0; i < nx; ++i) {
        const int qx = i <= nx / 2 ? i : i - nx;
        const double kx = 2.0 * M_PI * qx / (nx * grid.step_x());
        for (int j = 0; j < ny; ++j) {
            const int qy = j <= ny / 2 ? j : j - ny;
            const double ky = 2.0 * M_PI * qy / (ny * grid.step_y());
            out(i * ny + j) *= 0.5 * (kx * kx + ky * ky);
        }
    }
    fftw_execute_dft(backward, data, data);
    out /= static_cast<double>(n);
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
    return out;
}

Eigen::VectorXcd apply_egpe(
    const Eigen::VectorXcd& psi,
    const std::shared_ptr<const gpes::Grid<Dim>>& grid,
    const gpes::PhysConfig& input_phys,
    double mu)
{
    gpes::PhysConfig phys = input_phys;
    const double lz = 1.0 / std::sqrt(grid->omega_z());
    gpes::calc_inter_consts<Dim>(phys, lz);
    gpes::DipolarInteraction<Dim> ddi(grid, lz, phys.V_dd);
    Eigen::VectorXd phi;
    ddi.compute_DDI_term(psi, phi);
    Eigen::VectorXcd out = spectral_kinetic(psi, *grid);
    const int ny = static_cast<int>(grid->size_y());
    for (int i = 0; i < static_cast<int>(grid->size_x()); ++i) {
        for (int j = 0; j < ny; ++j) {
            const int k = i * ny + j;
            const double density = std::norm(psi(k));
            out(k) += (grid->potential()(i, j) - mu
                + phys.g_scat * density + phi(k)
                + phys.g_lhy * std::pow(density, 1.5)) * psi(k);
        }
    }
    return out;
}

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

} // namespace

int main() {
    try {
        std::cout << "GL64 constant integral: "
                  << gpes::detail::gl64::Integrate([](double) { return 1.0; }, 0.0, 1.0)
                  << '\n';
        const double q = 40.0;
        const double form_factor = 2.0 - 3.0 * std::sqrt(M_PI) * q
            * gpes::scaled_erfc_nonnegative(q);
        require(std::isfinite(form_factor), "high-k dipolar kernel is not finite");
        require(std::abs(form_factor + 1.0) < 1e-3, "high-k dipolar kernel has wrong limit");

        for (double epsilon_dd : {0.5, 2.0, 8.843}) {
            constexpr int samples = 200000;
            double midpoint_integral = 0.0;
            for (int i = 0; i < samples; ++i) {
                const double u = (static_cast<double>(i) + 0.5) / samples;
                const double argument = 1.0 - epsilon_dd + 3.0 * epsilon_dd * u * u;
                if (argument > 0.0) midpoint_integral += std::pow(argument, 2.5);
            }
            midpoint_integral /= samples;
            const double q5 = gpes::detail::ComputeQ5(epsilon_dd);
            std::cout << "Q5(" << epsilon_dd << "): quadrature=" << q5
                      << ", midpoint=" << midpoint_integral << '\n';
            require(std::abs(q5 - midpoint_integral) / q5 < 1e-9,
                "Q5 quadrature does not match its defining real-part integral");
        }

        auto mutable_grid = std::make_shared<gpes::Grid<Dim>>(8, 10, -4.0, -5.0, 0.7, 1.2, 6.0);
        std::shared_ptr<const gpes::Grid<Dim>> grid = mutable_grid;
        const int n = static_cast<int>(grid->size_x() * grid->size_y());
        Eigen::VectorXcd psi(n);
        Eigen::VectorXcd perturbation(n);
        for (int i = 0; i < n; ++i) {
            const double x = 0.17 * static_cast<double>(i + 1);
            psi(i) = Complex(0.3 + 0.04 * std::cos(x), 0.02 * std::sin(0.7 * x));
            perturbation(i) = Complex(std::sin(1.3 * x), 0.6 * std::cos(0.9 * x));
        }

        gpes::PhysConfig phys{0.006, 0.008843, 100.0};
        const double mu = 0.83;
        gpes::bdg::BdGOperator2D bdg(grid, psi, phys, mu);
        Eigen::VectorXcd input(2 * n), output(2 * n);
        input.head(n) = perturbation;
        input.tail(n) = perturbation.conjugate();
        bdg.perform_op(input.data(), output.data());

        const double epsilon = 1e-6;
        const Eigen::VectorXcd finite_difference =
            (apply_egpe(psi + epsilon * perturbation, grid, phys, mu)
             - apply_egpe(psi - epsilon * perturbation, grid, phys, mu))
            / (2.0 * epsilon);
        const double relative_error = (output.head(n) - finite_difference).norm()
            / std::max(finite_difference.norm(), 1e-14);
        require(relative_error < 2e-7, "BdG operator does not match the eGPE Jacobian");
        const double lower_relative_error =
            (output.tail(n) + finite_difference.conjugate()).norm()
            / std::max(finite_difference.norm(), 1e-14);
        require(lower_relative_error < 2e-7,
            "BdG lower block violates the particle-hole convention");

        std::cout << "BdG finite-difference relative error: " << relative_error << '\n';
        std::cout << "BdG lower-block relative error: " << lower_relative_error << '\n';
        std::cout << "High-k dipolar form factor: " << form_factor << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
