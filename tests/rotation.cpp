#include "Solvers/AngularMomentum.hpp"
#include "Solvers/SplitStep.hpp"
#include "Core/utility.hpp"

#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {
constexpr auto Dim = gpes::Dimension::Two;
using Grid = gpes::Grid<Dim>;
using State = gpes::WaveFunction<Dim>;
using Complex = std::complex<double>;
using Stepper = gpes::solvers::SplitStep<Dim>;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

Eigen::VectorXcd vortex(const Grid& grid, int charge, double particles = 1.0) {
    Eigen::VectorXcd result(grid.size_x() * grid.size_y());
    for (unsigned int i = 0; i < grid.size_x(); ++i) {
        const double x = grid.start_pos_x() + i * grid.step_x();
        for (unsigned int j = 0; j < grid.size_y(); ++j) {
            const double y = grid.start_pos_y() + j * grid.step_y();
            const Complex winding = charge == 0
                ? Complex(1.0, 0.0)
                : std::pow(Complex(x, charge < 0 ? -y : y), std::abs(charge));
            result(i * grid.size_y() + j) = winding * std::exp(-0.5 * (x*x + y*y));
        }
    }
    result *= std::sqrt(particles / (result.squaredNorm() * grid.step_x() * grid.step_y()));
    return result;
}

void check_operator(int nx, int ny) {
    auto grid = std::make_shared<Grid>(nx, ny, -8.0, -10.0, 1.0, 1.0, 8.0);
    gpes::operators::AngularMomentum2D lz(grid);
    for (const int charge : {0, 1, -1, 2}) {
        const auto psi = vortex(*grid, charge, 3.0);
        Eigen::VectorXcd result;
        lz.apply(psi, result);
        const double error = (result - charge * psi).norm() / psi.norm();
        std::cout << "Lz vortex nx=" << nx << " ny=" << ny << " charge=" << charge << " error=" << error << '\n';
        require(error < 2e-10, "Lz does not reproduce analytic circulation");
        require(std::abs(lz.expectation(psi) - 3.0 * charge) < 2e-10, "Total Lz is missing the integration measure or N");
    }
    Eigen::VectorXcd u(nx * ny), v(nx * ny), lu, lv, conjugate;
    for (int k = 0; k < nx * ny; ++k) {
        u(k) = Complex(std::sin(0.37 * k), std::cos(0.19 * k));
        v(k) = Complex(std::sin(0.11 * k), std::cos(0.41 * k));
    }
    lz.apply(u, lu);
    lz.apply(v, lv);
    require(std::abs(u.dot(lv) - lu.dot(v)) < 1e-12 * u.norm() * lv.norm(), "Discrete Lz is not Hermitian");
    lz.apply(u.conjugate().eval(), conjugate);
    require((conjugate + lu.conjugate()).norm() < 1e-12 * lu.norm(), "Discrete Lz* is not -Lz");
    auto aliased = u;
    lz.apply(aliased, aliased);
    require((aliased - lu).norm() < 1e-12 * lu.norm(), "Lz cannot operate in place");
    if (nx % 2 == 0 && ny % 2 == 0) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) u(i * ny + j) = (i + j) % 2 ? -1.0 : 1.0;
        }
        lz.apply(u, lu);
        require(lu.norm() < 1e-10, "Even-grid Nyquist derivative convention is inconsistent");
    }
}

void check_diagnostics() {
    auto grid = std::make_shared<Grid>(64, 80, -8.0, -10.0, 1.0, 1.0, 8.0);
    const gpes::PhysConfig physics{1e-12, 0.0, 3.0};
    const double omega = 0.3;
    gpes::SimConfig sim;
    for (const int charge : {0, 1, -1, 2}) {
        State psi(grid);
        psi.vec() = vortex(*grid, charge, physics.num_of_prt);
        Stepper stepper(grid, psi, physics, sim, omega);
        const auto d = stepper.diagnostics(psi);
        const double expected = std::abs(charge) + 1.0 - omega * charge;
        require(std::abs(d.energy / physics.num_of_prt - expected) < 1e-9, "Wrong rotating energy");
        require(std::abs(d.lab_energy / physics.num_of_prt - (std::abs(charge) + 1.0)) < 1e-9, "Wrong lab energy");
        require(std::abs(d.energy - (d.lab_energy - omega * d.angular_momentum)) < 1e-12, "E_Omega != E - Omega Lz");
        require(std::abs(d.chemical_potential - expected) < 1e-9, "Chemical potential ignores rotation");
        require(d.relative_residual < 1e-9 && d.scaled_residual < 1e-9, "Residual does not use the rotating Hamiltonian");
        require(std::abs(d.angular_momentum_per_particle - charge) < 1e-10, "Wrong angular momentum per particle");
    }
    State psi(grid);
    psi.vec() = vortex(*grid, 0);
    bool rejected = false;
    try { Stepper unconfined(grid, psi, physics, sim, 1.0); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "Unconfined rotating ground-state request was accepted");
}

double real_time_error(double dt, double omega, int charge) {
    auto grid = std::make_shared<Grid>(64, 80, -8.0, -10.0, 1.0, 1.0, 8.0);
    const gpes::PhysConfig physics{1e-12, 0.0, 1.0};
    State psi(grid);
    const auto initial = vortex(*grid, charge);
    psi.vec() = initial;
    gpes::SimConfig sim;
    sim.ground_state = false;
    sim.dt = dt;
    Stepper stepper(grid, psi, physics, sim, omega);
    const int steps = static_cast<int>(std::lround(0.2 / dt));
    for (int k = 0; k < steps; ++k) stepper.step(psi);
    require(std::abs(psi.vec().norm() / initial.norm() - 1.0) < 2e-12, "Rotating real-time step is not norm preserving");
    const double energy = std::abs(charge) + 1.0 - omega * charge;
    return (psi.vec() - std::exp(Complex(0.0, -energy * steps * dt)) * initial).norm() / initial.norm();
}

void check_imaginary_time() {
    auto grid = std::make_shared<Grid>(64, 80, -8.0, -10.0, 1.0, 1.0, 8.0);
    const gpes::PhysConfig physics{1e-12, 0.0, 1.0};
    const auto plus = vortex(*grid, 1), minus = vortex(*grid, -1);
    for (double omega : {-0.3, 0.3}) {
        State psi(grid);
        psi.vec() = (plus + minus) / std::sqrt(2.0);
        gpes::SimConfig sim;
        sim.dt = 0.004;
        Stepper stepper(grid, psi, physics, sim, omega);
        for (int k = 0; k < 25; ++k) stepper.step(psi);
        stepper.set_time_step(0.002);
        for (int k = 0; k < 50; ++k) stepper.step(psi);
        const Complex ratio = plus.dot(psi.vec()) / minus.dot(psi.vec());
        require(std::abs(ratio - std::exp(2.0 * omega * 0.2)) < 2e-5, "Imaginary-time propagation has the wrong rotation sign or time step");
        require(std::abs(psi.vec().squaredNorm() * grid->step_x() * grid->step_y() - 1.0) < 2e-12, "Imaginary-time normalization changed N");
    }
}

} // namespace

int main() {
    try {
        check_operator(64, 80);
        check_operator(65, 81);
        check_diagnostics();
        for (double omega : {-0.3, 0.0, 0.3}) {
            for (int charge : {-1, 1}) {
                const double coarse = real_time_error(0.02, omega, charge);
                const double fine = real_time_error(0.01, omega, charge);
                std::cout << "Omega=" << omega << ", m=" << charge << ": coarse=" << coarse << ", fine=" << fine << '\n';
                require(fine < 0.3 * coarse && fine < 1e-4, "Linear rotating split step fails second-order refinement");
            }
        }
        check_imaginary_time();
        std::cout << "Angular momentum, rotating diagnostics and propagation passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
