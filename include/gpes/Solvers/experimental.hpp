#pragma once

#include <cmath>
#include <complex>
#include <stdexcept>

#include <Eigen/Dense>
#include <fftw3.h>

#include "Core/utility.hpp"
#include "grid/grid.hpp"
#include "wavefunction/wavefunction.hpp"

namespace gpes::solvers::experimental {

struct SystemEnergy {
    double kinetic = 0.0;
    double potential = 0.0;
    double contact = 0.0;
    double dipolar = 0.0;
    double lhy = 0.0;

    double total() const {
        return kinetic + potential + contact + dipolar + lhy;
    }
};

inline double trap_weight_1d(int i, int size) {
    if (size <= 1) {
        return 1.0;
    }
    return (i == 0 || i == size - 1) ? 0.5 : 1.0;
}

inline SystemEnergy energy_1d(
    const Grid<Dimension::One>& grid,
    const Eigen::VectorXcd& psi,
    const PhysConfig& phys,
    const Eigen::VectorXd* ddi_potential = nullptr
) {
    const int n = static_cast<int>(grid.size());
    if (psi.size() != n) {
        throw std::invalid_argument("energy_1d: wavefunction and grid size mismatch");
    }
    if (ddi_potential && ddi_potential->size() != n) {
        throw std::invalid_argument("energy_1d: DDI potential size mismatch");
    }
    if (n <= 0 || !std::isfinite(grid.step()) || grid.step() == 0.0) {
        throw std::invalid_argument("energy_1d: invalid grid");
    }

    SystemEnergy out;
    const double dx = grid.step();
    const double L = static_cast<double>(n) * dx;
    const double dk = 2.0 * M_PI / L;

    Eigen::VectorXcd psi_k = psi;
    auto* data = reinterpret_cast<fftw_complex*>(psi_k.data());
    fftw_plan plan = fftw_plan_dft_1d(n, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (int m = 0; m < n; ++m) {
        const int q = (m <= n / 2) ? m : (m - n);
        const double k = dk * static_cast<double>(q);
        out.kinetic += 0.5 * (k * k) * std::norm(psi_k(m));
    }
    out.kinetic *= dx / static_cast<double>(n);

    for (int i = 0; i < n; ++i) {
        const double w = trap_weight_1d(i, n);
        const double density = std::norm(psi(i));
        const double measure = dx * w;

        out.potential += measure * grid.potential()(i) * density;
        out.contact += measure * 0.5 * phys.g_scat * density * density;
        out.lhy += measure * 0.4 * phys.g_lhy * std::pow(density, 2.5);
        if (ddi_potential) {
            out.dipolar += measure * 0.5 * (*ddi_potential)(i) * density;
        }
    }

    return out;
}

inline SystemEnergy energy_1d(
    const WaveFunction<Dimension::One>& psi,
    const PhysConfig& phys,
    const Eigen::VectorXd* ddi_potential = nullptr
) {
    if (!psi.grid()) {
        throw std::invalid_argument("energy_1d: wavefunction has no grid");
    }
    return energy_1d(*psi.grid(), psi.vec(), phys, ddi_potential);
}

inline SystemEnergy energy_2d(
    const Grid<Dimension::Two>& grid,
    const Eigen::VectorXcd& psi,
    const PhysConfig& phys,
    const Eigen::VectorXd* ddi_potential = nullptr
) {
    const int nx = static_cast<int>(grid.size_x());
    const int ny = static_cast<int>(grid.size_y());
    const int n = nx * ny;

    if (psi.size() != n) {
        throw std::invalid_argument("energy_2d: wavefunction and grid size mismatch");
    }
    if (ddi_potential && ddi_potential->size() != n) {
        throw std::invalid_argument("energy_2d: DDI potential size mismatch");
    }
    if (nx <= 0 || ny <= 0 ||
        !std::isfinite(grid.step_x()) || grid.step_x() == 0.0 ||
        !std::isfinite(grid.step_y()) || grid.step_y() == 0.0) {
        throw std::invalid_argument("energy_2d: invalid grid");
    }

    SystemEnergy out;
    const double dx = grid.step_x();
    const double dy = grid.step_y();
    const double Lx = static_cast<double>(nx) * dx;
    const double Ly = static_cast<double>(ny) * dy;
    const double dkx = 2.0 * M_PI / Lx;
    const double dky = 2.0 * M_PI / Ly;

    Eigen::VectorXcd psi_k = psi;
    auto* data = reinterpret_cast<fftw_complex*>(psi_k.data());
    fftw_plan plan = fftw_plan_dft_2d(nx, ny, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (int i = 0; i < nx; ++i) {
        const int qx = (i <= nx / 2) ? i : (i - nx);
        const double kx = dkx * static_cast<double>(qx);
        for (int j = 0; j < ny; ++j) {
            const int qy = (j <= ny / 2) ? j : (j - ny);
            const double ky = dky * static_cast<double>(qy);
            const int idx = i * ny + j;
            out.kinetic += 0.5 * (kx * kx + ky * ky) * std::norm(psi_k(idx));
        }
    }
    out.kinetic *= (dx * dy) / static_cast<double>(n);

    for (int i = 0; i < nx; ++i) {
        const double wx = trap_weight_1d(i, nx);
        for (int j = 0; j < ny; ++j) {
            const double wy = trap_weight_1d(j, ny);
            const int idx = i * ny + j;
            const double density = std::norm(psi(idx));
            const double measure = dx * dy * wx * wy;

            out.potential += measure * grid.potential()(i, j) * density;
            out.contact += measure * 0.5 * phys.g_scat * density * density;
            out.lhy += measure * 0.4 * phys.g_lhy * std::pow(density, 2.5);
            if (ddi_potential) {
                out.dipolar += measure * 0.5 * (*ddi_potential)(idx) * density;
            }
        }
    }

    return out;
}

inline SystemEnergy energy_2d(
    const WaveFunction<Dimension::Two>& psi,
    const PhysConfig& phys,
    const Eigen::VectorXd* ddi_potential = nullptr
) {
    if (!psi.grid()) {
        throw std::invalid_argument("energy_2d: wavefunction has no grid");
    }
    return energy_2d(*psi.grid(), psi.vec(), phys, ddi_potential);
}

} // namespace gpes::solvers::experimental

