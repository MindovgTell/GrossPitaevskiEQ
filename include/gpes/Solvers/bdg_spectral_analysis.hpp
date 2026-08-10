#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <fftw3.h>

#include "Solvers/BogoliubovdeGennes.hpp"

// Companion utilities for gpes::bdg.
// The BdG header is included explicitly so this file is self-contained.

namespace gpes::bdg::analysis {

using Complex = std::complex<double>;

struct RadialSpectrumConfig {
    int n_k_bins = 128;
    int n_omega_bins = 256;
    double k_max = -1.0;            // <= 0: use full FFT-grid maximum
    double omega_min = 0.0;
    double omega_max = -1.0;        // <= 0: infer from selected modes
    double eta = 0.02;              // Lorentzian HWHM
    int j_max = -1;                 // <= 0: all supplied modes
    bool positive_frequency_only = true;
    bool positive_norm_only = true;
    double max_abs_imag_omega = 1e-7;
    double density_floor_fraction = 1e-12;
};

struct ModeMomentumSummary {
    std::size_t mode = 0;
    Complex omega{};
    double bdg_norm = 0.0;
    double density_weight = 0.0;
    double k_peak = 0.0;
    double k_mean = 0.0;
    double k_rms = 0.0;
    std::vector<double> radial_weight;
};

struct DynamicStructureFactor {
    std::vector<double> k;
    std::vector<double> omega;
    // Row-major: S[ik * omega.size() + iw]
    std::vector<double> values;
    std::vector<ModeMomentumSummary> modes;
};

struct BdGValidationReport {
    double stationary_relative_residual = std::numeric_limits<double>::quiet_NaN();
    double goldstone_residual = std::numeric_limits<double>::quiet_NaN();
    double max_pairing_error = std::numeric_limits<double>::quiet_NaN();
    int unstable_mode_count = 0;
    int near_zero_mode_count = 0;
};

inline int flat_index(int i, int j, int ny) { return i * ny + j; }

inline Eigen::VectorXcd density_fluctuation(
    const Eigen::VectorXcd& psi0,
    const Eigen::VectorXcd& u,
    const Eigen::VectorXcd& v)
{
    if (psi0.size() != u.size() || u.size() != v.size())
        throw std::invalid_argument("density_fluctuation: vector size mismatch");

    // General complex-background convention. For real psi0 this is psi0 * (u + v).
    return (psi0.conjugate().array() * u.array()
          + psi0.array() * v.array()).matrix();
}

inline Eigen::VectorXcd phase_fluctuation(
    const Eigen::VectorXcd& psi0,
    const Eigen::VectorXcd& u,
    const Eigen::VectorXcd& v,
    double density_floor_fraction = 1e-12)
{
    if (psi0.size() != u.size() || u.size() != v.size())
        throw std::invalid_argument("phase_fluctuation: vector size mismatch");

    Eigen::VectorXcd out(psi0.size());
    double max_density = 0.0;
    for (Eigen::Index i = 0; i < psi0.size(); ++i)
        max_density = std::max(max_density, std::norm(psi0(i)));
    const double floor = density_floor_fraction * max_density;

    for (Eigen::Index i = 0; i < psi0.size(); ++i) {
        const double n = std::norm(psi0(i));
        if (n <= floor) {
            out(i) = 0.0;
        } else {
            out(i) = (std::conj(psi0(i)) * u(i) - psi0(i) * v(i)) /
                     Complex(0.0, 2.0 * n);
        }
    }
    return out;
}

inline Eigen::VectorXcd fft2_forward(
    const Eigen::VectorXcd& input,
    int nx,
    int ny,
    double dx,
    double dy)
{
    if (input.size() != nx * ny)
        throw std::invalid_argument("fft2_forward: input size mismatch");

    Eigen::VectorXcd output = input;
    auto* ptr = reinterpret_cast<fftw_complex*>(output.data());
    fftw_plan plan = fftw_plan_dft_2d(nx, ny, ptr, ptr, FFTW_FORWARD, FFTW_ESTIMATE);
    if (!plan)
        throw std::runtime_error("fft2_forward: FFTW plan creation failed");
    fftw_execute_dft(plan, ptr, ptr);
    fftw_destroy_plan(plan);

    // Continuum-transform convention: integral d^2r f(r) exp(-i k.r).
    output *= dx * dy;
    return output;
}

inline std::vector<double> fft_k_magnitudes(
    int nx, int ny, double dx, double dy)
{
    const double lx = nx * dx;
    const double ly = ny * dy;
    const double dkx = 2.0 * M_PI / lx;
    const double dky = 2.0 * M_PI / ly;
    std::vector<double> kmag(static_cast<std::size_t>(nx * ny));

    for (int i = 0; i < nx; ++i) {
        const int qx = (i <= nx / 2) ? i : i - nx;
        const double kx = dkx * qx;
        for (int j = 0; j < ny; ++j) {
            const int qy = (j <= ny / 2) ? j : j - ny;
            const double ky = dky * qy;
            kmag[static_cast<std::size_t>(flat_index(i, j, ny))] = std::hypot(kx, ky);
        }
    }
    return kmag;
}

inline double lorentzian(double x, double eta)
{
    return eta / (M_PI * (x * x + eta * eta));
}

inline DynamicStructureFactor compute_dynamic_structure_factor(
    const BdGResult& result,
    const Eigen::VectorXcd& psi0,
    const Grid<Dimension::Two>& grid,
    RadialSpectrumConfig cfg = {})
{
    const int nx = static_cast<int>(grid.size_x());
    const int ny = static_cast<int>(grid.size_y());
    const int n = nx * ny;
    const double dx = grid.step_x();
    const double dy = grid.step_y();

    if (psi0.size() != n)
        throw std::invalid_argument("compute_dynamic_structure_factor: psi0/grid mismatch");
    if (cfg.n_k_bins < 2 || cfg.n_omega_bins < 2 || cfg.eta <= 0.0)
        throw std::invalid_argument("compute_dynamic_structure_factor: invalid binning/broadening");

    const auto kmag = fft_k_magnitudes(nx, ny, dx, dy);
    const double full_kmax = *std::max_element(kmag.begin(), kmag.end());
    const double kmax = cfg.k_max > 0.0 ? std::min(cfg.k_max, full_kmax) : full_kmax;
    const double dk = kmax / static_cast<double>(cfg.n_k_bins);

    const int mode_limit = cfg.j_max > 0
        ? std::min<int>(cfg.j_max, static_cast<int>(result.modes.size()))
        : static_cast<int>(result.modes.size());

    std::vector<ModeMomentumSummary> selected;
    selected.reserve(static_cast<std::size_t>(mode_limit));
    double inferred_omega_max = cfg.omega_min;

    for (int m = 0; m < mode_limit; ++m) {
        const auto& mode = result.modes[static_cast<std::size_t>(m)];
        if (mode.u.size() != n || mode.v.size() != n)
            throw std::invalid_argument("compute_dynamic_structure_factor: mode/grid mismatch");

        const double qnorm = bdg_mode_norm(mode, dx, dy);
        if (cfg.positive_frequency_only && mode.omega.real() <= 0.0) continue;
        if (cfg.positive_norm_only && qnorm <= 0.0) continue;
        if (std::abs(mode.omega.imag()) > cfg.max_abs_imag_omega) continue;

        const Eigen::VectorXcd dn = density_fluctuation(psi0, mode.u, mode.v);
        const Eigen::VectorXcd dnk = fft2_forward(dn, nx, ny, dx, dy);

        ModeMomentumSummary summary;
        summary.mode = static_cast<std::size_t>(m);
        summary.omega = mode.omega;
        summary.bdg_norm = qnorm;
        summary.radial_weight.assign(static_cast<std::size_t>(cfg.n_k_bins), 0.0);
        std::vector<int> counts(static_cast<std::size_t>(cfg.n_k_bins), 0);

        // |delta n(k)|^2. Radial bin sum approximates angular integration.
        for (int p = 0; p < n; ++p) {
            const double k = kmag[static_cast<std::size_t>(p)];
            if (k > kmax) continue;
            int ik = static_cast<int>(std::floor(k / dk));
            if (ik == cfg.n_k_bins) ik = cfg.n_k_bins - 1;
            if (ik < 0 || ik >= cfg.n_k_bins) continue;
            summary.radial_weight[static_cast<std::size_t>(ik)] += std::norm(dnk(p));
            ++counts[static_cast<std::size_t>(ik)];
        }

        // Divide by number of lattice points in an annulus to avoid the trivial 2*pi*k DOS factor.
        for (int ik = 0; ik < cfg.n_k_bins; ++ik) {
            if (counts[static_cast<std::size_t>(ik)] > 0)
                summary.radial_weight[static_cast<std::size_t>(ik)] /=
                    static_cast<double>(counts[static_cast<std::size_t>(ik)]);
        }

        double wsum = 0.0, kwsum = 0.0, k2wsum = 0.0;
        int peak_bin = 0;
        for (int ik = 0; ik < cfg.n_k_bins; ++ik) {
            const double k = (ik + 0.5) * dk;
            const double w = summary.radial_weight[static_cast<std::size_t>(ik)];
            if (w > summary.radial_weight[static_cast<std::size_t>(peak_bin)]) peak_bin = ik;
            wsum += w;
            kwsum += k * w;
            k2wsum += k * k * w;
        }
        summary.density_weight = wsum;
        summary.k_peak = (peak_bin + 0.5) * dk;
        if (wsum > 0.0) {
            summary.k_mean = kwsum / wsum;
            summary.k_rms = std::sqrt(k2wsum / wsum);
        }

        inferred_omega_max = std::max(inferred_omega_max, mode.omega.real());
        selected.push_back(std::move(summary));
    }

    const double omega_max = cfg.omega_max > cfg.omega_min
        ? cfg.omega_max
        : std::max(cfg.omega_min + 10.0 * cfg.eta, 1.05 * inferred_omega_max);

    DynamicStructureFactor dsf;
    dsf.k.resize(static_cast<std::size_t>(cfg.n_k_bins));
    dsf.omega.resize(static_cast<std::size_t>(cfg.n_omega_bins));
    dsf.values.assign(static_cast<std::size_t>(cfg.n_k_bins * cfg.n_omega_bins), 0.0);
    dsf.modes = selected;

    for (int ik = 0; ik < cfg.n_k_bins; ++ik)
        dsf.k[static_cast<std::size_t>(ik)] = (ik + 0.5) * dk;
    const double domega = (omega_max - cfg.omega_min) /
                          static_cast<double>(cfg.n_omega_bins - 1);
    for (int iw = 0; iw < cfg.n_omega_bins; ++iw)
        dsf.omega[static_cast<std::size_t>(iw)] = cfg.omega_min + iw * domega;

    for (const auto& mode : dsf.modes) {
        for (int ik = 0; ik < cfg.n_k_bins; ++ik) {
            const double wk = mode.radial_weight[static_cast<std::size_t>(ik)];
            for (int iw = 0; iw < cfg.n_omega_bins; ++iw) {
                const double broaden = lorentzian(
                    dsf.omega[static_cast<std::size_t>(iw)] - mode.omega.real(), cfg.eta);
                dsf.values[static_cast<std::size_t>(ik * cfg.n_omega_bins + iw)] += wk * broaden;
            }
        }
    }
    return dsf;
}

inline void save_dynamic_structure_factor(
    const DynamicStructureFactor& dsf,
    const std::filesystem::path& folder)
{
    std::filesystem::create_directories(folder);

    std::ofstream modes(folder / "mode_momentum_summary.csv");
    if (!modes) throw std::runtime_error("cannot write mode_momentum_summary.csv");
    modes << "mode,omega_real,omega_imag,bdg_norm,density_weight,k_peak,k_mean,k_rms\n";
    modes << std::setprecision(17);
    for (const auto& m : dsf.modes) {
        modes << m.mode << ',' << m.omega.real() << ',' << m.omega.imag() << ','
              << m.bdg_norm << ',' << m.density_weight << ',' << m.k_peak << ','
              << m.k_mean << ',' << m.k_rms << '\n';
    }

    std::ofstream out(folder / "dynamic_structure_factor.csv");
    if (!out) throw std::runtime_error("cannot write dynamic_structure_factor.csv");
    out << "k,omega,S\n" << std::setprecision(17);
    for (std::size_t ik = 0; ik < dsf.k.size(); ++ik)
        for (std::size_t iw = 0; iw < dsf.omega.size(); ++iw)
            out << dsf.k[ik] << ',' << dsf.omega[iw] << ','
                << dsf.values[ik * dsf.omega.size() + iw] << '\n';
}

inline double quartet_pairing_error(const std::vector<BdGMode>& modes)
{
    if (modes.empty()) return 0.0;
    double max_error = 0.0;
    for (const auto& a : modes) {
        const Complex targets[3] = {-a.omega, std::conj(a.omega), -std::conj(a.omega)};
        for (Complex target : targets) {
            double best = std::numeric_limits<double>::infinity();
            for (const auto& b : modes) best = std::min(best, std::abs(b.omega - target));
            max_error = std::max(max_error, best / std::max(1.0, std::abs(a.omega)));
        }
    }
    return max_error;
}

inline BdGValidationReport validate_spectrum(
    const BdGResult& result,
    double zero_tolerance = 1e-7,
    double instability_tolerance = 1e-7)
{
    BdGValidationReport report;
    report.max_pairing_error = quartet_pairing_error(result.modes);
    for (const auto& mode : result.modes) {
        if (std::abs(mode.omega) < zero_tolerance) ++report.near_zero_mode_count;
        if (std::abs(mode.omega.imag()) > instability_tolerance) ++report.unstable_mode_count;
    }
    return report;
}

// Apply any eGPE stationary operator supplied by the caller and calculate
// ||H[psi]psi - mu psi|| / max(||H[psi]psi||, |mu| ||psi||).
template<class ApplyStationaryOperator>
inline double stationary_relative_residual(
    const Eigen::VectorXcd& psi0,
    double mu,
    ApplyStationaryOperator&& apply_h)
{
    Eigen::VectorXcd hpsi(psi0.size());
    apply_h(psi0, hpsi);
    const double numerator = (hpsi - mu * psi0).norm();
    const double denominator = std::max({hpsi.norm(), std::abs(mu) * psi0.norm(), 1e-300});
    return numerator / denominator;
}

// Matrix-free Goldstone check: q=(psi0,-psi0*) should satisfy L q approximately zero.
template<class BdGOperatorType>
inline double goldstone_relative_residual(
    const BdGOperatorType& op,
    const Eigen::VectorXcd& psi0)
{
    Eigen::VectorXcd q(2 * psi0.size()), lq(2 * psi0.size());
    q.head(psi0.size()) = psi0;
    q.tail(psi0.size()) = -psi0.conjugate();
    op.perform_op(q.data(), lq.data());
    return lq.norm() / std::max(q.norm(), 1e-300);
}

// Symplectic Gram-Schmidt against the phase mode. Useful for an initial vector,
// but a fully projected Arnoldi iteration requires projection inside perform_op().
inline void project_out_phase_mode(
    Eigen::VectorXcd& x,
    const Eigen::VectorXcd& psi0,
    double dx,
    double dy)
{
    const Eigen::Index n = psi0.size();
    if (x.size() != 2 * n) throw std::invalid_argument("project_out_phase_mode: size mismatch");
    Eigen::VectorXcd q(2 * n);
    q.head(n) = psi0;
    q.tail(n) = -psi0.conjugate();

    auto j_inner = [n, dx, dy](const Eigen::VectorXcd& a, const Eigen::VectorXcd& b) {
        return dx * dy * (a.head(n).dot(b.head(n)) - a.tail(n).dot(b.tail(n)));
    };

    const Complex denom = j_inner(q, q);
    if (std::abs(denom) > 1e-14) x -= q * (j_inner(q, x) / denom);
}

struct TimingEstimate {
    double seconds_per_operator_application = 0.0;
    double estimated_seconds = 0.0;
    std::size_t assumed_operator_applications = 0;
};

template<class BdGOperatorType>
inline TimingEstimate estimate_arnoldi_runtime(
    const BdGOperatorType& op,
    int requested_modes,
    int ncv,
    int repetitions = 30,
    double operation_factor = 4.0)
{
    if (repetitions <= 0 || requested_modes <= 0 || ncv <= requested_modes)
        throw std::invalid_argument("estimate_arnoldi_runtime: invalid arguments");

    Eigen::VectorXcd x = Eigen::VectorXcd::Random(op.cols());
    Eigen::VectorXcd y(op.rows());
    op.perform_op(x.data(), y.data()); // warm-up

    const auto t0 = std::chrono::steady_clock::now();
    for (int r = 0; r < repetitions; ++r) {
        op.perform_op(x.data(), y.data());
        x.swap(y);
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double per_op = std::chrono::duration<double>(t1 - t0).count() / repetitions;

    // Arnoldi operation count is problem-dependent. This gives a planning estimate,
    // not a guarantee. Replace with eigs.num_operations() after a smaller pilot run.
    const std::size_t assumed = static_cast<std::size_t>(
        std::ceil(operation_factor * ncv * std::max(1, requested_modes / 8)));
    return {per_op, per_op * static_cast<double>(assumed), assumed};
}

} // namespace gpes::bdg::analysis
