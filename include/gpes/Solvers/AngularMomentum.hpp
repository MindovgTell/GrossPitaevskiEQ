#pragma once

#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <numbers>
#include <stdexcept>

#include <Eigen/Core>
#include <fftw3.h>

#include "grid/grid.hpp"

namespace gpes::operators {

namespace detail {

inline double fourier_wavenumber(int index, int size, double spacing) {
    const int mode = index <= size / 2 ? index : index - size;
    return 2.0 * std::numbers::pi * mode / (size * spacing);
}

inline double first_derivative_wavenumber(int index, int size, double spacing) {
    // The even-grid Nyquist mode has no unique signed derivative. Zeroing it
    // preserves a real differentiation matrix and Lz* = -Lz on the grid.
    return size % 2 == 0 && index == size / 2
        ? 0.0 : fourier_wavenumber(index, size, spacing);
}

} // namespace detail

// Lz = -i (x d/dy - y d/dx), in units of hbar, on the solver's periodic grid.
// FFT work buffers are reused; one instance must not be used concurrently.
class AngularMomentum2D {
public:
    explicit AngularMomentum2D(std::shared_ptr<const Grid<Dimension::Two>> grid)
        : grid_(std::move(grid)) {
        if (!grid_) throw std::invalid_argument("AngularMomentum2D: missing grid");
        nx_ = static_cast<int>(grid_->size_x());
        ny_ = static_cast<int>(grid_->size_y());
        if (nx_ <= 0 || ny_ <= 0 ||
            static_cast<long long>(nx_) * ny_ > std::numeric_limits<int>::max() ||
            !(grid_->step_x() > 0.0) || !std::isfinite(grid_->step_x()) ||
            !(grid_->step_y() > 0.0) || !std::isfinite(grid_->step_y())) {
            throw std::invalid_argument("AngularMomentum2D: invalid Fourier grid");
        }
        fourier_.resize(nx_ * ny_);
        derivative_.resize(nx_ * ny_);
        auto* spectrum = reinterpret_cast<fftw_complex*>(fourier_.data());
        auto* derivative = reinterpret_cast<fftw_complex*>(derivative_.data());
        forward_ = fftw_plan_dft_2d(nx_, ny_, spectrum, spectrum, FFTW_FORWARD, FFTW_ESTIMATE);
        backward_ = fftw_plan_dft_2d(nx_, ny_, derivative, derivative, FFTW_BACKWARD, FFTW_ESTIMATE);
        if (!forward_ || !backward_) {
            if (forward_) fftw_destroy_plan(forward_);
            if (backward_) fftw_destroy_plan(backward_);
            throw std::runtime_error("AngularMomentum2D: FFTW plan creation failed");
        }
    }

    ~AngularMomentum2D() {
        fftw_destroy_plan(forward_);
        fftw_destroy_plan(backward_);
    }
    AngularMomentum2D(const AngularMomentum2D&) = delete;
    AngularMomentum2D& operator=(const AngularMomentum2D&) = delete;

    void apply(const Eigen::VectorXcd& psi, Eigen::VectorXcd& result) {
        if (psi.size() != fourier_.size()) {
            throw std::invalid_argument("AngularMomentum2D: state/grid size mismatch");
        }
        fourier_ = psi;
        fftw_execute(forward_);
        const std::complex<double> imaginary(0.0, 1.0);
        const double inverse_size = 1.0 / static_cast<double>(fourier_.size());
        result.resize(psi.size());

        for (int i = 0; i < nx_; ++i) {
            const double kx = detail::first_derivative_wavenumber(i, nx_, grid_->step_x());
            for (int j = 0; j < ny_; ++j) {
                const int k = i * ny_ + j;
                derivative_(k) = imaginary * kx * fourier_(k);
            }
        }
        fftw_execute(backward_);
        for (int i = 0; i < nx_; ++i) {
            for (int j = 0; j < ny_; ++j) {
                const int k = i * ny_ + j;
                const double y = grid_->start_pos_y() + j * grid_->step_y();
                result(k) = imaginary * y * derivative_(k) * inverse_size;
            }
        }

        for (int i = 0; i < nx_; ++i) {
            for (int j = 0; j < ny_; ++j) {
                const int k = i * ny_ + j;
                const double ky = detail::first_derivative_wavenumber(j, ny_, grid_->step_y());
                derivative_(k) = imaginary * ky * fourier_(k);
            }
        }
        fftw_execute(backward_);
        for (int i = 0; i < nx_; ++i) {
            const double x = grid_->start_pos_x() + i * grid_->step_x();
            for (int j = 0; j < ny_; ++j) {
                const int k = i * ny_ + j;
                result(k) -= imaginary * x * derivative_(k) * inverse_size;
            }
        }
    }

    double expectation(const Eigen::VectorXcd& psi) {
        Eigen::VectorXcd lz_psi;
        apply(psi, lz_psi);
        return std::real(psi.dot(lz_psi)) * grid_->step_x() * grid_->step_y();
    }

private:
    std::shared_ptr<const Grid<Dimension::Two>> grid_;
    int nx_ = 0, ny_ = 0;
    Eigen::VectorXcd fourier_, derivative_;
    fftw_plan forward_ = nullptr, backward_ = nullptr;
};

} // namespace gpes::operators
