#pragma once

#include "CLI11.hpp"
#include "Core/utility.hpp"

#include <chrono>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>

namespace workflow {

namespace fs = std::filesystem;
constexpr auto Dim = gpes::Dimension::Two;
using Grid = gpes::Grid<Dim>;
using Wavefunction = gpes::WaveFunction<Dim>;

struct Parameters {
    int nx = 1024;
    int ny = 1024;
    double x_half_length = 25.0;
    double y_half_length = 25.0;
    double a_s = 0.00265;
    double a_dd = 0.008843;
    double num_particles = 80000.0;
    double omega_x = 1.0;
    double omega_y = 5.0;
    double omega_z = 8.0;
    double rotation_frequency = 0.0;

    gpes::PhysConfig physics() const { return {a_s, a_dd, num_particles}; }
    std::shared_ptr<Grid> grid() const {
        return std::make_shared<Grid>(nx, ny, -x_half_length, -y_half_length,
                                      omega_x, omega_y, omega_z);
    }
};

inline void positive(double value, const std::string& name, bool allow_zero = false) {
    if (!std::isfinite(value) || (allow_zero ? value < 0.0 : value <= 0.0)) {
        throw std::invalid_argument(name + " must be finite and " + (allow_zero ? "nonnegative" : "positive"));
    }
}

inline void validate(const Parameters& p) {
    if (p.nx < 8 || p.ny < 8 || static_cast<long long>(p.nx) * p.ny > std::numeric_limits<int>::max() / 4) {
        throw std::invalid_argument("Grid dimensions must be at least 8 and fit the solver's integer indexing");
    }
    positive(p.x_half_length, "x_half_length");
    positive(p.y_half_length, "y_half_length");
    positive(p.a_s, "a_s");
    positive(p.a_dd, "a_dd", true);
    positive(p.num_particles, "num_particles");
    positive(p.omega_x, "omega_x");
    positive(p.omega_y, "omega_y");
    positive(p.omega_z, "omega_z");
    if (!std::isfinite(p.rotation_frequency) ||
        !(std::abs(p.rotation_frequency) < std::min(p.omega_x, p.omega_y))) {
        throw std::invalid_argument("A trapped rotating ground state requires finite |Omega| < min(omega_x, omega_y)");
    }
    if (p.num_particles != std::round(p.num_particles) ||
        p.num_particles > std::numeric_limits<unsigned int>::max()) {
        throw std::invalid_argument("num_particles must be an integer representable as unsigned int");
    }
}

inline Parameters read_parameters(const fs::path& run, bool allow_rotation = false) {
    const auto path = run / "wavefunction/parameters.txt";
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Missing parameter manifest: " + path.string());
    Parameters p;
    CLI::App manifest;
    manifest.config_formatter(std::make_shared<CLI::ConfigINI>());
    manifest.add_option("--nx", p.nx)->required();
    manifest.add_option("--ny", p.ny)->required();
    manifest.add_option("--x_half_length", p.x_half_length)->required();
    manifest.add_option("--y_half_length", p.y_half_length)->required();
    manifest.add_option("--a_s", p.a_s)->required();
    manifest.add_option("--a_dd", p.a_dd)->required();
    manifest.add_option("--num_particles", p.num_particles)->required();
    manifest.add_option("--omega_x", p.omega_x)->required();
    manifest.add_option("--omega_y", p.omega_y)->required();
    manifest.add_option("--omega_z", p.omega_z)->required();
    // Only the ground-state caller opts in. The unchanged BdG reader must not
    // silently discard a rotating Hamiltonian that it does not implement.
    if (allow_rotation) manifest.add_option("--rotation_frequency", p.rotation_frequency);
    manifest.parse_from_stream(in);
    validate(p);
    return p;
}

inline void write_text(const fs::path& path, const std::string& text) {
    fs::create_directories(path.parent_path());
    const fs::path temporary = path.string() + ".tmp";
    std::ofstream out(temporary);
    if (!out) throw std::runtime_error("Cannot write " + temporary.string());
    out << text;
    out.close();
    if (!out) throw std::runtime_error("Write failed: " + temporary.string());
    fs::rename(temporary, path);
}

inline void write_parameters(const fs::path& run, const Parameters& p) {
    std::ostringstream out;
    out << std::setprecision(17)
        << "nx=" << p.nx << '\n' << "ny=" << p.ny << '\n'
        << "x_half_length=" << p.x_half_length << '\n'
        << "y_half_length=" << p.y_half_length << '\n'
        << "a_s=" << p.a_s << '\n' << "a_dd=" << p.a_dd << '\n'
        << "num_particles=" << p.num_particles << '\n'
        << "omega_x=" << p.omega_x << '\n' << "omega_y=" << p.omega_y << '\n'
        << "omega_z=" << p.omega_z << '\n';
    // Absence means Omega=0, preserving compatibility with older saved states.
    if (p.rotation_frequency != 0.0) out << "rotation_frequency=" << p.rotation_frequency << '\n';
    write_text(run / "wavefunction/parameters.txt", out.str());
}

inline fs::path default_root() {
    const char* scratch = std::getenv("SCRATCH");
    return scratch && *scratch ? fs::path(scratch) / "GPES" : fs::path("runs");
}

inline fs::path new_output(const std::string& requested, const fs::path& parent, const std::string& prefix) {
    if (!requested.empty()) {
        const auto path = fs::absolute(requested).lexically_normal();
        fs::create_directories(path.parent_path());
        if (!fs::create_directory(path)) {
            throw std::runtime_error("Output directory already exists; choose a new path: " + path.string());
        }
        return path;
    }
    fs::create_directories(parent);
    const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::ostringstream stamp;
    stamp << prefix << '_' << std::put_time(std::gmtime(&now), "%Y%m%d_%H%M%S");
    std::random_device random;
    for (int trial = 0; trial < 100; ++trial) {
        const auto path = parent / (stamp.str() + '_' + std::to_string(random()));
        if (fs::create_directory(path)) return fs::absolute(path).lexically_normal();
    }
    throw std::runtime_error("Cannot create a unique output directory");
}

inline double particle_number(const Wavefunction& psi) {
    return psi.vec().squaredNorm() * psi.grid()->step_x() * psi.grid()->step_y();
}

inline void check_normalization(const Wavefunction& psi, const Parameters& p) {
    const double norm = particle_number(psi);
    if (!psi.vec().allFinite() || !std::isfinite(norm) || std::abs(norm - p.num_particles) / p.num_particles > 1e-8) {
        throw std::runtime_error("Wavefunction normalization does not match the parameter manifest");
    }
}

inline double boundary_ratio(const Wavefunction& psi) {
    const int nx = static_cast<int>(psi.grid()->size_x());
    const int ny = static_cast<int>(psi.grid()->size_y());
    double peak = 0.0, boundary = 0.0;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            const double density = std::norm(psi.vec()(i * ny + j));
            peak = std::max(peak, density);
            if (i == 0 || j == 0 || i == nx - 1 || j == ny - 1) boundary = std::max(boundary, density);
        }
    }
    return peak > 0.0 ? boundary / peak : std::numeric_limits<double>::infinity();
}

inline void model_warning(const Parameters& p) {
    if (p.omega_z <= std::max(p.omega_x, p.omega_y)) {
        std::cerr << "Model warning: omega_z does not exceed both in-plane frequencies; "
                     "the frozen-z quasi-2D approximation is not controlled.\n";
    }
}

inline void record_failure(const fs::path& run, const std::string& error) noexcept {
    if (run.empty()) return;
    try {
        write_text(run / "FAILED", error + '\n');
        fs::remove(run / "RUNNING");
    }
    catch (const std::exception& e) { std::cerr << "Cannot save failure report: " << e.what() << '\n'; }
}

} // namespace workflow
