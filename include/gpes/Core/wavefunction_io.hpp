#pragma once

#include <complex>
#include <cctype>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "wavefunction/wavefunction.hpp"

namespace gpes::io {

namespace detail {

inline std::string trim_copy(const std::string& s) {
    std::size_t begin = 0;
    while (begin < s.size() && std::isspace(static_cast<unsigned char>(s[begin])) != 0) {
        ++begin;
    }

    std::size_t end = s.size();
    while (end > begin && std::isspace(static_cast<unsigned char>(s[end - 1])) != 0) {
        --end;
    }
    return s.substr(begin, end - begin);
}

inline std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, ',')) {
        fields.push_back(trim_copy(field));
    }
    return fields;
}

inline bool is_comment_or_empty(const std::string& line) {
    const std::string trimmed = trim_copy(line);
    return trimmed.empty() || trimmed[0] == '#';
}

} // namespace detail

inline std::vector<std::filesystem::path> findWavefunctionCsvFiles(const std::filesystem::path& folderPath)
{
    std::vector<std::filesystem::path> results;

    if (!std::filesystem::exists(folderPath)) {
        throw std::runtime_error("Folder does not exist: " + folderPath.string());
    }

    if (!std::filesystem::is_directory(folderPath)) {
        throw std::runtime_error("Path is not a directory: " + folderPath.string());
    }

    for (const auto& entry : std::filesystem::recursive_directory_iterator(folderPath)) {
        if (entry.is_regular_file() &&
            entry.path().filename() == "wavefunction.csv") {
            results.push_back(entry.path());
        }
    }

    return results;
}

inline void save_wavefunction_csv(
    const WaveFunction<Dimension::One>& psi,
    const std::filesystem::path& filepath
) {
    if (!psi.grid()) {
        throw std::invalid_argument("save_wavefunction_csv(1D): wavefunction has no grid");
    }

    std::ofstream out(filepath);
    if (!out) {
        throw std::runtime_error("save_wavefunction_csv(1D): cannot open file: " + filepath.string());
    }

    out << "index,x,real,imag,abs2\n";
    out << std::setprecision(17);

    const auto& grid = *psi.grid();
    const auto& vec = psi.vec();
    const int n = static_cast<int>(vec.size());
    const double start = grid.start();
    const double step = grid.step();

    for (int i = 0; i < n; ++i) {
        const double x = start + static_cast<double>(i) * step;
        const std::complex<double> v = vec(i);
        out << i << "," << x << "," << v.real() << "," << v.imag() << "," << std::norm(v) << "\n";
    }
}

inline void save_wavefunction_csv(
    const WaveFunction<Dimension::Two>& psi,
    const std::filesystem::path& filepath
) {
    if (!psi.grid()) {
        throw std::invalid_argument("save_wavefunction_csv(2D): wavefunction has no grid");
    }

    std::ofstream out(filepath);
    if (!out) {
        throw std::runtime_error("save_wavefunction_csv(2D): cannot open file: " + filepath.string());
    }

    out << "ix,iy,x,y,real,imag,abs2\n";
    out << std::setprecision(17);

    const auto& grid = *psi.grid();
    const auto& vec = psi.vec();
    const int nx = static_cast<int>(grid.size_x());
    const int ny = static_cast<int>(grid.size_y());
    const double start_x = grid.start_pos_x();
    const double start_y = grid.start_pos_y();
    const double step_x = grid.step_x();
    const double step_y = grid.step_y();

    for (int i = 0; i < nx; ++i) {
        const double x = start_x + static_cast<double>(i) * step_x;
        for (int j = 0; j < ny; ++j) {
            const double y = start_y + static_cast<double>(j) * step_y;
            const int idx = i * ny + j;
            const std::complex<double> v = vec(idx);
            out << i << "," << j << "," << x << "," << y << ","
                << v.real() << "," << v.imag() << "," << std::norm(v) << "\n";
        }
    }
}

inline std::filesystem::path save_state_energy_history_csv(
    const std::vector<double>& energies,
    const std::filesystem::path& wavefunction_csv_path,
    const std::string& filename = "state_energy_history.csv"
) {
    const std::filesystem::path output_path = wavefunction_csv_path.has_parent_path()
        ? (wavefunction_csv_path.parent_path() / filename)
        : std::filesystem::path(filename);

    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error(
            "save_state_energy_history_csv: cannot open file: " + output_path.string()
        );
    }

    out << "step,energy\n";
    out << std::setprecision(17);
    for (std::size_t i = 0; i < energies.size(); ++i) {
        out << i << "," << energies[i] << "\n";
    }

    return output_path;
}

inline void load_wavefunction_csv(
    WaveFunction<Dimension::One>& psi,
    const std::filesystem::path& filepath
) {
    if (!psi.grid()) {
        throw std::invalid_argument("load_wavefunction_csv(1D): wavefunction has no grid");
    }

    std::ifstream in(filepath);
    if (!in) {
        throw std::runtime_error("load_wavefunction_csv(1D): cannot open file: " + filepath.string());
    }

    Eigen::VectorXcd loaded = Eigen::VectorXcd::Zero(psi.size());
    const int n = static_cast<int>(loaded.size());
    int sequential_index = 0;

    std::string line;
    while (std::getline(in, line)) {
        if (detail::is_comment_or_empty(line)) {
            continue;
        }

        const auto fields = detail::split_csv_line(line);
        if (fields.empty()) {
            continue;
        }

        try {
            if (fields.size() >= 4) {
                const int idx = std::stoi(fields[0]);
                if (idx < 0 || idx >= n) {
                    throw std::out_of_range("index out of range");
                }
                loaded(idx) = std::complex<double>(std::stod(fields[2]), std::stod(fields[3]));
            } else if (fields.size() >= 2) {
                if (sequential_index >= n) {
                    break;
                }
                loaded(sequential_index++) = std::complex<double>(std::stod(fields[0]), std::stod(fields[1]));
            }
        } catch (const std::exception&) {
            // Skip non-data/header lines.
        }
    }

    psi.set_vec(loaded);
}

inline void load_wavefunction_csv(
    WaveFunction<Dimension::Two>& psi,
    const std::filesystem::path& filepath
) {
    if (!psi.grid()) {
        throw std::invalid_argument("load_wavefunction_csv(2D): wavefunction has no grid");
    }

    std::ifstream in(filepath);
    if (!in) {
        throw std::runtime_error("load_wavefunction_csv(2D): cannot open file: " + filepath.string());
    }

    const auto& grid = *psi.grid();
    const int nx = static_cast<int>(grid.size_x());
    const int ny = static_cast<int>(grid.size_y());
    const int n = nx * ny;

    Eigen::VectorXcd loaded = Eigen::VectorXcd::Zero(n);
    int sequential_index = 0;

    std::string line;
    while (std::getline(in, line)) {
        if (detail::is_comment_or_empty(line)) {
            continue;
        }

        const auto fields = detail::split_csv_line(line);
        if (fields.empty()) {
            continue;
        }

        try {
            if (fields.size() >= 6) {
                const int ix = std::stoi(fields[0]);
                const int iy = std::stoi(fields[1]);
                if (ix < 0 || ix >= nx || iy < 0 || iy >= ny) {
                    throw std::out_of_range("index out of range");
                }
                const int idx = ix * ny + iy;
                loaded(idx) = std::complex<double>(std::stod(fields[4]), std::stod(fields[5]));
            } else if (fields.size() >= 4) {
                const int ix = std::stoi(fields[0]);
                const int iy = std::stoi(fields[1]);
                if (ix < 0 || ix >= nx || iy < 0 || iy >= ny) {
                    throw std::out_of_range("index out of range");
                }
                const int idx = ix * ny + iy;
                loaded(idx) = std::complex<double>(std::stod(fields[2]), std::stod(fields[3]));
            } else if (fields.size() >= 2) {
                if (sequential_index >= n) {
                    break;
                }
                loaded(sequential_index++) = std::complex<double>(std::stod(fields[0]), std::stod(fields[1]));
            }
        } catch (const std::exception&) {
            // Skip non-data/header lines.
        }
    }

    psi.set_vec(loaded);
}

} // namespace gpes::io
