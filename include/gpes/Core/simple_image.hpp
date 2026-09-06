#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>

namespace gpes::image {

namespace detail {

inline void write_u16(std::ostream& out, std::uint16_t value) {
    const char bytes[2] = {
        static_cast<char>(value & 0xffu),
        static_cast<char>((value >> 8u) & 0xffu)};
    out.write(bytes, 2);
}

inline void write_u32(std::ostream& out, std::uint32_t value) {
    const char bytes[4] = {
        static_cast<char>(value & 0xffu),
        static_cast<char>((value >> 8u) & 0xffu),
        static_cast<char>((value >> 16u) & 0xffu),
        static_cast<char>((value >> 24u) & 0xffu)};
    out.write(bytes, 4);
}

inline std::array<unsigned char, 3> inferno_like(double x) {
    x = std::clamp(x, 0.0, 1.0);
    const std::array<std::array<double, 3>, 5> colors{{
        {0.0, 0.0, 4.0}, {87.0, 15.0, 109.0}, {187.0, 55.0, 84.0},
        {249.0, 142.0, 9.0}, {252.0, 255.0, 164.0}}};
    const double scaled = x * static_cast<double>(colors.size() - 1);
    const std::size_t lo = std::min<std::size_t>(
        static_cast<std::size_t>(scaled), colors.size() - 2);
    const double t = scaled - static_cast<double>(lo);
    std::array<unsigned char, 3> rgb{};
    for (int c = 0; c < 3; ++c) {
        rgb[static_cast<std::size_t>(c)] = static_cast<unsigned char>(std::lround(
            (1.0 - t) * colors[lo][static_cast<std::size_t>(c)]
            + t * colors[lo + 1][static_cast<std::size_t>(c)]));
    }
    return rgb;
}

inline std::string xml_escape(const std::string& value) {
    std::string escaped;
    for (char c : value) {
        switch (c) {
            case '&': escaped += "&amp;"; break;
            case '<': escaped += "&lt;"; break;
            case '>': escaped += "&gt;"; break;
            case '"': escaped += "&quot;"; break;
            default: escaped += c; break;
        }
    }
    return escaped;
}

} // namespace detail

inline void save_scalar_bmp(
    const Eigen::VectorXd& values,
    int nx,
    int ny,
    const std::filesystem::path& path,
    bool logarithmic = false)
{
    if (nx <= 0 || ny <= 0 || values.size() != nx * ny) {
        throw std::invalid_argument("save_scalar_bmp: invalid image dimensions");
    }
    if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());

    double maximum = 0.0;
    for (Eigen::Index i = 0; i < values.size(); ++i) {
        if (std::isfinite(values(i))) maximum = std::max(maximum, values(i));
    }
    const std::uint32_t row_size = static_cast<std::uint32_t>((3 * ny + 3) & ~3);
    const std::uint32_t pixel_bytes = row_size * static_cast<std::uint32_t>(nx);

    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("save_scalar_bmp: cannot open " + path.string());
    out.write("BM", 2);
    detail::write_u32(out, 54u + pixel_bytes);
    detail::write_u32(out, 0u);
    detail::write_u32(out, 54u);
    detail::write_u32(out, 40u);
    detail::write_u32(out, static_cast<std::uint32_t>(ny));
    detail::write_u32(out, static_cast<std::uint32_t>(nx));
    detail::write_u16(out, 1u);
    detail::write_u16(out, 24u);
    detail::write_u32(out, 0u);
    detail::write_u32(out, pixel_bytes);
    detail::write_u32(out, 2835u);
    detail::write_u32(out, 2835u);
    detail::write_u32(out, 0u);
    detail::write_u32(out, 0u);

    std::vector<char> row(row_size, 0);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double normalized = maximum > 0.0 ? std::max(0.0, values(i * ny + j)) / maximum : 0.0;
            if (logarithmic) normalized = std::log1p(99.0 * normalized) / std::log(100.0);
            const auto rgb = detail::inferno_like(normalized);
            row[3 * j] = static_cast<char>(rgb[2]);
            row[3 * j + 1] = static_cast<char>(rgb[1]);
            row[3 * j + 2] = static_cast<char>(rgb[0]);
        }
        out.write(row.data(), static_cast<std::streamsize>(row.size()));
    }
}

inline void save_line_svg(
    const std::vector<double>& values,
    const std::filesystem::path& path,
    const std::string& title,
    const std::string& y_label)
{
    if (values.empty()) return;
    if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
    const auto [min_it, max_it] = std::minmax_element(values.begin(), values.end());
    const double ymin = *min_it;
    const double ymax = *max_it;
    const double span = std::max(ymax - ymin, 1e-15);
    std::ofstream out(path);
    if (!out) throw std::runtime_error("save_line_svg: cannot open " + path.string());
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"650\" viewBox=\"0 0 1000 650\">\n"
        << "<rect width=\"1000\" height=\"650\" fill=\"white\"/><g stroke=\"#222\" fill=\"none\">"
        << "<path d=\"M90 40V570H970\"/></g>\n"
        << "<text x=\"500\" y=\"28\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"20\">"
        << detail::xml_escape(title) << "</text>\n"
        << "<text x=\"25\" y=\"305\" transform=\"rotate(-90 25 305)\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"16\">"
        << detail::xml_escape(y_label) << "</text><polyline fill=\"none\" stroke=\"#16697a\" stroke-width=\"2\" points=\"";
    for (std::size_t i = 0; i < values.size(); ++i) {
        const double x = values.size() == 1 ? 90.0 : 90.0 + 880.0 * i / (values.size() - 1);
        const double y = 570.0 - 530.0 * (values[i] - ymin) / span;
        out << x << ',' << y << ' ';
    }
    out << "\"/><text x=\"90\" y=\"610\" font-family=\"sans-serif\" font-size=\"14\">check index</text>"
        << "<text x=\"95\" y=\"60\" font-family=\"monospace\" font-size=\"13\">max=" << std::setprecision(8) << ymax
        << "</text><text x=\"95\" y=\"565\" font-family=\"monospace\" font-size=\"13\">min=" << ymin
        << "</text></svg>\n";
}

inline void save_spectrum_svg(
    const std::vector<std::complex<double>>& omega,
    const std::filesystem::path& path)
{
    if (omega.empty()) return;
    if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
    double bound = 1e-12;
    for (const auto& value : omega) bound = std::max({bound, std::abs(value.real()), std::abs(value.imag())});
    std::ofstream out(path);
    if (!out) throw std::runtime_error("save_spectrum_svg: cannot open " + path.string());
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"650\" viewBox=\"0 0 1000 650\">"
        << "<rect width=\"1000\" height=\"650\" fill=\"white\"/><path d=\"M90 305H970M90 40V570\" stroke=\"#333\"/>"
        << "<text x=\"500\" y=\"28\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"20\">BdG excitation spectrum</text>"
        << "<text x=\"500\" y=\"625\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"16\">mode index</text>";
    for (std::size_t i = 0; i < omega.size(); ++i) {
        const double x = omega.size() == 1 ? 90.0 : 90.0 + 880.0 * i / (omega.size() - 1);
        const double yr = 305.0 - 265.0 * omega[i].real() / bound;
        const double yi = 305.0 - 265.0 * omega[i].imag() / bound;
        out << "<circle cx=\"" << x << "\" cy=\"" << yr << "\" r=\"4\" fill=\"#16697a\"/>"
            << "<circle cx=\"" << x << "\" cy=\"" << yi << "\" r=\"3\" fill=\"#c33c54\"/>";
    }
    out << "<text x=\"110\" y=\"65\" font-family=\"sans-serif\" font-size=\"14\" fill=\"#16697a\">Re omega</text>"
        << "<text x=\"110\" y=\"85\" font-family=\"sans-serif\" font-size=\"14\" fill=\"#c33c54\">Im omega</text></svg>\n";
}

} // namespace gpes::image
