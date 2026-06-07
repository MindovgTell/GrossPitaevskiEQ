#pragma once

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <memory>
#include <numeric>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "grid/grid.hpp"
#include "wavefunction/wavefunction.hpp"
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wsign-compare"
#endif
#include "third_party/matplotlibcpp.h"
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
#include "definitions.hpp"

namespace plt = matplotlibcpp;
namespace gpes::plot {

    inline void draw(Grid<Dimension::One>& grid, std::string xlabel = "x",  std::string ylabel = "V");

    inline void draw(std::shared_ptr<const Grid<Dimension::One>> grid, std::string xlabel = "x",  std::string ylabel = "V"){
        if (!grid) {
            std::cerr << "draw: grid is null\n";
            return;
        }
        try {
            int size = grid->size();
            if (size <= 0) {
                std::cerr << "draw: grid size must be positive\n";
                return;
            }
            double start = grid->start();
            Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(size, start, -1 * start);
            const Eigen::VectorXd& Psi = grid->potential();
            std::vector<double> x_vec(x.data(), x.data() + x.size());
            std::vector<double> y_vec(Psi.data(), Psi.data() + Psi.size());
            plt::figure();
            plt::plot(x_vec,y_vec, std::string("b-"),{{"label", "data trend"}});
            plt::xlabel(xlabel);
            plt::ylabel(ylabel);
            plt::legend();
            plt::grid();
            plt::show();
        } catch (const std::exception& e) {
            std::cerr << "draw: plotting failed: " << e.what() << "\n";
        }
    }

    void draw(WaveFunction<Dimension::One>& wave, std::string xlabel = "x",  std::string ylabel = "V"){
        int size = wave.size();
        double start = wave.grid()->start();
        Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(size, start, -1 * start);
        Eigen::VectorXd Psi = wave.prob();
        std::vector<double> x_vec(x.data(), x.data() + x.size());
        std::vector<double> y_vec(Psi.data(), Psi.data() + Psi.size());
        plt::figure();
        plt::plot(x_vec,y_vec, std::string("b-"),{{"label", "Initial state"}});
        plt::xlabel("");
        plt::ylabel("");
        plt::legend();
        plt::grid();
        plt::show();
    }

    void draw(Grid<Dimension::One>& grid,  WaveFunction<Dimension::One>& wave){
        int size = grid.size();
        double start = grid.start();

        Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(size, start, -1 * start);
        const Eigen::VectorXd& potential = grid.potential();
        Eigen::VectorXd Psi = wave.prob();


        std::vector<double> x_vec(x.data(), x.data() + x.size());
        std::vector<double> y_vec(potential.data(), potential.data() + potential.size());
        std::vector<double> z_vec(Psi.data(), Psi.data() + Psi.size());

        plt::figure();
        plt::plot(x_vec,y_vec, std::string("b-"),{{"label", "Initial state"}});
        plt::plot(x_vec,z_vec, std::string("r-"),{{"label", "Final state"}});
        plt::xlabel("");
        plt::ylabel("");
        plt::legend();
        plt::grid();
        plt::show();
    }

void draw(WaveFunction<Dimension::One>& first_wave,  WaveFunction<Dimension::One>& second_wave){
    //Add support for exceptions with different size of vectors

    int size = first_wave.grid()->size();
    double start = first_wave.grid()->start();

    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(size, start, -1 * start);

    Eigen::VectorXd Psi = first_wave.prob();
    Eigen::VectorXd Fin = second_wave.prob();


    std::vector<double> x_vec(x.data(), x.data() + x.size());
    std::vector<double> y_vec(Psi.data(), Psi.data() + Psi.size());
    std::vector<double> z_vec(Fin.data(), Fin.data() + Fin.size());

    plt::figure();
    plt::plot(x_vec,y_vec, std::string("b-"),{{"label", "First state"}});
    plt::plot(x_vec,z_vec, std::string("r-"),{{"label", "Second state"}});
    plt::xlabel("");
    plt::ylabel("");
    plt::legend();
    plt::grid();
    plt::show();
}



void draw_energy(const std::vector<double>& vec_of_energies){

    std::vector<double> en_len(vec_of_energies.size());
    std::iota(en_len.begin(), en_len.end(), 1);

    std::vector<double> x_vec(en_len.data(), en_len.data() + en_len.size());

    plt::figure();
    plt::plot(x_vec, vec_of_energies, std::string("b-"),{{"label", "data trend"}});
    plt::xlabel("time [s]");
    plt::ylabel("observation [m]");
    plt::legend();
    plt::grid();
    plt::show();
}

inline void save_energy_plot(
    const std::vector<double>& vec_of_energies,
    const std::filesystem::path& output_path
) {
    if (vec_of_energies.empty()) {
        std::cerr << "save_energy_plot: empty energy vector\n";
        return;
    }

    std::vector<double> step_idx(vec_of_energies.size());
    std::iota(step_idx.begin(), step_idx.end(), 0.0);

    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    try {
        plt::figure();
        plt::plot(step_idx, vec_of_energies, std::string("b-"), {{"label", "state energy"}});
        plt::xlabel("step");
        plt::ylabel("energy");
        plt::legend();
        plt::grid();
        plt::savefig(output_path.string());
        plt::close();
    } catch (const std::exception& e) {
        std::cerr << "save_energy_plot: plotting failed: " << e.what() << "\n";
        try {
            plt::close();
        } catch (...) {
        }
    }
}

void draw_energy_relative_fluctuations(
    const std::vector<double>& vec_of_energies,
    bool in_percent = true
){
    if (vec_of_energies.empty()) {
        std::cerr << "draw_energy_relative_fluctuations: empty energy vector\n";
        return;
    }

    std::vector<double> step_idx(vec_of_energies.size());
    std::iota(step_idx.begin(), step_idx.end(), 1.0);

    const double e_ref = vec_of_energies.front();
    const double ref_abs = std::max(std::abs(e_ref), 1e-14);

    std::vector<double> rel_fluct(vec_of_energies.size(), 0.0);
    for (std::size_t i = 0; i < vec_of_energies.size(); ++i) {
        rel_fluct[i] = (vec_of_energies[i] - e_ref) / ref_abs;
        if (in_percent) {
            rel_fluct[i] *= 100.0;
        }
    }

    double max_abs = 0.0;
    for (double v : rel_fluct) {
        max_abs = std::max(max_abs, std::abs(v));
    }

    plt::figure();
    plt::plot(step_idx, rel_fluct, std::string("r-"), {{"label", "relative energy fluctuation"}});
    plt::xlabel("step");
    plt::ylabel(in_percent ? "dE/E0 [%]" : "dE/E0");
    if (max_abs > 0.0) {
        plt::ylim(-1.1 * max_abs, 1.1 * max_abs);
    }
    plt::legend();
    plt::grid();
    plt::show();
}


// Two dimensional functions
void draw(Grid<Dimension::Two>& grid);
void draw(std::shared_ptr<const Grid<Dimension::Two>> grid) {
    int size_x = grid->size_x();
    int size_y = grid->size_y();
    double start_x = grid->start_pos_x();
    double start_y = grid->start_pos_y();

    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(size_x, start_x, -1 * start_x);
    const Eigen::MatrixXd& Mat = grid->potential();

    Eigen::VectorXd V_x = Mat.row(size_y / 2);
    Eigen::VectorXd V_y = Mat.col(size_x / 2);

    std::vector<double> x_vec(x.data(), x.data() + x.size());

    std::vector<double> y_vec(V_x.data(), V_x.data() + V_x.size());
    std::vector<double> z_vec(V_y.data(), V_y.data() + V_y.size());

    try {
        // plt::backend("Agg");

        // First subplot (left)
        plt::figure();  // (rows, cols, index)
        plt::plot(x_vec, y_vec, "r-", {{"label", "Quadratic"}});
        // plt::title("Quadratic Growth");
        plt::xlabel("X Values");
        plt::ylabel("Y Values");
        plt::grid();
        plt::legend();

        // Second subplot (right)
        plt::figure();
        plt::plot(x_vec, z_vec, "b-", {{"label", "Custom Pattern"}});
        // plt::title("Custom Growth Pattern");
        plt::xlabel("X Values");
        plt::ylabel("Y Values");
        plt::grid();
        plt::legend();

        // // Adjust spacing between subplots
        // plt::tight_layout();

        // Display or save
        plt::show();
    }
    catch (const std::exception& e) {
        std::cout << "Caught exception: " << e.what() << std::endl;
    }
}


std::vector<std::vector<double>> Eigen_to_vector2D(std::shared_ptr<const Grid<Dimension::Two>> mat) {
    std::vector<std::vector<double>> result(mat->size_x(), std::vector<double>(mat->size_y()));
    for (int i = 0; i < mat->size_x(); ++i)
        for (int j = 0; j < mat->size_y(); ++j)
            result[i][j] = mat->potential()(i, j);
    return result;
}

void heatmap(std::shared_ptr<const Grid<Dimension::Two>> mat){

    try {
        std::vector< std::vector<double> > heatmap_data = Eigen_to_vector2D(mat);

        // Optional: verify data
        if (heatmap_data.empty() || heatmap_data[0].empty()) {
            throw std::runtime_error("Data is empty.");
        }

        // Plot
        plt::imshow(heatmap_data); //  , {{"cmap", "'plasma'"}, {"interpolation", "'nearest'"}}
        plt::colorbar();
        plt::title("Eigen Heatmap");
        // plt::xlim(mat.get_start_position_x(), -1. * mat.get_start_position_x());
        // plt::ylim(mat.get_start_position_y(), -1. * mat.get_start_position_y());

        plt::xlabel("X");
        plt::ylabel("Y");
        plt::legend();
        plt::show();
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }

}





std::vector<std::vector<double>> Eigen_to_vector2D(WaveFunction<Dimension::Two>& wave) {
    std::vector<std::vector<double>> result(wave.grid()->size_x(), std::vector<double>(wave.grid()->size_y()));
    for (int i = 0; i < wave.grid()->size_x(); ++i)
        for (int j = 0; j < wave.grid()->size_y(); ++j){
            int index = wave.get_index(i, j);
            result[i][j] = wave.prob(index);
        }
    return result;
}

void heatmap(WaveFunction<Dimension::Two>& wave, std::string name = "Wavefunction"){
    try {
        std::vector< std::vector<double> > heatmap_data = Eigen_to_vector2D(wave);

        // Optional: verify data
        if (heatmap_data.empty() || heatmap_data[0].empty()) {
            throw std::runtime_error("Data is empty.");
        }

        // Plot
        plt::imshow(heatmap_data); //  , {{"cmap", "'plasma'"}, {"interpolation", "'nearest'"}}
        plt::colorbar();
        plt::title(name);
        // plt::xlim(wave.get_start_position_x(), -1. * wave.get_start_position_x());
        // plt::ylim(wave.get_start_position_y(), -1. * wave.get_start_position_y());
        plt::xlabel("X", {{"fontsize", "18"}});
        plt::ylabel("Y", {{"fontsize", "18"}});
        plt::show();
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
}

inline void draw_slices(WaveFunction<Dimension::Two>& wave, const std::string& name = "Wavefunction") {
    try {
        const int size_x = wave.grid()->size_x();
        const int size_y = wave.grid()->size_y();
        const double start_x = wave.grid()->start_pos_x();
        const double start_y = wave.grid()->start_pos_y();
        const double step_x = wave.grid()->step_x();
        const double step_y = wave.grid()->step_y();

        std::vector<double> x_axis(size_x);
        std::vector<double> y_axis(size_y);
        for (int i = 0; i < size_x; ++i) {
            x_axis[i] = start_x + i * step_x;
        }
        for (int j = 0; j < size_y; ++j) {
            y_axis[j] = start_y + j * step_y;
        }

        const Eigen::VectorXd x_slice = wave.get_x_slice();
        const Eigen::VectorXd y_slice = wave.get_y_slice();
        std::vector<double> x_slice_vec(x_slice.data(), x_slice.data() + x_slice.size());
        std::vector<double> y_slice_vec(y_slice.data(), y_slice.data() + y_slice.size());

        plt::figure();
        plt::plot(x_axis, x_slice_vec, "b-", {{"label", "center slice along x"}});
        plt::plot(y_axis, y_slice_vec, "r-", {{"label", "center slice along y"}});
        plt::title(name + " slices");
        plt::xlabel("position", {{"fontsize", "18"}});
        plt::ylabel("|Psi|^2", {{"fontsize", "18"}});
        plt::grid();
        plt::legend();
        plt::show();
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
}

} // namespace gpes
