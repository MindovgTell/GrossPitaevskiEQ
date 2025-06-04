#ifndef VISUALISATION_HPP
#define VISUALISATION_HPP


#include <vector>
#include <Eigen/Dense>
#include "grid.hpp"
#include "wavefunciton.hpp"
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;


namespace GPES {


void draw(std::vector<double>& vec){
    int size = vec.size();

    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(size, 0, size);

    std::vector<double> x_vec(x.data(), x.data() + x.size());


    plt::figure();
    plt::plot(x_vec, vec, std::string("b-"),{{"label", "data trend"}});
    plt::xlabel("time [s]");
    plt::ylabel("observation [m]");
    plt::legend();
    plt::grid();
    plt::show(); 
}

void draw(Grid<Dimension::One>& grid){
    int size = grid.get_size_of_grid();
    double start = grid.get_start_position();

    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(size, start, -1 * start);
    Eigen::VectorXd& Psi = grid.get_potential();

    std::vector<double> x_vec(x.data(), x.data() + x.size());
    std::vector<double> y_vec(Psi.data(), Psi.data() + Psi.size());

    plt::figure();
    plt::plot(x_vec,y_vec, std::string("b-"),{{"label", "data trend"}});
    plt::xlabel("time [s]");
    plt::ylabel("observation [m]");
    plt::legend();
    plt::grid();
    plt::show(); 
}

void draw(WaveFunction<Dimension::One>& wave){

    int size = wave.get_size_of_grid();
    double start = wave.get_start_position();

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

    int size = grid.get_size_of_grid();
    double start = grid.get_start_position();

    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(size, start, -1 * start);
    Eigen::VectorXd& potential = grid.get_potential();
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

    int size = first_wave.get_size_of_grid();
    double start = first_wave.get_start_position();

    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(size, start, -1 * start);

    Eigen::VectorXd Psi = first_wave.prob();
    Eigen::VectorXd Fin = second_wave.prob();


    std::vector<double> x_vec(x.data(), x.data() + x.size());
    std::vector<double> y_vec(Psi.data(), Psi.data() + Psi.size());
    std::vector<double> z_vec(Fin.data(), Fin.data() + Fin.size());

    plt::figure();
    plt::plot(x_vec,y_vec, std::string("b-"),{{"label", "Initial state"}});
    plt::plot(x_vec,z_vec, std::string("r-"),{{"label", "Final state"}});
    plt::xlabel("");
    plt::ylabel("");
    plt::legend();
    plt::grid();
    plt::show(); 
}

void draw(WaveFunction<Dimension::One>& first_wave,  WaveFunction<Dimension::One>& second_wave, WaveFunction<Dimension::One>& third_wave ){
    //Add support for exceptions with different size of vectors

    int size = first_wave.get_size_of_grid();
    double start = first_wave.get_start_position();

    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(size, start, -1 * start);

    Eigen::VectorXd Psi = first_wave.prob();
    Eigen::VectorXd Fin = second_wave.prob();
    Eigen::VectorXd TF  = third_wave.prob();

    std::vector<double> x_vec(x.data(), x.data() + x.size());
    std::vector<double> y_vec(Psi.data(), Psi.data() + Psi.size());
    std::vector<double> z_vec(Fin.data(), Fin.data() + Fin.size());
    std::vector<double> t_vec(TF.data(), TF.data() + TF.size());

    plt::figure();
    plt::plot(x_vec,y_vec, std::string("b-"),{{"label", "Initial state"}});
    plt::plot(x_vec,z_vec, std::string("r-"),{{"label", "Final state"}});
    plt::plot(x_vec,t_vec, std::string("g-"),{{"label", "Thomas-Fermi limit"}});
    plt::xlabel("x / l$_{\\perp}$",{{"fontsize", "20"}});
    plt::ylabel("$|\\Psi(x)|^2$", {{"fontsize", "20"}});
    plt::legend();
    plt::grid();
    plt::show(); 
}

void draw_energy(std::vector<double>& x, std::vector<double>& vec, std::vector<double>& TM_en){
    std::vector<double> x_vec(x.data(), x.data() + x.size());

    plt::figure();
    plt::plot(x_vec, vec, std::string("b-"),{{"label", "data trend"}});
    plt::plot(x_vec, TM_en, std::string("g--"),{{"label", "TM state energy"}});
    plt::xlabel("time [s]");
    plt::ylabel("observation [m]");
    plt::legend();
    plt::grid();
    plt::show(); 
}

void draw_energy(std::vector<double>& vec_of_energies, std::vector<double>& TM_en){

    std::vector<double> en_len(vec_of_energies.size());
    std::iota(en_len.begin(), en_len.end(), 1);

    std::vector<double> x_vec(en_len.data(), en_len.data() + en_len.size());

    plt::figure();
    plt::plot(x_vec, vec_of_energies, std::string("b-"),{{"label", "data trend"}});
    plt::plot(x_vec, TM_en, std::string("g--"),{{"label", "TM state energy"}});
    plt::xlabel("simulation step");
    plt::ylabel("final state energy, $\\hbar \\omega$");
    plt::legend();
    plt::grid();
    plt::show(); 
}

void draw_energy(std::vector<double>& vec_of_energies){

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


void draw3(Eigen::VectorXd& x, Eigen::VectorXd& Psi,  Eigen::VectorXd& Fin, Eigen::VectorXd& V){
    std::vector<double> x_vec(x.data(), x.data() + x.size());
    std::vector<double> y_vec(Psi.data(), Psi.data() + Psi.size());
    std::vector<double> z_vec(Fin.data(), Fin.data() + Fin.size());
    std::vector<double> t_vec(V.data(), V.data() + V.size());

    plt::figure();
    plt::plot(x_vec,y_vec, std::string("b-"),{{"label", "Initial state"}});
    plt::plot(x_vec,z_vec, std::string("r-"),{{"label", "Final state"}});
    plt::plot(x_vec,t_vec, std::string("g-"),{{"label", "Potential well"}});
    plt::xlabel("");
    plt::ylabel("");
    plt::legend();
    plt::grid();
    plt::show(); 
}



//--------//-----------------------------------//--------//
//--------//Additional drawing functions for 2D//--------//
//--------//-----------------------------------//--------//


void draw(Grid<Dimension::Two>& grid) {
    int size_x = grid.get_size_of_grid_x();
    int size_y = grid.get_size_of_grid_y();
    double start_x = grid.get_start_position_x();
    double start_y = grid.get_start_position_y();

    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(size_x, start_x, -1 * start_x);
    Eigen::MatrixXd& Mat = grid.get_potential();

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

std::vector<std::vector<double>> Eigen_to_vector2D(Grid<Dimension::Two>& mat) {
    std::vector<std::vector<double>> result(mat.get_size_of_grid_x(), std::vector<double>(mat.get_size_of_grid_y()));
    for (int i = 0; i < mat.get_size_of_grid_x(); ++i)
        for (int j = 0; j < mat.get_size_of_grid_y(); ++j)
            result[i][j] = mat(i, j);
    return result;
}

void heatmap(Grid<Dimension::Two>& mat){

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
    std::vector<std::vector<double>> result(wave.get_size_of_grid_x(), std::vector<double>(wave.get_size_of_grid_y()));
    for (int i = 0; i < wave.get_size_of_grid_x(); ++i)
        for (int j = 0; j < wave.get_size_of_grid_y(); ++j){
            int index = wave.get_index(i, j);
            result[i][j] = wave.prob(index);
        }
    return result;
}

void heatmap(WaveFunction<Dimension::Two>& wave, std::string name){

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



std::vector<std::vector<double>> Eigen_to_vector2D(const Eigen::MatrixXd& mat) {
    std::vector<std::vector<double>> result(mat.rows(), std::vector<double>(mat.cols()));
    for (int i = 0; i < mat.rows(); ++i)
        for (int j = 0; j < mat.cols(); ++j)
            result[i][j] = mat(i, j);
    return result;
}

void heatmap(const Eigen::MatrixXd& mat){

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

        plt::show();
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }

}

}

#endif