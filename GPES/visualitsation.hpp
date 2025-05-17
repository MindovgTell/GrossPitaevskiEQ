#ifndef VISUALISATION_HPP
#define VISUALISATION_HPP


#include <Eigen/Dense>
#include "grid.hpp"
#include "matplotlibcpp.h"


namespace GPES {


void draw(Grid<Dimension::One>& grid){
    int size = grid.get_size_of_grid();
    double start = grid.get_

    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(size, start, -1 * start);


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


void draw1(Eigen::VectorXd& x, Eigen::VectorXd& Psi){
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

void draw2(Eigen::VectorXd& x, Eigen::VectorXd& Psi,  Eigen::VectorXd& Fin){
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