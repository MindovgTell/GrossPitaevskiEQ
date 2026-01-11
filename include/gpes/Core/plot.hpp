#pragma once

#include <vector>
#include <Eigen/Dense>
#include <memory>
#include "grid/grid.hpp"
#include "wavefunction/wavefunction.hpp"
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wsign-compare"
#endif
#include "matplotlibcpp.h"
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


} // namespace gpes
