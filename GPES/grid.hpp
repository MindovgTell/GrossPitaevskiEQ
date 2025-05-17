#ifndef GRID_HPP
#define GRID_HPP

#include "definitions.hpp"
#include <Eigen/Dense>


namespace GPES {

//********************************/***********/********************************//
//                                                                             //
//********************************/1DFunctions/********************************//
//                                                                             //
//********************************/***********/********************************//

template <Dimension dim>
class Grid;


template<>
class Grid<Dimension::One> {

private: 
    unsigned int _size;
    double _step, _start, _omega;

    Eigen::VectorXd _potential;
    
    //initialize grid size
    void init_grid();

public:
    //Construction
    Grid(unsigned int number_of_nodes, double step);
    // Delete copy and move constructors
    Grid(const GPES::Grid<Dimension::One>&) = delete;
    Grid(GPES::Grid<Dimension::One> &&)     = delete;


    //Getters
    unsigned int get_size_of_grid() {return _size; }
    double get_start_position() { return _start; }
    double get_step_size() { return _step; }
    double get_omega() { return _omega; }

    Eigen::VectorXd get_potential(){ return _potential;}



    void set_harmonic_potential(double omega);


    // //Harmonic potential terms
    // Eigen::VectorXd create_harmonic_potential();
    // Eigen::VectorXd create_potential();



    //operator() overloading

    double& operator()(int i) {
        return _potential(i);
    }

    double operator()(int i) const {
        return _potential(i);
    }
};

inline Grid<Dimension::One>::Grid(unsigned int number_of_nodes, double start): _size(number_of_nodes), _start(start) {
    //Setting step
    _step = (-2 * _start) / _size;
    _omega = 0;
    init_grid();
}

void inline Grid<Dimension::One>::init_grid(){
    _potential = Eigen::VectorXd::Zero(_size);
}

void inline Grid<Dimension::One>::set_harmonic_potential(double omega){
    Eigen::VectorXd V(_size);
    V.setZero();
    for(int i = 0; i != _size - 0; ++i){
        double x = _start + i * _step;
        V(i) = 0.5 * (std::pow(omega * x,2.));
    }

    _omega = omega;

    _potential = V;
}



// Eigen::VectorXd CrankNicolson::create_harmonic_potential_1D(){
//     Eigen::VectorXd V(this->m_size);
//     V.setZero();
//     double x = this->_start;
//     for(int i = 0; i < this->m_size; ++i){
//         V(i) = 0.5 * std::pow(x,2.);
//         x += this->step;
//     }
//     return V;
// }

// Eigen::VectorXd CrankNicolson::create_potential_1D(){
//     Eigen::VectorXd V(this->m_size);
//     V.setZero();
//     V(0) = V_0;
//     V(this->m_size - 1) = V_0;
//     return V;
// }


//********************************/***********/********************************//
//                                                                             //
//********************************/2DFunctions/********************************//
//                                                                             //
//********************************/***********/********************************//

template<>
class Grid<Dimension::Two> {

private: 
    int _size_x, _size_y;
    double _step_x, _step_y, _start_x, _start_y;
    
    //harmonic potential frequencies
    double _omega_x, _omega_y, _lambda;

    Eigen::MatrixXd _potential;

    //initialize grid
    void init_grid();

public:
    //Constructor 
    Grid(int N_x, int N_y, double start_x, double start_y);
    // Delete copy and move constructors
    Grid(const GPES::Grid<Dimension::Two>&) = delete;
    Grid(GPES::Grid<Dimension::Two> &&)     = delete;

    //Getters
    unsigned int get_size_of_grid_x(){ return _size_x; }
    unsigned int get_size_of_grid_y(){ return _size_y; }

    double get_start_position_x(){ return _start_x; }
    double get_start_position_y(){ return _start_y; }

    double get_step_size_x(){ return _step_x; }
    double get_step_size_y(){ return _step_y; }

    Eigen::MatrixXd get_potential() {return _potential; }


    void set_harmonic_potential(double omega_x, double omega_y);

    //operator() overloading

    double& operator()(int i, int j) {
        return _potential(i, j);
    }

    double operator()(int i, int j) const {
        return _potential(i, j);
    }

};

inline Grid<Dimension::Two>::Grid(int N_x, int N_y, double start_x, double start_y):
    _size_x(N_x), _size_y(N_y), _start_x(start_x), _start_y(start_y) {
        //Setting steps
        _step_x = (-2 * _start_x) / _size_x;
        _step_y = (-2 * _start_y) / _size_y;

        _omega_x = 0;
        _omega_y = 0;
        _lambda = 0;

        this->init_grid();
}

void inline Grid<Dimension::Two>::init_grid(){
    // Eigen::MatrixXd U(_size_x, _size_y);
    // U.setZero();
    // this->_potential = U;

    _potential = Eigen::MatrixXd::Zero(_size_x,_size_y);
}

void inline Grid<Dimension::Two>::set_harmonic_potential(double omega_x, double omega_y){
    Eigen::MatrixXd V(_size_x, _size_y);
    V.setZero();

    this->_omega_x = omega_x;
    this->_omega_y = omega_y;

    for(int i = 1; i != _size_x - 1; ++i){
        double x = this->_start_x + i * this->_step_x;
        for(int j = 1; j != _size_y-1; ++j){
            double y = this->_start_y + j * this->_step_y;
            V(i,j) = 0.5 * (std::pow(omega_x * x,2.) + std::pow(omega_y * y,2.));
        }
    }

    this->_lambda = this->_omega_x / this->_omega_y;

    this->_potential = V;
}


} // namespace GPES

#endif


// //Code example
// GPES::Grid<Dimension::Two> grid(100,100,-20.,-20.);
// grid.set_harmonic_potential(100,100);