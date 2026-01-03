#ifndef GRID_HPP
#define GRID_HPP

#include <iostream>
#include <Eigen/Dense>
#include "definitions.hpp"




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
    double _step, _start, _omega, _omega_t;

    Eigen::VectorXd _potential;
    
    //initialize grid size
    void init_grid();

public:
    //Construction
    Grid(unsigned int number_of_nodes, double start);
    // Delete copy and move constructors
    Grid(const GPES::Grid<Dimension::One>&) = delete;
    Grid(GPES::Grid<Dimension::One> &&)     = delete;


    //Getters
    unsigned int size() {return _size; }
    const unsigned int size() const {return _size; }

    double start_pos() { return _start; }
    const double start_pos() const { return _start; }

    double step() { return _step; }
    const double step() const { return _step; }

    double get_omega() { return _omega; }
    double get_omega_t() { return _omega_t; }
    const Eigen::VectorXd& potential() const { return _potential;}
    
    //Setters
    void set_harmonic_potential(double omega);
    void set_transverse(double omega) { _omega_t = omega;}

    //operator() overloading
    double& operator()(int i) { return _potential(i); }
    double operator()(int i) const { return _potential(i); }
};

Grid<Dimension::One>::Grid(unsigned int number_of_nodes, double start): _size(number_of_nodes), _start(start) {
    //Setting step
    _step = (-2 * _start) / _size;
    _omega = 0;
    _omega_t = 0;
    init_grid();
}

void Grid<Dimension::One>::init_grid(){
    _potential = Eigen::VectorXd::Zero(_size);
}

void Grid<Dimension::One>::set_harmonic_potential(double omega){
    Eigen::VectorXd V(_size);
    V.setZero();
    for(int i = 0; i != _size - 0; ++i){
        double x = _start + i * _step;
        V(i) = 0.5 * (std::pow(omega * x,2.));
    }
    _omega = omega;
    _potential = V;
}


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
    double _omega_x, _omega_y, _omega_z;

    Eigen::MatrixXd _potential;

    //initialize grid
    void init_grid();

public:
    //Constructor 
    Grid(int number_of_nodes_x, int number_of_nodes_y, double start_x, double start_y);
    Grid();
    // Delete copy and move constructors
    // Grid(const GPES::Grid<Dimension::Two>&) = delete;
    // Grid(GPES::Grid<Dimension::Two> &&) 

    //Getters
    unsigned int size_x(){ return _size_x; }
    unsigned int size_y(){ return _size_y; }
    double start_pos_x(){ return _start_x; }
    double start_pos_y(){ return _start_y; }
    double step_x(){ return _step_x; }
    double step_y(){ return _step_y; }
    double omega_x() {return _omega_x; }
    double omega_y() {return _omega_y; }
    double omega_z() {return _omega_z; }
    const Eigen::MatrixXd& potential() const {return _potential;}

    void set_harmonic_potential(double omega_x, double omega_y);

    //operator() overloading
    double& operator()(int i, int j) { return _potential(i, j); }
    double operator()(int i, int j) const { return _potential(i, j); }

    void set_z_freq(double omega_z) {_omega_z = omega_z;}

};

Grid<Dimension::Two>::Grid(int number_of_nodes_x, int number_of_nodes_y, double start_x, double start_y):
    _size_x(number_of_nodes_x), _size_y(number_of_nodes_y), _start_x(start_x), _start_y(start_y) {
        //Setting steps
        _step_x = std::abs((2 * _start_x) / _size_x);
        _step_y = std::abs((2 * _start_y) / _size_y);

        _omega_x = 0;
        _omega_y = 0;
        _omega_z = 0;

        init_grid();
}

void Grid<Dimension::Two>::init_grid(){
    _potential = Eigen::MatrixXd::Zero(_size_x,_size_y);
}

void Grid<Dimension::Two>::set_harmonic_potential(double omega_x, double omega_y){
    Eigen::MatrixXd V(_size_x, _size_y);
    V.setZero();

    _omega_x = omega_x;
    _omega_y = omega_y;

    for(int i = 0; i != _size_x; ++i){
        double x = _start_x + i * _step_x;
        for(int j = 0; j != _size_y; ++j){
            double y = _start_y + j * _step_y;
            V(i,j) = 0.5 * (std::pow(omega_x * x,2.) + std::pow(omega_y * y,2.));
        }
    }
    _potential = V;
}



} // namespace GPES

#endif


// //Code example
// GPES::Grid<Dimension::Two> grid(100,100,-20.,-20.);
// grid.set_harmonic_potential(100,100);