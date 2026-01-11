#ifndef GRID_HPP
#define GRID_HPP

#include <iostream>
#include <Eigen/Dense>
#include "Core/definitions.hpp"




namespace gpes {

//********************************/***********/********************************//
//                                                                             //
//********************************/1DFunctions/********************************//
//                                                                             //
//********************************/***********/********************************//

template <Dimension dim>
class Grid;


template<>
class Grid<Dimension::One> {
public:
    using VecType = Eigen::VectorXd; // Type for potential, neccesery to each grid object
    static constexpr Dimension Dim = Dimension::One;

private: 
    size_t size_;
    double step_, start_, omega_, omega_t_;

    VecType potential_;

    //initialize grid size
    void init_grid() { potential_ = VecType::Zero(size_); }

public:
    //Construction
    Grid(unsigned int number_of_nodes, double start): size_(number_of_nodes), start_(start) {
        //Setting step
        step_ = (-2 * start_) / size_;
        omega_ = 0;
        omega_t_ = 0;
        init_grid();
    }
    Grid(unsigned int number_of_nodes, double start, double omega, double omega_t): 
        size_(number_of_nodes), 
        start_(start),
        omega_(omega),
        omega_t_(omega_t)
    {
        //Setting step
        step_ = (-2 * start_) / size_;
        set_harmonic_potential(5);
    }
    // Delete copy and move constructors
    Grid(const Grid&) = delete;
    Grid& operator=(Grid&) = delete;

    Grid(Grid&&)     = delete;
    Grid& operator=(Grid&&) = delete;

    //Getters
    size_t size() {return size_; }
    size_t size() const {return size_; }

    double start() { return start_; }
    double start() const { return start_; }

    double step() { return step_; }
    double step() const { return step_; }

    double omega() { return omega_; }
    double omega() const { return omega_; }

    double omega_t() { return omega_t_; }
    double omega_t() const { return omega_t_; }

    const VecType& potential() const { return potential_;}
    
    //Setters
    void set_harmonic_potential(double omega) {
        VecType V(size_);
        V.setZero();
        for(size_t i = 0; i != size_; ++i){
            double x = start_ + i * step_;
            V(i) = 0.5 * (std::pow(omega * x,2.));
        }
        omega_ = omega;
        potential_ = V;
    }

    void set_transverse(double omega) { omega_t_ = omega;}

    //operator() overloading
    double& operator()(int i) { return potential_(i); }
    double operator()(int i) const { return potential_(i); }
};

//********************************/***********/********************************//
//                                                                             //
//********************************/2DFunctions/********************************//
//                                                                             //
//********************************/***********/********************************//

template<>
class Grid<Dimension::Two> {
public:

    using VecType = Eigen::VectorXd; // Type for potential, neccesery to each grid object
    static constexpr Dimension Dim = Dimension::Two;

private: 
    int size_x_, size_y_;
    double step_x_, step_y_, start_x_, start_y_;
    
    //harmonic potential frequencies
    double omega_x_, omega_y_, omega_z_;

    Eigen::MatrixXd potential_;

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
    unsigned int size_x(){ return size_x_; }
    unsigned int size_x() const { return size_x_; }
    unsigned int size_y(){ return size_y_; }
    unsigned int size_y() const { return size_y_; }

    double start_pos_x(){ return start_x_; }
    double start_pos_x() const { return start_x_; }
    double start_pos_y(){ return start_y_; }
    double start_pos_y() const { return start_y_; }

    double step_x(){ return step_x_; }
    double step_x() const { return step_x_; }
    double step_y(){ return step_y_; }
    double step_y() const { return step_y_; }

    double omega_x() {return omega_x_; }
    double omega_x() const {return omega_x_; }
    double omega_y() {return omega_y_; }
    double omega_y() const {return omega_y_; }
    double omega_z() {return omega_z_; }
    double omega_z() const {return omega_z_; }
    const Eigen::MatrixXd& potential() const {return potential_;}

    void set_harmonic_potential(double omega_x, double omega_y);

    //operator() overloading
    double& operator()(int i, int j) { return potential_(i, j); }
    double operator()(int i, int j) const { return potential_(i, j); }

    void set_z_freq(double omega_z) {omega_z_ = omega_z;}

};

inline Grid<Dimension::Two>::Grid(int number_of_nodes_x, int number_of_nodes_y, double start_x, double start_y):
    size_x_(number_of_nodes_x), size_y_(number_of_nodes_y), start_x_(start_x), start_y_(start_y) {
        //Setting steps
        step_x_ = std::abs((2 * start_x_) / size_x_);
        step_y_ = std::abs((2 * start_y_) / size_y_);

        omega_x_ = 0;
        omega_y_ = 0;
        omega_z_ = 0;

        init_grid();
}

inline void Grid<Dimension::Two>::init_grid(){
    potential_ = Eigen::MatrixXd::Zero(size_x_,size_y_);
}

inline void Grid<Dimension::Two>::set_harmonic_potential(double omega_x, double omega_y){
    Eigen::MatrixXd V(size_x_, size_y_);
    V.setZero();

    omega_x_ = omega_x;
    omega_y_ = omega_y;

    for(int i = 0; i != size_x_; ++i){
        double x = start_x_ + i * step_x_;
        for(int j = 0; j != size_y_; ++j){
            double y = start_y_ + j * step_y_;
            V(i,j) = 0.5 * (std::pow(omega_x * x,2.) + std::pow(omega_y * y,2.));
        }
    }
    potential_ = V;
}



} // namespace GPES

#endif
