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
        unsigned int _N;
        double _step, _start;
        double _omega;

        Eigen::VectorXd _potential;

    public:
        //Construction
        Grid(unsigned int number_of_nodes, double step);
        // Delete copy and move constructors
        Grid(const GPES::Grid<Dimension::One>&) = delete;
        Grid(GPES::Grid<Dimension::One> &&)     = delete;


        //Getters
        unsigned int size_of_grid() {return this->_N;}
        double start_position() { return this->_start;}
        double step_size() {return this->_step;}

        Eigen::VectorXd get_potential(){ return this->_potential;}

        //initialize grid size
        void init_grid();

        void set_harmonic_potential(double omega);
    };

    inline Grid<Dimension::One>::Grid(unsigned int number_of_nodes, double start): _N(number_of_nodes), _start(start) {
        //Setting step
        this->_step = (-2 * _start) / this->_N;
        this->_omega = 0;
        this->init_grid();
    }

    void inline Grid<Dimension::One>::init_grid(){
        Eigen::VectorXd U(_N);
        U.setZero();
        this->_potential = U;
    }

    void inline Grid<Dimension::One>::set_harmonic_potential(double omega){
        Eigen::VectorXd V(_N);
        V.setZero();
        for(int i = 1; i != _N - 1; ++i){
            double x = this->_start + i * this->_step;
            V(i) = 0.5 * (std::pow(omega * x,2.));
        }

        this->_omega = omega;

        this->_potential = V;
    }



    //********************************/***********/********************************//
    //                                                                             //
    //********************************/2DFunctions/********************************//
    //                                                                             //
    //********************************/***********/********************************//

    template<>
    class Grid<Dimension::Two> {

    private: 
        int _N_x, _N_y;
        double _step_x, _step_y, _start_x, _start_y;
        
        //harmonic potential frequencies
        double _omega_x, _omega_y, _lambda;

        Eigen::MatrixXd _potential;

    public:
        //Constructor 
        Grid(int N_x, int N_y, double start_x, double start_y);
        // Delete copy and move constructors
        Grid(const GPES::Grid<Dimension::Two>&) = delete;
        Grid(GPES::Grid<Dimension::Two> &&)     = delete;

        //Getters
        unsigned int size_of_grid_x(){ return this->_N_x; }
        unsigned int size_of_grid_y(){ return this->_N_y; }

        double start_position_x(){ return this->_start_x; }
        double start_position_x(){ return this->_start_y; }

        double step_size_x(){ return this->_step_x; }
        double step_size_y(){ return this->_step_y; }

        Eigen::MatrixXd get_potential() {return this->_potential; }


        //initialize grid
        void init_grid();

        void set_harmonic_potential(double omega_x, double omega_y);

    };

    inline Grid<Dimension::Two>::Grid(int N_x, int N_y, double start_x, double start_y):
        _N_x(N_x), _N_y(N_y), _start_x(start_x), _start_y(start_y) {
            //Setting steps
            this->_step_x = (-2 * _start_x) / this->_N_x;
            this->_step_y = (-2 * _start_y) / this->_N_y;

            _omega_x = 0;
            _omega_y = 0;
            _lambda = 0;

            this->init_grid();
    }

    void inline Grid<Dimension::Two>::init_grid(){
        Eigen::MatrixXd U(_N_x, _N_y);
        U.setZero();
        this->_potential = U;
    }

    void inline Grid<Dimension::Two>::set_harmonic_potential(double omega_x, double omega_y){
        Eigen::MatrixXd V(_N_x, _N_y);
        V.setZero();

        this->_omega_x = omega_x;
        this->_omega_y = omega_y;

        for(int i = 1; i != _N_x - 1; ++i){
            double x = this->_start_x + i * this->_step_x;
            for(int j = 1; j != _N_y-1; ++j){
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