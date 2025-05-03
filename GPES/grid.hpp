#include "definitions.hpp"
// #include <Eigen/Dense>
// #include <Eigen/Sparse>

namespace GPES {
    
    template <Dimension dim>
    class Grid;


    template<>
    class Grid<Dimension::One> {

    private: 
        unsigned int _N;
        double _step, _start;

    public:
        Grid(unsigned int number_of_nodes, double step);
        unsigned int size_of_grid() {return this->_N;}
        double start_location() { return this->_start;}
        double step_size() {return this->_step;}
    };

    inline Grid<Dimension::One>::Grid(unsigned int number_of_nodes, double start): _N(number_of_nodes), _start(start) {
        this->_step = (-2 * _start) / this->_N;
    }


    template<>
    class Grid<Dimension::Two> {

    private: 
        int _N_x, _N_y;
        double _step_x, _step_y, _start_x, _start_y;

    public:
        Grid();
        unsigned int size();
        double start_location();
        double step_size();
    };

    inline Grid<Dimension::Two>::Grid()  {}


} // namespace GPES