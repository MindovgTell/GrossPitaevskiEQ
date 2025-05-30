#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP

#include <fstream>
#include <Eigen/Dense>
#include <complex>
#include "definitions.hpp"
#include "grid.hpp"


namespace GPES
{

template<Dimension dim>
class WaveFunction;

//********************************/***********/********************************//
//                                                                             //
//************************/One dimensional wavefunction/***********************//
//                                                                             //
//********************************/***********/********************************//

template<>
class WaveFunction<Dimension::One> {
private:
    double _a_s, _a_dd, _omega, _g_scattering;
    unsigned int _Num, _size;
    Eigen::VectorXcd _Psi;

    double _step, _start;

public:
    WaveFunction(): _a_s(0), _a_dd(0), _Num(0), _size(0), _step(0), _start(0), _g_scattering(0), _Psi(Eigen::VectorXcd(0)) {}

    WaveFunction(Grid<Dimension::One>& grid, double a_s, double a_dd, int number_of_particles): _a_s(a_s), _a_dd(a_dd), _Num(number_of_particles) {
        _omega = grid.get_omega();
        _size = grid.get_size_of_grid();
        _step = grid.get_step_size();
        _start = grid.get_start_position();
        _g_scattering = -2. / _a_s;
        _Psi = Eigen::VectorXcd::Zero(_size);
    }

    WaveFunction(const Eigen::VectorXcd& vec): _Psi(vec){
        _size = vec.size();
        _step = 0;
        _start = 0;
    }

    // WaveFunction(const WaveFunction<Dimension::One>& vec): {
    //     _Psi = vec.get_wavefunction();
    //     _size = vec.get_size_of_grid();
    //     _step = 0;
    //     _start = 0;
    // }


    //Thomas-Fermi ansatz
    double thomas_fermi_state(double x);
    Eigen::VectorXcd TM_state();

    //Gauss wave function
    std::complex<double> gauss_wave_packet(double sigma_x, double x, double x_c);
    double square_func(double x);

    //Initialization of starting 1D state
    void set_state_TF(double x_c);
    void set_state_Square(double x_c);
    void set_state_Gauss(double x_c, double sigma_x);

    void set_vec(Eigen::VectorXcd& vec) {_Psi = vec;}
    void set_Num(unsigned int Num) {_Num = Num;}
    void set_size_of_grid(unsigned int size) {_size = size;}
    void set_start_position(double start) {_start = start;}
    void set_step_size(double step) {_step = step;}

    //Getters
    Eigen::VectorXcd& get_wavefunction() {return _Psi;}
    int get_Num() { return _Num; }
    unsigned int get_size_of_grid() {return _size; }
    double get_start_position() { return _start; }
    double get_step_size() { return _step; }
    double get_g_scat() { return _g_scattering; }
    double get_a_s() {return _a_s; }
    double get_a_dd() {return _a_dd; }
    // Acceptors function

    int size() { return _Psi.size(); }

    void setZero() { _Psi.setZero(); }

    //Overloading operators
    std::complex<double>& operator()(int i) { return _Psi(i); }
    std::complex<double> operator()(int i) const{ return _Psi(i); }

    //Functions for saving vectors into csv tables
    void save_vector_to_csv(std::string filename, const Eigen::VectorXd& v);
    void save_state_to_csv(std::string filename);

    void normalize(Eigen::VectorXcd& vec);
    Eigen::VectorXd prob(Eigen::VectorXcd &vec);
    Eigen::VectorXd prob();

};


double WaveFunction<Dimension::One>::square_func(double x){
    double R_tf = std::cbrt(1.5 * _Num * _g_scattering);
    double high = this->_Num / R_tf;
    double potential = R_tf * R_tf * 0.5;
    double out = potential * (1 - std::pow(x/R_tf,2.));

    if(out > 0)
        return high;
    else 
        return 0;
    
}

std::complex<double> WaveFunction<Dimension::One>::gauss_wave_packet(double sigma_x,  double x,  double x_c){
    std::complex<double> i(0, 1); // Define the imaginary unit
    double exponent = -(pow(x-x_c,2) / (2 * pow(sigma_x,2))); //this->m_size/2 //-x_c
    //std::complex<double> phase = i * p_x * (x - x_c); 
    return std::exp(exponent); //+ phase
}

double WaveFunction<Dimension::One>::thomas_fermi_state(double x){
    double R_tf = std::cbrt(1.5 * _Num * _g_scattering); // mult times _g_scattering
    double potential = R_tf * R_tf * 0.5;
    double out = potential * (1 - std::pow(x/R_tf,2.));
    if(out > 0)
        return std::sqrt(out);
    else 
        return 0;
}

Eigen::VectorXcd WaveFunction<Dimension::One>::TM_state(){
    Eigen::VectorXcd U(_size);

    double x = this->_start;
    std::complex<double> psum = 0;

    for(int i = 0; i < _size; ++i){
        x += _step;
        std::complex<double> c = thomas_fermi_state(x); // Add parameters for Thomas-Fermi function
        U(i) = c;
        psum += std::norm(c);
    }
    std::complex<double> normalization_factor = std::sqrt(_Num) / std::sqrt(psum * _step);

    U *= std::abs(normalization_factor);
    return U;
}

void WaveFunction<Dimension::One>::set_state_TF(double x_c){
    Eigen::VectorXcd U(_size);
    double psum = 0;
    double x = _start;
    for(int i = 0; i < _size; ++i){
        std::complex<double> c = thomas_fermi_state(x - x_c); // Add parameters for Thomas-Fermi function
        U(i) = c;
        psum += std::norm(c);

        x += _step;
    }

    std::complex<double> normalization_factor = std::sqrt(_Num) / std::sqrt(psum * _step);
    _Psi = U * std::abs(normalization_factor);
}

void WaveFunction<Dimension::One>::set_state_Square(double x_c){
    Eigen::VectorXcd U(_size);
    double psum = 0;
    double x = _start;
    for(int i = 0; i < _size; ++i){
        std::complex<double> c = square_func(x);
        U(i) = c;
        psum += std::norm(c);

        x += _step;
    }

    std::complex<double> normalization_factor = std::sqrt(_Num) / std::sqrt(psum * _step);
    _Psi = U * std::abs(normalization_factor);
}

void WaveFunction<Dimension::One>::set_state_Gauss(double x_c, double sigma_x){
    Eigen::VectorXcd U(_size);
    double psum = 0;
    double x = _start;
    for(int i = 0; i < _size; ++i){
        std::complex<double> c = gauss_wave_packet(sigma_x, x, x_c);
        U(i) = c;
        psum += std::norm(c);
        x += _step;
    }

    std::complex<double> normalization_factor = std::sqrt(_Num) / std::sqrt(psum * _step);
    _Psi = U * std::abs(normalization_factor);
} 

// Returning the probability density vector
Eigen::VectorXd WaveFunction<Dimension::One>::prob(Eigen::VectorXcd &vec)
{
    int size = vec.size();
    Eigen::VectorXd pr(size);
    for(int i = 0; i < size; ++i){
        pr(i) = std::norm(vec(i));
    }
    return pr;
}

Eigen::VectorXd WaveFunction<Dimension::One>::prob()
{
    Eigen::VectorXd pr(_size);
    for(int i = 0; i < _size; ++i){
        pr(i) = std::norm(_Psi(i));
    }
    return pr;
}

void WaveFunction<Dimension::One>::normalize(Eigen::VectorXcd &vec){
        int size = vec.size();
        std::complex<double> psum = 0;
        for(int i = 0; i < size; ++i){
            psum += std::norm(vec(i));
        }
        std::complex<double> normalization_factor = std::sqrt(_Num) / std::sqrt(psum * _step); // 
        vec *= std::abs(normalization_factor);
    }

void WaveFunction<Dimension::One>::save_vector_to_csv(std::string filename, const Eigen::VectorXd& v){
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ",","\n");
    int size = v.size();
    Eigen::VectorXd vector(size);

    std::ofstream file(filename);
    if(file.is_open()){
        //file << vec.format(CSVFormat) << '\n';
        file << vector.format(CSVFormat);
        file.close();
    }
}

void WaveFunction<Dimension::One>::save_state_to_csv(std::string filename){
    Eigen::VectorXd pr = prob();
    save_vector_to_csv(filename, pr);
}




//********************************/***********/********************************//
//                                                                             //
//************************/Two dimensional wavefunction/***********************//
//                                                                             //
//********************************/***********/********************************//

template<>
class WaveFunction<Dimension::Two> {
private:
    double _a_s, _a_dd, _omega_x, _omega_y, _g_scattering;
    unsigned int _Num, _size_x, _size_y;
    double _step_x, _step_y, _start_x, _start_y;

    Eigen::VectorXcd _Psi;

public:
    WaveFunction(): _a_s(0), _g_scattering(0), _a_dd(0), _Num(0), _size_x(0), _step_x(0), _start_x(0), _size_y(0), _step_y(0), _start_y(0), _Psi(Eigen::VectorXcd(0)) {}

    WaveFunction(Grid<Dimension::Two>& grid, double a_s, double a_dd, int number_of_particles): _a_s(a_s), _a_dd(a_dd), _Num(number_of_particles) {
        _omega_x = grid.get_omega_x();
        _omega_y = grid.get_omega_y();
        _size_x = grid.get_size_of_grid_x();
        _size_y = grid.get_size_of_grid_y();
        _step_x = grid.get_step_size_x();
        _step_y = grid.get_step_size_y();
        _start_x = grid.get_start_position_x();
        _start_y = grid.get_start_position_y();
        _g_scattering = std::sqrt(8. * M_PI) * a_s;
        _Psi = Eigen::VectorXcd::Zero(_size_x * _size_y);
    }


    int get_index(int i, int j){ return i*_size_x + j; } // M is 
    
    //Thomas-Fermi anzats
    double thomas_fermi_radius_x();
    double thomas_fermi_radius_y();
    double chem_potential();
    double thomas_fermi_state(double x, double y);
    Eigen::VectorXcd TF_state();

    std::complex<double> gauss_wave_packet(double x, double y, double x_c, double y_c, double sigma_x, double sigma_y);

    //Initialization of starting 1D state
    void set_state_TF(double x_c, double y_c);
    void set_state_Gauss(double x_c, double y_c, double sigma_x, double sigma_y);

    void set_vec(Eigen::VectorXcd& vec) {_Psi = vec;}
    void set_Num(unsigned int Num) {_Num = Num;}
    void set_size_of_grid_x(unsigned int size) {_size_x = size;}
    void set_size_of_grid_y(unsigned int size) {_size_y = size;}
    void set_start_position_x(double start) {_start_x = start;}
    void set_start_position_y(double start) {_start_y = start;}
    void set_step_size_x(double step) {_step_x = step;}
    void set_step_size_y(double step) {_step_y = step;}

    //Getters
    Eigen::VectorXcd& get_wavefunction() {return _Psi;}
    unsigned int get_Num() { return _Num; }
    unsigned int get_size_of_grid_x() {return _size_x; }
    unsigned int get_size_of_grid_y() {return _size_y; }
    double get_start_position_x() { return _start_x; }
    double get_start_position_y() { return _start_y; }
    double get_step_size_x() { return _step_x; }
    double get_step_size_y() { return _step_y; }

    double get_g_scattering() { return _g_scattering;}
    double get_a_s() { return _a_s;}
    double get_a_dd() { return _a_dd;}
    // Acceptors function

    unsigned int size() { return _Psi.size(); }

    void setZero() { _Psi.setZero(); }

    //Overloading operators
    std::complex<double>& operator()(int i) { return _Psi(i); }
    std::complex<double> operator()(int i) const{ return _Psi(i); }


    //Probability density funcitons
    Eigen::VectorXd prob(Eigen::VectorXcd &vec);
    Eigen::VectorXd prob();
    double prob(int index);
};

double WaveFunction<Dimension::Two>::thomas_fermi_radius_x(){
    double numerator = 4.0 * _omega_y * _Num;
    double denumerator = std::pow(_omega_x, 3.) * M_PI;
    return std::pow(( numerator / denumerator), 0.25);
}

double WaveFunction<Dimension::Two>::thomas_fermi_radius_y(){
    double numerator = 4.0 * _omega_x * _Num;
    double denumerator = std::pow(_omega_y, 3.) * M_PI;
    return std::pow(( numerator / denumerator), 0.25);
}

double WaveFunction<Dimension::Two>::chem_potential(){
    double numerator = 1 * _omega_x * _omega_y * _Num; //add g_scattering
    return std::sqrt(numerator / M_PI);
}

double WaveFunction<Dimension::Two>::thomas_fermi_state(double x, double y){
    double out, R_x, R_y;
    R_x = thomas_fermi_radius_x();
    R_y = thomas_fermi_radius_y();

    out = chem_potential() * (1 - std::pow(x/R_x, 2.) - std::pow(y/R_y, 2.)); //  /,  we need to add definition for 

    if (out > 0)
        return std::sqrt(out);
    else 
        return 0;
}

std::complex<double> WaveFunction<Dimension::Two>::gauss_wave_packet(double x, double y, double x_c, double y_c, double sigma_x, double sigma_y){ //, double p_x, double p_y
    //std::complex<double> i(0, 1); // Define the imaginary unit
    double exponent = -1.0 * (std::pow(x - x_c ,2.) / (2 * std::pow(sigma_x,2.))) - (std::pow(y - y_c, 2.) / (2*std::pow(sigma_y, 2.)));

    //std::complex<double> phase = i * (p_x * (x - x_c) + p_y * (y - y_c));
    return std::exp(exponent); // + phase
}

Eigen::VectorXcd WaveFunction<Dimension::Two>::TF_state(){
    int size = _size_x * _size_y;
    Eigen::VectorXcd U(size);
    std::complex<double> psum = 0;
    for(int i = 0; i < _size_x; ++i){
        double x = _start_x + i * _step_x;
        for(int j = 0; j < _size_y; ++j){
            double y = _start_y + j * _step_y;
            //Initial state function
            std::complex<double> c = thomas_fermi_state(x,y);

            int index = get_index(i,j);
            U(index) = c;
            psum += std::norm(c);
        }
    }
    std::complex<double> normalization_factor = std::sqrt(_Num) / std::sqrt(psum * _step_x * _step_y);
    U *= std::abs(normalization_factor);
    return U;
}

// Function for initializing wave function in Thomas-Fermi limit without ddi 
void WaveFunction<Dimension::Two>::set_state_TF(double x_c, double y_c){
    int size = _size_x * _size_y;
    Eigen::VectorXcd U(size);
    std::complex<double> psum = 0;
    
    for(int i = 0; i < _size_x; ++i){
        double x = _start_x + i * _step_x;
        
        for(int j = 0; j < _size_y; ++j){
            double y = _start_y + j * _step_y;
            //Initial state function
            std::complex<double> c = thomas_fermi_state(x,y);
            int index = get_index(i,j);
            U(index) = c;
            psum += std::norm(c);
        }
    }

    std::complex<double> normalization_factor = std::sqrt(_Num) / std::sqrt(psum * _step_x * _step_y);
    U *= std::abs(normalization_factor);
    _Psi = U;
}

// Function for initializing wave function in Thomas-Fermi limit without ddi 
void WaveFunction<Dimension::Two>::set_state_Gauss(double x_c, double y_c, double sigma_x, double sigma_y){
    int size = _size_x * _size_y;
    Eigen::VectorXcd U(size);
    std::complex<double> psum = 0;
    for(int i = 0; i < _size_x; ++i){
        double x = _start_x + i * _step_x;
        for(int j = 0; j < _size_y; ++j){
            double y = _start_y + j * _step_y;
            //Initial state function
            std::complex<double> c = gauss_wave_packet(x, y, x_c, y_c, sigma_x, sigma_y); // ,p_x, p_y
            int index = get_index(i,j);
            U(index) = c;
            psum += std::norm(c);
        }
    }

    std::complex<double> normalization_factor = std::sqrt(_Num) / std::sqrt(psum * _step_x * _step_y);
    U *= std::abs(normalization_factor);
    _Psi = U;
}

            

// Returning the probability density vector
Eigen::VectorXd WaveFunction<Dimension::Two>::prob(Eigen::VectorXcd &vec){
    int size = vec.size();
    Eigen::VectorXd pr(size);
    for(int i = 0; i < size; ++i){
        pr(i) = std::norm(vec(i));
    }
    return pr;
}

Eigen::VectorXd WaveFunction<Dimension::Two>::prob(){
    int size = _size_x * _size_y;
    Eigen::VectorXd pr(size);
    for(int i = 0; i < size; ++i){
        pr(i) = std::norm(_Psi(i));
    }
    return pr;
}

double WaveFunction<Dimension::Two>::prob(int index){
    int size = _size_x * _size_y;
    if(index > 0 && index < size)
        return std::norm(_Psi(index));
    else 
        return 0;
}



    
}

#endif
