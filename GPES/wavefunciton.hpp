#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP

#include <Eigen/Dense>
#include <complex>
#include "definitions.hpp"
#include "grid.hpp"


namespace GPES
{

template<Dimension dim>
class WaveFunction;

template<>
class WaveFunction<Dimension::One>{
private:
    double _a_s, _a_dd, _omega;
    uint32_t _Num;
    Eigen::VectorXcd _Psi;

    double _size, _step, _start;


public:
    WaveFunction(Grid<Dimension::One>& grid, double a_s, double a_dd, uint32_t N): _a_s(a_s), _a_dd(a_dd), _Num(N) {
        _omega = grid.get_omega();
        _size = grid.get_size_of_grid();
        _step = grid.get_step_size();
        _start = grid.get_start_position();
        _Psi = Eigen::VectorXd::Zero(_size);
    }


    //Thomas-Fermi ansatz
    double thomas_fermi_state(double x);
    Eigen::VectorXcd TM_state();

    //Gauss wave function
    std::complex<double> gauss_wave_packet(double sigma_x, double x, double x_c, double p_x);
    double square_func(double x);

    //Initialization of starting 1D state
    void set_start_state_TM(double x_c);
    void set_start_state_Square(double x_c);
    void set_start_state_Gauss(double x_c, double sigma_x, double p_x);


    Eigen::VectorXcd get_vector() {return _Psi;}
    uint32_t get_Num() { return _Num; }

    // Acceptors function

    int size() { return _Psi.size(); }

    void setZero() { _Psi.setZero(); }

    //Overloading operators
    std::complex<double>& operator()(int i) { return _Psi(i); }

    std::complex<double> operator()(int i) const{ return _Psi(i); }


    void save_vector_to_csv(std::string filename, const Eigen::VectorXd& v);

    void normalize(Eigen::VectorXcd& vec);

    Eigen::VectorXd prob(Eigen::VectorXcd &vec);

};


double WaveFunction<Dimension::One>::square_func(double x){
    double R_tf = std::cbrt(1.5 * _Num);
    double high = this->_Num / R_tf;
    double potential = R_tf * R_tf * 0.5;
    double out = potential * (1 - std::pow(x/R_tf,2.));

    if(out > 0)
        return high;
    else 
        return 0;
    
}

std::complex<double> WaveFunction<Dimension::One>::gauss_wave_packet(double sigma_x,  double x,  double x_c, double p_x){
    std::complex<double> i(0, 1); // Define the imaginary unit
    double exponent = -(pow(x-x_c,2) / (2 * pow(sigma_x,2))); //this->m_size/2 //-x_c
    std::complex<double> phase = i * p_x * (x - x_c); 
    return std::exp(exponent); //+ phase
}

double WaveFunction<Dimension::One>::thomas_fermi_state(double x){
    double R_tf = std::cbrt(1.5 * _Num); // mult times _g_scattering
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

void WaveFunction<Dimension::One>::set_start_state_TM(double x_c){
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
    this->_Psi = U * std::abs(normalization_factor);
}

void WaveFunction<Dimension::One>::set_start_state_Square(double x_c){
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
    this->_Psi = U * std::abs(normalization_factor);
}

void WaveFunction<Dimension::One>::set_start_state_Gauss(double x_c, double sigma_x, double p_x){
    Eigen::VectorXcd U(_size);
    double psum = 0;
    double x = _start;
    for(int i = 0; i < _size; ++i){
        std::complex<double> c = gauss_wave_packet(sigma_x, x, x_c, p_x);
        U(i) = c;
        psum += std::norm(c);

        x += _step;
    }

    std::complex<double> normalization_factor = std::sqrt(_Num) / std::sqrt(psum * _step);
    this->_Psi = U * std::abs(normalization_factor);
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


template<>
class WaveFunction<Dimension::Two>{

};
    
}

#endif
