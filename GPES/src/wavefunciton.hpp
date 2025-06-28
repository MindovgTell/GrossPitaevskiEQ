#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <iomanip>
#include <map>
#include <fftw3.h>
#include "definitions.hpp"
#include "grid.hpp"


namespace GPES
{

using ParamList = std::map<std::string,double>;
using MatR = Eigen::Matrix<std::complex<double>,
                           Eigen::Dynamic,
                           Eigen::Dynamic,
                           Eigen::RowMajor>;

template<Dimension dim>
class WaveFunction;


static std::vector<std::string> split(const std::string& str, char delim = ',')
{
    std::vector<std::string> tokens;
    std::stringstream        ss(str);
    std::string              item;

    while (std::getline(ss, item, delim))
        tokens.emplace_back(item);

    return tokens;
}


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
        if (_a_s != 0)
            _g_scattering =  _a_s;
        else
            _g_scattering = 0;
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

    void set_vec(Eigen::VectorXcd& vec)         { _Psi = vec; }
    void set_Num(unsigned int Num)              { _Num = Num; }
    void set_size_of_grid(unsigned int size)    { _size = size; }
    void set_start_position(double start)       { _start = start; }
    void set_step_size(double step)             { _step = step; }
    void set_a_s(double a_s)                    { _a_s = a_s; }
    void set_a_dd(double a_dd)                  { _a_dd = a_dd; }
    void set_omega(double omega)                { _omega = omega; }

    //Getters
    Eigen::VectorXcd get_wavefunction()         { return _Psi; }
    int get_Num()                               { return _Num; }
    unsigned int get_size_of_grid()             { return _size; }
    double get_start_position()                 { return _start; }
    double get_step_size()                      { return _step; }
    double get_g_scat()                         { return _g_scattering; }
    double get_a_s()                            { return _a_s; }
    double get_a_dd()                           { return _a_dd; }
    // Acceptors function

    int size() { return _Psi.size(); }

    void setZero() { _Psi.setZero(); }

    //Overloading operators
    std::complex<double>& operator()(int i) { return _Psi(i); }
    std::complex<double> operator()(int i) const{ return _Psi(i); }

    //Functions for saving vectors into csv tables
    void savecsv(const std::string filename, const Eigen::VectorXd& v);
    void savecsv(const std::string filename, const Eigen::VectorXcd& v);
    void savecsv_prob(const std::string filename);
    void savecsv_state(const std::string filename);

    void readcsv(const std::string filename);
    // WaveFunction<Dimension::One> read_from_csv(
    //                                            const std::string&               filename,
    //                                            Eigen::VectorXcd&                v
    //                                           );

    void print_params();


    void normalize(Eigen::VectorXcd& vec);
    Eigen::VectorXd prob(Eigen::VectorXcd &vec);
    Eigen::VectorXd prob();


    Eigen::VectorXcd momentum_space_transform() const;
};


double WaveFunction<Dimension::One>::square_func(double x){
    double R_tf = std::cbrt(1.5 * _Num * _g_scattering);
    double high = _Num / R_tf;
    double potential = R_tf * R_tf * 0.5;
    double out = potential * (1 - std::pow(x/R_tf,2.));

    if(out > 0)
        return high;
    else 
        return 0;
    
}

std::complex<double> WaveFunction<Dimension::One>::gauss_wave_packet(double sigma_x,  double x,  double x_c){
    // std::complex<double> i(0, 1); // Define the imaginary unit
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

void WaveFunction<Dimension::One>::savecsv(std::string filename, const Eigen::VectorXd& v){
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ",","\n");
    std::ofstream file(filename);
    if(file.is_open()){    
        //file << vec.format(CSVFormat) << '\n';
        file << v.format(CSVFormat);
        file.close();
    }
}

// void WaveFunction<Dimension::One>::savecsv(std::string filename, const Eigen::VectorXcd& v){
//     const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ",","\n");
//     std::ofstream file(filename);
//     if(file.is_open()){    
//         //file << vec.format(CSVFormat) << '\n';
//         file << v.format(CSVFormat);
//         file.close();
//     }
// }

//test
void WaveFunction<Dimension::One>::savecsv(const std::string filename, const Eigen::VectorXcd& v){


    std::ofstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Could not open file: " + filename);


    ParamList params{
        {"a_s",              _a_s},
        {"a_dd",             _a_dd},
        {"omega",            _omega},
        {"Num_of_particle",  static_cast<double>(_Num)},
        {"size_of_grid",     static_cast<double>(_size)},
        {"step_size",        _step},
        {"grid_start_point", _start}
    };

    /* ---------- 1) metadata names ---------- */
    bool first = true;
    for (const auto& kv : params) {
        if (!first) file << ',';
        file << kv.first;
        first = false;
    }
    file << '\n';

    /* ---------- 2) metadata values ---------- */
    first = true;
    file << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const auto& kv : params) {
        if (!first) file << ',';
        file << kv.second;
        first = false;
    }
    file << "\n\n";                           // blank line before the data block

    /* ---------- 3) the actual wave-function ---------- */
    file << "index,real,imag\n";
    for (Eigen::Index i = 0; i < v.size(); ++i) {
        file << i << ',' << v[i].real() << ',' << v[i].imag() << '\n';
    }
}

void WaveFunction<Dimension::One>::savecsv_prob(const std::string filename){
    Eigen::VectorXd pr = prob();
    savecsv(filename, pr);
}


void WaveFunction<Dimension::One>::savecsv_state(const std::string filename){
    savecsv(filename, _Psi);
}



void WaveFunction<Dimension::One>::readcsv(const std::string filename){
    // Open and check file competability
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Cannot open file: " + filename);

    /* Physical parameters keys */
    std::string line;
    if (!std::getline(file, line))
        throw std::runtime_error("File is empty: " + filename);
    std::vector<std::string> keys = split(line);

    /* Physical parameters values */
    if (!std::getline(file, line))
        throw std::runtime_error("Missing value line in " + filename);
    std::vector<std::string> valStr = split(line);

    if (keys.size() != valStr.size())
        throw std::runtime_error("Key/value count mismatch in " + filename);

    ParamList params;
    for (std::size_t i = 0; i < keys.size(); ++i)
        params[keys[i]] = std::stod(valStr[i]);

    /* ---------- 3) skip blank lines + verify data header ------------------ */
    while (std::getline(file, line) && line.empty()) /*skip*/ ;
    if (line != "index,real,imag")
        throw std::runtime_error("Unexpected data header: " + line);

    /* ---------- 4) read the complex vector -------------------------------- */
    std::vector<std::complex<double>> data;
    while (std::getline(file, line))
    {
        if (line.empty()) continue;
        auto cols = split(line);
        if (cols.size() != 3)
            throw std::runtime_error("Bad data row: " + line);
        data.emplace_back(std::stod(cols[1]), std::stod(cols[2]));
    }

    Eigen::VectorXcd vec(data.size());
    for (Eigen::Index i = 0; i < vec.size(); ++i) vec[i] = data[i];

    /*Initialization of wavefunciton parameters*/
    _Psi    =   vec;
    _a_dd   =   params.at("a_dd");
    _a_s    =   params.at("a_s");
    _omega  =   params.at("omega");
    _Num    =   static_cast<std::size_t>(params.at("Num_of_particle"));
    _size   =   static_cast<std::size_t>(params.at("size_of_grid"));
    _step   =   params.at("step_size");
    _start  =   params.at("grid_start_point");
}


// WaveFunction<Dimension::One> WaveFunction<Dimension::One>::read_from_csv(
//                                                                          const std::string&               filename,
//                                                                          Eigen::VectorXcd&                v
//                                                                         ) {
//     std::ifstream file(filename);
//     if (!file.is_open())
//         throw std::runtime_error("Cannot open file: " + filename);
//     /* ---------- 1) metadata keys ---------- */
//     std::string line;
//     if (!std::getline(file, line))
//         throw std::runtime_error("File is empty: " + filename);
//     std::vector<std::string> keys = split(line);
//     /* ---------- 2) metadata values -------- */
//     if (!std::getline(file, line))
//         throw std::runtime_error("Missing value line in " + filename);
//     std::vector<std::string> valStr = split(line);
//     if (keys.size() != valStr.size())
//         throw std::runtime_error("Key/value count mismatch in " + filename);
//     ParamList params;
//     params.reserve(keys.size());
//     for (std::size_t i = 0; i < keys.size(); ++i)
//         params.emplace_back(keys[i], std::stod(valStr[i]));
//     /* ---------- 3) skip blank lines + verify data header ------------------ */
//     while (std::getline(file, line) && line.empty()) /*skip*/ ;
//     if (line != "index,real,imag")
//         throw std::runtime_error("Unexpected data header: " + line);
//     /* ---------- 4) read the complex vector -------------------------------- */
//     std::vector<std::complex<double>> data;
//     while (std::getline(file, line))
//     {
//         if (line.empty()) continue;
//         auto cols = split(line);
//         if (cols.size() != 3)
//             throw std::runtime_error("Bad data row: " + line);
//         data.emplace_back(std::stod(cols[1]), std::stod(cols[2]));
//     }
//     Eigen::VectorXcd vec(data.size());
//     for (Eigen::Index i = 0; i < vec.size(); ++i) vec[i] = data[i];
//     /* ---------- 5) build and return the object ---------------------------- */
//     return WaveFunction<Dimension::One>(std::move(vec), std::move(params));
// }


void WaveFunction<Dimension::One>::print_params(){
    int width = 10;
    std::cout << std::setw(width) << "grid_size" << std::setw(width) << "start" << std::setw(width) << "N_of_mol"  << std::setw(width) << "a_s" 
    << std::setw(width) << "a_dd" << std::setw(width) << std::endl;

    std::cout << std::setw(width) << _size << std::setw(width) << _start << std::setw(width)
    << _Num << std::setw(width) << _a_s << std::setw(width) << _a_dd << std::setw(width) << std::endl;
}



Eigen::VectorXcd WaveFunction<Dimension::One>::momentum_space_transform() const {
    const int N = _size;
    Eigen::VectorXcd result(N);

    // Allocate aligned memory
    fftw_complex* in = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * N));
    fftw_complex* out = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * N));

    // Fill input
    for (int i = 0; i < N; ++i) {
        in[i][0] = _Psi(i).real();
        in[i][1] = _Psi(i).imag();
    }

    // Create and execute plan
    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_MEASURE);
    fftw_execute(plan);

    // Normalize and shift in one pass
    const int halfN = N / 2;
    for (int i = 0; i < N; ++i) {
        result(i) = std::complex<double>(out[i][0], out[i][1]) /  static_cast<double>(N);
    }

    // Cleanup
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return result;
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
    WaveFunction(): _a_s(0), _g_scattering(0), _a_dd(0), _Num(0), _size_x(0), _step_x(0), _start_x(0), _size_y(0), _step_y(0), _start_y(0), _omega_x(0), _omega_y(0),_Psi(Eigen::VectorXcd(0)) {}

    WaveFunction(Grid<Dimension::Two>& grid, double a_s, double a_dd, int number_of_particles): _a_s(a_s), _a_dd(a_dd), _Num(number_of_particles) {
        _omega_x        =   grid.get_omega_x();
        _omega_y        =   grid.get_omega_y();
        _size_x         =   grid.get_size_of_grid_x();
        _size_y         =   grid.get_size_of_grid_y();
        _step_x         =   grid.get_step_size_x();
        _step_y         =   grid.get_step_size_y();
        _start_x        =   grid.get_start_position_x();
        _start_y        =   grid.get_start_position_y();
        _g_scattering   =   std::sqrt(8. * M_PI) * a_s;
        _Psi            =   Eigen::VectorXcd::Zero(_size_x * _size_y);
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

    void set_vec(Eigen::VectorXcd& vec) {if(vec.size() == _size_x*_size_y) _Psi = vec;}
    void set_Num(unsigned int Num) {_Num = Num;}
    void set_size_of_grid_x(unsigned int size) {_size_x = size;}
    void set_size_of_grid_y(unsigned int size) {_size_y = size;}
    void set_start_position_x(double start) {_start_x = start;}
    void set_start_position_y(double start) {_start_y = start;}
    void set_step_size_x(double step) {_step_x = step;}
    void set_step_size_y(double step) {_step_y = step;}
    void set_a_s(double a_s) {_a_s = a_s;}
    void set_a_dd(double a_dd) {_a_dd = a_dd;}
    void set_omega_x(double omega) {_omega_x = omega;}
    void set_omega_y(double omega) {_omega_y = omega;}

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

    Eigen::VectorXd get_x_slice();
    Eigen::VectorXd get_y_slice();
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

    void savecsv(const std::string filename, const Eigen::VectorXcd& v);
    void savecsv_state(const std::string filename);
    void savecsv_prob(const std::string filename);

    void readcsv(const std::string filename);

    double calc_state_energy();

    void print_params();
    Eigen::VectorXcd momentum_space_transform() const;
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

Eigen::VectorXd WaveFunction<Dimension::Two>::get_x_slice(){
    int y_center = _size_y / 2;
    Eigen::VectorXd slice(_size_x);
    for(int i = 0; i < _size_x; ++i){
        slice(i) = std::norm(_Psi(y_center * _size_x + i));
    }
    return slice;
}


Eigen::VectorXd WaveFunction<Dimension::Two>::get_y_slice(){
    int x_center = _size_x / 2;
    Eigen::VectorXd slice(_size_y);
    for(int j = 0; j < _size_y; ++j){
        slice(j) = std::norm(_Psi(j * _size_x + x_center));
    }
    return slice;
}
    

static std::string parent_dir(const std::string& path) {
    auto pos = path.find_last_of("/\\");
    return (pos == std::string::npos ? "" : path.substr(0, pos));
}

void WaveFunction<Dimension::Two>::savecsv(const std::string file_path, const Eigen::VectorXcd& v){
    // 1) Ensure parent directory exists
    std::string dir = parent_dir(file_path);
    if (!dir.empty()) {
        // mkdir -p dir
        std::string cmd = "mkdir -p '" + dir + "'";
        if (std::system(cmd.c_str()) != 0) {
            throw std::runtime_error("Failed to create directory: " + dir);
        }
    }

    std::ofstream file(file_path, std::ios::out /*| std::ios::trunc is implicit*/);
    if (!file) {
        std::cerr << "Error: cannot open " << file_path << " for writing\n";
        return;
    }

    ParamList params{
        {"a_s",                _a_s                        },
        {"a_dd",               _a_dd                       },
        {"omega_x",            _omega_x                    },
        {"omega_y",            _omega_y                    },
        {"Num_of_particle",    static_cast<double>(_Num)   },
        {"size_of_grid_x",     static_cast<double>(_size_x)},
        {"size_of_grid_y",     static_cast<double>(_size_y)},
        {"step_size_x",        _step_x                     },
        {"step_size_y",        _step_y                     },
        {"grid_start_point_x", _start_x                    },
        {"grid_start_point_y", _start_y                    }
    };

    /* ---------- 1) metadata names ---------- */
    bool first = true;
    for (const auto& kv : params) {
        if (!first) file << ',';
        file << kv.first;
        first = false;
    }
    file << '\n';

    /* ---------- 2) metadata values ---------- */
    first = true;
    file << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const auto& kv : params) {
        if (!first) file << ',';
        file << kv.second;
        first = false;
    }
    file << "\n\n";                           // blank line before the data block

    /* ---------- 3) the actual wave-function ---------- */
    file << "index,real,imag\n";
    for (Eigen::Index i = 0; i < v.size(); ++i) {
        file << i << ',' << v[i].real() << ',' << v[i].imag() << '\n';
    }

    std::cout << "State have been saved" << std::endl;
}

void WaveFunction<Dimension::Two>::savecsv_prob(const std::string filename){
    Eigen::VectorXd pr = prob();
    savecsv(filename, pr);
}

void WaveFunction<Dimension::Two>::savecsv_state(const std::string filename){
    savecsv(filename, _Psi);
}


static double safe_parse_double(const std::string& str) {
    try {
        return std::stod(str);
    }
    catch (const std::invalid_argument&) {
        return 0.0;
    }
    catch (const std::out_of_range&) {
        return 0.0;
    }
}


void WaveFunction<Dimension::Two>::readcsv(const std::string filename){
    // Open and check file competability
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Cannot open file: " + filename);

    /* Physical parameters keys */
    std::string line;
    if (!std::getline(file, line))
        throw std::runtime_error("File is empty: " + filename);
    std::vector<std::string> keys = split(line);

    /* Physical parameters values */
    if (!std::getline(file, line))
        throw std::runtime_error("Missing value line in " + filename);
    std::vector<std::string> valStr = split(line);

    if (keys.size() != valStr.size())
        throw std::runtime_error("Key/value count mismatch in " + filename);

    ParamList params;
    for (std::size_t i = 0; i < keys.size(); ++i)
        params[keys[i]] = std::stod(valStr[i]);

    /* ---------- 3) skip blank lines + verify data header ------------------ */
    while (std::getline(file, line) && line.empty()) /*skip*/ ;
    if (line != "index,real,imag")
        throw std::runtime_error("Unexpected data header: " + line);

    /* ---------- 4) read the complex vector -------------------------------- */
    std::vector<std::complex<double>> data;
    data.reserve(10000);
    while (std::getline(file, line))
    {
        if (line.empty()) continue;
        auto cols = split(line);
        if (cols.size() != 3)
            throw std::runtime_error("Bad data row: " + line); 
        data.emplace_back(safe_parse_double(cols[1]), safe_parse_double(cols[2]));
    }

    Eigen::VectorXcd vec(data.size());
    for (Eigen::Index i = 0; i < vec.size(); ++i) vec[i] = data[i];

    /*Initialization of wavefunciton parameters*/
    _Psi        =   vec;
    _a_dd       =   params.at("a_dd");
    _a_s        =   params.at("a_s");
    _omega_x    =   params.at("omega_x");
    _omega_y    =   params.at("omega_y");
    _Num        =   static_cast<std::size_t>(params.at("Num_of_particle"));
    _size_x     =   static_cast<std::size_t>(params.at("size_of_grid_x"));
    _size_y     =   static_cast<std::size_t>(params.at("size_of_grid_y"));
    _step_x     =   params.at("step_size_x");
    _step_y     =   params.at("step_size_y");
    _start_x    =   params.at("grid_start_point_x");
    _start_y    =   params.at("grid_start_point_y");
}


void WaveFunction<Dimension::Two>::print_params(){
    int width = 10;
    std::cout << std::setw(width) << "grid_size_x"<< std::setw(width) << "grid_size_y" << std::setw(width) << "start_x" << std::setw(width) << "start_y"
    << std::setw(width) << "N_of_mol"  << std::setw(width) << "a_s"  << std::setw(width) << "a_dd"  << std::setw(width) << "omega_x" << std::setw(width) << "omega_y" << std::setw(width) << std::endl;

    std::cout << std::setw(width) << _size_x << std::setw(width) << _size_y  << std::setw(width) << _start_x << std::setw(width) << _start_y << std::setw(width)
    << _Num << std::setw(width) << _a_s << std::setw(width) << _a_dd << std::setw(width) << _omega_x << std::setw(width) << _omega_y << std::setw(width) << std::endl;
}


Eigen::VectorXcd WaveFunction<Dimension::Two>::momentum_space_transform() const {
    int Nx = _size_x;
    int Ny = _size_y;
    int N = Nx * Ny;

    Eigen::VectorXcd result(N);

    // Allocate FFTW input/output
    fftw_complex* in = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * N));
    fftw_complex* out = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * N));

    // Fill input array from _Psi (row-major assumed: y rows, x columns)
    for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
            int i = y * Nx + x; // make sure get_index is get_index(row, col)
            in[i][0] = _Psi(i).real();
            in[i][1] = _Psi(i).imag();

            // Optional: Apply a window function like Hann or Gaussian here
            // double wx = 0.5 * (1 - std::cos(2 * M_PI * x / (Nx - 1)));
            // double wy = 0.5 * (1 - std::cos(2 * M_PI * y / (Ny - 1)));
            // double window = wx * wy;
            // in[i][0] *= window;
            // in[i][1] *= window;
        }
    }

    // Plan and execute FFT
    fftw_plan plan = fftw_plan_dft_2d(Ny, Nx, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    // Normalize result by total number of points (preserve amplitude scale)
    for (int i = 0; i < N; ++i) {
        result(i) = std::complex<double>(out[i][0], out[i][1]) / static_cast<double>(N);
    }

    // Apply fftshift to center zero momentum
    Eigen::VectorXcd shifted_result(N);
    for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
            int orig_idx = y * Nx + x;
            int sx = (x + Nx / 2) % Nx;
            int sy = (y + Ny / 2) % Ny;
            int shift_idx = sy * Nx + sx;
            shifted_result(shift_idx) = result(orig_idx);
        }
    }

    // Cleanup
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return shifted_result;
}



} // end of the namespace

#endif
