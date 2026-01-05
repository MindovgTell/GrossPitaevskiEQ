#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP

#include <memory>

#include "grid.hpp"


namespace gpes
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

public: 
    using VectorType = Eigen::VectorXcd;
    using Scalar     = VectorType::Scalar;
    using ShrdPtrGridDim1 = std::shared_ptr<const Grid<Dimension::One>>;

    static constexpr Dimension Dim = Dimension::One;

private:
    Eigen::VectorXcd data_;
    ShrdPtrGridDim1 gridptr_;

public:
    WaveFunction() = default; 

    explicit WaveFunction(ShrdPtrGridDim1 grid) : gridptr_(std::move(grid)) {
        if (!gridptr_) throw std::runtime_error("WaveFunction needs a grid");
        int size = gridptr_->size();
        data_   = VectorType::Zero(size);
    }

    WaveFunction(ShrdPtrGridDim1 grid, const Eigen::VectorXcd& vec): data_(vec), gridptr_(std::move(grid)) {
        if (!gridptr_) throw std::runtime_error("WaveFunction needs a grid");
    }

    // Copy and Move constructors
    WaveFunction(const WaveFunction&)            = default;
    WaveFunction(WaveFunction&&) noexcept        = default;

    // Copy and Move assignment operators
    WaveFunction& operator=(const WaveFunction&) = default;
    WaveFunction& operator=(WaveFunction&&) noexcept = default;

    const ShrdPtrGridDim1& grid() const { return gridptr_; }

    // Destructor
    ~WaveFunction() = default; 

    //Overloaded operators
    Scalar& operator()(int i) { return data_(i); }
    Scalar operator()(int i) const { return data_(i); }

    // TODO: overload few more operators for using wavefunction almost as eigen Vector



    // Thomas-Fermi ansatz
    // x -- point on the grid
    // Num -- number of particles
    // g_scat -- scattering constant
    double thomas_fermi_state(double x, unsigned int Num, double g_scat);
    // Setting the wave function to Thomas-Fermi state
    // Num -- number of particles
    // g_scat -- scattering constant
    Eigen::VectorXcd TM_state(unsigned int Num, double g_scat);

    //Gauss wave function
    std::complex<double> gauss_wave_packet(double sigma_x, double x, double x_c);
    double square_func(double x, unsigned int Num, double g_scat);

    //Initialization of starting 1D state
    void set_state_TF(double x_c, unsigned int Num, double g_scat);
    void set_state_Square(double x_c, unsigned int Num, double g_scat);
    void set_state_Gauss(double x_c, double sigma_x, unsigned int Num);

    void set_vec(Eigen::VectorXcd& vec) { data_ = vec; }

    //Getters
    Eigen::VectorXcd& get_wavefunction() { return data_; }

    // Acceptors function

    void setZero() { data_.setZero(); }

    VectorType& vec() noexcept {return data_;}
    const VectorType& vec() const noexcept {return data_;}

    Scalar* data() noexcept { return data_.data(); }
    const Scalar* data() const noexcept { return data_.data(); }

    Eigen::Index size() const noexcept {return data_.size();}

    void normalize(Eigen::VectorXcd& vec, unsigned int Num);
    Eigen::VectorXd prob(Eigen::VectorXcd &vec);
    Eigen::VectorXd prob();
};


double WaveFunction<Dimension::One>::square_func(double x, unsigned int Num, double g_scat){
    double R_tf = std::cbrt(1.5 * Num * g_scat);
    double high = Num / R_tf;
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

double WaveFunction<Dimension::One>::thomas_fermi_state(double x,  unsigned int Num, double g_scat){
    double R_tf = std::cbrt(1.5 * Num * g_scat); // mult times _g_scattering
    double potential = R_tf * R_tf * 0.5;
    double out = potential * (1 - std::pow(x/R_tf,2.));
    if(out > 0)
        return std::sqrt(out);
    else 
        return 0;
}

Eigen::VectorXcd WaveFunction<Dimension::One>::TM_state(unsigned int Num, double g_scat){
    size_t size = gridptr_->size();
    double x = gridptr_->start();
    double step = gridptr_->step();

    Eigen::VectorXcd U(size);

    double x = gridptr_->start();
    double step = gridptr_->step();
    std::complex<double> psum = 0;

    for(int i = 0; i < size; ++i){
        x += step;
        std::complex<double> c = thomas_fermi_state(x, Num, g_scat); // Add parameters for Thomas-Fermi function
        U(i) = c;
        psum += std::norm(c);
    }
    std::complex<double> normalization_factor = std::sqrt(Num) / std::sqrt(psum * step);

    U *= std::abs(normalization_factor);
    return U;
}

void WaveFunction<Dimension::One>::set_state_TF(double x_c, unsigned int Num, double g_scat){
    size_t size = gridptr_->size();
    double x = gridptr_->start();
    double step = gridptr_->step();

    Eigen::VectorXcd U(size);
    double psum = 0;
    for(int i = 0; i < size; ++i){
        std::complex<double> c = thomas_fermi_state(x - x_c, Num, g_scat); // Add parameters for Thomas-Fermi function
        U(i) = c;
        psum += std::norm(c);

        x += step;
    }

    std::complex<double> normalization_factor = std::sqrt(Num) / std::sqrt(psum * step);
    data_ = U * std::abs(normalization_factor);
}

void WaveFunction<Dimension::One>::set_state_Square(double x_c, unsigned int Num, double g_scat){
    size_t size = gridptr_->size();
    double x = gridptr_->start();
    double step = gridptr_->step();

    Eigen::VectorXcd U(size);
    double psum = 0;
    for(int i = 0; i < size; ++i){
        std::complex<double> c = square_func(x, Num, g_scat);
        U(i) = c;
        psum += std::norm(c);

        x += step;
    }

    std::complex<double> normalization_factor = std::sqrt(Num) / std::sqrt(psum * step);
    data_ = U * std::abs(normalization_factor);
}

void WaveFunction<Dimension::One>::set_state_Gauss(double x_c, double sigma_x, unsigned int Num){
    size_t size = gridptr_->size();
    double x = gridptr_->start();
    double step = gridptr_->step();

    Eigen::VectorXcd U(size);
    double psum = 0;

    for(int i = 0; i < size; ++i){
        std::complex<double> c = gauss_wave_packet(sigma_x, x, x_c);
        U(i) = c;
        psum += std::norm(c);
        x += step;
    }

    std::complex<double> normalization_factor = std::sqrt(Num) / std::sqrt(psum * step);
    data_ = U * std::abs(normalization_factor);
} 

// Returning the probability density vector
Eigen::VectorXd WaveFunction<Dimension::One>::prob(Eigen::VectorXcd &vec) {
    int size = vec.size();
    Eigen::VectorXd pr(size);
    for(int i = 0; i < size; ++i){
        pr(i) = std::norm(vec(i));
    }
    return pr;
}

Eigen::VectorXd WaveFunction<Dimension::One>::prob() {
    size_t size = gridptr_->size();

    Eigen::VectorXd pr(size);
    for(int i = 0; i < size; ++i){
        pr(i) = std::norm(data_(i));
    }
    return pr;
}

void WaveFunction<Dimension::One>::normalize(Eigen::VectorXcd &vec, unsigned int Num){
    double step = gridptr_->step();
    int size = vec.size();
    std::complex<double> psum = 0;
    for(int i = 0; i < size; ++i){
        psum += std::norm(vec(i));
    }
    std::complex<double> normalization_factor = std::sqrt(Num) / std::sqrt(psum * step); // 
    vec *= std::abs(normalization_factor);
}



//********************************/***********/********************************//
//                                                                             //
//************************/Two dimensional wavefunction/***********************//
//                                                                             //
//********************************/***********/********************************//

template<>
class WaveFunction<Dimension::Two> {
public: 
    using VectorType = Eigen::VectorXcd;
    using Scalar = VectorType::Scalar;
    using ShrdPtrGridDim2 = std::shared_ptr<const Grid<Dimension::Two>>;

private:
    VectorType data_;
    ShrdPtrGridDim2 gridptr_;

public:

    // Constructors 
    WaveFunction() = default; 

    explicit WaveFunction(ShrdPtrGridDim2 grid): gridptr_(std::move(grid)) {
        if (!gridptr_) throw std::runtime_error("WaveFunction needs a grid");
        int sx = gridptr_->size_x();
        int sy = gridptr_->size_y();
        data_   = VectorType::Zero(sx * sy);
    }

    WaveFunction(ShrdPtrGridDim2 grid, VectorType& vec): gridptr_(std::move(grid)) {
        if (!gridptr_) throw std::runtime_error("WaveFunction needs a grid");
        int sx = gridptr_->size_x();
        int sy = gridptr_->size_y();
        data_   = VectorType::Zero(sx * sy);
    }

    // Copy and Move constructors
    WaveFunction(const WaveFunction&)            = default;
    WaveFunction(WaveFunction&&) noexcept        = default;

    // Copy and Move assignment operators
    WaveFunction& operator=(const WaveFunction&) = default;
    WaveFunction& operator=(WaveFunction&&) noexcept = default;

    const ShrdPtrGridDim2& grid() const { return gridptr_; }

    // Destructor
    ~WaveFunction() = default; 

    //Overloaded operators
    Scalar& operator()(int i) { return data_(i); }
    Scalar operator()(int i) const{ data_(i); }

    WaveFunction& operator*=(Scalar a) {
        data_ *= a;
        return *this;
    }

    WaveFunction& operator*=(double a) {
        data_ *= a;
        return *this;
    } 


    // TODO: overload few more operators for using wavefunction almost as eigen Vector



    int get_index(int i, int j){ 
        int sx = gridptr_->size_x();
        return i*sx + j; 
    }
    
    //Thomas-Fermi anzats
    double thomas_fermi_radius_x(const double &wx, const double &wy, unsigned int _Num);
    double thomas_fermi_radius_y(const double &wx, const double &wy, unsigned int _Num);
    double chem_potential(const double &wx, const double &wy, unsigned int _Num);
    double thomas_fermi_state(double x, double y, unsigned int number_of_prt);
    VectorType TF_state(unsigned int number_of_prt);

    std::complex<double> gauss_wave_packet(double x, double y, double x_c, double y_c, double sigma_x, double sigma_y);

    //Initialization of starting 1D state
    void set_state_TF(double x_c, double y_c, unsigned int number_of_prt);
    void set_state_Gauss(double x_c, double y_c, double sigma_x, double sigma_y, unsigned int number_of_prt);

    void set_vec(Eigen::VectorXcd& vec) {   
        int sx = gridptr_->size_x();
        int sy = gridptr_->size_y();
        if(vec.size() == sx * sy) 
            data_ = vec;
    }

    //Getters
    VectorType& vec() noexcept {return data_;}
    const VectorType& vec() const noexcept {return data_;}

    Scalar* data() noexcept { return data_.data(); }
    const Scalar* data() const noexcept { return data_.data(); }

    Eigen::Index size() const noexcept {return data_.size();}

    unsigned int grid_size_x() {return gridptr_->size_x(); }
    unsigned int grid_size_y() {return gridptr_->size_y(); }
    double start_x() { return gridptr_->start_pos_x(); }
    double start_y() { return gridptr_->start_pos_y(); }
    double step_x() { return gridptr_->step_x(); }
    double step_y() { return gridptr_->step_y(); }

    Eigen::VectorXd get_x_slice();
    Eigen::VectorXd get_y_slice();

    // Acceptors function

    void setZero() { data_.setZero(); }

    //Probability density funcitons
    Eigen::VectorXd prob(Eigen::VectorXcd &vec);
    Eigen::VectorXd prob();
    double prob(int index);
};

double WaveFunction<Dimension::Two>::thomas_fermi_radius_x(const double &wx, const double &wy, unsigned int _Num){
    double numerator = 4.0 * wy * _Num;
    double denumerator = std::pow(wx, 3.) * M_PI;
    return std::pow(( numerator / denumerator), 0.25);
}

double WaveFunction<Dimension::Two>::thomas_fermi_radius_y(const double &wx, const double &wy, unsigned int _Num){
    double numerator = 4.0 * wx * _Num;
    double denumerator = std::pow(wy, 3.) * M_PI;
    return std::pow(( numerator / denumerator), 0.25);
}

double WaveFunction<Dimension::Two>::chem_potential(const double &wx, const double &wy, unsigned int _Num){
    double numerator = 1 * wx * wy * _Num; //add g_scattering
    return std::sqrt(numerator / M_PI);
}

double WaveFunction<Dimension::Two>::thomas_fermi_state(double x, double y, unsigned int number_of_prt){

    double wx = gridptr_->omega_x();
    double wy = gridptr_->omega_y();

    double out, R_x, R_y;
    R_x = thomas_fermi_radius_x(wx, wy, number_of_prt);
    R_y = thomas_fermi_radius_y(wx, wy, number_of_prt);

    out = chem_potential(wx, wy, number_of_prt) * (1 - std::pow(x/R_x, 2.) - std::pow(y/R_y, 2.)); // 

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

Eigen::VectorXcd WaveFunction<Dimension::Two>::TF_state(unsigned int number_of_prt){

    int sizex = gridptr_->size_x();
    int sizey = gridptr_->size_y();

    double startx = gridptr_->start_pos_x();
    double starty = gridptr_->start_pos_y();

    double stepx = gridptr_->step_x();
    double stepy = gridptr_->step_y();

    int size = sizex * sizey;
    Eigen::VectorXcd U(size);
    std::complex<double> psum = 0;
    for(int i = 0; i < sizex; ++i){
        double x = startx + i * stepx;
        for(int j = 0; j < sizey; ++j){
            double y = starty + j * stepy;
            //Initial state function
            std::complex<double> c = thomas_fermi_state(x,y, number_of_prt);

            int index = get_index(i,j);
            U(index) = c;
            psum += std::norm(c);
        }
    }

    std::complex<double> normalization_factor = std::sqrt(number_of_prt) / std::sqrt(psum * stepx * stepy);
    U *= std::abs(normalization_factor);
    return U;
}

// Function for initializing wave function in Thomas-Fermi limit without ddi 
void WaveFunction<Dimension::Two>::set_state_TF(double x_c, double y_c, unsigned int number_of_prt){

    int sizex = gridptr_->size_x();
    int sizey = gridptr_->size_y();

    double startx = gridptr_->start_pos_x();
    double starty = gridptr_->start_pos_y();

    double stepx = gridptr_->step_x();
    double stepy = gridptr_->step_y();

    int size = sizex * sizey;
    Eigen::VectorXcd U(size);
    std::complex<double> psum = 0;
    
    for(int i = 0; i < sizex; ++i){
        double x = startx + i * stepx;
        
        for(int j = 0; j < sizey; ++j){
            double y = starty + j * stepy;
            //Initial state function
            std::complex<double> c = thomas_fermi_state(x,y, number_of_prt);
            int index = get_index(i,j);
            U(index) = c;
            psum += std::norm(c);
        }
    }

    std::complex<double> normalization_factor = std::sqrt(number_of_prt) / std::sqrt(psum * stepx * stepy);
    U *= std::abs(normalization_factor);
    data_ = U;
}

// Function for initializing wave function in Thomas-Fermi limit without ddi 
void WaveFunction<Dimension::Two>::set_state_Gauss(double x_c, double y_c, double sigma_x, double sigma_y, unsigned int number_of_prt){

    int sizex = gridptr_->size_x();
    int sizey = gridptr_->size_y();

    double startx = gridptr_->start_pos_x();
    double starty = gridptr_->start_pos_y();

    double stepx = gridptr_->step_x();
    double stepy = gridptr_->step_y();


    int size = sizex * sizey;
    Eigen::VectorXcd U(size);
    std::complex<double> psum = 0;
    for(int i = 0; i < sizex; ++i){
        double x = startx + i * stepx;
        for(int j = 0; j < sizey; ++j){
            double y = starty + j * stepy;
            //Initial state function
            std::complex<double> c = gauss_wave_packet(x, y, x_c, y_c, sigma_x, sigma_y); // ,p_x, p_y
            int index = get_index(i,j);
            U(index) = c;
            psum += std::norm(c);
        }
    }

    std::complex<double> normalization_factor = std::sqrt(number_of_prt) / std::sqrt(psum * stepx * stepy);
    U *= std::abs(normalization_factor);
    data_ = U;
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
    int sizex =  gridptr_->size_x();
    int sizey =  gridptr_->size_y();
    int size = sizex * sizey;
    Eigen::VectorXd pr(size);
    for(int i = 0; i < size; ++i){
        pr(i) = std::norm(data_(i));
    }
    return pr;
}

double WaveFunction<Dimension::Two>::prob(int index){
    int sizex =  gridptr_->size_x();
    int sizey =  gridptr_->size_y();
    int size = sizex * sizey;
    if(index > 0 && index < size)
        return std::norm(data_(index));
    else 
        return 0;
}

Eigen::VectorXd WaveFunction<Dimension::Two>::get_x_slice(){
    int sizex =  gridptr_->size_x();
    int sizey =  gridptr_->size_y();
    int y_center = sizey / 2;
    Eigen::VectorXd slice(sizex);
    for(int i = 0; i < sizex; ++i){
        slice(i) = std::norm(data_(y_center * sizex + i));
    }
    return slice;
}


Eigen::VectorXd WaveFunction<Dimension::Two>::get_y_slice(){
    int sizex = gridptr_->size_x();
    int sizey =  gridptr_->size_y();
    int x_center = sizex / 2;
    Eigen::VectorXd slice(sizey);
    for(int j = 0; j < sizey; ++j){
        slice(j) = std::norm(data_(j * sizex + x_center));
    }
    return slice;
}

} // end of the namespace

#endif
