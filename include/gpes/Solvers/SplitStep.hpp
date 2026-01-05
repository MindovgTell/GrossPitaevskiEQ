#ifndef SPLITSTEP_HPP
#define SPLITSTEP_HPP

#include "dipoleinteraction.hpp"
#include "wavefunction.hpp"
#include "grid.hpp"

namespace gpes {

    template <Dimension Dim>
    class SplitStep;

    //********************************/***********/********************************//
    //                                                                             //
    //**************************/One dimensional solver/***************************//
    //                                                                             //
    //********************************/***********/********************************//

    template<>
    class SplitStep<Dimension::One> {
    private:
        double dt, N_steps;
     
        Eigen::VectorXcd U_k;
        Eigen::VectorXcd U_v;

        const gpes::Grid<Dimension::One>& _Grid;
        gpes::WaveFunction<Dimension::One>& _Psi;
        
    public: 

        SplitStep(const Grid<Dimension::One>& grid, WaveFunction<Dimension::One>& psi, double delta_t, int Nsteps);

    private:
        void init_kinetic_op();
        void calc_potential_op(const Eigen::VectorXcd& psi);
        void sim();
        
    };


    SplitStep<Dimension::One>::SplitStep(const Grid<Dimension::One>& grid, WaveFunction<Dimension::One>& psi, double delta_t, int Nsteps) : _Grid(grid), _Psi(psi), dt(delta_t), N_steps(Nsteps) {}


    void SplitStep<Dimension::One>::init_kinetic_op() {
        int N = _Grid.size();
        std::vector<double> k(N);
        double dx = _Grid.step();
        double L = N * dx;
        const double dk = 2.0 * M_PI / L;

        U_k.resize(N);

        for (int m = 0; m < N; ++m) {
            int q = (m <= N/2) ? m : m - N;  // integer wave index (positive + negative)
            k[m] = dk * q;                   // physical k
        }

        for(int j = 0; j < N; j++){
            U_k(j) = std::exp(std::complex<double>(0.0, -1.0 * k[j]*k[j] * dt / 4));
        }    
    }

    void SplitStep<Dimension::One>::calc_potential_op(const Eigen::VectorXcd& vec){
        int N = vec.size();

        // calculate_DDI(vec);
        std::vector<double> U = GPES::calculate_DDI_not_FFT(_Grid, vec);

        U_v.resize(N);

        double g_scat = _Psi.get_g_scat();
        double g_lhy  = _Psi.get_g_lhy();

        for(int i = 0; i < N; ++i){

            double U_trap =  _Grid.potential()(i);
            double U_scattering =  g_scat*std::norm(vec(i));
            double U_lhy = g_lhy * std::pow(std::norm(vec(i)), 1.5);


            U_v(i) = std::exp(std::complex(0.0, -dt*(U_trap + U_scattering + U[i] + U_lhy)));
        }
    }

    void SplitStep<Dimension::One>::sim() {

        int N = _Psi.size();

        Eigen::VectorXcd space_wf(N);
        Eigen::VectorXcd momentum_wf(N);


        // Setting up Fourier transform
        fftw_plan   plan_fwd, plan_bwd;
        
        plan_fwd = fftw_plan_dft_1d(
            N,
            reinterpret_cast<fftw_complex*>(space_wf.data()), 
            reinterpret_cast<fftw_complex*>(momentum_wf.data()), 
            FFTW_FORWARD,
            FFTW_MEASURE
        );
        plan_bwd = fftw_plan_dft_1d(
            N,
            reinterpret_cast<fftw_complex*>(momentum_wf.data()), 
            reinterpret_cast<fftw_complex*>(space_wf.data()), 
            FFTW_FORWARD,
            FFTW_MEASURE
        );

        // Time stepping 
        for(int step = 0; step < N_steps; step++){

            // 1) Half kinetic step
            //      Transform wavefunction to momentum space
            fftw_execute(plan_fwd);
            //      Calculate the action of the kinetic operator on the wavefunction in Fourier space
            for(int i = 0; i < N; i++) momentum_wf(i) *= U_k(i);
            //      Return to the real space
            fftw_execute(plan_bwd);
            //      Normalization
            for(int i = 0; i < N; i++) space_wf(i) /= double(N);

            // 2) Full potential step
            calc_potential_op(space_wf);
            for(int i = 0; i < N; i++) space_wf(i) *= U_v(i);

            // 3) Half kinetic step
            //      Transform wavefunction to momentum space
            fftw_execute(plan_fwd);
            //      Calculate the action of the kinetic operator on the wavefunction in Fourier space
            for(int i = 0; i < N; i++) momentum_wf(i) *= U_k(i);
            //      Return to the real space
            fftw_execute(plan_bwd);
            //      Normalization
            for(int i = 0; i < N; i++) space_wf(i) /= double(N);
        }

        _Psi = space_wf;

        fftw_destroy_plan(plan_fwd);
        fftw_destroy_plan(plan_bwd);
    }



    //********************************/***********/********************************//
    //                                                                             //
    //**************************/Two dimensional solver/***************************//
    //                                                                             //
    //********************************/***********/********************************//

    class SplitStep<Dimension::Two> {
    private:
        double dt, N_steps;

        Eigen::VectorXcd U_k;
        Eigen::VectorXcd U_v;
        
        GPES::Grid<Dimension::Two>& _Grid;
        GPES::WaveFunction<Dimension::Two>& _Psi;
    public:
        void sim();

    private:
        void init_kinetic_op();
        void calc_potential_op(const Eigen::VectorXcd& vec);

    };


    void GPES::SplitStep<Dimension::Two>::sim() {
        // Implementation of the 2D split-step method

        int N_x = _Psi.grid_size_x();
        int N_y = _Psi.grid_size_y();


        // Setting up Fourier transform


        

    }



    void SplitStep<Dimension::Two>::calc_potential_op(const Eigen::VectorXcd& vec){
        int N = vec.size();

        // calculate_DDI(vec);
        std::vector<double> U = GPES::calculate_DDI_not_FFT(_Grid, vec);
        

        U_v.resize(N);

        double g_scat = _Psi.get_g_scat();
        double g_lhy  = _Psi.get_g_lhy();

        for(int i = 0; i < N; ++i){

            std::complex<double> U_trap = dt*0.5 * _Grid.potential()(i);
            std::complex<double> U_scattering =  (dt*0.5 * g_scat)*std::norm(vec(i));
            std::complex<double> U_dd = 1.0 * (dt*0.5) * U[i];  //;
            std::complex<double> U_lhy = 1.0 * (dt*0.5) * g_lhy * std::pow(std::norm(vec(i)), 1.5);


            U_v(i) = U_trap + U_scattering + U_dd + U_lhy;
        }
    }






}

#endif