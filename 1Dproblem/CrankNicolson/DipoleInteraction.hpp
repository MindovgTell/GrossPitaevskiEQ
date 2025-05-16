#include <fftw3.h>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <algorithm>

// // Complementary error function approximation (from Numerical Recipes)
// double erfc_approx(double x) {
//     double t, z, ans;
//     z = std::abs(x);
//     t = 1.0 / (1.0 + 0.5 * z);
//     ans = t * exp(-z*z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
//           t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
//           t * (-0.82215223 + t * 0.17087277)))))))));
//     return x >= 0.0 ? ans : 2.0 - ans;
// }

class DipolarInteraction2D {
private:
    int Nx, Ny;
    double Lx, Ly;
    double lz;
    double Cdd;
    double dx, dy;
    
    // FFTW plans and arrays
    fftw_plan plan_forward, plan_backward;
    double* density_real;
    fftw_complex* density_fourier;
    fftw_complex* potential_fourier;
    double* potential_real;
    
    // Precomputed U_tilde
    Eigen::MatrixXd U_tilde;
    
public:
    DipolarInteraction2D(int nx, int ny, double lx, double ly, double confinement_length, double interaction_strength)
        : Nx(nx), Ny(ny), Lx(lx), Ly(ly), lz(confinement_length), Cdd(interaction_strength) {
        
        dx = Lx / Nx;
        dy = Ly / Ny;
        
        // Allocate memory for FFTW
        density_real = (double*)fftw_malloc(sizeof(double) * Nx * Ny);
        density_fourier = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx * (Ny/2 + 1));
        potential_fourier = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx * (Ny/2 + 1));
        potential_real = (double*)fftw_malloc(sizeof(double) * Nx * Ny);
        
        // Create FFTW plans
        plan_forward = fftw_plan_dft_r2c_2d(Nx, Ny, density_real, density_fourier, FFTW_MEASURE);
        plan_backward = fftw_plan_dft_c2r_2d(Nx, Ny, potential_fourier, potential_real, FFTW_MEASURE);
        
        // Precompute U_tilde
        precompute_U_tilde();
    }
    
    ~DipolarInteraction2D() {
        fftw_destroy_plan(plan_forward);
        fftw_destroy_plan(plan_backward);
        fftw_free(density_real);
        fftw_free(density_fourier);
        fftw_free(potential_fourier);
        fftw_free(potential_real);
    }
    
    void precompute_U_tilde() {
        U_tilde.resize(Nx, Ny/2 + 1);
        
        for (int i = 0; i < Nx; ++i) {
            double kx = (i < Nx/2) ? (2.0 * M_PI * i / Lx) : (2.0 * M_PI * (i - Nx) / Lx);
            
            for (int j = 0; j < Ny/2 + 1; ++j) {
                double ky = 2.0 * M_PI * j / Ly;
                double k_perp = std::sqrt(kx*kx + ky*ky);
                double kappa = k_perp * lz / std::sqrt(2.0);
                
                if (k_perp == 0.0) {
                    U_tilde(i, j) = 2.0 * Cdd / 3.0;
                } else {
                    U_tilde(i, j) = (Cdd / 3.0) * (2.0 - 3.0 * M_PI * kappa * std::exp(kappa*kappa) * std::erfc(kappa)); //* erfc_approx(kappa));
                }
            }
        }
    }
    
    void compute_DDI_term(const Eigen::MatrixXcd& psi, Eigen::MatrixXd& Phi_DDI) {
        // Copy |psi|^2 to real array

        assert(psi.rows() == Nx && psi.cols() == Ny && "psi dimensions must match Nx x Ny");


        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                // density_real[i*Ny + j] = std::norm(psi(i, j));
            }
        }
        
        // Forward FFT
        fftw_execute(plan_forward);
        
        // Multiply by U_tilde in Fourier space
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny/2 + 1; ++j) {
                potential_fourier[i*(Ny/2 + 1) + j][0] = U_tilde(i, j) * density_fourier[i*(Ny/2 + 1) + j][0];
                potential_fourier[i*(Ny/2 + 1) + j][1] = U_tilde(i, j) * density_fourier[i*(Ny/2 + 1) + j][1];
            }
        }
        
        // Inverse FFT
        fftw_execute(plan_backward);
        
        // Copy result to Eigen matrix and normalize
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                Phi_DDI(i, j) = potential_real[i*Ny + j] / (Nx * Ny);
            }
        }
    }
};