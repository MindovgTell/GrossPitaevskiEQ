# Gross-Pitaevskii Equation Solver

This project provides tools for numerically solving the extended Gross-Pitaevskii equation, a nonlinear Schrödinger equation commonly used in Bose-Einstein condensate simulations. It implements both the Crank-Nicolson finite difference scheme and the Split Step method for efficient and accurate approximation of the ground state of the system and time evolution, respectively.

## Mathematical Formulation

The Gross-Pitaevskii equation (GPE) with dipole-dipole interaction and quantum fluctuations is given by:

$$
i\hbar \frac{\partial \psi(\mathbf{r}, t)}{\partial t} = \left[ -\frac{\hbar^2}{2m} \nabla^2 + V_{\text{ext}}(\mathbf{r}) + g |\psi(\mathbf{r}, t)|^2 + \Phi_{\text{dd}}(\mathbf{r}, t) + \gamma_{\text{QF}} |\psi(\mathbf{r}, t)|^3 \right] \psi(\mathbf{r}, t)
$$

where:
- $\psi(\mathbf{r}, t)$ is the condensate wavefunction,
- $V_{\text{ext}}(\mathbf{r})$ is the external potential,
- $g$ is the contact interaction strength,
- $\Phi_{\text{dd}}(\mathbf{r}, t)$ is the dipole-dipole interaction.
- $\gamma_{\text{QF}}$ is the quantum fluctuation coefficient.

## Installation and Build

## Installation and Build

### Requirements
- [Eigen](https://eigen.tuxfamily.org/) (C++ library for linear algebra)
- [FFTW](http://www.fftw.org/) (Fast Fourier Transform library)
- Python 3.8+ (for post-processing)
- [matplotlib](https://matplotlib.org/) (Python plotting library)

### Clone and Build
```bash
git clone https://github.com/yourusername/GrossPitaevski.git
cd GrossPitaevski
mkdir build && cd build
cmake ..
cmake --build .
```

## Example of usage

To perform a calculation, first initialize a `Grid` object to define the spatial domain. Next, create a `WaveFunction` object on this grid to represent the condensate state. Finally, choose and apply a solver—either the Crank-Nicolson or Split Step Fourier algorithm—to evolve the wavefunction in time.  

```cpp
#include "GPES/Grid.hpp"
#include "GPES/WaveFunction.hpp"
#include "GPES/CrankNicolson.hpp"

int main() {
    try {
        // Initialize grid 
        GPES::Grid<Dimension::Two> grid(300,300, -10, -10); 
        // Setting up harmonic potential
        grid.set_harmonic_potential(1,1); 
        // Tunning the frequency of the trap
        grid.set_z_freq(5); 
        // Defining the wavefunction
        GPES::WaveFunction<Dimension::Two> Phi(grid,0.00607, 0.00882, 10000);
        // Initialize solver for finding ground state
        GPES::CrankNicolson<Dimension::Two> solver(Phi,grid,0.001,0.1, 0.614);
        // Running simulation
        solver.simulation();
        // Save wavefunction
        GPES::save_wavefunction(Phi, "directory_for_saving");
    }
    catch (const std::exception& ex) {           // catches runtime_error and any std::exception
        std::cout << "Error while saving wave-function: " << ex.what() << '\n';
        return;                    // tell the OS something went wrong
    }
    catch (...)                                 // catches non-standard exceptions, just in case
    {
        std::cout << "Unknown error.\n";
        return;
    }
}