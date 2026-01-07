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

### Requirements
- [Eigen](https://eigen.tuxfamily.org/) (C++ library for linear algebra)
- [FFTW](http://www.fftw.org/) (Fast Fourier Transform library)
- Python 3.8+ (for post-processing)
- [matplotlib](https://matplotlib.org/) (Python plotting library)
- [spdlog](https://github.com/gabime/spdlog) (Logging information)
- [CLI11] (https://github.com/CLIUtils/CLI11) (CLI tools)

### Clone and Build
```bash
git clone https://github.com/MindowgTell/GrossPitaevskiEQ.git
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

int main(int argc, char **argv) {
    try {
        auto [phys, sim] = parse_cli(argc, argv);
        print_parsed_params(argc, argv);

        const gpes::Dimension Dim = gpes::Dimension::One;

        auto grid_ptr = std::make_shared<gpes::Grid<Dim>>(100,-10);
        grid_ptr->set_harmonic_potential(5);

        gpes::WaveFunction<Dim> Psi(grid_ptr);
        Psi.set_state_Gauss(0, 5, 1000);

        gpes::solvers::CrankNicolson<Dim> stepper(grid_ptr, Psi, phys, sim);

        gpes::TDSESolver<Dim, decltype(stepper)> solver(grid_ptr, stepper, Psi, sim);

    }
    catch (...) {                          
        std::cout << "Unknown error.\n";
        return;
    }
}