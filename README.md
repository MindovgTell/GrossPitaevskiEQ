# Gross-Pitaevskii Equation Solver

This project provides tools for numerically solving the Gross-Pitaevskii equation, a nonlinear Schrödinger equation commonly used in Bose-Einstein condensate simulations. It implements both the Crank-Nicolson finite difference scheme and the Split Step method for efficient and accurate time evolution.

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

1. **Clone the repository:**
    ```bash
    git clone https://github.com/yourusername/GrossPitaevski.git
    cd GrossPitaevski
    ```

2. **Install dependencies:**
    - Make sure you have Python 3.8+ installed.
    - Install required Python packages:
      ```bash
      pip install -r requirements.txt
      ```

3. **Build and run:**
    - Run the main solver script:
      ```bash
      python main.py
      ```

## Example of usage

To perform a calculation, first initialize a `Grid` object to define the spatial domain. Next, create a `WaveFunction` object on this grid to represent the condensate state. Finally, choose and apply a solver—either the Crank-Nicolson or Split Step Fourier algorithm—to evolve the wavefunction in time.

```cpp
#include "GPES/Grid.hpp"
#include "GPES/WaveFunction.hpp"
#include "GPES/CrankNicolson.hpp"
#include "GPES/SplitStep.hpp"

int main() {
    // Initialize grid
    Grid grid(-10.0, 10.0, 256);

    // Define initial wavefunction
    WaveFunction psi0(grid, "gaussian");

    // Choose a solver
    SplitStepSolver solver(grid, psi0, parameters);

    // Run the simulation
    solver.evolve(1.0, 0.001);

    return 0;
}```