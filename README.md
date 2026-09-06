# Quasi-2D dipolar eGPE and BdG solver

This code computes stationary states of a harmonically trapped dipolar Bose gas by normalized imaginary-time split-step Fourier propagation, then linearizes the same extended Gross-Pitaevskii equation (eGPE) and computes low-energy Bogoliubov-de Gennes (BdG) modes with a matrix-free shift-invert Arnoldi method.

## Model and units

The implementation is dimensionless. Length, frequency, and energy are measured in

`l_x = sqrt(hbar / (m omega_x))`, `omega_x`, and `hbar omega_x`.

Consequently, `a_s` and `a_dd` must be supplied in units of `l_x`, and all trap frequencies are ratios to the reference `omega_x`. The two-dimensional wavefunction is normalized as

`integral dx dy |psi|^2 = N`.

The frozen `z` wavefunction is a harmonic Gaussian with `l_z = 1/sqrt(omega_z)`. For dipoles polarized normal to the plane, the reduced interactions used by the code are

```text
g_2D       = sqrt(8 pi) a_s / l_z
g_dd,2D    = sqrt(8 pi) a_dd / l_z
U_dd(k)    = g_dd,2D [2 - 3 sqrt(pi) q exp(q^2) erfc(q)]
q          = |k| l_z / sqrt(2)
gamma_QF   = [128/(3 pi^(1/4))] sqrt(2/5) a_s^(5/2) Q5(a_dd/a_s) / l_z^(3/2)
```

The high-momentum kernel is evaluated through a stable scaled complementary error function. The real part of the dipolar LHY angular integral is used for `Q5`. The eGPE is

```text
i dpsi/dt = [-laplacian/2 + V + g_2D n + Phi_dd + gamma_QF n^(3/2)] psi,
n = |psi|^2,  Phi_dd = U_dd * n.
```

The local-density LHY prescription is the usual three-dimensional dipolar expression projected onto the frozen Gaussian. Strictly two-dimensional beyond-mean-field corrections are nonuniversal, so this model should be described as a frozen-axis quasi-2D eGPE, not a universal 2D theory.

## BdG definition

The perturbation convention is

```text
psi(t) = exp(-i mu t) [psi0 + u exp(-i omega t) + v* exp(+i omega* t)].
```

For `n0=|psi0|^2`, the implemented upper block is

```text
A u + B v + psi0 delta_Phi_dd,
A = T + V - mu + 2 g_2D n0 + Phi_dd + (5/2) gamma_QF n0^(3/2),
B = g_2D psi0^2 + (3/2) gamma_QF sqrt(n0) psi0^2,
delta_Phi_dd = U_dd * (psi0* u + psi0 v).
```

The lower block is its bosonic particle-hole partner `-B* u - A* v - psi0* delta_Phi_dd`. Modes are normalized with `integral (|u|^2-|v|^2) dx dy = +/-1`. The run is rejected rather than saved as valid when the inner shift-invert systems or final eigenpairs fail their residual checks. `run_metadata.json` also records the stationary and Goldstone residuals.

## Requested 40-point sweep

[parameters/sweep_40.tsv](parameters/sweep_40.tsv) contains 8 evenly spaced scattering lengths from `0.001` to `0.009`, crossed with 5 aspect ratios from `omega_y/omega_z=1` to `10`, for exactly 40 runs. It fixes `a_dd=0.008843`, `N=80000`, `omega_x=1`, and `omega_z=8`. The production Slurm script supplies a `1024x1024` grid on `[-25,25) x [-8,8)`.

Important: the implemented dimensional reduction integrates out `z`, which is controlled only when `omega_z` is substantially larger than both `omega_x` and `omega_y`. The literal requested range `omega_y/omega_z >= 1` does not satisfy that quasi-2D condition. If the intended ratio was `omega_z/omega_y`, or if a different physical axis is frozen, change the table before interpreting the sweep physically. The program records the literal ratio in every run's metadata to prevent an axis-label mistake from going unnoticed.

A single noisy initial state is used per table row to permit broken-symmetry droplet solutions. Near first-order transitions, one initial condition does not prove the global minimum. Publication-quality phase boundaries require several seeds and structured droplet initial states per parameter point, followed by selection of the lowest converged energy.

## Build and run on Slurm

Dependencies are CMake 3.24+, a C++20 compiler, Eigen 3.3+, and FFTW3. Spectra 1.2, spdlog, and fmt headers are vendored for network-independent compute-node builds. The build script loads the Helios module defaults currently encoded in the script; override `GPES_MODULES` if those module names differ.

```bash
bash scripts/submit_sweep.sh
```

The submission wrapper first submits the build, then submits the 40-task array with an `afterok` dependency. It limits the default concurrency to 8 tasks; each task requests 32 GiB and writes all numerical results and Slurm logs below:

```text
$SCRATCH/GPES/sweep_YYYYMMDD_HHMMSS_XXXXXX/
```

Every default submission creates a new timestamped sweep directory containing the 40 parameter-set directories and `slurm/` logs. The random suffix prevents collisions between submissions in the same second. The wrapper prints the exact directory before submitting jobs.

Useful submission overrides include:

```bash
MAX_CONCURRENT=4 BDG_MODES=32 bash scripts/submit_sweep.sh
RUN_ROOT=/alternate/path bash scripts/submit_sweep.sh parameters/sweep_40.tsv
```

To resume an existing sweep, explicitly reuse its directory:

```bash
RUN_ROOT="$SCRATCH/GPES/sweep_20260906_120000_ABC123" bash scripts/submit_sweep.sh
```

With that explicit `RUN_ROOT`, interrupted tasks resume from atomically replaced `wavefunction/wavefunction.csv` checkpoints. A task with a `COMPLETED` marker is skipped on resubmission. Do not resubmit into that directory while its previous tasks are still running; concurrent writers are not locked.
The adjacent checkpoint parameter manifest is checked before every resume, so a state cannot silently be reused with a changed grid or Hamiltonian.

## Outputs

Each run directory contains:

```text
wavefunction/wavefunction.csv       final state and restart checkpoint
wavefunction/parameters.txt         checked restart parameter manifest
wavefunction/energy_history.csv     convergence history
bdg/spectrum.csv                    complex frequencies, norms, and residuals
bdg/excited_states/mode_*.csv       full u,v fields
images/ground_state_density.bmp     dependency-free density image
images/energy_history.svg           convergence image
images/bdg_spectrum.svg             real and imaginary frequencies
images/bdg_modes/*.bmp              density-fluctuation images for every mode
run_metadata.json                   parameters and numerical quality diagnostics
COMPLETED                           successful-run marker
```

Embedded Python/matplotlib plotting is optional. The dependency-free BMP and SVG products are always generated, including in the default headless cluster build.

At `1024x1024`, one wavefunction CSV is about 0.1 GiB and 24 full BdG mode CSVs are typically several GiB per parameter point. Budget roughly 0.15-0.25 TiB of scratch space for all 40 runs, depending on numeric text widths. `--no-mode-fields` substantially reduces storage, but should only be used when eigenvalues and mode images are sufficient for the analysis.

## Local verification

```bash
cmake -S . -B build -DGPES_ENABLE_PLOTTING=OFF -DGPES_BUILD_TESTS=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The numerical regression test checks the high-`k` dipolar limit and compares the complete BdG action against a central finite-difference derivative of the nonlinear eGPE, including contact, LHY, and nonlocal dipolar terms.
