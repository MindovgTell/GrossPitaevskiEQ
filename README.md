# Quasi-2D dipolar eGPE and BdG solver

The ground-state program computes stationary states of a harmonically trapped dipolar Bose gas by normalized imaginary-time split-step Fourier propagation. A separate BdG program computes low-energy modes of nonrotating stationary states with a matrix-free shift-invert Arnoldi method.

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
i dpsi/dt = [-laplacian/2 + V + g_2D n + Phi_dd + gamma_QF n^(3/2) - Omega L_z] psi,
n = |psi|^2,  Phi_dd = U_dd * n.
L_z = -i (x d/dy - y d/dx),  E_Omega = E - Omega <L_z>.
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

The lower block is its bosonic particle-hole partner `-B* u - A* v - psi0* delta_Phi_dd`. Modes with a nonzero Bogoliubov norm are normalized with `integral (|u|^2-|v|^2) dx dy = +/-1`; zero-norm modes cannot use this normalization. Failed inner solves or eigenpair residual checks never produce a validated result. The separate `bdg` program retains returned candidates and failure diagnostics; the legacy combined executable rejects them without saving the candidate spectrum.

## Separate ground-state and BdG programs

`app/groundstate.cpp` and `app/bdg.cpp` link against `gpes::gpes`. Both are headless: they produce BMP/SVG images without requiring Python. The Slurm sweep uses `groundstate` only. `bdg` remains a separate, nonrotating program: a saved state with nonzero `rotation_frequency` is not a valid input for its current operator.

### Build

Build all programs inside the repository:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DGPES_ENABLE_PLOTTING=OFF
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

Alternatively, install the current library and build `app` as a separate CMake project:

```bash
cmake -S . -B build/library -DCMAKE_BUILD_TYPE=Release \
  -DGPES_BUILD_APP=OFF -DGPES_BUILD_TESTS=OFF -DGPES_ENABLE_PLOTTING=OFF
cmake --build build/library --parallel 2
cmake --install build/library --prefix "$PWD/build/install"
cmake -S app -B build/apps -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$PWD/build/install"
cmake --build build/apps --parallel 2
```

The first build places executables in `bin/`; the independent build uses `build/apps/bin/`. The standalone app build uses the installed GPES headers/library, and the repository's vendored CLI11 (or `-DGPES_CLI11_INCLUDE_DIR=/path/to/CLI11`). Reinstall GPES after changing its public headers before rebuilding external consumers.

### Calculate one state

For one local `200x200` ground-state example, build and run with:

```bash
bash scripts/build_local.sh
bash scripts/run_solver.sh
```

`build_local.sh` builds only `bin/groundstate` and disables optional Python plotting by default. `run_solver.sh --build` performs both commands. The example uses `a_s=0.009`, `a_dd=0.008843`, `N=80000`, `omega_y/omega_x=2`, `omega_z=8`, and `200000` maximum imaginary-time steps. Its unique result directory is printed and lies under `runs/`; no previous run is overwritten. A `200x200` result is for local workflow checks, not a resolution-converged physical result. For a short check, run `bash scripts/run_solver.sh -- --max-ground-steps 10`. A step-limited run normally exits with status 2 and saves `NOT_CONVERGED`, a checkpoint, and images; use `-- --resume /path/to/previous/run` for another attempt in a new directory. Other ground-state CLI options may be passed after `--`, overriding the example defaults.

```bash
GS="$SCRATCH/GPES/my_state_001"  # Must not already exist.
./bin/groundstate \
  --a-s 0.004 --a-dd 0.008843 --num-particles 80000 \
  --w-x 1 --lambda 0.8 --w-z 5 --rotation-frequency 0.1 \
  --modulation-file parameters/three_peak_seed.tsv \
  --nx 1024 --ny 1024 --x-half-length 25 --y-half-length 25 \
  --max-ground-steps 200000 --output-dir "$GS"
```

Without `--output-dir`, every invocation creates a unique `$SCRATCH/GPES/groundstate_<UTC-time>_<suffix>` directory, or `runs/groundstate_<...>` when `SCRATCH` is unset. No existing destination is overwritten. The grid defaults to `1024x1024`; `--nx 100 --ny 100` changes it for a workflow test, not a resolution-converged production calculation. `--lambda` means `omega_x/omega_y` and sets `omega_y=omega_x/lambda`; alternatively supply `--w-y`. Changing the trap does not automatically change the box size. The example interaction parameters and seed are illustrative, not a calibrated phase boundary.

The default initial state is a Gaussian with width parameters 2.5 in both directions, multiplied by complex noise (`--seed-noise 0.001`, `--random-seed 12345`). Set `--background-width-x` and `--background-width-y` for an anisotropic background. `--modulation-file` adds the proposal's Gaussian peaks *before* noise. Each noncomment row contains `x y sigma_x sigma_y relative_amplitude`; the amplitude is measured against the original Gaussian's peak, and the state is renormalized to `N`. The 2D library method is `WaveFunction<Dimension::Two>::add_gaussian_modulations(components)`, where `components` is a vector of `GaussianComponent` arrays in the same field order. The supplied [three-peak example](parameters/three_peak_seed.tsv) is a zigzag seed. All options are listed by `./bin/groundstate --help`.

Energy convergence now uses

```text
|E_new - E_old| <= absolute_energy_tol + relative_energy_tol * max(|E_new|, |E_old|)
```

Defaults are `--relative-energy-tol 1e-8` and `--absolute-energy-tol 1e-12`. This means approximately eight-significant-digit energy stability, not eight decimal places, literal equality of rounded digits, or eight-digit physical accuracy. Unlike the old `max(...,1)` scaling, small nonzero energies use their actual scale plus the explicit absolute floor. Larger tolerances are less strict. Energy and `--stationary-residual-tol 1e-4` must pass `--convergence-checks-required 3` consecutive checks, spaced `--convergence-check-interval 100` steps apart. The step starts at `--dt 0.001` and can refine down to `--min-ground-dt 1e-4` when energy stalls before stationarity is reached. Successful output also requires boundary/peak density no larger than `--boundary-density-tol 1e-8`.

The state is checkpointed atomically every `--snapshot-interval 5000` steps and at the end. A step-limit failure still saves the final field, diagnostics, and images. To continue it, with the same grid and physics but a new output directory:

```bash
./bin/groundstate --resume "$GS" --max-ground-steps 200000
```

Restart parameters are read from the saved manifest. Physics/grid overrides with `--resume` are rejected. Iteration counters, convergence history, time-step controls, and tolerances start anew for the new invocation; supply numerical overrides explicitly as needed. The original run stays untouched.

### Calculate its spectrum

```bash
./bin/bdg --input-dir "$GS" --modes 24 --ncv 64 --shift 0.05
```

`bdg` reads `wavefunction/parameters.txt` and `wavefunction/wavefunction.csv`, checks grid coordinates, field finiteness, and normalization, then recomputes the chemical potential and stationary residual using its nonrotating eGPE. It never evolves or renormalizes the input. Use it only for `Omega=0` states; manifests with nonzero rotation are rejected. A previous `COMPLETED` marker is not required: the loaded state is checked.

Every attempt creates `$GS/bdg_attempts/attempt_<UTC-time>_<suffix>/`, or a new explicit `--output-dir`. Retrying BdG therefore never repeats ground-state propagation or overwrites an earlier spectrum.

| Option | Default | Meaning |
| --- | --- | --- |
| `--modes` | 24 | Eigenpairs near the target shift, not necessarily 24 positive physical excitations. |
| `--ncv` | 64 | Arnoldi basis size; effective size is `min(2*nx*ny, max(ncv, 2*modes+1))`. Larger values cost memory. |
| `--max-iter`, `--tol` | 1500, `1e-7` | Outer Arnoldi iteration limit and tolerance. |
| `--shift`, `--shift-imag` | 0.05, 0 | Real and imaginary parts of the spectral target. Shift-invert selects modes near it. |
| `--linear-max-iter` | 1000 | Maximum inner GMRES iterations per inverse application. |
| `--linear-restart` | 48 | GMRES restart length; increasing it costs memory. |
| `--linear-tol` | `1e-8` | Inner shifted-system relative residual tolerance. |
| `--residual-tol` | `1e-4` | Independent direct residual threshold for returned eigenpairs. |
| `--stationary-residual-tol` | `1e-4` | Recomputed input eGPE residual threshold. |
| `--boundary-density-tol` | `1e-8` | Input boundary/peak density threshold. |
| `--dense-max-dim` | 0 | Tiny-grid dense inverse only if `2*nx*ny` is at most this value. Maximum allowed is 4096; leave disabled for production. |
| `--no-mode-fields` | off | Omit full `u,v` CSVs, but keep mode images and eigenvalues. |
| `--allow-unconverged` | off | Bypass input stationarity/boundary checks for diagnostics; such results remain incomplete. |

The implementation applies the BdG operator through FFTs and pointwise operations. Shift-invert applies `(L - sigma I)^(-1)` using restarted GMRES with the existing Fourier preconditioner; Spectra performs the outer Arnoldi iteration. This program separation does not change that eigensolver or implement phase-mode deflation. Inspect `linear_failures`, `max_linear_residual`, `max_eigen_residual`, and `stationary_relative_residual` in `bdg_status.txt` before tuning controls. Raising only the outer iteration limit cannot remedy inaccurate inner solves.

To process a saved sweep in one process:

```bash
./bin/bdg --input-dir "$SCRATCH/GPES/sweep_DIRECTORY" --recursive
```

This processes nonrotating states **sequentially**, continues after individual failures, and writes `state_0001/`, `state_0002/`, etc. under one new attempt directory. Each contains `source.txt`; the attempt root contains `batch_status.txt`. The Slurm scripts below submit ground-state jobs only.

### Outputs and exit status

```text
Ground-state run:
  wavefunction/parameters.txt       grid and Hamiltonian, full precision
  wavefunction/wavefunction.csv     complex field and restart checkpoint
  wavefunction/energy_history.csv   rotating-frame energy convergence
  wavefunction/diagnostics.csv      E_Omega, E, Lz, mu, residuals at each check
  wavefunction/seed.txt             initial background and random seed
  wavefunction/modulations.tsv      Gaussian peak seed, when supplied
  ground_state.txt                 convergence controls and final diagnostics
  images/initial_density.bmp
  images/initial_phase.bmp
  images/ground_state_density.bmp
  images/ground_state_phase.bmp
  images/energy_history.svg
  log/

BdG attempt:
  source.txt, parameters.txt        input directory and copied Hamiltonian manifest
  bdg_status.txt                   controls, input quality, and solver diagnostics
  spectrum.csv                     returned candidate eigenpairs and residuals
  validation.csv                   per-mode direct-residual pass/fail flags
  excited_states/mode_*.csv         full u,v fields, unless disabled
  images/bdg_spectrum.svg           available when finite eigenpairs were returned
  images/bdg_modes/mode_*.bmp        finite density-fluctuation magnitudes
  log/
```

Both programs return 0 and write `COMPLETED` only on success. Exit 2 means a saved but unconverged/diagnostic calculation (`NOT_CONVERGED` or `INCOMPLETE`). Exit 1 means a run/input/output error, recorded in `FAILED` when a destination exists; CLI syntax errors may use other nonzero codes. `RUNNING` without a terminal marker can indicate a killed job. A BdG convergence failure retains any returned candidate modes, even zero-mode CSV headers, but **these are not a validated spectrum**. An exception before the eigensolver returns can leave only failure diagnostics, not mode files. For recursive processing any unsuccessful state makes the overall exit status nonzero.

## Default 80-point sweep

[parameters/sweep_80.tsv](parameters/sweep_80.tsv) contains 16 evenly spaced scattering lengths from `0.001` to `0.009`, crossed with the 5 in-plane aspect ratios `omega_y/omega_x = 1, 2, 3, 4, 5`, for exactly 80 ground-state attempts. It fixes `a_dd=0.008843`, `N=80000`, `omega_x=1`, and `omega_z=8`. The Slurm script supplies a `1024x1024` grid on `[-25,25) x [-25,25)`. Every default row has `Omega=0` and a noisy Gaussian seed. The previous [40-point table](parameters/sweep_40.tsv) remains available for explicit submission. Neither table is the proposal's `omega_x/omega_y` scan with several rotation frequencies and seed families.

The implemented dimensional reduction integrates out `z`. The existing `omega_z=8` exceeds the new maximum `omega_y=5`, but the ratio is only `8/5` at that end of the sweep; verify that frozen-axis quasi-2D assumptions hold for the intended physical setup. The saved frequencies allow either ratio to be reconstructed; `ground_state.txt` explicitly records `lambda_omega_x_over_omega_y`.

The fixed box still requires a boundary-density check; at fixed 1024 points, enlarging it reduces spatial resolution. Old `test_100.tsv` files are not regenerated automatically: recreate any one-row test file from the current sweep table.

One initial state is used per table row. Near first-order transitions, one converged state does not prove the global minimum. For the proposal's branch search, repeat each physical parameter set with several Gaussian peak patterns and vortex seeds, then compare converged `E_Omega` values. The Gaussian peak method does not create vortices by itself.

## Build and run on Slurm

Dependencies are CMake 3.24+, a C++20 compiler, Eigen 3.3+, and FFTW3. Spectra 1.2, spdlog, and fmt headers are vendored for network-independent compute-node builds. The build script loads the Helios module defaults currently encoded in the script; override `GPES_MODULES` if those module names differ.

```bash
bash scripts/submit_sweep.sh
```

The wrapper submits one ground-state build/test job and then a dependent array of 80 ground-state tasks. Up to 8 tasks run concurrently by default; each requests one CPU and 32 GiB. The submission command returns after queuing jobs. The Slurm build does not build or test BdG, and no BdG job is submitted. Results and logs are stored below:

```text
$SCRATCH/GPES/sweep_YYYYMMDD_HHMMSS_XXXXXX/
```

Every default submission creates a new timestamped sweep directory containing the 80 parameter-set directories, `slurm/` logs, and `parameter_table.tsv`, a copy of the submitted inputs. It also records grid and seed settings and copies requested modulation files to `seeds/` before queuing, so later edits to the original files do not change queued jobs. The random suffix prevents collisions between submissions in the same second. The wrapper prints the directory and Slurm job IDs.

The table accepts one row per run with the following whitespace-separated columns:

```text
a_s a_dd N dt legacy_duration omega_x omega_y omega_z run_name [Omega] [modulation_file]
```

`legacy_duration` is retained for compatibility with the sweep tables, but the ground-state program stops by convergence or `MAX_GROUND_STEPS`, not by this duration. With no `Omega`, rotation is zero. A modulation path may be absolute or relative to the table's directory; it must follow an explicit `run_name` and `Omega`. For example:

```text
0.004 0.008843 80000 0.0005 0.1 1 1.25 5 lambda_0p8_omega_0p1 0.1 three_peak_seed.tsv
```

This row has `lambda=omega_x/omega_y=0.8`, one three-peak seed, and one Slurm task. To scan several seeds at the same physical parameters, give each row a distinct `run_name` and seed file. Such a sweep is not generated automatically.

Useful submission overrides include:

```bash
MAX_CONCURRENT=4 bash scripts/submit_sweep.sh
NX=100 NY=100 MAX_GROUND_STEPS=1000 bash scripts/submit_sweep.sh
RUN_ROOT=/alternate/path bash scripts/submit_sweep.sh parameters/sweep_80.tsv
```

To resume an existing sweep, explicitly reuse its directory:

```bash
RUN_ROOT="$SCRATCH/GPES/sweep_20260906_120000_ABC123" bash scripts/submit_sweep.sh
```

With explicit `RUN_ROOT`, the wrapper requires the same parameter table and seed files. Completed tasks are skipped. An unfinished task resumes its latest saved checkpoint into a new sibling directory named `_retry_...`; previous outputs are not overwritten. Do not submit the same `RUN_ROOT` while its earlier tasks are still running. Use a new `RUN_ROOT` for changed physics, grid, or seeds.

Ground-state convergence uses the relative-plus-absolute energy criterion described above. Slurm supplies a relative energy tolerance of `1e-6` (the standalone CLI still defaults to `1e-8`) and an absolute floor of `1e-12`. The stationary residual must still be at most `1e-4` for three consecutive checks, spaced 100 steps apart. The Slurm step limit is `200000`; override it with `MAX_GROUND_STEPS`. Set `RELATIVE_ENERGY_TOL` or `--relative-energy-tol` to change the energy threshold. Relaxing energy alone does not bypass the stationary residual requirement.

## Outputs

Each Slurm ground-state run directory contains:

```text
wavefunction/wavefunction.csv       final state and restart checkpoint
wavefunction/parameters.txt         checked restart parameter manifest
wavefunction/energy_history.csv     rotating-frame energy history
wavefunction/diagnostics.csv        energy, angular momentum and residual history
wavefunction/seed.txt               initial Gaussian width and noise seed
wavefunction/modulations.tsv        peak seed, when supplied
images/ground_state_density.bmp     dependency-free density image
images/ground_state_phase.bmp       cyclic phase image
images/energy_history.svg           convergence image
ground_state.txt                    convergence and final diagnostics
COMPLETED / NOT_CONVERGED / FAILED  run status
```

Embedded Python/matplotlib plotting is optional. The dependency-free BMP and SVG products are always generated, including in the default headless cluster build.

At `1024x1024`, the complex wavefunction CSV is much larger than the in-memory array. Check scratch quotas before launching many runs; this Slurm workflow saves no BdG mode fields.

## Local verification

The library can also be built and installed independently of the CLI:

```bash
cmake -S . -B build/library \
  -DCMAKE_BUILD_TYPE=Release \
  -DGPES_BUILD_APP=OFF -DGPES_BUILD_TESTS=OFF \
  -DGPES_ENABLE_PLOTTING=OFF -DGPES_ALLOW_FETCHCONTENT=OFF
cmake --build build/library --parallel
cmake --install build/library --prefix "$PWD/build/install"
```

A separate CMake project then consumes the installed library:

```cmake
cmake_minimum_required(VERSION 3.24)
project(my_simulation LANGUAGES CXX)
find_package(gpes 1 CONFIG REQUIRED)
add_executable(my_simulation main.cpp)
target_link_libraries(my_simulation PRIVATE gpes::gpes)
```

Configure that consumer with `-DCMAKE_PREFIX_PATH=/absolute/path/to/install`. Include `gpes.hpp` in its C++ sources; the target supplies the include directories, C++20 requirement, and transitive dependencies. The installed package carries the Spectra headers selected at build time and, for the default vendored build, spdlog/fmt headers. Eigen >= 3.3 and FFTW remain external dependencies found on the consumer machine, including through cluster modules or explicit `EIGEN3_INCLUDE_DIR`, `FFTW3_INCLUDE_DIR`, and `FFTW3_LIBRARY` settings. Builds using external spdlog require its CMake package or the matching dependency target in the consumer.

The core package does not require Python. If GPES was built with plotting enabled, request `find_package(gpes 1 CONFIG REQUIRED COMPONENTS plot)` and link `gpes::plot`; that component also requires Python development files and NumPy. Runtime plotting requires matplotlib. The dependency-free BMP/SVG helpers remain usable with `gpes::gpes`.

`tests/package_consumer` is a standalone consumer fixture: configure it against the install prefix, build it, and run CTest there to exercise the installed split-step and BdG APIs. This tests package integration on a tiny grid, not scientific convergence. Installations are relocatable, provided their external dependencies and compiler ABI remain compatible.

The existing executable and numerical regression test are built as before:

```bash
cmake -S . -B build -DGPES_ENABLE_PLOTTING=OFF -DGPES_BUILD_TESTS=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The numerical regression test checks the high-`k` dipolar limit and compares the complete BdG action against a central finite-difference derivative of the nonlinear eGPE, including contact, LHY, and nonlocal dipolar terms.
