# Shallow Water Equations


![Shallow Water Example](Shallow_water.png)
---

## Overview

This project aims to simulate fluids flow governed by the SWE using C++ for heavy computation and python for visualizations.

---

## Theory

The governing PDEs for this simulation are the so called "Shallow Water Equations". They are a solution of the full Navier-Stokes equations under the approximations of 
- constant water density
- newtonian fluid
- shallow water domain
    - vertical velocities are negligible
    - hydrostatic pressure balance
- vertical depth average

It`s big advantage compared to the full Navier-Stokes equation lies in it's variety of possible applications (although a lot of approximations were used, we can find an abondon amount of cases, where those approximations are valid) and in it's time scale to be solved. It is numerically much faster to solve the SWE instead of full Navier-Stokes. This is useful (and even mandatory) for short term predictions of tsunamis for example.

---

$$
h = \text{water depth},\qquad
u = \text{\(x\)-velocity},\qquad
v = \text{\(y\)-velocity},\qquad
$$

$$
g = \text{gravitational acceleration},\qquad
b = \text{bottom elevation}
$$

---

The PDEs are, without any addition of bottom topography or friction or any kind of source term:

>Mass conservation:

$$
\partial _t h + \partial _x (hu) + \partial _y (hv) = 0
$$

>Momenntum conservation:

$$
\partial _t (hu) + \partial _x (hu^2 + \frac{1}{2}gh^2) + \partial _y (huv) = 0
$$

$$
\partial _t (hv) + \partial _y (hv^2 + \frac{1}{2}gh^2) + \partial _x (huv) = 0
$$

The closed form PDE can be written as:

$$
\partial _t U + \partial _x F(U) + \partial _y G(U) = 0
$$

with 
$$
 U = \begin{bmatrix} h \\ hu \\ hv \end{bmatrix} \qquad F(U) = \begin{bmatrix} hu \\ hu^2 +  \frac{1}{2}gh^2\\ huv \end{bmatrix}  \qquad G(U) = \begin{bmatrix} hv \\ huv \\ hv^2 +  \frac{1}{2}gh^2 \end{bmatrix} 
$$

For a more general purpose library, an arbitrary bottom potpography is available. In the PDE, this is portraid as a source term in the following form:

$$
\partial _t U + \partial _x F(U) + \partial _y G(U) = S(U,b)
$$

with

$$
S(U,b) = \begin{bmatrix} 0 \\ -gh\partial_x b \\ -gh\partial_y b \end{bmatrix} 
$$

## Numerical Model

Finite-volume methods are known to be well suited to hyperbolic conservation laws because they preserve the conservative form of the governing equations and can robustly capture shocks and discontinuities. In addition, several well-balanced and positivity-preserving finite-volume schemes have been developed for source terms such as bottom topography. Therefore, the numerical model is based on a two-dimensional finite-volume method and implements a [well-balanced hydrostatic reconstruction](docs/bathymetry.md) scheme.



### Finite Volume

The Finite Volume method is based on a structured rectangular grid, containing $N_x \times N_y$ cells. It's central part are the fluxes through the boarders of every cell. These fluxes depend on height and velocity values of it's boardering cells. Since we only have a finite discretization of our physical domain, this boardering values are in general discontinous and the calculation of the flux therefore constitutes a "Riemann Problem".

### Riemann Problem

For the Riemann Problem, there exist numerous approximative ways to solve. They are well known and all have advantages and disatvantageous. In order to compare and define the best suited Riemann solver to optimize the library for, three different Riemann Solvers were implemented:
- Rusanov
- HLL
- ROE

They're respective properties can be found in [SWE_second.pptx](SWE_second.pptx) (only a qualitative analysis, gained from testing them).

### Reconstruction

The Riemann solvers all need values from their boardering cells. First order Finite Volume works with a constant cell value (a cell average). In order to get values that correspond better to the actual physical value at the boundary one can use a reconstruction method to approximate this physical value better. In this library the generic first order (Piecewise constant) and second order (MUSCL) reconstruction is provided.

## Structure

![](docs/Solver_structure.svg)

The library consists of a highly modular, but not optimized 
- SerialSolver
    - good for comparison of different Riemann and Reconstruction methods
    - easily extendable 
    - not very fast

Furthermore four different optimized versions are provided, that all use only the HLL Riemann Integrator and the MUSCL Reconstruction method. Those two (MUSCL+HLL) were chosen due to their good stability and second order accuracy in smooth enough regions.

- `HPC_solver`
    - serial solver optimized
- `omp_solver`
    - OpenMP accelerated
- `cuda_solver`
    - GPU accelerated
- `cuda_solver_4gpu`
    - runs on 4 GPUs

## Requirements

To build and run the project, you need:

- **CMake** at least 3.15
- **C++** compiler supporting C++17
- **Python 3.11** or newer
- **NetCDF-C**
- **NetCDF-C++4**

### Recommended
- **OpenMP:** OpenMP-capable compiler and runtime library
- **CUDA:** Toolkit and a supported NVIDIA GPU for CUDA backends

## Compile and Run

Configure the project:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release 
```

Build it:

```
cmake --build build --parallel
```

Run the simulation:

```
./build/swe_solver
```

---

Alternatively, use the provided helper script:
```
./scripts/run.sh
```



### MacOS with MacPorts

For MacOS the path of netCDF is often not findable for CMake or pkgconfig. In that case the manual prefix has to be added. For example:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/opt/local
```

or 

```
CMAKE_PREFIX_PATH=/opt/local ./scripts/run.sh
```

#### OpenMP

The OpenMP and CUDA backends are conditionally compiled. If either dependency is unavailable or explicitly disabled, the corresponding accelerated backend is excluded while the serial solver remains available.

CMake automatically searches for OpenMP and if it is found the OpenMP accelerated backend gets compiled to enable a toml change only backend choosing.

If no OpenMP is found, the solver compiles only the serial solver.

To require OpenMP and fail configuration if no OpenMP found, run:
```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DENABLE_OPENMP=ON -DREQUIRE_OPENMP=ON
```

## Cluster execution

The submission scripts under `scripts/` are specific to the clusters on which they were developed. Module names, queues, partitions, GPU resource requests, compiler versions, and library paths may need to be adapted for other systems.

Run all submission commands from the repository root because the scripts configure the project with `cmake -S .`

### Serial on LEO3E

The serial script disables all accelerated backends:

```bash
qsub scripts/leo3e/run_serial.sh
```

The selected simulation configuration must use a serial backend.

### OpenMP on LEO3E

```bash
qsub scripts/leo3e/run_openmp.sh
```

The script requests four OpenMP slots and sets `OMP_NUM_THREADS` from the scheduler-provided `NSLOTS` value. It configures the project with OpenMP enabled and required.

The selected simulation configuration can use the OpenMP backend or the serial ones.

### Single-GPU CUDA job

```bash
sbatch scripts/gpu/run_cuda.sh
```

The script requests one CUDA-capable node, disables OpenMP, builds the CUDA backend, and launches one solver process using `srun`.

The selected configuration must use the CUDA backend.

### Four-GPU CUDA job

```bash
sbatch scripts/gpu/run_cuda_4gpu.sh
```

The script requests one CUDA-capable node and starts one process that manages the four devices.

The selected configuration must use the four-GPU CUDA backend.

---

### Numerical Verification

The numerical correctness of the solver implementations can be monitored using optional sanity checks. These checks evaluate important physical and numerical properties during a simulation and can help detect unstable configurations, implementation errors, or invalid states.

The checks are intended primarily for debugging and verification. They can increase the runtime of large simulations and are therfore modularly choosable in the [simulation_config](configs/simulation_config.toml).

#### Mass Conservation

For a closed domain without mass source terms, the total water mass should remain constant:

$$
M(t) = \int_{\Omega} h(x,y,t)\,\mathrm{d}\Omega.
$$

In the finite-volume discretization, the total mass is approximated by

$$
M^n = \sum_{i,j} h_{i,j}^n\,\Delta x\,\Delta y.
$$

The sanity check compares the current mass with the initial mass and reports the relative deviation:

$$
\varepsilon^n
=
\frac{\left|M^n-M^0\right|}
     {M^0}.
$$

For refelecting boundary conditions, mass should be conserved up to floating-point roundoff or up to the numerical tolerance manually set in [simulation_config](configs/simulation_config.toml).

#### Energy Monitoring

For the inviscid shallow water equations without external forcing, friction, or energy flux through the boundary, the total mechanical energy is

$$
E(t)
=
\int_{\Omega}
\left[
\frac{1}{2}h\left(u^2+v^2\right)
+
\frac{1}{2}gh^2
\right]
\,\mathrm{d}\Omega.
$$

Using the conserved momenta $hu$ and $hv$, the discrete energy can be evaluated as

$$
E^n
=
\sum_{i,j}
\left[
\frac{(hu)_{i,j}^2+(hv)_{i,j}^2}
     {2h_{i,j}}
+
\frac{1}{2}g h_{i,j}^2
\right]
\Delta x\,\Delta y.
$$

If bottom topography is included, the gravitational potential-energy contribution associated with the bottom elevation must also be taken into account:

$$
E^n
=
\sum_{i,j}
\left[
\frac{(hu)_{i,j}^2+(hv)_{i,j}^2}
     {2h_{i,j}}
+
\frac{1}{2}g h_{i,j}^2
+
g h_{i,j} b_{i,j}
\right]
\Delta x\,\Delta y.
$$

Unlike mass, total energy is not necessarily conserved exactly by the numerical method. Approximate Riemann solvers and slope limiters introduce numerical dissipation, so a gradual decrease in energy can be expected. A sudden increase, a non-finite value, or a large unexplained deviation may indicate numerical instability.


#### Positivity

A physically valid shallow-water state requires a non-negative water depth:

$$
h_{i,j} \geq 0
\qquad
\text{for all cells }(i,j).
$$

The positivity check monitors the minimum water depth in the computational domain and detects negative states and throws runtime errors if detected.

#### Scaling Analysis

For a detailed scaling analysis and some plots, take a look at [SWE_final](docs/SWE_final.pptx).

---

### Visualization
Create and activate the virtual environment needed to run the visualization files:

```
python3 -m venv venv
source venv/bin/activate
python3 -m pip install -r configs/requirements.txt
```

For a 3D animation, run:

```
python3 python/visualize_3D_animated.py
```

For a 2D animation, run:
```
python3 python/visualize_2D_animated.py
```

The stated visualization scripts use an interpolation between different cell values. For a fine enough grid that is accurate enough. For sparse grids and to inspect the exact call values, the following scripts are provided: `python/visualize_correct_grid2D.py` and `python/visualize_correct_grid3D.py` 

#### Sanity Checks

If the `debug` option is set in the [simulation_config](configs/simulation_config.toml), the monitored sanity checks may be inspected with

```
python3 python/sanity_checks/viz_mass.py
python3 python/sanity_checks/viz_energy.py
```

---

## Remarks

Every possible solver combination conserves mass up to floating point roundoff and depending on the artificial diffusivity of the solver combination also conserves energy in a physically reasonable and expectable range.

Depending on the solver combination and the chosen initial conditions, the solver preserves positivity, but for very difficult initial conditions (steep gradients in very shallow water) the numerical simulation might fail. This behaviour is known to happen and might need further limiters and special treatment.

## Additional Material

To inspect the numerical scheme without inclusion of a bottom topography, take a look at [numerics_without_bathymetry](docs/numerics_without_bathymetry.md).


With bottom topography and therefor a well balanced hydrostatic reconstruction, [bathymetry](docs/bathymetry.md).

The presentation slides used to present [initial outlook](docs/SWE_first.pptx), [progress](docs/SWE_second.pptx) and [final report](docs/SWE_final.pptx) of the solver library are provided in `docs/*.pptx`.



## License
Copyright (c) 2026 Lukas Blanck

MIT License —  [LICENSE](LICENSE).