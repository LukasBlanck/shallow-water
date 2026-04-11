# leo3e

## System
- **Kernel:** `Linux 3.10.0-1160.105.1.el7.x86_64`
- **OS:** `CentOS Linux 7 (Core)`
- **Architecture:** `x86_64`

## CPU
- **Model:** `Intel(R) Xeon(R) CPU E5-2650 v3 @ 2.30GHz`
- **Sockets:** `2`
- **Cores per socket:** `10`
- **Threads per core:** `1`
- **Total CPUs:** `20`
- **NUMA nodes:** `2`
- **Instruction support:** includes `AVX`, `AVX2`, `FMA`

## Memory
- **RAM:** `188G`
- **Available:** `140G`
- **Swap:** `191G`

## Current module state
- **Loaded module:** `python/3.11.0-numpy-conda-2023.03`

### Compilers
- `gcc/10.3.0`
- `gcc/8.2.0`
- `gcc/6.5.0`
- `intel/15.0(default)`
- `intel/18.0u1`

### MPI
- `openmpi/4.0.3`
- `openmpi/3.1.4`
- `openmpi/3.1.1`

### Build tools
- `cmake/3.15.3`

### Scientific libraries
- `netcdf-4/4.6.1`
- `netcdf-4/4.3.3.1`
- `netcdf/4.3.3.1`
- `netcdf-4-fortran/4.4.4`
- `hdf5/1.8.20`
- `hdf5/1.8.15`
- `boost/1.58.0`
- `eigen/3.3.0`
- `fftw/3.3.8`
- `gsl/1.16`
- `curl/7.66.0`
- `zlib/1.2.8`
- `intel-mkl/2018u1`

## Current toolchain in PATH
- `gcc`: `/usr/bin/gcc` → version `4.8.5`
- `g++`: `/usr/bin/g++`
- `gfortran`: `/usr/bin/gfortran` → version `4.8.5`
- `cmake`: `/usr/bin/cmake` → version `2.8.12.2`
- `make`: `/usr/bin/make`
- `python3`: Anaconda 2023.03
- `mpicc/mpicxx/mpifort`: **not currently available in PATH**

## Important observation
The current default environment is **old** for modern C++ development:
- GCC `4.8.5`
- CMake `2.8.12.2`
- no MPI compiler wrappers loaded

For a C++ + netCDF project, a better starting point is likely:
```bash
module load gcc/10.3.0
module load cmake/3.15.3
module load openmpi/4.0.3
module load netcdf-4/4.6.1
module load hdf5/1.8.20
```



----------------------------------------------
----------------------------------------------


----------------------------------------------
----------------------------------------------

# GPU cluster environment summary

## System
- **Host:** `mp-gpu3-login`
- **Kernel:** `Linux 3.10.0-957.1.3.el7.x86_64`
- **OS:** `CentOS Linux 7 (Core)`
- **Architecture:** `x86_64`

##  GPU execution
A CUDA test job ran successfully on a compute node and detected **4 × NVIDIA Tesla V100-SXM2-16GB** GPUs.

## Important distinction
- The login node `mp-gpu3-login` has CPU `Intel Xeon W-2102`
- The compute node used for the CUDA job had **4 V100 GPUs**
- Therefore, `lscpu` on the login node does not describe the GPU hardware used in the job

## CUDA environment
- CUDA was initially loaded as `cuda/11.2`
- During the job, modules were reloaded and CUDA changed to `cuda/12.4`
- The successful compilation/execution therefore likely used **CUDA 12.4**

## GPU properties observed
- **GPU model:** `Tesla V100-SXM2-16GB`
- **GPU count:** `4`
- **Max threads per block:** `1024`
- **Unified Addressing:** supported

## CPU
- **Model:** `Intel(R) Xeon(R) W-2102 CPU @ 2.90GHz`
- **Sockets:** `1`
- **Cores per socket:** `4`
- **Threads per core:** `1`
- **Total CPUs:** `4`
- **NUMA nodes:** `1`
- **Instruction support:** includes `AVX`, `AVX2`, `AVX-512`, `FMA`

## Memory
- **RAM:** `30G`
- **Available:** `23G`
- **Swap:** none

## Current module state
Loaded modules:
- `autotools`
- `prun/1.2`
- `gnu7/7.3.0`
- `openmpi3/3.1.0`
- `ohpc`
- `cuda/11.2`

### Compilers
- `gnu7/7.3.0`
- `gnu8/8.2.0`
- `gnu11/11.3.0`
- `pgi/19.7`
- `pgi/16.5`

### MPI
- `openmpi3/3.1.0`
- `mpich/3.2.1`
- `mvapich2/2.2`

### Build tools
- `cmake/3.12.2`
- `cmake/3.24.3`
- `cmake/3.26.3`

### GPU / CUDA
- `cuda/11.2`
- `cuda/12.4`
- `cuda/10.1`
- `cuda/10.0`
- `cuDNN/7.4.2`
- `nccl/2.4.2`

### Scientific libraries
- `netcdf/4.6.1`
- `netcdf-cxx/4.3.0`
- `pnetcdf/1.9.0`
- `hdf5/1.10.2`
- `phdf5/1.10.2`
- `boost/1.67.0`
- `fftw/3.3.7`
- `gsl/2.4`
- `openblas/0.2.20`
- `scalapack/2.0.2`

## Current toolchain in PATH
- `gcc`: `/opt/ohpc/pub/compiler/gcc/7.3.0/bin/gcc` → version `7.3.0`
- `g++`: `/opt/ohpc/pub/compiler/gcc/7.3.0/bin/g++`
- `gfortran`: `/opt/ohpc/pub/compiler/gcc/7.3.0/bin/gfortran` → version `7.3.0`
- `make`: `/usr/bin/make`
- `mpicc`: `/opt/ohpc/pub/mpi/openmpi3-gnu7/3.1.0/bin/mpicc`
- `mpicxx`: `/opt/ohpc/pub/mpi/openmpi3-gnu7/3.1.0/bin/mpicxx`
- `mpifort`: `/opt/ohpc/pub/mpi/openmpi3-gnu7/3.1.0/bin/mpifort`
- `nvcc`: available → CUDA `11.2.152`
- `cmake`: **not currently in PATH**
- `python3`: **not currently in PATH**


For a modern C++ project, it would still be sensible to additionally load a newer CMake:
```bash
module load cmake/3.26.3
module load netcdf/4.6.1
module load netcdf-cxx/4.3.0
module load hdf5/1.10.2
```