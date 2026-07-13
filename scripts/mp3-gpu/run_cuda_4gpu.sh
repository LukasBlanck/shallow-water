#!/bin/bash
#SBATCH -p nvltv
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -t 00:20:00
#SBATCH -o /dev/null

set -euo pipefail

CONFIG_FILE="${1:-configs/simulation_config.toml}"
BUILD_TYPE="${BUILD_TYPE:-Release}"
BUILD_DIR="${BUILD_DIR:-build-cuda-4gpu}"
EXECUTABLE="${EXECUTABLE:-${BUILD_DIR}/swe_solver}"
OUTFILE="swe_cuda_4gpu.out"
ARCHIVE="swe_cuda_4gpu.outs_old"

section() {
    echo
    echo "============================================================"
    echo "$1"
    echo "============================================================"
}

archive_old_output() {
    if [[ -f "${OUTFILE}" ]]; then
        {
            echo "===== old ${OUTFILE} from $(date '+%Y-%m-%d %H:%M:%S') ====="
            cat "${OUTFILE}"
            echo
        } >> "${ARCHIVE}"
        rm -f "${OUTFILE}"
    fi
    exec > "${OUTFILE}" 2>&1
}

capture_netcdf_paths() {
    section "Discovering NetCDF paths from gnu7/openmpi3 stack"
    module purge
    module load gnu7/7.3.0
    module load openmpi3/3.1.0
    module load netcdf/4.6.1
    module load netcdf-cxx/4.3.0

    NETCDF_INCLUDE_DIR="${NETCDF_INC}"
    NETCDF_LIBRARY_DIR="${NETCDF_LIB}"
    NETCDF_CXX_INCLUDE_DIR="${NETCDF_CXX_INC}"
    NETCDF_CXX_LIBRARY_DIR="${NETCDF_CXX_LIB}"
    HDF5_INCLUDE_DIR="${HDF5_INC}"
    HDF5_LIBRARY_DIR="${HDF5_LIB}"
    MPI_LIBRARY_DIR="${MPI_DIR}/lib"

    test -f "${NETCDF_INCLUDE_DIR}/netcdf.h"
    test -f "${NETCDF_LIBRARY_DIR}/libnetcdf.so"
    test -f "${NETCDF_CXX_LIBRARY_DIR}/libnetcdf_c++4.so"
}

load_build_toolchain() {
    section "Loading CUDA build toolchain"
    module purge
    module load gnu11/11.3.0
    module load cuda/12.4
    module load cmake/3.26.3

    GCC_LIBDIR="$(dirname "$(g++ -print-file-name=libstdc++.so.6)")"

    export CC="$(which gcc)"
    export CXX="$(which g++)"
    export CUDAHOSTCXX="${CXX}"
    export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

    export CPATH="${NETCDF_CXX_INCLUDE_DIR}:${NETCDF_INCLUDE_DIR}:${HDF5_INCLUDE_DIR}:${CPATH:-}"
    export LIBRARY_PATH="${NETCDF_CXX_LIBRARY_DIR}:${NETCDF_LIBRARY_DIR}:${HDF5_LIBRARY_DIR}:${MPI_LIBRARY_DIR}:${LIBRARY_PATH:-}"
    export LD_LIBRARY_PATH="${GCC_LIBDIR}:${NETCDF_CXX_LIBRARY_DIR}:${NETCDF_LIBRARY_DIR}:${HDF5_LIBRARY_DIR}:${MPI_LIBRARY_DIR}:/opt/sw/cuda-12.4/lib64:${LD_LIBRARY_PATH:-}"
}

print_diagnostics() {
    section "CUDA SWE 4-GPU job"
    echo "JOB-ID:      ${SLURM_JOB_ID:-unknown}"
    echo "Node:        $(hostname)"
    echo "Workdir:     $(pwd)"
    echo "Config:      ${CONFIG_FILE}"
    echo "Build dir:   ${BUILD_DIR}"
    echo "Executable:  ${EXECUTABLE}"
    echo "Start time:  $(date)"

    section "Toolchain"
    module list || true
    echo "gcc:   $(which gcc)"; gcc --version | head -n 1
    echo "g++:   $(which g++)"; g++ --version | head -n 1
    echo "nvcc:  $(which nvcc)"; nvcc --version
    echo "cmake: $(which cmake)"; cmake --version | head -n 1

    section "GPU diagnostics"
    nvidia-smi || true
    nvidia-smi -L || true
    echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-not set}"
}

configure_build() {
    section "Configuring CUDA 4-GPU build"
    rm -rf "${BUILD_DIR}"
    cmake -S . -B "${BUILD_DIR}" \
        -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
        -DCMAKE_C_COMPILER="${CC}" \
        -DCMAKE_CXX_COMPILER="${CXX}" \
        -DCMAKE_CUDA_HOST_COMPILER="${CUDAHOSTCXX}" \
        -DENABLE_CUDA_4=ON \
        -DREQUIRE_CUDA=ON \
        -DENABLE_OPENMP=OFF \
        -DREQUIRE_OPENMP=OFF \
        -DCMAKE_CUDA_ARCHITECTURES=70 \
        -DNETCDF_INCLUDE_DIR="${NETCDF_INCLUDE_DIR}" \
        -DNETCDF_LIBRARY_DIR="${NETCDF_LIBRARY_DIR}" \
        -DNETCDF_CXX_INCLUDE_DIR="${NETCDF_CXX_INCLUDE_DIR}" \
        -DNETCDF_CXX_LIBRARY_DIR="${NETCDF_CXX_LIBRARY_DIR}"
}

build_and_run() {
    section "Building"
    cmake --build "${BUILD_DIR}" -j "${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-4}}"

    section "Running CUDA 4-GPU SWE solver"
    [[ -x "${EXECUTABLE}" ]] || { echo "ERROR: executable not found: ${EXECUTABLE}"; find "${BUILD_DIR}" -maxdepth 4 -type f -executable || true; exit 1; }
    srun -n 1 -c "${SLURM_CPUS_PER_TASK:-1}" "${EXECUTABLE}" "${CONFIG_FILE}"

    section "CUDA 4-GPU SWE solver finished successfully"
    echo "End time: $(date)"
}

archive_old_output
[[ -f "${CONFIG_FILE}" ]] || { echo "ERROR: config file not found: ${CONFIG_FILE}"; exit 1; }
capture_netcdf_paths
load_build_toolchain
print_diagnostics
configure_build
build_and_run
