#!/bin/bash
#$ -N swe_serial
#$ -q short.q,std.q
#$ -cwd
#$ -j yes
#$ -l h_rt=0:20:0

set -euo pipefail

CONFIG_FILE="${CONFIG_FILE:-configs/simulation_config.toml}"
BUILD_DIR="${BUILD_DIR:-build-serial}"
EXECUTABLE="${BUILD_DIR}/swe_solver"
BUILD_TYPE="${BUILD_TYPE:-Release}"

echo "============================================================"
echo "Serial SWE job"
echo "Started:           $(date)"
echo "Host:              $(hostname)"
echo "Working directory: $(pwd)"
echo "Configuration:     ${CONFIG_FILE}"
echo "Build directory:   ${BUILD_DIR}"
echo "Build type:        ${BUILD_TYPE}"
echo "============================================================"

# ============================================================
# Modules
# ============================================================

echo "Loading modules..."

module purge
module load gcc/10.3.0
module load cmake/3.15.3
module load netcdf-4/4.3.3.1

# ============================================================
# NetCDF C from cluster module
# ============================================================

if [[ -n "${NETCDF:-}" ]]; then
    export NETCDF_ROOT="${NETCDF}"
fi

if [[ -n "${UIBK_NETCDF_4_INC:-}" ]]; then
    export NETCDF_INCLUDE_DIR="${UIBK_NETCDF_4_INC}"
    export CPATH="${UIBK_NETCDF_4_INC}:${CPATH:-}"
fi

if [[ -n "${UIBK_NETCDF_4_LIB:-}" ]]; then
    export NETCDF_LIBRARY_DIR="${UIBK_NETCDF_4_LIB}"
    export LIBRARY_PATH="${UIBK_NETCDF_4_LIB}:${LIBRARY_PATH:-}"
    export LD_LIBRARY_PATH="${UIBK_NETCDF_4_LIB}:${LD_LIBRARY_PATH:-}"
fi

# ============================================================
# HDF5 dependency
# ============================================================

if [[ -n "${UIBK_HDF5_INC:-}" ]]; then
    export CPATH="${UIBK_HDF5_INC}:${CPATH:-}"
fi

if [[ -n "${UIBK_HDF5_LIB:-}" ]]; then
    export LIBRARY_PATH="${UIBK_HDF5_LIB}:${LIBRARY_PATH:-}"
    export LD_LIBRARY_PATH="${UIBK_HDF5_LIB}:${LD_LIBRARY_PATH:-}"
fi

if [[ -n "${UIBK_HDF5_BIN:-}" ]]; then
    export PATH="${UIBK_HDF5_BIN}:${PATH}"
fi

# ============================================================
# NetCDF C++4 from user installation
# ============================================================

export NETCDF_CXX_ROOT="${NETCDF_CXX_ROOT:-${HOME}/local/netcdf-cxx4}"
export NETCDF_CXX_INCLUDE_DIR="${NETCDF_CXX_ROOT}/include"
export NETCDF_CXX_LIBRARY_DIR="${NETCDF_CXX_ROOT}/lib"

export CPATH="${NETCDF_CXX_INCLUDE_DIR}:${CPATH:-}"
export LIBRARY_PATH="${NETCDF_CXX_LIBRARY_DIR}:${LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="${NETCDF_CXX_LIBRARY_DIR}:${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="${NETCDF_CXX_LIBRARY_DIR}/pkgconfig:${PKG_CONFIG_PATH:-}"

# ============================================================
# Diagnostics
# ============================================================

echo "============================================================"
echo "Loaded modules"
echo "============================================================"

module list || true

echo "============================================================"
echo "Toolchain"
echo "============================================================"

echo "Compiler:"
command -v g++
g++ --version

echo
echo "CMake:"
command -v cmake
cmake --version

echo
echo "NetCDF paths:"
echo "NETCDF_ROOT=${NETCDF_ROOT:-not set}"
echo "NETCDF_INCLUDE_DIR=${NETCDF_INCLUDE_DIR:-not set}"
echo "NETCDF_LIBRARY_DIR=${NETCDF_LIBRARY_DIR:-not set}"
echo "NETCDF_CXX_ROOT=${NETCDF_CXX_ROOT}"
echo "NETCDF_CXX_INCLUDE_DIR=${NETCDF_CXX_INCLUDE_DIR}"
echo "NETCDF_CXX_LIBRARY_DIR=${NETCDF_CXX_LIBRARY_DIR}"

# ============================================================
# Build
# ============================================================

echo "============================================================"
echo "Configuring serial build"
echo "============================================================"

rm -rf "${BUILD_DIR}"

cmake -S . -B "${BUILD_DIR}" \
    -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
    -DENABLE_OPENMP=OFF \
    -DENABLE_CUDA=OFF \
    -DENABLE_CUDA_4=OFF \
    -DNETCDF_ROOT="${NETCDF_ROOT:-}" \
    -DNETCDF_INCLUDE_DIR="${NETCDF_INCLUDE_DIR:-}" \
    -DNETCDF_LIBRARY_DIR="${NETCDF_LIBRARY_DIR:-}" \
    -DNETCDF_CXX_ROOT="${NETCDF_CXX_ROOT}" \
    -DNETCDF_CXX_INCLUDE_DIR="${NETCDF_CXX_INCLUDE_DIR}" \
    -DNETCDF_CXX_LIBRARY_DIR="${NETCDF_CXX_LIBRARY_DIR}"

echo "============================================================"
echo "Compiling serial executable"
echo "============================================================"

cmake --build "${BUILD_DIR}"

if [[ ! -x "${EXECUTABLE}" ]]; then
    echo "ERROR: executable '${EXECUTABLE}' was not created."
    echo "Executable files found under '${BUILD_DIR}':"
    find "${BUILD_DIR}" -maxdepth 3 -type f -executable || true
    exit 1
fi

# ============================================================
# Run
# ============================================================

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "ERROR: configuration file '${CONFIG_FILE}' does not exist."
    exit 1
fi

echo "============================================================"
echo "Running serial simulation"
echo "============================================================"

"${EXECUTABLE}" "${CONFIG_FILE}"

echo "============================================================"
echo "Job finished: $(date)"
echo "============================================================"