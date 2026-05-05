#!/bin/bash
#$ -N swe_serial
#$ -q short.q,std.q
#$ -cwd
#$ -j yes
#$ -l h_rt=0:60:0
#$ -pe openmp 1

set -euo pipefail

CONFIG_FILE="${1:-simulation_config.toml}"
EXECUTABLE="./build/my_project"

echo "============================================================"
echo "Job started on: $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Config file: ${CONFIG_FILE}"
echo "NSLOTS: ${NSLOTS:-1}"
echo "============================================================"

echo "Loading modules..."
module purge
module load gcc/10.3.0
module load cmake/3.15.3

# Important:
# This is the NetCDF C module against which your local netcdf-cxx4 was built.
# Do not switch this to netcdf-4/4.6.1.
module load netcdf-4/4.3.3.1

# ============================================================
# Normalize NetCDF C from module
# ============================================================

if [[ -n "${NETCDF:-}" ]]; then
    export NETCDF_ROOT="${NETCDF}"
    export NETCDF_DIR="${NETCDF}"
    export NetCDF_ROOT="${NETCDF}"
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
# Normalize HDF5 from module dependency
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
# NetCDF C++ from local user installation
# ============================================================

export NETCDF_CXX_ROOT="${HOME}/local/netcdf-cxx4"
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
which g++
g++ --version

echo
echo "CMake:"
which cmake
cmake --version

echo
echo "Pkg-config:"
which pkg-config || true
pkg-config --version || true

echo "============================================================"
echo "NetCDF diagnostics"
echo "============================================================"

echo "NETCDF_ROOT=${NETCDF_ROOT:-not set}"
echo "NETCDF_INCLUDE_DIR=${NETCDF_INCLUDE_DIR:-not set}"
echo "NETCDF_LIBRARY_DIR=${NETCDF_LIBRARY_DIR:-not set}"
echo "NETCDF_CXX_ROOT=${NETCDF_CXX_ROOT:-not set}"
echo "NETCDF_CXX_INCLUDE_DIR=${NETCDF_CXX_INCLUDE_DIR:-not set}"
echo "NETCDF_CXX_LIBRARY_DIR=${NETCDF_CXX_LIBRARY_DIR:-not set}"

echo
echo "CPATH=${CPATH:-not set}"
echo "LIBRARY_PATH=${LIBRARY_PATH:-not set}"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-not set}"
echo "PKG_CONFIG_PATH=${PKG_CONFIG_PATH:-not set}"

echo
echo "NetCDF C headers:"
find "${NETCDF_INCLUDE_DIR:-/does/not/exist}" -maxdepth 1 -type f \
    | grep -Ei 'netcdf' || true

echo
echo "NetCDF C libraries:"
find "${NETCDF_LIBRARY_DIR:-/does/not/exist}" -maxdepth 1 \( -type f -o -type l \) \
    | grep -Ei 'libnetcdf' || true

echo
echo "NetCDF C++ headers:"
find "${NETCDF_CXX_INCLUDE_DIR:-/does/not/exist}" -maxdepth 1 -type f \
    | grep -Ei 'ncfile|ncgroup|ncvar|netcdf' || true

echo
echo "NetCDF C++ libraries:"
find "${NETCDF_CXX_LIBRARY_DIR:-/does/not/exist}" -maxdepth 1 \( -type f -o -type l \) \
    | grep -Ei 'libnetcdf.*c\+\+|libnetcdf_c' || true

echo
echo "ncxx4-config:"
which ncxx4-config || true
ncxx4-config --all 2>/dev/null || true

echo "============================================================"
echo "Building project inside the job"
echo "============================================================"

rm -rf build

cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DNETCDF_ROOT="${NETCDF_ROOT:-}" \
    -DNETCDF_INCLUDE_DIR="${NETCDF_INCLUDE_DIR:-}" \
    -DNETCDF_LIBRARY_DIR="${NETCDF_LIBRARY_DIR:-}" \
    -DNETCDF_CXX_ROOT="${NETCDF_CXX_ROOT:-}" \
    -DNETCDF_CXX_INCLUDE_DIR="${NETCDF_CXX_INCLUDE_DIR:-}" \
    -DNETCDF_CXX_LIBRARY_DIR="${NETCDF_CXX_LIBRARY_DIR:-}"

cmake --build build -j "${NSLOTS:-1}"

if [[ ! -x "${EXECUTABLE}" ]]; then
    echo "ERROR: executable ${EXECUTABLE} not found or not executable."
    echo "Available executables in build/:"
    find build -maxdepth 3 -type f -executable || true
    exit 1
fi

"${EXECUTABLE}" "${CONFIG_FILE}"

echo "============================================================"
echo "Job finished on: $(date)"
echo "============================================================"