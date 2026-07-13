#!/usr/bin/env bash
set -euo pipefail

BUILD_DIR="${BUILD_DIR:-build}"
BUILD_TYPE="${BUILD_TYPE:-Release}"
CONFIG_FILE="${1:-simulation_config.toml}"

echo "Configuring ${BUILD_TYPE} build..."
cmake -S . -B "${BUILD_DIR}" \
    -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"

echo "Compiling..."
cmake --build "${BUILD_DIR}" --parallel

echo "Running simulation..."
"${BUILD_DIR}/swe_solver" "${CONFIG_FILE}"

echo "Done."