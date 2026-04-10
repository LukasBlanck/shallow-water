#!/usr/bin/env bash
set -e  # bash stops if a fail occurs

echo "Configuring..."
cmake -S . -B build -G Ninja -DCMAKE_PREFIX_PATH=/opt/local

echo "Compiling..."
cmake --build build

echo "Running project..."
./build/my_project

echo "Running visualization..."
python3 python/visualize.py

echo "Done."