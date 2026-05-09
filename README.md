# Shallow Water Equations


> **Status:** This project is currently under active development.

---

## Overview

This project aims to simulate fluids flow governed by the SWE using C++ for heavy computation and python3 for visualization.

---

## Features

Current / planned project scope:

- [x] CMake-based build for local setup 
- [x] Modular parameters/solvers via .toml input
- [x] Python visualization
- [x] SanityChecks::Positivity
- [x] SanityChecks::Mass_conservation
- [x] Debug output
---
- [x] Reconstruction::PiecewiseConst 
- [x] Reconstruction::MUSCL
- [x] Riemann::Rusanov
- [x] Riemann::HLL
- [x] Riemann::ROE

---
- [x] Non-flat bathymetry
- [x] CMake-based build for laptop setup
- [ ] CMake-based build for cluster setup 
- [ ] Highly optimized HLL-MUSCL solver:
    - [x] Serial
    - [x] OpenMP (Cluster)
    - [ ] CUDA (Cluster)

---
- [ ] Timing capabilities
- [ ] Weak scaling analysis
- [ ] Strong scaling analysis

## Requirements

To build and run the project, you currently need:

- **CMake**
- **Ninja**
- **C/C++ compiler**
- **Python 3**
- NetCDF library

#### Laptop

Compile and run on laptop:
```
./run.sh
```


#### LEO3E
Compile and run on leo3e:

```
./run_leo3e.sh
tail -f swe_serial.o<JOB-ID>
```

---
Visualization:
```
python3 -m venv venv
venv source/bin/activate
python3 pip install configs/requirements.txt
```
```
python3 python/visualize_3D_animated.py
```

## Solver Structure

![](docs/Solver_structure.svg)

## License
Copyright (c) 2026 Lukas Blanck

MIT License —  [LICENSE](LICENSE).