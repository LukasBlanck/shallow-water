#pragma once
#include "configs/config.hpp"
#include <string>

inline std::string backend_name_from_cfg(const SimulationConfig &cfg) {
    switch (cfg.backend.type) {
    case BackendType::Serial:
        return "Serial";
    case BackendType::OptimizedSerial:
        return "OptimizedSerial";
    case BackendType::OpenMP:
        return "OpenMP";
    }

    return "UnknownBackend";
}

inline std::string bathymetry_name_from_cfg(const SimulationConfig &cfg) {
    switch (cfg.bathymetry.type) {
    case BathymetryType::None:
        return "None";
    case BathymetryType::Flat:
        return "Flat";
    case BathymetryType::GaussHill:
        return "Gaussian Hill";
    }

    return "UnknownBathymetry";
}