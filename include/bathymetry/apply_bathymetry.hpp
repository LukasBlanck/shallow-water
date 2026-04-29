#pragma once

#include <stdexcept>

#include "configs/config.hpp"
#include "include/bathymetry/flat_bathymetry.hpp"
#include "include/bathymetry/gaussian_hill.hpp"
#include "include/bathymetry/no_bathymetry.hpp"
#include "include/core/array2D.hpp"
#include "include/core/grid.hpp"

inline void apply_bathymetry(const SimulationConfig &cfg, const Grid &grid, Array2D &B) {
    switch (cfg.bathymetry.type) {
    case BathymetryType::Flat: {
        Flat bathy(cfg.bathymetry.b0);
        bathy.apply(grid, B);
        return;
    }

    case BathymetryType::None: {
        None bathy;
        bathy.apply(grid, B);
        return;
    }

    case BathymetryType::GaussHill: {
        // Temporary fallback until GaussHill exists.
        GaussHill bathy(cfg.bathymetry.bathy_peak_height, cfg.bathymetry.bathy_sigma_x,
                        cfg.bathymetry.bathy_sigma_y, cfg.bathymetry.bathy_x0,
                        cfg.bathymetry.bathy_y0, cfg.bathymetry.b0);
        bathy.apply(grid, B);
        return;
    }

    default:
        throw std::runtime_error("Unsupported Bathymetry type");
    }
}