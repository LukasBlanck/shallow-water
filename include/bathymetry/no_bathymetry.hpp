#pragma once

#include "include/bathymetry/bathymetry.hpp"
#include "include/core/array2D.hpp"

class None : public Bathymetry {
  public:
    None() {};
    static constexpr bool enabled = false;

    void apply(const Grid &grid, Array2D &B) const override {
        const int nG = grid.nG();
        for (int i = 0; i < grid.Nx_total(); i++) {
            for (int j = 0; j < grid.Ny_total(); j++) {
                B(i, j) = 0.0;
            }
        }
    }
};