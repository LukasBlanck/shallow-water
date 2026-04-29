#pragma once

#include "include/core/array2D.hpp"
#include "include/core/grid.hpp"

class Bathymetry {
  public:
    virtual ~Bathymetry() = default;
    virtual void apply(const Grid &grid, Array2D &U) const = 0;
};