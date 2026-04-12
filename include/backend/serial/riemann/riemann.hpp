#pragma once

#include "include/core/cell_state.hpp"
class RiemannSolver {
  public:
    virtual ~RiemannSolver() = default;

    virtual CellState x_flux(const CellState &U_L, const CellState &U_R) const = 0;
    virtual CellState y_flux(const CellState &U_B, const CellState &U_T) const = 0;
};