#pragma once

// this class builds the cell values of the neighbouring states

#include "include/core/cell_state.hpp"
#include "include/core/state.hpp"
class Reconstruction {
  public:
    Reconstruction();

    // for x fluxes
    CellState virtual U_L(const State &U, int i, int j) const = 0;
    CellState virtual U_R(const State &U, int i, int j) const = 0;

    // for y fluxes
    CellState virtual U_B(const State &U, int i, int j) const = 0;
    CellState virtual U_T(const State &U, int i, int j) const = 0;
};