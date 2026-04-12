#pragma once

// this class builds the cell values of the neighbouring states

#include "include/core/cell_state.hpp"
#include "include/core/state.hpp"
class Reconstruction {
  public:
    Reconstruction();

    // for x fluxes
    CellState virtual U_L(State &U, int i, int j) const = 0;
    CellState virtual U_R(State &U, int i, int j) const = 0;

    // for y fluxes
    CellState virtual U_B(State &U, int i, int j) const = 0;
    CellState virtual U_T(State &U, int i, int j) const = 0;
};