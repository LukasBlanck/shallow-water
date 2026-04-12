// U_L = U_bar_i,j
// U_R = U_bar_i+1,j
// U_T = U_bar_i,j+1
// U_B = U_bar_i,j

#pragma once

#include "include/backend/serial/reconstruction/reconstruction.hpp"
#include "include/core/cell_state.hpp"

class PiecewiseConst : public Reconstruction {
    // These functions return the values of the cells (i,j) neighbouring the face
  public:
    PiecewiseConst() = default;

    // x neighbours   (left and right)
    CellState U_L(const State &U, int i, int j) const override {
        return CellState(U.h()(i, j), U.hu()(i, j), U.hv()(i, j));
    }
    CellState U_R(const State &U, int i, int j) const override {
        return CellState(U.h()(i + 1, j), U.hu()(i + 1, j), U.hv()(i + 1, j));
    }

    // y neighbours   (bottom and top)
    CellState U_B(const State &U, int i, int j) const override {
        return CellState(U.h()(i, j), U.hu()(i, j), U.hv()(i, j));
    }

    CellState U_T(const State &U, int i, int j) const override {
        return CellState(U.h()(i, j + 1), U.hu()(i, j + 1), U.hv()(i, j + 1));
    }
};