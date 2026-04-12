// U_L = U_bar_i,j
// U_R = U_bar_i+1,j
// U_T = U_bar_i,j+1
// U_B = U_bar_i,j

#include "include/backend/serial/reconstruction/reconstruction.hpp"
#include "include/core/cell_state.hpp"

class PiecewiseConst : public Reconstruction {
    // These functions return the values of the cells (i,j) neighbouring the face
  public:
    // x neighbours   (left and right)
    CellState U_L(State &U, int i, int j) const override {
        return CellState(U.h()(i, j), U.hu()(i, j), U.hv()(i, j));
    }
    CellState U_R(State &U, int i, int j) const override {
        return CellState(U.h()(i + 1, j), U.hu()(i + 1, j), U.hv()(i + 1, j));
    }

    // y neighbours   (bottom and top)
    CellState U_B(State &U, int i, int j) const override {
        return CellState(U.h()(i, j), U.hu()(i, j), U.hv()(i, j));
    }

    CellState U_T(State &U, int i, int j) const override {
        return CellState(U.h()(i, j + 1), U.hu()(i, j + 1), U.hv()(i, j + 1));
    }
};