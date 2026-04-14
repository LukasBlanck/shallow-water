#pragma once

#include "include/core/cell_state.hpp"
class RiemannSolver {
  public:
    virtual ~RiemannSolver() = default;

    virtual CellState x_flux(const CellState &U_L, const CellState &U_R) const = 0;
    virtual CellState y_flux(const CellState &U_B, const CellState &U_T) const = 0;

  protected:
    CellState F(const CellState &U) const {
        return CellState(U.hu(), U.u() * U.u() * U.h() + 0.5 * constants::g * U.h() * U.h(),
                         U.hu() * U.v());
    }

    CellState G(const CellState &U) const {
        return CellState(U.hv(), U.hu() * U.v(),
                         U.h() * U.v() * U.v() + 0.5 * constants::g * U.h() * U.h());
    }
};