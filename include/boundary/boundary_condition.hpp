#pragma once

#include "include/core/array2D.hpp"
#include "include/core/grid.hpp"
#include "include/core/state.hpp"

class BoundaryCondition {
  public:
    BoundaryCondition(Grid &grid)
        : Nx_total_(grid.Nx_total()), Ny_total_(grid.Ny_total()), nG_(grid.nG()) {}

    virtual ~BoundaryCondition() = default;

    virtual void apply_BC(State &U) { (void)U; }

    virtual void apply_BC(Array2D &B) { (void)B; }

  protected:
    int Nx_total_;
    int Ny_total_;
    int nG_;
};