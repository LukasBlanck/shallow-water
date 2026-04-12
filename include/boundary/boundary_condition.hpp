#pragma once

#include "include/core/state.hpp"
class BoundaryCondition {
  public:
    BoundaryCondition(Grid &grid)
        : Nx_total_(grid.Nx_total()), Ny_total_(grid.Ny_total()), nG_(grid.nG()) {};

    virtual void apply_BC(State &U) {};

  protected:
    int Nx_total_;
    int Ny_total_;

    int nG_;
};