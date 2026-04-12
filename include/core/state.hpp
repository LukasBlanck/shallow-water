#pragma once
#include "include/core/array2D.hpp"
#include "include/core/grid.hpp"

class State {
    // This class contains double h, hu and hv for ALL cells (inclusively ghost cells)
  public:
    State(const Grid &grid)
        : h_(grid.Nx_total(), grid.Ny_total()), hu_(grid.Nx_total(), grid.Ny_total()),
          hv_(grid.Nx_total(), grid.Ny_total()) {};

    Array2D &h() { return h_; } // get and set to reference of array
    Array2D &hu() { return hu_; }
    Array2D &hv() { return hv_; }

    const Array2D &h() const { return h_; }
    const Array2D &hu() const { return hu_; }
    const Array2D &hv() const { return hv_; }

  private:
    Array2D h_; // stores the average vals of all cells
    Array2D hu_;
    Array2D hv_;
};