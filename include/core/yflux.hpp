#pragma once

#include "include/core/grid.hpp"
#include <stdexcept>
#include <vector>

class YFlux {
  public:
    YFlux(const Grid &grid, double init = 0.0)
        : Nx_total_(grid.Nx_total()), Ny_total_(grid.Ny_total()),
          values_(Nx_total_ * (Ny_total_ + 1), init) {};

    double &operator()(int i, int j) { // returns reference to element (get and set)
        return values_[index(i, j)];
    }
    const double &operator()(int i, int j) const { return values_[index(i, j)]; }

  private:
    int Nx_total_;
    int Ny_total_;

    std::vector<double> values_;

    int index(int i, int j) const {
        if (i < 0 || j < 0 || i >= Nx_total_ || j >= (Ny_total_ + 1)) {
            throw std::out_of_range("Array2D index out of range!");
        }
        return i * (Ny_total_ + 1) + j; // row major
    }
};