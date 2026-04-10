#pragma once
#include <vector> // size is not known at compile time -> vector

class Array2D {
  public:
    Array2D(int Nx, int Ny, double init = 0.0)
        : Nx_(Nx), Ny_(Ny), values_(Nx * Ny, init) {}; // initialize with zeros

    double &operator()(int i, int j) { // returns reference to element (get and set)
        return values_[index(i, j)];
    }
    const double &operator()(int i, int j) const { return values_[index(i, j)]; }

    int Nx() const noexcept { return Nx_; } // size of array
    int Ny() const noexcept { return Ny_; }

  private:
    int Nx_;
    int Ny_;
    std::vector<double> values_;

    int index(int i, int j) const {
        if (i < 0 || j < 0 || i >= Nx_ || j >= Ny_) {
            throw std::out_of_range("Array2D index out of range!");
        }
        return i * Ny_ + j; // row major
    }
};
