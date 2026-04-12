#pragma once
#include "include/core/array2D.hpp"
#include "include/core/grid.hpp"

class State {
    // This class contains arrays of double h, hu and hv for ALL cells (inclusively ghost cells)
  public:
    State(const Grid &grid)
        : h_(grid.Nx_total(), grid.Ny_total()), hu_(grid.Nx_total(), grid.Ny_total()),
          hv_(grid.Nx_total(), grid.Ny_total()) {};

    State(int Nx, int Ny) : h_(Nx, Ny), hu_(Nx, Ny), hv_(Nx, Ny) {}

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

inline void check_same_shape(const State &a, const State &b, const char *op) {
    if (a.h().Nx() != b.h().Nx() || a.h().Ny() != b.h().Ny()) {
        throw std::invalid_argument(std::string("State sizes do not match for ") + op);
    }
}

inline State &operator+=(State &lhs, const State &rhs) {
    check_same_shape(lhs, rhs, "addition");
    lhs.h() += rhs.h();
    lhs.hu() += rhs.hu();
    lhs.hv() += rhs.hv();
    return lhs;
}

inline State &operator-=(State &lhs, const State &rhs) {
    check_same_shape(lhs, rhs, "subtraction");
    lhs.h() -= rhs.h();
    lhs.hu() -= rhs.hu();
    lhs.hv() -= rhs.hv();
    return lhs;
}

inline State &operator*=(State &lhs, double s) noexcept {
    lhs.h() *= s;
    lhs.hu() *= s;
    lhs.hv() *= s;
    return lhs;
}

inline State operator+(State lhs, const State &rhs) {
    lhs += rhs;
    return lhs;
}

inline State operator-(State lhs, const State &rhs) {
    lhs -= rhs;
    return lhs;
}

inline State operator*(State lhs, double s) noexcept {
    lhs *= s;
    return lhs;
}

inline State operator*(double s, State rhs) noexcept {
    rhs *= s;
    return rhs;
}