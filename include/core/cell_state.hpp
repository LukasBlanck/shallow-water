#pragma once

#include "include/constants.hpp"

class CellState {
  public:
    CellState(double h, double hu, double hv) : h_(h), hu_(hu), hv_(hv) {};
    CellState() : h_(0.0), hu_(0.0), hv_(0.0) {}; // generic constructor

    double h() const noexcept { return h_; }
    double hu() const noexcept { return hu_; }
    double hv() const noexcept { return hv_; }

    double &h() noexcept { return h_; }
    double &hu() noexcept { return hu_; }
    double &hv() noexcept { return hv_; }

    double u() const noexcept {
        const double u = (h_ > constants::eps) ? (hu_ / h_) : 0.0;
        return u;
    }
    double v() const noexcept {
        const double v = (h_ > constants::eps) ? (hv_ / h_) : 0.0;
        return v;
    }

  private:
    double h_;
    double hu_;
    double hv_;
};

inline CellState operator+(const CellState &a, const CellState &b) {
    return CellState(a.h() + b.h(), a.hu() + b.hu(), a.hv() + b.hv());
}

inline CellState operator-(const CellState &a, const CellState &b) {
    return CellState(a.h() - b.h(), a.hu() - b.hu(), a.hv() - b.hv());
}

inline CellState operator*(double s, const CellState &a) {
    return CellState(s * a.h(), s * a.hu(), s * a.hv());
}

inline CellState operator*(const CellState &a, double s) { return s * a; }