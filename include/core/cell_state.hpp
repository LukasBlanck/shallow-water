#pragma once

#include "include/constants.hpp"
#include <array>

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

inline CellState operator*(const CellState &a, double s) {
    return s * a;
}

inline CellState operator/(const CellState &a, double s) {
    return CellState(a.h() / s, a.hu() / s, a.hv() / s);
}

// vector operators for ROE
inline std::array<double, 3> operator*(double a, const std::array<double, 3> &s) {
    return std::array<double, 3>{s[0] * a, s[1] * a, s[2] * a};
}

inline std::array<double, 3> operator+(const std::array<double, 3> &a,
                                       const std::array<double, 3> &s) {
    return std::array<double, 3>{s[0] + a[0], s[1] + a[1], s[2] + a[2]};
}

inline CellState operator-(const CellState &a, const std::array<double, 3> &s) {
    return CellState(a.h() - s[0], a.hu() - s[1], a.hv() - s[2]);
}