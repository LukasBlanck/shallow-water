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

    // for fast operators:
    double *data() noexcept { return values_.data(); } // without bound checks!
    const double *data() const noexcept { return values_.data(); }
    std::size_t size() const noexcept { return values_.size(); }

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
    void check_same_shape(const Array2D &rhs, const char *op) const {
        if (Nx_ != rhs.Nx_ || Ny_ != rhs.Ny_) {
            throw std::invalid_argument(std::string("Array2D sizes do not match for ") + op);
        }
    }
};

inline Array2D &operator+=(Array2D &lhs, const Array2D &rhs) {
    if (lhs.Nx() != rhs.Nx() || lhs.Ny() != rhs.Ny()) {
        throw std::invalid_argument("Array2D sizes do not match for addition");
    }

    double *lhs_ptr = lhs.data();
    const double *rhs_ptr = rhs.data();
    const std::size_t n = lhs.size();

    for (std::size_t k = 0; k < n; ++k) {
        lhs_ptr[k] += rhs_ptr[k];
    }
    return lhs;
}

inline Array2D &operator-=(Array2D &lhs, const Array2D &rhs) {
    if (lhs.Nx() != rhs.Nx() || lhs.Ny() != rhs.Ny()) {
        throw std::invalid_argument("Array2D sizes do not match for subtraction");
    }

    double *lhs_ptr = lhs.data();
    const double *rhs_ptr = rhs.data();
    const std::size_t n = lhs.size();

    for (std::size_t k = 0; k < n; ++k) {
        lhs_ptr[k] -= rhs_ptr[k];
    }
    return lhs;
}

inline Array2D &operator*=(Array2D &lhs, double s) noexcept {
    double *ptr = lhs.data();
    const std::size_t n = lhs.size();

    for (std::size_t k = 0; k < n; ++k) {
        ptr[k] *= s;
    }
    return lhs;
}

inline Array2D operator+(Array2D lhs, const Array2D &rhs) {
    lhs += rhs;
    return lhs;
}

inline Array2D operator-(Array2D lhs, const Array2D &rhs) {
    lhs -= rhs;
    return lhs;
}

inline Array2D operator*(Array2D lhs, double s) noexcept {
    lhs *= s;
    return lhs;
}

inline Array2D operator*(double s, Array2D rhs) noexcept {
    rhs *= s;
    return rhs;
}