#pragma once

#include "include/core/state.hpp"
#include <stdexcept>

class SSPRK3 {
  public:
    explicit SSPRK3(double dt) : dt_(dt) {}

    void compute_U1(State &U1, const State &Un, const State &L_of_U) const {
        check_same_shape(U1, Un, "SSPRK3::compute_U1");
        check_same_shape(U1, L_of_U, "SSPRK3::compute_U1");

        saxpy_state(U1, Un, L_of_U, 1.0, dt_);
    }

    void compute_U2(State &U2, const State &U1, const State &Un, const State &L_of_U1) const {
        check_same_shape(U2, U1, "SSPRK3::compute_U2");
        check_same_shape(U2, Un, "SSPRK3::compute_U2");
        check_same_shape(U2, L_of_U1, "SSPRK3::compute_U2");

        combine_state(U2, Un, U1, L_of_U1,
                      3.0 / 4.0,  // a * Un
                      1.0 / 4.0,  // b * U1
                      dt_ / 4.0); // c * L_of_U1
    }

    void compute_U_next(State &U_next, const State &U2, const State &Un,
                        const State &L_of_U2) const {
        check_same_shape(U_next, U2, "SSPRK3::compute_U_next");
        check_same_shape(U_next, Un, "SSPRK3::compute_U_next");
        check_same_shape(U_next, L_of_U2, "SSPRK3::compute_U_next");

        combine_state(U_next, Un, U2, L_of_U2,
                      1.0 / 3.0,        // a * Un
                      2.0 / 3.0,        // b * U2
                      2.0 * dt_ / 3.0); // c * L_of_U2
    }

  private:
    double dt_;

    static void saxpy_array(Array2D &out, const Array2D &A, const Array2D &B, double alpha,
                            double beta) {
        if (out.Nx() != A.Nx() || out.Ny() != A.Ny() || out.Nx() != B.Nx() || out.Ny() != B.Ny()) {
            throw std::invalid_argument("Array sizes do not match in SSPRK3");
        }

        double *out_ptr = out.data();
        const double *a_ptr = A.data();
        const double *b_ptr = B.data();
        const std::size_t n = out.size();

        for (std::size_t k = 0; k < n; ++k) {
            out_ptr[k] = alpha * a_ptr[k] + beta * b_ptr[k];
        }
    }

    static void combine_array(Array2D &out, const Array2D &A, const Array2D &B, const Array2D &C,
                              double a, double b, double c) {
        if (out.Nx() != A.Nx() || out.Ny() != A.Ny() || out.Nx() != B.Nx() || out.Ny() != B.Ny() ||
            out.Nx() != C.Nx() || out.Ny() != C.Ny()) {
            throw std::invalid_argument("Array sizes do not match in SSPRK3");
        }

        double *out_ptr = out.data();
        const double *a_ptr = A.data();
        const double *b_ptr = B.data();
        const double *c_ptr = C.data();
        const std::size_t n = out.size();

        for (std::size_t k = 0; k < n; ++k) {
            out_ptr[k] = a * a_ptr[k] + b * b_ptr[k] + c * c_ptr[k];
        }
    }

    static void saxpy_state(State &out, const State &A, const State &B, double alpha, double beta) {
        saxpy_array(out.h(), A.h(), B.h(), alpha, beta);
        saxpy_array(out.hu(), A.hu(), B.hu(), alpha, beta);
        saxpy_array(out.hv(), A.hv(), B.hv(), alpha, beta);
    }

    static void combine_state(State &out, const State &A, const State &B, const State &C, double a,
                              double b, double c) {
        combine_array(out.h(), A.h(), B.h(), C.h(), a, b, c);
        combine_array(out.hu(), A.hu(), B.hu(), C.hu(), a, b, c);
        combine_array(out.hv(), A.hv(), B.hv(), C.hv(), a, b, c);
    }
};