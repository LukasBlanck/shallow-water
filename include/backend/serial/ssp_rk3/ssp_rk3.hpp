#pragma once

#include "include/core/state.hpp"

class SSPRK3 {
    SSPRK3(double dt) : dt_(dt) {};

    State compute_U1(const State &Un, const State &L_of_U) { return (Un + dt_ * L_of_U); }

    State compute_U2(const State &U1, const State &Un, const State &L_of_U1) {
        return ((3.0 / 4.0) * Un) + ((1.0 / 4.0) * (U1 + dt_ * L_of_U1));
    }

    State compute_U_next(const State &U2, const State &Un, const State &L_of_U2) {
        return ((1.0 / 3.0) * Un) + ((2.0 / 3.0) * (U2 + dt_ * L_of_U2));
    }

  private:
    double dt_;
};