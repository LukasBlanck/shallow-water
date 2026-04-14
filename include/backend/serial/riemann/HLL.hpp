#pragma once

#include "include/backend/serial/riemann/riemann.hpp"
#include "include/constants.hpp"
#include "include/core/cell_state.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
class HLL : public RiemannSolver {
  public:
    CellState x_flux(const CellState &U_L, const CellState &U_R) const override {
        double sR = dissipation_speed_RT(U_L.u(), U_L.h(), U_R.u(), U_R.h());
        double sL = dissipation_speed_LB(U_L.u(), U_L.h(), U_R.u(), U_R.h());

        if (sL >= 0.0) {
            return F(U_L);
        } else if (sL < 0.0 && sR > 0.0) {
            return (sR * F(U_L) - sL * F(U_R) + sL * sR * (U_R - U_L)) / (sR - sL);
        } else if (sR <= 0.0) {
            return F(U_R);
        } else {
            throw std::runtime_error("HLL could not provide dissipation speed.");
        }
    }
    CellState y_flux(const CellState &U_B, const CellState &U_T) const override {
        double sT = dissipation_speed_RT(U_B.u(), U_B.h(), U_T.u(), U_T.h());
        double sB = dissipation_speed_LB(U_B.u(), U_B.h(), U_T.u(), U_T.h());

        if (sB >= 0.0) {
            return G(U_B);
        } else if (sB < 0.0 && sT > 0.0) {
            return (sT * G(U_B) - sB * G(U_T) + sB * sT * (U_T - U_B)) / (sT - sB);
        } else if (sT <= 0.0) {
            return G(U_T);
        } else {
            throw std::runtime_error("HLL could not provide dissipation speed.");
        }
    };

  private:
    double dissipation_speed_LB(double uLB, double hLB, double uRT, double hRT) const {
        return std::min(uLB - std::sqrt(constants::g * hLB), uRT - std::sqrt(constants::g * hRT));
    }
    double dissipation_speed_RT(double uLB, double hLB, double uRT, double hRT) const {
        return std::min(uLB + std::sqrt(constants::g * hLB), uRT + std::sqrt(constants::g * hRT));
    }
};