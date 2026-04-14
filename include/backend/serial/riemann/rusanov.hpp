#pragma once

#include "include/backend/serial/riemann/riemann.hpp"
#include "include/constants.hpp"
#include "include/core/cell_state.hpp"
#include <algorithm>

class Rusanov : public RiemannSolver {
  public:
    Rusanov() = default;

    CellState x_flux(const CellState &U_L, const CellState &U_R) const override {
        double a = get_dissipation_speed(U_L.u(), U_L.h(), U_R.u(), U_R.h());
        CellState lhs = 0.5 * (F(U_L) + F(U_R)) - 0.5 * a * (U_R - U_L);
        return lhs;
    }

    CellState y_flux(const CellState &U_B, const CellState &U_T) const override {
        double a = get_dissipation_speed(U_B.v(), U_B.h(), U_T.v(), U_T.h());
        CellState lhs = 0.5 * (G(U_B) + G(U_T)) - 0.5 * a * (U_T - U_B);
        return lhs;
    }

  private:
    double get_dissipation_speed(double u_L_B, double h_L_B, double u_R_T, double h_R_T) const {
        h_R_T = std::max(h_R_T, 0.0); // to avoid negative square root
        h_L_B = std::max(h_L_B, 0.0); // to avoid negative square root

        return std::max((std::abs(u_L_B) + std::sqrt(constants::g * h_L_B)),
                        (std::abs(u_R_T) + std::sqrt(constants::g * h_R_T)));
    }
};