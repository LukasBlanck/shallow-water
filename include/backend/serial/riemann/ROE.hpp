#pragma once

#include "include/backend/serial/riemann/riemann.hpp"
#include "include/core/cell_state.hpp"
#include <array>
#include <cmath>
#include <cstdlib>
class ROE : public RiemannSolver {
  public:
    CellState x_flux(const CellState &U_L, const CellState &U_R) const {
        double h = (U_L.h() + U_R.h()) / 2.0;
        double c = std::sqrt(constants::g * h);

        double u_roe = compute_roe_avg(U_L.u(), U_L.h(), U_R.u(), U_R.h());
        double v_roe = compute_roe_avg(U_L.v(), U_L.h(), U_R.v(), U_R.h());

        double lamda1 = u_roe - c;
        double lamda2 = u_roe;
        double lamda3 = u_roe + c;

        std::array<double, 3> r1{1, lamda1, v_roe};
        std::array<double, 3> r2{0, 0, 1};
        std::array<double, 3> r3{1, lamda3, v_roe};

        CellState dU = U_R - U_L;
        const double dh = dU.h();
        const double dhu = dU.hu();
        const double dhv = dU.hv();

        double alpha1 = 0.5 * (dh - ((dhu - u_roe * dh) / c));
        double alpha3 = 0.5 * (dh + ((dhu - u_roe * dh) / c));
        double alpha2 = dhv - v_roe * dh;

        return 0.5 * (F(U_L) + F(U_R)) -
               0.5 * (std::abs(lamda1) * alpha1 * r1 + std::abs(lamda2) * alpha2 * r2 +
                      std::abs(lamda3) * alpha3 * r3);
    };
    CellState y_flux(const CellState &U_B, const CellState &U_T) const {
        double h = (U_B.h() + U_T.h()) / 2.0;
        double c = std::sqrt(constants::g * h);

        double u_roe = compute_roe_avg(U_B.u(), U_B.h(), U_T.u(), U_T.h());
        double v_roe = compute_roe_avg(U_B.v(), U_B.h(), U_T.v(), U_T.h());

        double lamda1 = v_roe - c;
        double lamda2 = v_roe;
        double lamda3 = v_roe + c;

        std::array<double, 3> r1{1, u_roe, lamda1};
        std::array<double, 3> r2{0, 1, 0};
        std::array<double, 3> r3{1, u_roe, lamda3};

        CellState dU = U_T - U_B;
        const double dh = dU.h();
        const double dhu = dU.hu();
        const double dhv = dU.hv();

        double alpha1 = 0.5 * (dh - ((dhv - v_roe * dh) / c));
        double alpha3 = 0.5 * (dh + ((dhv - v_roe * dh) / c));
        double alpha2 = dhu - u_roe * dh;

        return 0.5 * (G(U_B) + G(U_T)) -
               0.5 * (std::abs(lamda1) * alpha1 * r1 + std::abs(lamda2) * alpha2 * r2 +
                      std::abs(lamda3) * alpha3 * r3);
    };

  private:
    double compute_roe_avg(double u_one, double h_one, double u_two, double h_two) const {
        return (((std::sqrt(h_one) * u_one) + (std::sqrt(h_two) * u_two)) /
                (std::sqrt(h_one) + std::sqrt(h_two)));
    }
};