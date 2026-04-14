#pragma once

#include "include/backend/serial/reconstruction/reconstruction.hpp"
#include "include/core/cell_state.hpp"
#include <cmath>
#include <cstdlib>

inline double minmod_scalar(double a, double b) {
    if (a * b <= 0.0) return 0.0;
    return (a > 0.0 ? 1.0 : -1.0) * std::min(std::abs(a), std::abs(b));
}

class MUSCL : public Reconstruction {
  public:
    // for x fluxes
    CellState U_L(const State &U, int i, int j) const override {
        CellState U_ij(U.h()(i, j), U.hu()(i, j), U.hv()(i, j));
        CellState sigma_x = minmod_x(U, i, j);
        return U_ij + 0.5 * sigma_x;
    };
    CellState U_R(const State &U, int i, int j) const override {
        CellState U_ip1j(U.h()(i + 1, j), U.hu()(i + 1, j), U.hv()(i + 1, j));
        CellState sigma_xp1 = minmod_x(U, i + 1, j);
        return U_ip1j - 0.5 * sigma_xp1;
    };

    // for y fluxes
    CellState U_B(const State &U, int i, int j) const override {
        CellState U_ij(U.h()(i, j), U.hu()(i, j), U.hv()(i, j));
        CellState sigma_y = minmod_y(U, i, j);
        return U_ij + 0.5 * sigma_y;
    };
    CellState U_T(const State &U, int i, int j) const override {
        CellState U_ijp1(U.h()(i, j + 1), U.hu()(i, j + 1), U.hv()(i, j + 1));
        CellState sigma_y = minmod_y(U, i, j + 1);
        return U_ijp1 - 0.5 * sigma_y;
    };

  private:
    CellState minmod_x(const State &U, int i, int j) const {
        double dhL = U.h()(i, j) - U.h()(i - 1, j);
        double dhR = U.h()(i + 1, j) - U.h()(i, j);

        double dhuL = U.hu()(i, j) - U.hu()(i - 1, j);
        double dhuR = U.hu()(i + 1, j) - U.hu()(i, j);

        double dhvL = U.hv()(i, j) - U.hv()(i - 1, j);
        double dhvR = U.hv()(i + 1, j) - U.hv()(i, j);

        return CellState(minmod_scalar(dhL, dhR), minmod_scalar(dhuL, dhuR),
                         minmod_scalar(dhvL, dhvR));
    }

    CellState minmod_y(const State &U, int i, int j) const {
        double dhB = U.h()(i, j) - U.h()(i, j - 1);
        double dhT = U.h()(i, j + 1) - U.h()(i, j);

        double dhuB = U.hu()(i, j) - U.hu()(i, j - 1);
        double dhuT = U.hu()(i, j + 1) - U.hu()(i, j);

        double dhvB = U.hv()(i, j) - U.hv()(i, j - 1);
        double dhvT = U.hv()(i, j + 1) - U.hv()(i, j);

        return CellState(minmod_scalar(dhB, dhT), minmod_scalar(dhuB, dhuT),
                         minmod_scalar(dhvB, dhvT));
    }
};
