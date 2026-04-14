#pragma once

#include "include/initial_condition/initial_condition.hpp"
class DamBreak : public InitialCondition {

  public:
    DamBreak(const double dam_height, const double h0, const double dam_x)
        : dam_height_(dam_height), h0_(h0), dam_x_(dam_x) {}
    void apply(const Grid &grid, State &U) const {

        const int nG = grid.nG();

        for (int i = 0; i < grid.Nx(); i++) {
            for (int j = 0; j < grid.Ny(); j++) {
                const double x = grid.x_center(i);
                const double y = grid.y_center(j);

                if (x < dam_x_) {
                    U.h()(i + nG, j + nG) = h0_;
                    U.hu()(i + nG, j + nG) = 0.0;
                    U.hv()(i + nG, j + nG) = 0.0;
                } else {
                    U.h()(i + nG, j + nG) = dam_height_;
                    U.hu()(i + nG, j + nG) = 0.0;
                    U.hv()(i + nG, j + nG) = 0.0;
                }
            }
        }
    }

  private:
    const double dam_height_;
    const double h0_;
    const double dam_x_;
};