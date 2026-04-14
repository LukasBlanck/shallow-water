#pragma once

#include "include/initial_condition/initial_condition.hpp"

class DamBreakRadial : public InitialCondition {
  public:
    DamBreakRadial(double inner_height, double outer_height, double center_x, double center_y,
                   double radius)
        : inner_height_(inner_height), outer_height_(outer_height), center_x_(center_x),
          center_y_(center_y), radius_(radius) {}

    void apply(const Grid &grid, State &U) const override {
        const int nG = grid.nG();
        const double r2 = radius_ * radius_;

        for (int i = 0; i < grid.Nx(); ++i) {
            for (int j = 0; j < grid.Ny(); ++j) {
                const double x = grid.x_center(i);
                const double y = grid.y_center(j);

                const double dx = x - center_x_;
                const double dy = y - center_y_;
                const bool inside = (dx * dx + dy * dy <= r2);

                U.h()(i + nG, j + nG) = inside ? inner_height_ : outer_height_;
                U.hu()(i + nG, j + nG) = 0.0;
                U.hv()(i + nG, j + nG) = 0.0;
            }
        }
    }

  private:
    const double inner_height_;
    const double outer_height_;
    const double center_x_;
    const double center_y_;
    const double radius_;
};