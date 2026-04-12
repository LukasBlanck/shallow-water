#pragma once

#include "include/initial_condition/initial_condition.hpp"
#include <cmath>

class GaussInitial : public InitialCondition {
  public:
    GaussInitial(double peak_height, double sigma_x, double sigma_y, double x0, double y0,
                 double h0)
        : peak_height_(peak_height), sigma_x_(sigma_x), sigma_y_(sigma_y), x0_(x0), y0_(y0),
          h0_(h0) {};

    void apply(const Grid &grid, State &U) const override {
        const int nG = grid.nG();

        for (int i = 0; i < grid.Nx(); i++) {
            for (int j = 0; j < grid.Ny(); j++) {
                const double x = grid.x_center(i);
                const double y = grid.y_center(j);
                U.h()(i + nG, j + nG) =
                    h0_ +
                    peak_height_ * std::exp(-((x - x0_) * (x - x0_)) / (2 * sigma_x_ * sigma_x_) -
                                            ((y - y0_) * (y - y0_)) / (2 * sigma_y_ * sigma_y_));
                U.hu()(i + nG, j + nG) = 0.0;
                U.hv()(i + nG, j + nG) = 0.0;
            }
        }
    }

  private:
    double peak_height_;
    double sigma_x_;
    double sigma_y_;
    double x0_;
    double y0_;
    double h0_;
};