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
        for (int i = 0; i < grid.Nx_total(); i++) {
            for (int j = 0; j < grid.Ny_total(); j++) {
                U.h()(i, j) =
                    h0_ + peak_height_ * std::exp(-(i * grid.dx() - x0_) * (i * grid.dx() - x0_) /
                                                      (2 * sigma_x_ * sigma_x_) -
                                                  (j * grid.dy() - y0_) * (j * grid.dy() - y0_) /
                                                      (2 * sigma_y_ * sigma_y_));
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