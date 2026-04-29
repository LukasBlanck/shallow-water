#pragma once

#include "include/bathymetry/bathymetry.hpp"
#include "include/core/array2D.hpp"

#include <cmath>

class GaussHill : public Bathymetry {
  public:
    GaussHill(double peak_height, double sigma_x, double sigma_y, double x0, double y0, double b0)
        : peak_height_(peak_height), sigma_x_(sigma_x), sigma_y_(sigma_y), x0_(x0), y0_(y0),
          b0_(b0) {};

    static constexpr bool enabled = true;

    void apply(const Grid &grid, Array2D &B) const override {
        const int nG = grid.nG();

        for (int i = 0; i < grid.Nx(); i++) {
            for (int j = 0; j < grid.Ny(); j++) {
                const double x = grid.x_center(i);
                const double y = grid.y_center(j);

                B(i + nG, j + nG) =
                    b0_ +
                    peak_height_ * std::exp(-((x - x0_) * (x - x0_)) / (2.0 * sigma_x_ * sigma_x_) -
                                            ((y - y0_) * (y - y0_)) / (2.0 * sigma_y_ * sigma_y_));
            }
        }
    }

  private:
    double peak_height_;
    double sigma_x_;
    double sigma_y_;
    double x0_;
    double y0_;
    double b0_;
};