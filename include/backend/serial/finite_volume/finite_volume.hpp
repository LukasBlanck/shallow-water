#pragma once

#include "include/core/grid.hpp"
#include "include/core/state.hpp"
#include "include/core/xflux_field.hpp"
#include "include/core/yflux_field.hpp"

class FiniteVolume {
  public:
    FiniteVolume(const Grid &grid)
        : Nx_total_(grid.Nx_total()), Ny_total_(grid.Ny_total()), dx_(grid.dx()), dy_(grid.dy()),
          nG_(grid.nG()) {};

    void apply_spatial_operator(State &lhs, const XFluxField &Fx, const YFluxField &Fy) const {
        for (int i = nG_; i < (Nx_total_ - nG_); i++) {
            for (int j = nG_; j < (Ny_total_ - nG_); j++) {
                lhs.h()(i, j) = (-1 / dx_) * (Fx.h()(i, j) - Fx.h()(i - 1, j)) -
                                (1 / dy_) * (Fy.h()(i, j) - Fy.h()(i, j - 1));

                lhs.hu()(i, j) = (-1 / dx_) * (Fx.hu()(i, j) - Fx.hu()(i - 1, j)) -
                                 (1 / dy_) * (Fy.hu()(i, j) - Fy.hu()(i, j - 1));

                lhs.hv()(i, j) = (-1 / dx_) * (Fx.hv()(i, j) - Fx.hv()(i - 1, j)) -
                                 (1 / dy_) * (Fy.hv()(i, j) - Fy.hv()(i, j - 1));
            }
        }
    };

  private:
    int Nx_total_;
    int Ny_total_;

    double dx_;
    double dy_;

    int nG_;
};