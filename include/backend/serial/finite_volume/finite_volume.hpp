#pragma once

#include "include/core/grid.hpp"
#include "include/core/state.hpp"
#include "include/core/xflux_field.hpp"
#include "include/core/yflux_field.hpp"

class FiniteVolume {
  public:
    FiniteVolume(const Grid &grid, XFluxField &Fx, YFluxField &Fy)
        : Nx_total_(grid.Nx_total()), Ny_total_(grid.Ny_total()), dx_(grid.dx()), dy_(grid.dy()),
          Fx_(Fx), Fy_(Fy) {};

    void apply_spatial_operator(State &U, State &lhs) const {
        for (int i = nG_; i < (Nx_total_ - nG_); i++) {
            for (int j = nG_; j < (Ny_total_ - nG_); j++) {
                lhs.h()(i, j) = (-1 / dx_) * (Fx_.h()(i + 1, j) - Fx_.h()(i, j)) -
                                (1 / dy_) * (Fy_.h()(i, j + 1) - Fy_.h()(i, j));

                lhs.hu()(i, j) = (-1 / dx_) * (Fx_.hu()(i + 1, j) - Fx_.hu()(i, j)) -
                                 (1 / dy_) * (Fy_.hu()(i, j + 1) - Fy_.hu()(i, j));

                lhs.hv()(i, j) = (-1 / dx_) * (Fx_.hv()(i + 1, j) - Fx_.hv()(i, j)) -
                                 (1 / dy_) * (Fy_.hv()(i, j + 1) - Fy_.hv()(i, j));
            }
        }
    };

  private:
    int Nx_total_;
    int Ny_total_;

    double dx_;
    double dy_;

    double nG_;
    // those states come from piece.const.rec. or MUSCL
    XFluxField &Fx_;
    YFluxField &Fy_;
};