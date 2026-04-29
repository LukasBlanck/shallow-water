#pragma once

#include "include/core/grid.hpp"
#include "include/core/state.hpp"
#include "include/core/xflux_field.hpp"
#include "include/core/yflux_field.hpp"

class FiniteVolumeBathy {
  public:
    FiniteVolumeBathy(const Grid &grid)
        : Nx_total_(grid.Nx_total()), Ny_total_(grid.Ny_total()), dx_(grid.dx()), dy_(grid.dy()),
          nG_(grid.nG()) {};

    void apply_spatial_operator(State &lhs, const XFluxField &Fx_minus, const XFluxField &Fx_plus,
                                const YFluxField &Fy_minus, const YFluxField &Fy_plus) const {
        for (int i = nG_; i < (Nx_total_ - nG_); i++) {
            for (int j = nG_; j < (Ny_total_ - nG_); j++) {

                lhs.h()(i, j) = (-1.0 / dx_) * (Fx_minus.h()(i, j) - Fx_plus.h()(i - 1, j)) -
                                (1.0 / dy_) * (Fy_minus.h()(i, j) - Fy_plus.h()(i, j - 1));

                lhs.hu()(i, j) = (-1.0 / dx_) * (Fx_minus.hu()(i, j) - Fx_plus.hu()(i - 1, j)) -
                                 (1.0 / dy_) * (Fy_minus.hu()(i, j) - Fy_plus.hu()(i, j - 1));

                lhs.hv()(i, j) = (-1.0 / dx_) * (Fx_minus.hv()(i, j) - Fx_plus.hv()(i - 1, j)) -
                                 (1.0 / dy_) * (Fy_minus.hv()(i, j) - Fy_plus.hv()(i, j - 1));
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