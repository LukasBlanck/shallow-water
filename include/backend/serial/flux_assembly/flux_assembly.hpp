#pragma once

#include "include/backend/serial/reconstruction/reconstruction.hpp"
#include "include/backend/serial/riemann/riemann.hpp"
#include "include/core/state.hpp"
#include "include/core/xflux_field.hpp"
#include "include/core/yflux_field.hpp"

class FluxAssembly {
  public:
    FluxAssembly() = default;
    void compute_x_fluxes(const State &U, const Reconstruction &recon,
                          const RiemannSolver &riemann_solver, XFluxField &Fx, Grid &grid) {
        const int nG = grid.nG();
        for (int i = nG - 1; i < nG + grid.Nx(); i++) {
            for (int j = nG; j < nG + grid.Ny(); j++) {
                CellState UL = recon.U_L(U, i, j);
                CellState UR = recon.U_R(U, i, j);

                CellState F_hat = riemann_solver.x_flux(UL, UR);

                Fx.h()(i, j) = F_hat.h();
                Fx.hu()(i, j) = F_hat.hu();
                Fx.hv()(i, j) = F_hat.hv();
            }
        }
    }

    void compute_y_fluxes(const State &U, const Reconstruction &recon,
                          const RiemannSolver &riemann_solver, YFluxField &Fy, Grid &grid) {
        const int nG = grid.nG();
        for (int i = nG; i < nG + grid.Nx(); i++) {
            for (int j = nG - 1; j < nG + grid.Ny(); j++) {

                CellState UB = recon.U_B(U, i, j);
                CellState UT = recon.U_T(U, i, j);

                CellState F_hat = riemann_solver.y_flux(UB, UT);

                Fy.h()(i, j) = F_hat.h();
                Fy.hu()(i, j) = F_hat.hu();
                Fy.hv()(i, j) = F_hat.hv();
            }
        }
    }
};