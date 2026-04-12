#pragma once

#include "include/backend/serial/reconstruction/reconstruction.hpp"
#include "include/backend/serial/riemann/riemann.hpp"
#include "include/core/state.hpp"
#include "include/core/xflux_field.hpp"

class FluxAssembly {
    void compute_x_fluxes(const State &U, const Reconstruction &recon,
                          const RiemannSolver &riemann_solver, XFluxField &Fx, Grid &grid) {
        for (int i = 0; i < grid.Nx_total() + 1; ++i) {
            for (int j = 0; j < grid.Ny_total(); ++j) {
                CellState UL = recon.U_L(U, i, j);
                CellState UR = recon.U_R(U, i, j);

                CellState Fhat = riemann_solver.x_flux(UL, UR);

                Fx.h()(i, j) = Fhat.h();
                Fx.hu()(i, j) = Fhat.hu();
                Fx.hv()(i, j) = Fhat.hv();
            }
        }
    }
};