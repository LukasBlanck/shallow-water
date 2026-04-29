#pragma once

#include <algorithm>
#include <cmath>

#include "include/constants.hpp"
#include "include/core/array2D.hpp"
#include "include/core/cell_state.hpp"
#include "include/core/grid.hpp"
#include "include/core/state.hpp"
#include "include/core/xflux_field.hpp"
#include "include/core/yflux_field.hpp"

class FluxAssemblyBathy {
  public:
    FluxAssemblyBathy() = default;

    template <class Recon, class Riemann>
    void compute_x_fluxes(const State &U, const Array2D &B, const Recon &recon,
                          const Riemann &riemann_solver, XFluxField &Fx_minus, XFluxField &Fx_plus,
                          const Grid &grid) const {

        const int nG = grid.nG();

        // Build W = (eta, hu, hv), where eta = h + b.
        State W = U;
        W.h() += B;

        for (int i = nG - 1; i < nG + grid.Nx(); ++i) {
            for (int j = nG; j < nG + grid.Ny(); ++j) {

                // Reconstruct W at x-interface (i+1/2, j).
                // UL.h() and UR.h() are eta^- and eta^+.
                CellState WL = recon.U_L(W, i, j);
                CellState WR = recon.U_R(W, i, j);

                const double eta_minus = WL.h();
                const double eta_plus = WR.h();

                const double hu_minus = WL.hu();
                const double hv_minus = WL.hv();

                const double hu_plus = WR.hu();
                const double hv_plus = WR.hv();

                // Reconstruct bottom at x-interface. -> manually here with PiecewiseConst instead
                // of MUSCL const double b_minus = bottom_x_minus(B, i, j); const double b_plus =
                // bottom_x_plus(B, i, j);

                // PiecewiseConst Reconstruction of bottom
                const double b_minus = B(i, j);
                const double b_plus = B(i + 1, j);

                const double b_star = std::max(b_minus, b_plus);

                // Provisional water heights.
                const double h_minus = positive_height(eta_minus - b_minus);
                const double h_plus = positive_height(eta_plus - b_plus);

                // Safe velocities from provisional reconstructed states.
                const double u_minus = safe_velocity(hu_minus, h_minus);
                const double v_minus = safe_velocity(hv_minus, h_minus);

                const double u_plus = safe_velocity(hu_plus, h_plus);
                const double v_plus = safe_velocity(hv_plus, h_plus);

                // Hydrostatic reconstructed water heights.
                const double h_star_minus = positive_height(eta_minus - b_star);
                const double h_star_plus = positive_height(eta_plus - b_star);

                // Corrected left/right states for HLL.
                CellState U_minus_star(h_star_minus, h_star_minus * u_minus,
                                       h_star_minus * v_minus);

                CellState U_plus_star(h_star_plus, h_star_plus * u_plus, h_star_plus * v_plus);

                // HLL flux from hydrostatically corrected states.
                CellState F_hll = riemann_solver.x_flux(U_minus_star, U_plus_star);

                // Hydrostatic pressure corrections.
                const double corr_minus =
                    0.5 * constants::g * (h_minus * h_minus - h_star_minus * h_star_minus);

                const double corr_plus =
                    0.5 * constants::g * (h_plus * h_plus - h_star_plus * h_star_plus);

                CellState F_minus(F_hll.h(), F_hll.hu() + corr_minus, F_hll.hv());

                CellState F_plus(F_hll.h(), F_hll.hu() + corr_plus, F_hll.hv());

                Fx_minus.h()(i, j) = F_minus.h();
                Fx_minus.hu()(i, j) = F_minus.hu();
                Fx_minus.hv()(i, j) = F_minus.hv();

                Fx_plus.h()(i, j) = F_plus.h();
                Fx_plus.hu()(i, j) = F_plus.hu();
                Fx_plus.hv()(i, j) = F_plus.hv();
            }
        }
    }

    template <class Recon, class Riemann>
    void compute_y_fluxes(const State &U, const Array2D &B, const Recon &recon,
                          const Riemann &riemann_solver, YFluxField &Fy_minus, YFluxField &Fy_plus,
                          const Grid &grid) const {

        const int nG = grid.nG();

        // Build W = (eta, hu, hv), where eta = h + b.
        State W = U;
        W.h() += B;

        for (int i = nG; i < nG + grid.Nx(); ++i) {
            for (int j = nG - 1; j < nG + grid.Ny(); ++j) {

                // Reconstruct W at y-interface (i, j+1/2).
                // WB.h() and WT.h() are eta^- and eta^+.
                CellState WB = recon.U_B(W, i, j);
                CellState WT = recon.U_T(W, i, j);

                const double eta_minus = WB.h();
                const double eta_plus = WT.h();

                const double hu_minus = WB.hu();
                const double hv_minus = WB.hv();

                const double hu_plus = WT.hu();
                const double hv_plus = WT.hv();

                // Reconstruct bottom at y-interface.
                // const double b_minus = bottom_y_minus(B, i, j);
                // const double b_plus = bottom_y_plus(B, i, j);

                const double b_minus = B(i, j);
                const double b_plus = B(i, j + 1);

                const double b_star = std::max(b_minus, b_plus);

                // Provisional water heights.
                const double h_minus = positive_height(eta_minus - b_minus);
                const double h_plus = positive_height(eta_plus - b_plus);

                // Safe velocities from provisional reconstructed states.
                const double u_minus = safe_velocity(hu_minus, h_minus);
                const double v_minus = safe_velocity(hv_minus, h_minus);

                const double u_plus = safe_velocity(hu_plus, h_plus);
                const double v_plus = safe_velocity(hv_plus, h_plus);

                // Hydrostatic reconstructed water heights.
                const double h_star_minus = positive_height(eta_minus - b_star);
                const double h_star_plus = positive_height(eta_plus - b_star);

                // Corrected bottom/top states for HLL.
                CellState U_minus_star(h_star_minus, h_star_minus * u_minus,
                                       h_star_minus * v_minus);

                CellState U_plus_star(h_star_plus, h_star_plus * u_plus, h_star_plus * v_plus);

                // HLL flux from hydrostatically corrected states.
                CellState G_hll = riemann_solver.y_flux(U_minus_star, U_plus_star);

                // Hydrostatic pressure corrections.
                const double corr_minus =
                    0.5 * constants::g * (h_minus * h_minus - h_star_minus * h_star_minus);

                const double corr_plus =
                    0.5 * constants::g * (h_plus * h_plus - h_star_plus * h_star_plus);

                CellState G_minus(G_hll.h(), G_hll.hu(), G_hll.hv() + corr_minus);

                CellState G_plus(G_hll.h(), G_hll.hu(), G_hll.hv() + corr_plus);

                Fy_minus.h()(i, j) = G_minus.h();
                Fy_minus.hu()(i, j) = G_minus.hu();
                Fy_minus.hv()(i, j) = G_minus.hv();

                Fy_plus.h()(i, j) = G_plus.h();
                Fy_plus.hu()(i, j) = G_plus.hu();
                Fy_plus.hv()(i, j) = G_plus.hv();
            }
        }
    }

  private:
    static double positive_height(double h) { return std::max(0.0, h); }

    static double safe_velocity(double momentum, double h) {
        if (h < constants::eps) {
            return 0.0;
        }
        return momentum / h;
    }

    static double minmod(double a, double b) {
        if (a * b <= 0.0) {
            return 0.0;
        }

        const double sign = (a > 0.0) ? 1.0 : -1.0;
        return sign * std::min(std::abs(a), std::abs(b));
    }

    static double bottom_slope_x(const Array2D &B, int i, int j) {
        const double db_left = B(i, j) - B(i - 1, j);
        const double db_right = B(i + 1, j) - B(i, j);

        return minmod(db_left, db_right);
    }

    static double bottom_slope_y(const Array2D &B, int i, int j) {
        const double db_bottom = B(i, j) - B(i, j - 1);
        const double db_top = B(i, j + 1) - B(i, j);

        return minmod(db_bottom, db_top);
    }

    // b^- at x-interface (i+1/2, j), reconstructed from cell (i, j).
    static double bottom_x_minus(const Array2D &B, int i, int j) {
        return B(i, j) + 0.5 * bottom_slope_x(B, i, j);
    }

    // b^+ at x-interface (i+1/2, j), reconstructed from cell (i+1, j).
    static double bottom_x_plus(const Array2D &B, int i, int j) {
        return B(i + 1, j) - 0.5 * bottom_slope_x(B, i + 1, j);
    }

    // b^- at y-interface (i, j+1/2), reconstructed from cell (i, j).
    static double bottom_y_minus(const Array2D &B, int i, int j) {
        return B(i, j) + 0.5 * bottom_slope_y(B, i, j);
    }

    // b^+ at y-interface (i, j+1/2), reconstructed from cell (i, j+1).
    static double bottom_y_plus(const Array2D &B, int i, int j) {
        return B(i, j + 1) - 0.5 * bottom_slope_y(B, i, j + 1);
    }
};