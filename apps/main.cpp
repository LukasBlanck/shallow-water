#include "include/backend/serial/finite_volume/finite_volume.hpp"
#include "include/backend/serial/flux_assembly/flux_assembly.hpp"
#include "include/backend/serial/reconstruction/piecewise_const.hpp"
#include "include/backend/serial/riemann/rusanov.hpp"
#include "include/backend/serial/ssp_rk3/ssp_rk3.hpp"
#include "include/boundary/reflecting_walls.hpp"
#include "include/core/grid.hpp"
#include "include/core/state.hpp"
#include "include/core/xflux_field.hpp"
#include "include/core/yflux_field.hpp"
#include "include/initial_condition/gauss_initial.hpp"
// #include "include/initial_condition/still_water.hpp"
#include "include/io/netCDF_writer.hpp"

#include <cstddef>
#include <iostream>
double min_h(const State &U, const Grid &grid) {
    double m = 1e100;
    for (int i = grid.nG(); i < grid.nG() + grid.Nx(); ++i) {
        for (int j = grid.nG(); j < grid.nG() + grid.Ny(); ++j) {
            m = std::min(m, U.h()(i, j));
        }
    }
    return m;
}

int main() {
    // ---------------- Problem setup ----------------
    Grid grid(50, 50, 10.0, 10.0, 1);

    const double end_time = 5.0;
    const std::size_t steps = 2000;
    const double dt = end_time / static_cast<double>(steps);

    // write every Nth step
    const std::size_t save_every = 2;

    // ---------------- Solution / stage states ----------------
    State U(grid);      // current solution
    State U1(grid);     // RK stage 1
    State U2(grid);     // RK stage 2
    State U_next(grid); // next time level
    State rhs(grid);    // spatial operator L(U)

    // ---------------- Initial condition ----------------
    GaussInitial(10.0, 0.7, 0.7, 3.0, 4.0, 2.0).apply(grid, U);
    // StillWater(3.0).apply(grid, U);

    // ---------------- Numerical ingredients ----------------
    ReflectingWalls bc(grid);
    Rusanov rusanov;
    PiecewiseConst recon;
    FluxAssembly flux_assembly;
    FiniteVolume fv(grid);
    SSPRK3 rk3(dt);

    XFluxField x_flux_field(grid);
    YFluxField y_flux_field(grid);

    // ---------------- Output ----------------
    NetCDFWriter writer("output/simulation.nc", grid);

    double time = 0.0;

    // fill ghosts once before first output / first RHS
    bc.apply_BC(U);

    // always write initial condition
    writer.write_snapshot(U, time, dt, "Rusanov", "PiecewiseConst", "SSPRK3");

    // ---------------- Helper: compute RHS = L(U) ----------------
    auto compute_rhs = [&](State &U_in, State &L_out) {
        bc.apply_BC(U_in);

        flux_assembly.compute_x_fluxes(U_in, recon, rusanov, x_flux_field, grid);
        flux_assembly.compute_y_fluxes(U_in, recon, rusanov, y_flux_field, grid);

        fv.apply_spatial_operator(L_out, x_flux_field, y_flux_field);
    };

    // ---------------- Time loop ----------------
    for (std::size_t step = 0; step < steps; ++step) {
        // Stage 1
        compute_rhs(U, rhs);
        rk3.compute_U1(U1, U, rhs);

        // Stage 2
        compute_rhs(U1, rhs);
        rk3.compute_U2(U2, U1, U, rhs);
        double min_h2 = 1e100;

        // Stage 3
        compute_rhs(U2, rhs);
        rk3.compute_U_next(U_next, U2, U, rhs);
        double min_h3 = 1e100;

        std::cout << "min_h(U1)     = " << min_h(U1, grid) << "\n";
        std::cout << "min_h(U2)     = " << min_h(U2, grid) << "\n";
        std::cout << "min_h(U_next) = " << min_h(U_next, grid) << "\n";
        // Advance solution
        U = U_next;
        time += dt;

        // Keep ghost cells valid for output / next operations
        bc.apply_BC(U);

        const std::size_t step_number = step + 1;

        // write every save_every-th step, and always write the last step
        if (step_number % save_every == 0 || step_number == steps) {
            writer.write_snapshot(U, time);
        }

        std::cout << "Step " << step_number << "/" << steps << " completed, t = " << time << '\n';
    }

    return 0;
}