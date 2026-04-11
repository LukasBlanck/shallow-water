#include "include/core/grid.hpp"
#include "include/core/state.hpp"
#include "include/initial_condition/gauss_initial.hpp"
// #include "include/initial_condition/still_water.hpp"
#include "include/io/netCDF_writer.hpp"

int main() {

    Grid grid(4, 4, 20, 20, 1); // build Grid: int Nx, int Ny, double Lx, double Ly, int nG
    State U0(grid);             // constructor

    // ------- Initial Condition --------
    // StillWater(3.0).apply(grid, U0); // initial state
    GaussInitial(3.0, 3.0, 3.0, 12.0, 7.0, 2.0)
        .apply(grid, U0); // double peak_height, double sigma_x, double sigma_y, double x0, double
                          // y0, double h0

    NetCDFWriter("output/simulation.nc", grid).write_snapshot(U0, 0.0);

    return 0;
}