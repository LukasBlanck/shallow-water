#include "include/core/grid.hpp"
#include "include/core/state.hpp"
#include "include/initial_condition/still_water.hpp"
#include "include/io/netCDF_writer.hpp"

int main() {

    Grid grid(4, 4, 2, 2, 1);         // build Grid: int Nx, int Ny, double Lx, double Ly, int nG
    State U0(grid);                   // constructor
    still_water(3.0).apply(grid, U0); // initial state

    NetCDFWriter("output/simulation.nc", grid).write_snapshot(U0, 0.0);

    return 0;
}