#include "configs/config.hpp"
#include "configs/simulation_config.hpp"
#include "include/solver_assembly/solver_assembly.hpp"
#include <cstdio>

using namespace std;

int main() {
    SimulationConfig cfg_loader;
    printf("\n==============================================================\n");
    printf("Loading Configuration from configs/simulation_config.toml\n");
    printf("==============================================================\n");
    load_config("configs/simulation_config.toml");
    printf("Succesfully loaded configurations.\n\n");

    printf("\n==============================================================");
    printf("Assembling Solver...");
    printf("==============================================================\n");
    SolverAssembly SolverAssembly("configs/simulation_config.toml");
    printf("Solver succesfully assembled.\n\n");

    printf("\n==============================================================");
    printf("Running simulation...");
    printf("==============================================================\n");
    SolverAssembly.run(); // calls include/backend/serial/serial_solver.hpp
    printf("Done.\n\n");

    return 0;
}