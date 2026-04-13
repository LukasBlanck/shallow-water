#include "include/solver_assembly/solver_assembly.hpp"
#include <cstdlib>
#include <exception>

using namespace std;

int main() {

    printf("\n==============================================================");
    printf("Assembling Solver...");
    printf("==============================================================\n");
    try {
        SolverAssembly SolverAssembly(
            "configs/simulation_config.toml"); // calls configs/simulation_config.cpp to load toml

        printf("Solver succesfully assembled.\n\n");

        printf("\n==============================================================");
        printf("Running simulation...");
        printf("==============================================================\n");
        SolverAssembly.run(); // calls include/backend/serial/serial_solver.hpp
    } catch (const std::exception &e) {
        std::cerr << "\nSimulation failed: " << e.what() << "\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
    printf("Done.\n\n");

    return 0;
}