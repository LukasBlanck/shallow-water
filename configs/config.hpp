#pragma once

#include <cstddef>
#include <filesystem>
#include <string>

enum class BackendType { Serial };

enum class BoundaryType { ReflectingWalls };

enum class BathymetryType { Flat };

enum class InitialConditionType { GaussInitial, StillWater };

enum class ReconstructionType { PiecewiseConst };

enum class RiemannType { Rusanov };

enum class TimeIntegratorType { SSPRK3 };

struct MeshConfig {
    double Lx{};
    double Ly{};
    int nG{};
    int Nx{};
    int Ny{};
    std::string spatial_unit_x{};
    std::string spatial_unit_y{};
    std::string spatial_unit_h{};
};

struct BathymetryConfig {
    BathymetryType type{};
};

struct BoundaryConfig {
    BoundaryType type{};
};

struct InitialConditionConfig {
    InitialConditionType type{};

    double peak_height{};
    double x0{};
    double y0{};
    double sigma_x{};
    double sigma_y{};
    double h0{};
};

struct SolverConfig {
    ReconstructionType reconstruction{};
    RiemannType riemann{};
    TimeIntegratorType time{};

    std::string limiter{};
    bool positivity_preserving{};
};

struct TimeConfig {
    double end_time{};
    std::size_t time_steps{};
    double cfl{};
    int save_every{};
    std::string time_unit{};
};

struct BackendConfig {
    BackendType type{};
    int threads{};
};

struct OutputConfig {
    std::filesystem::path path;
};

struct SanityCheckConfig {
    bool mass_conservation{false};
    double mass_threshold{};
    bool positivity{false};
    bool debug{false};
    std::filesystem::path output_path;
};

struct SimulationConfig {
    MeshConfig mesh;
    BathymetryConfig bathymetry;
    BoundaryConfig boundary;
    InitialConditionConfig initial_condition;
    SolverConfig solver;
    TimeConfig time;
    BackendConfig backend;
    OutputConfig output;
    SanityCheckConfig sanity_checks;
};