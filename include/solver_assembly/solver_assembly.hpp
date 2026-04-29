#pragma once

#include "configs/config.hpp"
#include "configs/simulation_config.hpp"
#include "include/backend/serial/reconstruction/muscl.hpp"
#include "include/backend/serial/reconstruction/piecewise_const.hpp"
#include "include/backend/serial/riemann/ROE.hpp"
#include "include/backend/serial/serial_solver.hpp"
#include "include/bathymetry/flat_bathymetry.hpp"
#include "include/bathymetry/gaussian_hill.hpp"
#include "include/bathymetry/no_bathymetry.hpp"
#include "include/boundary/reflecting_walls.hpp"

#include <filesystem>

class SolverAssembly {
  public:
    SolverAssembly(const std::filesystem::path &path) : cfg_(load_config(path)) {};

    void run() {
        validate_config();

        switch (cfg_.backend.type) {
        case BackendType::Serial:
            if (cfg_.boundary.type == BoundaryType::ReflectingWalls &&
                cfg_.solver.reconstruction == ReconstructionType::PiecewiseConst &&
                cfg_.solver.riemann == RiemannType::Rusanov &&
                cfg_.solver.time == TimeIntegratorType::SSPRK3 &&
                cfg_.bathymetry.type == BathymetryType::None) {

                SerialSolver<ReflectingWalls, PiecewiseConst, Rusanov, SSPRK3, None> solver(cfg_);
                solver.run();
                return;
            }

            if (cfg_.boundary.type == BoundaryType::ReflectingWalls &&
                cfg_.solver.reconstruction == ReconstructionType::MUSCL &&
                cfg_.solver.riemann == RiemannType::Rusanov &&
                cfg_.solver.time == TimeIntegratorType::SSPRK3 &&
                cfg_.bathymetry.type == BathymetryType::None) {

                SerialSolver<ReflectingWalls, MUSCL, Rusanov, SSPRK3, None> solver(cfg_);
                solver.run();
                return;
            }

            if (cfg_.boundary.type == BoundaryType::ReflectingWalls &&
                cfg_.solver.reconstruction == ReconstructionType::MUSCL &&
                cfg_.solver.riemann == RiemannType::HLL &&
                cfg_.solver.time == TimeIntegratorType::SSPRK3 &&
                cfg_.bathymetry.type == BathymetryType::None) {

                SerialSolver<ReflectingWalls, MUSCL, HLL, SSPRK3, None> solver(cfg_);
                solver.run();
                return;
            }

            if (cfg_.boundary.type == BoundaryType::ReflectingWalls &&
                cfg_.solver.reconstruction == ReconstructionType::PiecewiseConst &&
                cfg_.solver.riemann == RiemannType::HLL &&
                cfg_.solver.time == TimeIntegratorType::SSPRK3 &&
                cfg_.bathymetry.type == BathymetryType::None) {

                SerialSolver<ReflectingWalls, PiecewiseConst, HLL, SSPRK3, None> solver(cfg_);
                solver.run();
                return;
            }
            // ROE + PiecewiseConst
            if (cfg_.boundary.type == BoundaryType::ReflectingWalls &&
                cfg_.solver.reconstruction == ReconstructionType::PiecewiseConst &&
                cfg_.solver.riemann == RiemannType::ROE &&
                cfg_.solver.time == TimeIntegratorType::SSPRK3 &&
                cfg_.bathymetry.type == BathymetryType::None) {

                SerialSolver<ReflectingWalls, PiecewiseConst, ROE, SSPRK3, None> solver(cfg_);
                solver.run();
                return;
            }
            // ROE + MUSCL
            if (cfg_.boundary.type == BoundaryType::ReflectingWalls &&
                cfg_.solver.reconstruction == ReconstructionType::MUSCL &&
                cfg_.solver.riemann == RiemannType::ROE &&
                cfg_.solver.time == TimeIntegratorType::SSPRK3 &&
                cfg_.bathymetry.type == BathymetryType::None) {

                SerialSolver<ReflectingWalls, MUSCL, ROE, SSPRK3, None> solver(cfg_);
                solver.run();
                return;
            }
            // HLL + MUSCL + FLAT
            if (cfg_.boundary.type == BoundaryType::ReflectingWalls &&
                cfg_.solver.reconstruction == ReconstructionType::MUSCL &&
                cfg_.solver.riemann == RiemannType::HLL &&
                cfg_.solver.time == TimeIntegratorType::SSPRK3 &&
                cfg_.bathymetry.type == BathymetryType::Flat) {

                SerialSolver<ReflectingWalls, MUSCL, HLL, SSPRK3, Flat> solver(cfg_);
                solver.run();
                return;
            }
            // HLL + MUSCL + GaussHill
            if (cfg_.boundary.type == BoundaryType::ReflectingWalls &&
                cfg_.solver.reconstruction == ReconstructionType::MUSCL &&
                cfg_.solver.riemann == RiemannType::HLL &&
                cfg_.solver.time == TimeIntegratorType::SSPRK3 &&
                cfg_.bathymetry.type == BathymetryType::GaussHill) {

                SerialSolver<ReflectingWalls, MUSCL, HLL, SSPRK3, GaussHill> solver(cfg_);
                solver.run();
                return;
            }
            // HLL + PiecewiseConst + GaussHill
            if (cfg_.boundary.type == BoundaryType::ReflectingWalls &&
                cfg_.solver.reconstruction == ReconstructionType::PiecewiseConst &&
                cfg_.solver.riemann == RiemannType::HLL &&
                cfg_.solver.time == TimeIntegratorType::SSPRK3 &&
                cfg_.bathymetry.type == BathymetryType::GaussHill) {

                SerialSolver<ReflectingWalls, PiecewiseConst, HLL, SSPRK3, GaussHill> solver(cfg_);
                solver.run();
                return;
            }
            // HLL + PiecewiseConst + Flat
            if (cfg_.boundary.type == BoundaryType::ReflectingWalls &&
                cfg_.solver.reconstruction == ReconstructionType::PiecewiseConst &&
                cfg_.solver.riemann == RiemannType::HLL &&
                cfg_.solver.time == TimeIntegratorType::SSPRK3 &&
                cfg_.bathymetry.type == BathymetryType::Flat) {

                SerialSolver<ReflectingWalls, PiecewiseConst, HLL, SSPRK3, Flat> solver(cfg_);
                solver.run();
                return;
            }

            throw std::runtime_error("Unsupported serial solver combination");
        }

        throw std::runtime_error("Unsupported backend");
    }

  private:
    void validate_config() const {
        if (cfg_.mesh.Nx <= 0 || cfg_.mesh.Ny <= 0) {
            throw std::runtime_error("Nx and Ny must be positive");
        }

        if (cfg_.mesh.nG < 1) {
            throw std::runtime_error("nG must be at least 1");
        }

        if (cfg_.time.end_time <= 0.0) {
            throw std::runtime_error("end_time must be positive");
        }

        if (cfg_.time.time_steps == 0) {
            throw std::runtime_error("time_steps must be > 0");
        }

        if (cfg_.initial_condition.type == InitialConditionType::GaussInitial) {
            if (cfg_.initial_condition.sigma_x <= 0.0 || cfg_.initial_condition.sigma_y <= 0.0) {
                throw std::runtime_error("Gaussian sigmas must be positive");
            }
        }

        if (cfg_.initial_condition.type == InitialConditionType::GaussInitial) {
            if (cfg_.initial_condition.h0 < 0.0) {
                throw std::runtime_error("h0 must be positive");
            }
        }

        if (cfg_.initial_condition.type == InitialConditionType::StillWater) {
            if (cfg_.initial_condition.h0 < 0.0) {
                throw std::runtime_error("h0 must be positive");
            }
        }

        if (cfg_.solver.reconstruction == ReconstructionType::MUSCL && cfg_.mesh.nG < 2) {
            throw std::runtime_error("MUSCL reconstruction requires nG >= 2");
        }
        if (cfg_.mesh.Lx <= 0.0 || cfg_.mesh.Ly <= 0.0) {
            throw std::runtime_error("Lx and Ly must be positive");
        }

        if (cfg_.time.cfl <= 0.0) {
            throw std::runtime_error("cfl must be positive");
        }

        if (cfg_.time.cfl > 1.0) {
            throw std::runtime_error("cfl should usually be <= 1.0");
        }

        if (cfg_.time.save_every <= 0) {
            throw std::runtime_error("save_every must be > 0");
        }

        if (cfg_.time.save_every > cfg_.time.time_steps) {
            throw std::runtime_error("save_every must not be larger than time_steps");
        }

        // Bathymetry
        if (cfg_.bathymetry.type == BathymetryType::Flat ||
            cfg_.bathymetry.type == BathymetryType::GaussHill) {
            if (!std::isfinite(cfg_.bathymetry.b0)) {
                throw std::runtime_error("bathymetry.b0 must be finite");
            }
        }

        if (cfg_.bathymetry.type == BathymetryType::GaussHill) {
            if (cfg_.bathymetry.bathy_sigma_x <= 0.0 || cfg_.bathymetry.bathy_sigma_y <= 0.0) {
                throw std::runtime_error("Gaussian bathymetry sigmas must be positive");
            }

            if (cfg_.bathymetry.bathy_peak_height < 0.0) {
                throw std::runtime_error("bathymetry.bathy_peak_height must be non-negative");
            }

            if (cfg_.bathymetry.bathy_x0 < 0.0 || cfg_.bathymetry.bathy_x0 > cfg_.mesh.Lx ||
                cfg_.bathymetry.bathy_y0 < 0.0 || cfg_.bathymetry.bathy_y0 > cfg_.mesh.Ly) {
                throw std::runtime_error("Gaussian bathymetry center must lie inside the domain");
            }
        }
        // Initial Conditions
        if (cfg_.initial_condition.h0 <= 0.0) {
            throw std::runtime_error("initial_condition.h0 must be positive");
        }
        if (cfg_.initial_condition.type == InitialConditionType::GaussInitial) {
            if (cfg_.initial_condition.sigma_x <= 0.0 || cfg_.initial_condition.sigma_y <= 0.0) {
                throw std::runtime_error("Gaussian initial-condition sigmas must be positive");
            }

            if (cfg_.initial_condition.x0 < 0.0 || cfg_.initial_condition.x0 > cfg_.mesh.Lx ||
                cfg_.initial_condition.y0 < 0.0 || cfg_.initial_condition.y0 > cfg_.mesh.Ly) {
                throw std::runtime_error(
                    "Gaussian initial-condition center must lie inside the domain");
            }

            if (cfg_.initial_condition.h0 + cfg_.initial_condition.peak_height <= 0.0) {
                throw std::runtime_error("Gaussian initial condition can produce non-positive "
                                         "water height: h0 + peak_height <= 0");
            }
        }
        if (cfg_.initial_condition.type == InitialConditionType::DamBreak) {
            if (cfg_.initial_condition.dam_height <= 0.0) {
                throw std::runtime_error("dam_height must be positive");
            }

            if (cfg_.initial_condition.dam_x < 0.0 || cfg_.initial_condition.dam_x > cfg_.mesh.Lx) {
                throw std::runtime_error("dam_x must lie inside the domain");
            }
        }
        if (cfg_.initial_condition.type == InitialConditionType::DamBreakRadial) {
            if (cfg_.initial_condition.dam_height <= 0.0) {
                throw std::runtime_error("dam_height must be positive");
            }

            if (cfg_.initial_condition.dam_radius <= 0.0) {
                throw std::runtime_error("dam_radius must be positive");
            }

            if (cfg_.initial_condition.dam_x0 < 0.0 ||
                cfg_.initial_condition.dam_x0 > cfg_.mesh.Lx ||
                cfg_.initial_condition.dam_y0 < 0.0 ||
                cfg_.initial_condition.dam_y0 > cfg_.mesh.Ly) {
                throw std::runtime_error("radial dam-break center must lie inside the domain");
            }
        }
        if (cfg_.backend.type == BackendType::Serial && cfg_.backend.threads != 1) {
            // ignore, warn, or throw depending your design
        }
    }
    SimulationConfig cfg_;
};