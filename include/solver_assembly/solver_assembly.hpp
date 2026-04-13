#pragma once

#include "configs/simulation_config.hpp"
#include "include/backend/serial/reconstruction/piecewise_const.hpp"
#include "include/backend/serial/riemann/rusanov.hpp"
#include "include/backend/serial/serial_solver.hpp"
#include "include/backend/serial/ssp_rk3/ssp_rk3.hpp"
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
                cfg_.solver.time == TimeIntegratorType::SSPRK3) {

                SerialSolver<ReflectingWalls, PiecewiseConst, Rusanov, SSPRK3> solver(cfg_);
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

        if (cfg_.backend.type == BackendType::Serial && cfg_.backend.threads != 1) {
            // ignore, warn, or throw depending your design
        }
    }
    SimulationConfig cfg_;
};