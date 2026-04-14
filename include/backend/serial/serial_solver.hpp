#pragma once

#include "configs/config.hpp"

#include "include/backend/serial/finite_volume/finite_volume.hpp"
#include "include/backend/serial/flux_assembly/flux_assembly.hpp"
#include "include/backend/serial/reconstruction/piecewise_const.hpp"
#include "include/backend/serial/riemann/rusanov.hpp"
#include "include/backend/serial/ssp_rk3/ssp_rk3.hpp"
#include "include/core/grid.hpp"
#include "include/core/state.hpp"
#include "include/core/xflux_field.hpp"
#include "include/core/yflux_field.hpp"
#include "include/initial_condition/apply_initial_condition.hpp"
#include "include/io/netCDF_writer.hpp"
#include "include/io/sanity_checks_netdcdf_writer.hpp"
#include "include/solver_assembly/time_utils.hpp"
#include "tests/sanity_checks/sanity_checks.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

template <class Boundary, class Reconstruction, class Riemann, class TimeIntegrator>
class SerialSolver {
  public:
    explicit SerialSolver(const SimulationConfig &cfg)
        : cfg_(cfg), grid_(cfg.mesh.Nx, cfg.mesh.Ny, cfg.mesh.Lx, cfg.mesh.Ly, cfg.mesh.nG),
          dt_(cfg.time.end_time / static_cast<double>(cfg.time.time_steps)),
          end_time_(cfg.time.end_time), steps_(cfg.time.time_steps),
          save_every_(static_cast<std::size_t>(cfg.time.save_every)), U_(grid_), U1_(grid_),
          U2_(grid_), U_next_(grid_), rhs_(grid_), bc_(grid_), recon_(), riemann_(),
          flux_assembly_(), fv_(grid_), time_integrator_(dt_), x_flux_field_(grid_),
          y_flux_field_(grid_),
          writer_(cfg.output.path.string(), grid_, cfg.mesh.spatial_unit_x, cfg.mesh.spatial_unit_y,
                  cfg.mesh.spatial_unit_h, cfg.time.time_unit, cfg.time.save_every),
          sanity_checks_(make_sanity_checks(cfg)),
          sanity_writer_(make_sanity_writer(cfg, dt_, riemann_name_static(),
                                            reconstruction_name_static(),
                                            time_integrator_name_static())),
          riemann_name_(riemann_name_static()), reconstruction_name_(reconstruction_name_static()),
          time_integrator_name_(time_integrator_name_static()) {
        apply_initial_condition(cfg_, grid_, U_);
        bc_.apply_BC(U_);

        for (auto &check : sanity_checks_) {
            check->initialize(U_, grid_);
        }
    }

    void run() {
        // check if dt is stable
        const double dt_stable = compute_stable_dt(U_, grid_, cfg_.time.cfl);
        if (dt_ > dt_stable) {
            throw std::runtime_error("Time step too large: dt = " + std::to_string(dt_) +
                                     ", but CFL requires dt <= " + std::to_string(dt_stable));
        }

        double time = 0.0;

        // ETA
        const bool estimate_eta = cfg_.output.compute_eta;
        const std::size_t eta_probe_step = std::min<std::size_t>(100, steps_);
        bool eta_printed = false;
        const auto start_wall = std::chrono::steady_clock::now();

        // output and sanity check for initial setup
        writer_.write_snapshot(U_, time, dt_, riemann_name_, reconstruction_name_,
                               time_integrator_name_);
        run_sanity_checks(time, 0);

        // SSP-RK3
        for (std::size_t step = 0; step < steps_; ++step) {
            compute_rhs(U_, rhs_);
            time_integrator_.compute_U1(U1_, U_, rhs_);

            compute_rhs(U1_, rhs_);
            time_integrator_.compute_U2(U2_, U1_, U_, rhs_);

            compute_rhs(U2_, rhs_);
            time_integrator_.compute_U_next(U_next_, U2_, U_, rhs_);

            U_ = U_next_;
            time += dt_;

            bc_.apply_BC(U_);

            const std::size_t step_number = step + 1;

            // ETA
            const double elapsed_sec =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start_wall)
                    .count();

            const double eta_sec = estimate_eta_seconds(step_number, steps_, elapsed_sec);
            std::cout << "ETA = " << format_duration(eta_sec) << '\n';

            // output & sanity checks
            if (step_number % save_every_ == 0 || step_number == steps_) {
                writer_.write_snapshot(U_, time);
                run_sanity_checks(time, step_number);
            }
        }
    }

  private:
    // serial helper functions
    void compute_rhs(State &U_in, State &L_out) {
        bc_.apply_BC(U_in);
        flux_assembly_.compute_x_fluxes(U_in, recon_, riemann_, x_flux_field_, grid_);
        flux_assembly_.compute_y_fluxes(U_in, recon_, riemann_, y_flux_field_, grid_);
        fv_.apply_spatial_operator(L_out, x_flux_field_, y_flux_field_);
    }

    void run_sanity_checks(double time, std::size_t step) {
        for (auto &check : sanity_checks_) {
            check->evaluate(U_, grid_, time, step, sanity_writer_.get());
        }
    }

    static std::unique_ptr<SanityCheckNetCDFWriter>
    make_sanity_writer(const SimulationConfig &cfg, double dt, const std::string &riemann_solver,
                       const std::string &reconstruction, const std::string &time_integrator) {

        if (!cfg.sanity_checks.debug) {
            return nullptr;
        }
        if (cfg.sanity_checks.output_path.empty()) {
            throw std::runtime_error(
                "sanity_checks.output_path must be set if debug mode is enabled");
        }

        return std::make_unique<SanityCheckNetCDFWriter>(
            cfg.sanity_checks.output_path.string(), cfg.time.time_unit, cfg.mesh.spatial_unit_h,
            cfg.time.save_every, dt, riemann_solver, reconstruction, time_integrator);
    }

    static std::string riemann_name_static() {
        if constexpr (std::is_same_v<Riemann, Rusanov>) {
            return "Rusanov";
        }
        return "UnknownRiemann";
    }

    static std::string reconstruction_name_static() {
        if constexpr (std::is_same_v<Reconstruction, PiecewiseConst>) {
            return "PiecewiseConst";
        }
        return "UnknownReconstruction";
    }

    static std::string time_integrator_name_static() {
        if constexpr (std::is_same_v<TimeIntegrator, SSPRK3>) {
            return "SSPRK3";
        }
        return "UnknownTimeIntegrator";
    }

  private:
    const SimulationConfig cfg_;

    Grid grid_;
    double dt_{};
    double end_time_{};
    std::size_t steps_{};
    std::size_t save_every_{};

    State U_;
    State U1_;
    State U2_;
    State U_next_;
    State rhs_;

    Boundary bc_;
    Reconstruction recon_;
    Riemann riemann_;
    FluxAssembly flux_assembly_;
    FiniteVolume fv_;
    TimeIntegrator time_integrator_;

    XFluxField x_flux_field_;
    YFluxField y_flux_field_;

    NetCDFWriter writer_;
    std::vector<std::unique_ptr<SanityCheck>> sanity_checks_;
    std::unique_ptr<SanityCheckNetCDFWriter> sanity_writer_;

    const std::string riemann_name_;
    const std::string reconstruction_name_;
    const std::string time_integrator_name_;
};