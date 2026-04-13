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
#include "include/initial_condition/gauss_initial.hpp"
#include "include/io/netCDF_writer.hpp"

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>

template <class Boundary, class Reconstruction, class Riemann, class TimeIntegrator>
class SerialSolver {
  public:
    explicit SerialSolver(const SimulationConfig &cfg)
        : cfg_(cfg), grid_(cfg.mesh.Nx, cfg.mesh.Ny, cfg.mesh.Lx, cfg.mesh.Ly, cfg.mesh.nG),
          dt_(cfg.time.end_time / static_cast<double>(cfg.time.time_steps)),
          end_time_(cfg.time.end_time), steps_(cfg.time.time_steps),
          save_every_(static_cast<std::size_t>(cfg.time.save_every)),

          U_(grid_), U1_(grid_), U2_(grid_), U_next_(grid_), rhs_(grid_),

          bc_(grid_), recon_(), riemann_(), flux_assembly_(), fv_(grid_), time_integrator_(dt_),

          x_flux_field_(grid_), y_flux_field_(grid_),

          writer_(cfg.output.path.string(), grid_, cfg.mesh.spatial_unit_x, cfg.mesh.spatial_unit_y,
                  cfg.mesh.spatial_unit_h, cfg.time.time_unit, cfg.time.save_every) {
        apply_initial_condition();

        // keep ghost cells valid before first write / rhs
        bc_.apply_BC(U_);
    }

    void run() {
        const double dt_stable = compute_stable_dt(U_);

        if (dt_ > dt_stable) {
            throw std::runtime_error("Time step too large: dt = " + std::to_string(dt_) +
                                     ", but CFL requires dt <= " + std::to_string(dt_stable));
        }

        double time = 0.0;

        writer_.write_snapshot(U_, time, dt_, riemann_name(), reconstruction_name(),
                               time_integrator_name());

        for (std::size_t step = 0; step < steps_; ++step) {
            // Stage 1
            compute_rhs(U_, rhs_);
            time_integrator_.compute_U1(U1_, U_, rhs_);

            // Stage 2
            compute_rhs(U1_, rhs_);
            time_integrator_.compute_U2(U2_, U1_, U_, rhs_);

            // Stage 3
            compute_rhs(U2_, rhs_);
            time_integrator_.compute_U_next(U_next_, U2_, U_, rhs_);

            std::cout << "min_h(U1)     = " << min_h(U1_) << "\n";
            std::cout << "min_h(U2)     = " << min_h(U2_) << "\n";
            std::cout << "min_h(U_next) = " << min_h(U_next_) << "\n";

            U_ = U_next_;
            time += dt_;

            bc_.apply_BC(U_);

            const std::size_t step_number = step + 1;
            if (step_number % save_every_ == 0 || step_number == steps_) {
                writer_.write_snapshot(U_, time);
            }

            std::cout << "Step " << step_number << "/" << steps_ << " completed, t = " << time
                      << '\n';
        }
    }

  private:
    void compute_rhs(State &U_in, State &L_out) {
        bc_.apply_BC(U_in);

        flux_assembly_.compute_x_fluxes(U_in, recon_, riemann_, x_flux_field_, grid_);
        flux_assembly_.compute_y_fluxes(U_in, recon_, riemann_, y_flux_field_, grid_);

        fv_.apply_spatial_operator(L_out, x_flux_field_, y_flux_field_);
    }

    double min_h(const State &U) const {
        double m = 1e100;
        for (int i = grid_.nG(); i < grid_.nG() + grid_.Nx(); ++i) {
            for (int j = grid_.nG(); j < grid_.nG() + grid_.Ny(); ++j) {
                m = std::min(m, U.h()(i, j));
            }
        }
        return m;
    }

    void apply_initial_condition() {
        switch (cfg_.initial_condition.type) {
        case InitialConditionType::GaussInitial: {
            GaussInitial ic(cfg_.initial_condition.peak_height, cfg_.initial_condition.sigma_x,
                            cfg_.initial_condition.sigma_y, cfg_.initial_condition.x0,
                            cfg_.initial_condition.y0, cfg_.initial_condition.h0);
            ic.apply(grid_, U_);
            return;
        }
        }

        throw std::runtime_error("Unsupported initial condition");
    }

    std::string riemann_name() const {
        if constexpr (std::is_same_v<Riemann, Rusanov>) {
            return "Rusanov";
        } else {
            return "UnknownRiemann";
        }
    }

    std::string reconstruction_name() const {
        if constexpr (std::is_same_v<Reconstruction, PiecewiseConst>) {
            return "PiecewiseConst";
        } else {
            return "UnknownReconstruction";
        }
    }

    std::string time_integrator_name() const {
        if constexpr (std::is_same_v<TimeIntegrator, SSPRK3>) {
            return "SSPRK3";
        } else {
            return "UnknownTimeIntegrator";
        }
    }

    double compute_stable_dt(const State &U) const {
        double dt_limit = std::numeric_limits<double>::max();

        for (int i = grid_.nG(); i < grid_.nG() + grid_.Nx(); i++) {
            for (int j = grid_.nG(); j < grid_.nG() + grid_.Ny(); j++) {
                const double h = U.h()(i, j);
                const double hu = U.hu()(i, j);
                const double hv = U.hv()(i, j);

                if (h <= 0.0) {
                    throw std::runtime_error(
                        "Could not perform stability check for dt becuase dt <= 0!");
                }

                const double u = hu / h;
                const double v = hv / h;

                const double dt_cell =
                    cfg_.time.cfl *
                    (1 / ((std::abs(u) + std::sqrt(constants::g * h)) / grid_.dx() +
                          (std::abs(v) + std::sqrt(constants::g * h)) / grid_.dy()));

                dt_limit = std::min(dt_cell, dt_limit);
            }
        }
        return dt_limit;
    }

  private:
    const SimulationConfig cfg_;

    Grid grid_;
    double dt_{};
    double end_time_{};
    std::size_t steps_{};
    std::size_t save_every_{2};

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
};