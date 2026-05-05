#pragma once

#include "configs/config.hpp"

#include "include/backend/serial/finite_volume/finite_volume.hpp"
#include "include/backend/serial/finite_volume/finite_volume_bathy.hpp"
#include "include/backend/serial/flux_assembly/flux_assembly.hpp"
#include "include/backend/serial/flux_assembly/flux_assembly_bathy.hpp"
#include "include/backend/serial/reconstruction/muscl.hpp"
#include "include/backend/serial/reconstruction/piecewise_const.hpp"
#include "include/backend/serial/riemann/HLL.hpp"
#include "include/backend/serial/riemann/ROE.hpp"
#include "include/backend/serial/riemann/rusanov.hpp"
#include "include/backend/serial/ssp_rk3/ssp_rk3.hpp"
#include "include/bathymetry/apply_bathymetry.hpp"
#include "include/boundary/reflecting_walls.hpp"
#include "include/core/array2D.hpp"
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
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

struct SolverTimingStats {
    using Clock = std::chrono::steady_clock;

    double cfl_check = 0.0;
    double boundary_conditions = 0.0;
    double x_fluxes = 0.0;
    double y_fluxes = 0.0;
    double finite_volume = 0.0;
    double time_integration = 0.0;
    double positivity_correction = 0.0;
    double output = 0.0;
    double sanity_checks = 0.0;

    std::size_t rhs_calls = 0;
    std::size_t time_steps = 0;

    static Clock::time_point now() { return Clock::now(); }

    static double seconds_since(const Clock::time_point start) {
        return std::chrono::duration<double>(Clock::now() - start).count();
    }

    double measured_solver_time() const {
        return boundary_conditions + x_fluxes + y_fluxes + finite_volume + time_integration +
               positivity_correction + output + sanity_checks;
    }

    void print_summary() const {
        const double total = measured_solver_time();

        std::cout << "\n========== Solver timing summary ==========" << '\n';
        std::cout << "time steps           : " << time_steps << '\n';
        std::cout << "rhs calls            : " << rhs_calls << '\n';
        std::cout << "CFL check            : " << cfl_check << " s" << '\n';
        std::cout << "measured solver time : " << total << " s" << '\n';

        if (total <= 0.0) {
            std::cout << "No positive measured solver time.\n";
            std::cout << "===========================================\n";
            return;
        }

        const auto print_row = [total](const char *name, const double seconds) {
            const double percent = 100.0 * seconds / total;
            std::cout << std::left << std::setw(24) << name << std::right << std::setw(12)
                      << seconds << " s  " << std::setw(8) << percent << " %\n";
        };

        std::cout << std::fixed << std::setprecision(6);
        print_row("boundary conditions", boundary_conditions);
        print_row("x fluxes", x_fluxes);
        print_row("y fluxes", y_fluxes);
        print_row("finite volume", finite_volume);
        print_row("time integration", time_integration);
        print_row("positivity correction", positivity_correction);
        print_row("output", output);
        print_row("sanity checks", sanity_checks);
        std::cout << "===========================================\n";
    }
};

template <class Boundary, class Reconstruction, class Riemann, class TimeIntegrator,
          class Bathymetry>
class SerialSolver {
  public:
    explicit SerialSolver(const SimulationConfig &cfg)
        : cfg_(cfg), grid_(cfg.mesh.Nx, cfg.mesh.Ny, cfg.mesh.Lx, cfg.mesh.Ly, cfg.mesh.nG),
          dt_(cfg.time.end_time / static_cast<double>(cfg.time.time_steps)),
          end_time_(cfg.time.end_time), steps_(cfg.time.time_steps),
          save_every_(static_cast<std::size_t>(cfg.time.save_every)),
          B_(grid_.Nx_total(), grid_.Ny_total()), U_(grid_), U1_(grid_), U2_(grid_), U_next_(grid_),
          rhs_(grid_), bc_(grid_), recon_(), riemann_(), flux_assembly_(), fv_(grid_),
          fv_bathy_(grid_), time_integrator_(dt_), x_flux_field_(grid_), y_flux_field_(grid_),
          x_flux_minus_(grid_), x_flux_plus_(grid_), y_flux_minus_(grid_), y_flux_plus_(grid_),
          writer_(cfg.output.path.string(), grid_, cfg.mesh.spatial_unit_x, cfg.mesh.spatial_unit_y,
                  cfg.mesh.spatial_unit_h, cfg.time.time_unit, cfg.time.save_every),
          sanity_checks_(make_sanity_checks(cfg)),
          sanity_writer_(make_sanity_writer(
              cfg, dt_, riemann_name_static(), reconstruction_name_static(),
              time_integrator_name_static(), boundary_name_static(), bathymetry_name_static())),
          riemann_name_(riemann_name_static()), reconstruction_name_(reconstruction_name_static()),
          time_integrator_name_(time_integrator_name_static()),
          bathymetry_name_(bathymetry_name_static()), boundary_name_(boundary_name_static()) {

        apply_bathymetry(cfg_, grid_, B_);
        bc_.apply_BC(B_);

        writer_.write_bathymetry(B_);

        apply_initial_condition(cfg_, grid_, U_);
        bc_.apply_BC(U_);

        for (auto &check : sanity_checks_) {
            check->initialize(U_, grid_);
        }
    }

    void run() {
        std::cout << std::unitbuf;      // writes terminal output on cluster immediately
        std::cerr << std::unitbuf;

        // Check if dt is stable. This is reported separately and is not included in the
        // solver-time percentages, because it is a one-time pre-loop check.
        double dt_stable = 0.0;
        {
            const auto start = SolverTimingStats::now();
            dt_stable = compute_stable_dt(U_, grid_, cfg_.time.cfl);
            timing_.cfl_check += SolverTimingStats::seconds_since(start);
        }

        if (dt_ > dt_stable) {
            throw std::runtime_error("Time step too large: dt = " + std::to_string(dt_) +
                                     ", but CFL requires dt <= " + std::to_string(dt_stable));
        }

        double time = 0.0;

        // ETA uses wall time and therefore still includes output/sanity checks, just like before.
        const bool estimate_eta = cfg_.output.compute_eta;
        const std::size_t eta_probe_step = std::min<std::size_t>(100, steps_);
        bool eta_printed = false;
        const auto start_wall = std::chrono::steady_clock::now();

        {
            const auto start = SolverTimingStats::now();
            writer_.write_snapshot(U_, time, dt_, riemann_name_, reconstruction_name_,
                                   time_integrator_name_, boundary_name_, bathymetry_name_);
            timing_.output += SolverTimingStats::seconds_since(start);
        }

        {
            const auto start = SolverTimingStats::now();
            run_sanity_checks(time, 0);
            timing_.sanity_checks += SolverTimingStats::seconds_since(start);
        }

        if constexpr (Bathymetry::enabled) {
            printf("\nINFO: Bathymetry enabled.\n\n");

            for (std::size_t step = 0; step < steps_; ++step) {
                compute_rhs_bathy(U_, rhs_);
                {
                    const auto start = SolverTimingStats::now();
                    time_integrator_.compute_U1(U1_, U_, rhs_);
                    timing_.time_integration += SolverTimingStats::seconds_since(start);
                }
                enforce_positivity_timed(U1_, "U1", step);

                compute_rhs_bathy(U1_, rhs_);
                {
                    const auto start = SolverTimingStats::now();
                    time_integrator_.compute_U2(U2_, U1_, U_, rhs_);
                    timing_.time_integration += SolverTimingStats::seconds_since(start);
                }
                enforce_positivity_timed(U2_, "U2", step);

                compute_rhs_bathy(U2_, rhs_);
                {
                    const auto start = SolverTimingStats::now();
                    time_integrator_.compute_U_next(U_next_, U2_, U_, rhs_);
                    timing_.time_integration += SolverTimingStats::seconds_since(start);
                }
                enforce_positivity_timed(U_next_, "U3", step);

                std::swap(U_, U_next_);
                time += dt_;

                {
                    const auto start = SolverTimingStats::now();
                    bc_.apply_BC(U_);
                    timing_.boundary_conditions += SolverTimingStats::seconds_since(start);
                }

                ++timing_.time_steps;
                const std::size_t step_number = step + 1;

                if (estimate_eta && !eta_printed && step_number == eta_probe_step) {
                    const auto now = std::chrono::steady_clock::now();
                    const double elapsed_sec =
                        std::chrono::duration<double>(now - start_wall).count();

                    const double eta_sec = estimate_eta_seconds(step_number, steps_, elapsed_sec);

                    std::cout << "ETA = " << format_duration(eta_sec) << '\n';
                    eta_printed = true;
                }

                if (step_number % save_every_ == 0 || step_number == steps_) {
                    {
                        const auto start = SolverTimingStats::now();
                        writer_.write_snapshot(U_, time);
                        timing_.output += SolverTimingStats::seconds_since(start);
                    }

                    {
                        const auto start = SolverTimingStats::now();
                        run_sanity_checks(time, step_number);
                        timing_.sanity_checks += SolverTimingStats::seconds_since(start);
                    }
                }
            }

        } else {
            for (std::size_t step = 0; step < steps_; ++step) {
                compute_rhs(U_, rhs_);
                {
                    const auto start = SolverTimingStats::now();
                    time_integrator_.compute_U1(U1_, U_, rhs_);
                    timing_.time_integration += SolverTimingStats::seconds_since(start);
                }

                compute_rhs(U1_, rhs_);
                {
                    const auto start = SolverTimingStats::now();
                    time_integrator_.compute_U2(U2_, U1_, U_, rhs_);
                    timing_.time_integration += SolverTimingStats::seconds_since(start);
                }

                compute_rhs(U2_, rhs_);
                {
                    const auto start = SolverTimingStats::now();
                    time_integrator_.compute_U_next(U_next_, U2_, U_, rhs_);
                    timing_.time_integration += SolverTimingStats::seconds_since(start);
                }

                std::swap(U_, U_next_);
                time += dt_;

                {
                    const auto start = SolverTimingStats::now();
                    bc_.apply_BC(U_);
                    timing_.boundary_conditions += SolverTimingStats::seconds_since(start);
                }

                ++timing_.time_steps;
                const std::size_t step_number = step + 1;

                if (estimate_eta && !eta_printed && step_number == eta_probe_step) {
                    const auto now = std::chrono::steady_clock::now();
                    const double elapsed_sec =
                        std::chrono::duration<double>(now - start_wall).count();

                    const double eta_sec = estimate_eta_seconds(step_number, steps_, elapsed_sec);

                    std::cout << "ETA = " << format_duration(eta_sec) << '\n';
                    eta_printed = true;
                }

                if (step_number % save_every_ == 0 || step_number == steps_) {
                    {
                        const auto start = SolverTimingStats::now();
                        writer_.write_snapshot(U_, time);
                        timing_.output += SolverTimingStats::seconds_since(start);
                    }

                    {
                        const auto start = SolverTimingStats::now();
                        run_sanity_checks(time, step_number);
                        timing_.sanity_checks += SolverTimingStats::seconds_since(start);
                    }
                }
            }
        }

        timing_.print_summary();
    }

  private:
    void compute_rhs(State &U_in, State &L_out) {
        ++timing_.rhs_calls;

        {
            const auto start = SolverTimingStats::now();
            bc_.apply_BC(U_in);
            timing_.boundary_conditions += SolverTimingStats::seconds_since(start);
        }

        {
            const auto start = SolverTimingStats::now();
            flux_assembly_.compute_x_fluxes(U_in, recon_, riemann_, x_flux_field_, grid_);
            timing_.x_fluxes += SolverTimingStats::seconds_since(start);
        }

        {
            const auto start = SolverTimingStats::now();
            flux_assembly_.compute_y_fluxes(U_in, recon_, riemann_, y_flux_field_, grid_);
            timing_.y_fluxes += SolverTimingStats::seconds_since(start);
        }

        {
            const auto start = SolverTimingStats::now();
            fv_.apply_spatial_operator(L_out, x_flux_field_, y_flux_field_);
            timing_.finite_volume += SolverTimingStats::seconds_since(start);
        }
    }

    void compute_rhs_bathy(State &U_in, State &L_out) {
        ++timing_.rhs_calls;

        {
            const auto start = SolverTimingStats::now();
            bc_.apply_BC(U_in);
            timing_.boundary_conditions += SolverTimingStats::seconds_since(start);
        }

        {
            const auto start = SolverTimingStats::now();
            flux_assembly_bathy_.compute_x_fluxes(U_in, B_, recon_, riemann_, x_flux_minus_,
                                                  x_flux_plus_, grid_);
            timing_.x_fluxes += SolverTimingStats::seconds_since(start);
        }

        {
            const auto start = SolverTimingStats::now();
            flux_assembly_bathy_.compute_y_fluxes(U_in, B_, recon_, riemann_, y_flux_minus_,
                                                  y_flux_plus_, grid_);
            timing_.y_fluxes += SolverTimingStats::seconds_since(start);
        }

        {
            const auto start = SolverTimingStats::now();
            fv_bathy_.apply_spatial_operator(L_out, x_flux_minus_, x_flux_plus_, y_flux_minus_,
                                             y_flux_plus_);
            timing_.finite_volume += SolverTimingStats::seconds_since(start);
        }
    }

    void run_sanity_checks(double time, std::size_t step) {
        for (auto &check : sanity_checks_) {
            check->evaluate(U_, grid_, time, step, sanity_writer_.get());
        }
    }

    static std::unique_ptr<SanityCheckNetCDFWriter>
    make_sanity_writer(const SimulationConfig &cfg, double dt, const std::string &riemann_solver,
                       const std::string &reconstruction, const std::string &time_integrator,
                       const std::string &boundary_condition, const std::string &bathymetry) {

        if (!cfg.sanity_checks.debug) {
            return nullptr;
        }
        if (cfg.sanity_checks.output_path.empty()) {
            throw std::runtime_error(
                "sanity_checks.output_path must be set if debug mode is enabled");
        }

        return std::make_unique<SanityCheckNetCDFWriter>(
            cfg.sanity_checks.output_path.string(), cfg.time.time_unit, cfg.mesh.spatial_unit_h,
            cfg.time.save_every, dt, riemann_solver, reconstruction, time_integrator,
            boundary_condition, bathymetry);
    }

    static std::string riemann_name_static() {
        if constexpr (std::is_same_v<Riemann, Rusanov>) {
            return "Rusanov";
        }
        if constexpr (std::is_same_v<Riemann, HLL>) {
            return "HLL";
        }
        if constexpr (std::is_same_v<Riemann, ROE>) {
            return "ROE";
        }
        return "UnknownRiemann";
    }

    static std::string reconstruction_name_static() {
        if constexpr (std::is_same_v<Reconstruction, PiecewiseConst>) {
            return "PiecewiseConst";
        }
        if constexpr (std::is_same_v<Reconstruction, MUSCL>) {
            return "MUSCL";
        }
        return "UnknownReconstruction";
    }

    static std::string time_integrator_name_static() {
        if constexpr (std::is_same_v<TimeIntegrator, SSPRK3>) {
            return "SSPRK3";
        }
        return "UnknownTimeIntegrator";
    }

    static std::string boundary_name_static() {
        if constexpr (std::is_same_v<Boundary, ReflectingWalls>) {
            return "Reflecting Walls";
        }
        return "UnknownBoundary";
    }

    static std::string bathymetry_name_static() {
        if constexpr (std::is_same_v<Bathymetry, None>) {
            return "None";
        }
        if constexpr (std::is_same_v<Bathymetry, Flat>) {
            return "Flat";
        }
        if constexpr (std::is_same_v<Bathymetry, GaussHill>) {
            return "Gaussian Hill";
        }
        return "UnknownBathymetry";
    }

    void enforce_positivity_timed(State &U, const char *stage, std::size_t step) {
        const auto start = SolverTimingStats::now();
        enforce_positivity(U, stage, step);
        timing_.positivity_correction += SolverTimingStats::seconds_since(start);
    }

    void enforce_positivity(State &U, const char *stage, std::size_t step) {
        constexpr double h_floor = 1e-8;
        std::size_t corrected = 0;

        const int nG = grid_.nG();

        for (int i = nG; i < nG + grid_.Nx(); ++i) {
            for (int j = nG; j < nG + grid_.Ny(); ++j) {
                if (U.h()(i, j) < h_floor) {
                    U.h()(i, j) = h_floor;
                    U.hu()(i, j) = 0.0;
                    U.hv()(i, j) = 0.0;
                    ++corrected;
                }
                if (!std::isfinite(U.h()(i, j)) || !std::isfinite(U.hu()(i, j)) ||
                    !std::isfinite(U.hv()(i, j))) {
                    throw std::runtime_error("non-finite state detected in positivity correction");
                }
            }
        }

        if (corrected > 0) {
            std::cerr << "WARNING: positivity correction in " << stage << " at step " << step
                      << " corrected " << corrected << " cells\n";
        }
    }

  private:
    const SimulationConfig cfg_;

    Grid grid_;
    double dt_{};
    double end_time_{};
    std::size_t steps_{};
    std::size_t save_every_{};

    Array2D B_; // Bathymetry

    State U_;
    State U1_;
    State U2_;
    State U_next_;
    State rhs_;

    Boundary bc_;
    Reconstruction recon_;
    Riemann riemann_;
    FluxAssembly flux_assembly_;
    FluxAssemblyBathy flux_assembly_bathy_; // Bathymetry
    FiniteVolume fv_;
    FiniteVolumeBathy fv_bathy_; // Bathymetry
    TimeIntegrator time_integrator_;

    XFluxField x_flux_field_;
    YFluxField y_flux_field_;

    XFluxField x_flux_minus_, x_flux_plus_; // Bathymetry
    YFluxField y_flux_minus_, y_flux_plus_; // Bathymetry

    NetCDFWriter writer_;
    std::vector<std::unique_ptr<SanityCheck>> sanity_checks_;
    std::unique_ptr<SanityCheckNetCDFWriter> sanity_writer_;

    SolverTimingStats timing_;

    const std::string riemann_name_;
    const std::string reconstruction_name_;
    const std::string time_integrator_name_;
    const std::string bathymetry_name_;
    const std::string boundary_name_;
};
