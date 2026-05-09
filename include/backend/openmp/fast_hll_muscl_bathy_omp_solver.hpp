#pragma once

#include "configs/config.hpp"
#include "configs/parse_config_helper.hpp"

#include "include/bathymetry/apply_bathymetry.hpp"
#include "include/core/array2D.hpp"
#include "include/core/grid.hpp"
#include "include/core/state.hpp"
#include "include/initial_condition/apply_initial_condition.hpp"
#include "include/io/netCDF_writer.hpp"
#include "include/io/sanity_checks_netdcdf_writer.hpp"
#include "include/solver_assembly/time_utils.hpp"
#include "tests/sanity_checks/sanity_checks.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#if USE_OPENMP
#include "include/backend/openmp/fast_hll_muscl_bathy_omp_kernels.hpp"
#include <omp.h>

#endif

struct FastOpenMPSolverTimingStats {
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

        std::cout << "\n========== Optimized OpenMP solver timing summary ==========" << '\n';
        std::cout << "time steps           : " << time_steps << '\n';
        std::cout << "rhs calls            : " << rhs_calls << '\n';
        std::cout << "CFL check            : " << cfl_check << " s" << '\n';
        std::cout << "measured solver time : " << total << " s" << '\n';

        if (total <= 0.0) {
            std::cout << "No positive measured solver time.\n";
            std::cout << "=====================================================\n";
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
        std::cout << "=====================================================\n";
    }
};

class FastHLLMUSCLBathyOpenMPSolver {
  public:
    explicit FastHLLMUSCLBathyOpenMPSolver(const SimulationConfig &cfg)
        : cfg_(cfg), grid_(cfg.mesh.Nx, cfg.mesh.Ny, cfg.mesh.Lx, cfg.mesh.Ly, cfg.mesh.nG),
          gv_{grid_.Nx(),       grid_.Ny(),       grid_.nG(), grid_.Nx_total(),
              grid_.Ny_total(), grid_.Ny_total(), grid_.dx(), grid_.dy()},
          dt_(cfg.time.end_time / static_cast<double>(cfg.time.time_steps)),
          end_time_(cfg.time.end_time), steps_(cfg.time.time_steps),
          save_every_(static_cast<std::size_t>(cfg.time.save_every)),
          n_total_(static_cast<std::size_t>(gv_.Nx_total) * static_cast<std::size_t>(gv_.Ny_total)),
          h_(n_total_), hu_(n_total_), hv_(n_total_), h1_(n_total_), hu1_(n_total_), hv1_(n_total_),
          h2_(n_total_), hu2_(n_total_), hv2_(n_total_), rhs_h_(n_total_), rhs_hu_(n_total_),
          rhs_hv_(n_total_), B_(n_total_), fxm_h_(n_total_), fxm_hu_(n_total_), fxm_hv_(n_total_),
          fxp_h_(n_total_), fxp_hu_(n_total_), fxp_hv_(n_total_), fym_h_(n_total_),
          fym_hu_(n_total_), fym_hv_(n_total_), fyp_h_(n_total_), fyp_hu_(n_total_),
          fyp_hv_(n_total_), io_state_(grid_), io_bathy_(grid_.Nx_total(), grid_.Ny_total()),
          writer_(cfg.output.path.string(), grid_, cfg.mesh.spatial_unit_x, cfg.mesh.spatial_unit_y,
                  cfg.mesh.spatial_unit_h, cfg.time.time_unit, cfg.time.save_every),
          sanity_checks_(make_sanity_checks(cfg)), sanity_writer_(make_sanity_writer(cfg, dt_)) {

        if (grid_.nG() < 2) {
            throw std::runtime_error("FastHLLMUSCLBathyOpenMPSolver requires nG >= 2 for MUSCL.");
        }

        Array2D B_tmp(grid_.Nx_total(), grid_.Ny_total());
        apply_bathymetry(cfg_, grid_, B_tmp);
        apply_scalar_reflecting_bc(B_tmp.data());
        copy_from_array2d(B_tmp, B_);
        copy_to_array2d(B_, io_bathy_);
        writer_.write_bathymetry(io_bathy_);

        State U_tmp(grid_);
        apply_initial_condition(cfg_, grid_, U_tmp);
        copy_from_state(U_tmp, h_, hu_, hv_);
        apply_reflecting_bc(h_.data(), hu_.data(), hv_.data());

        for (auto &check : sanity_checks_) {
            fill_io_state();
            check->initialize(io_state_, grid_);
        }
    }

    void run() {
        std::cout << std::unitbuf;
        std::cerr << std::unitbuf;

        std::cout << "INFO: OpenMP environment active\n";
#pragma omp parallel
        {
#pragma omp single
            {
                std::cout << "INFO: OpenMP actual threads in parallel region = "
                          << omp_get_num_threads() << ", max threads = " << omp_get_max_threads()
                          << '\n'
                          << std::flush;
            }
        }

        double dt_stable = 0.0;
        {
            const auto start = FastOpenMPSolverTimingStats::now();
            dt_stable = fast_hll_muscl_bathy_omp::compute_stable_dt(gv_, h_.data(), hu_.data(),
                                                                    hv_.data(), cfg_.time.cfl);
            timing_.cfl_check += FastOpenMPSolverTimingStats::seconds_since(start);
        }

        if (dt_ > dt_stable) {
            throw std::runtime_error("Time step too large: dt = " + std::to_string(dt_) +
                                     ", but CFL requires dt <= " + std::to_string(dt_stable));
        }

        double time = 0.0;
        const bool estimate_eta = cfg_.output.compute_eta;
        const std::size_t eta_probe_step = std::min<std::size_t>(100, steps_);
        bool eta_printed = false;
        const auto start_wall = std::chrono::steady_clock::now();

        {
            const auto start = FastOpenMPSolverTimingStats::now();
            fill_io_state();
            writer_.write_snapshot(io_state_, time, dt_, backend_name_from_cfg(cfg_), "HLL",
                                   "MUSCL", "SSPRK3", "Reflecting Walls",
                                   bathymetry_name_from_cfg(cfg_));
            timing_.output += FastOpenMPSolverTimingStats::seconds_since(start);
        }

        {
            const auto start = FastOpenMPSolverTimingStats::now();
            run_sanity_checks(time, 0);
            timing_.sanity_checks += FastOpenMPSolverTimingStats::seconds_since(start);
        }

        for (std::size_t step = 0; step < steps_; ++step) {
            rk3_step(step);
            time += dt_;

            const std::size_t step_number = step + 1;
            ++timing_.time_steps;

            if (estimate_eta && !eta_printed && step_number == eta_probe_step) {
                const auto now = std::chrono::steady_clock::now();
                const double elapsed_sec = std::chrono::duration<double>(now - start_wall).count();
                const double eta_sec = estimate_eta_seconds(step_number, steps_, elapsed_sec);
                std::cout << "ETA = " << format_duration(eta_sec) << '\n';
                eta_printed = true;
            }

            if (step_number % save_every_ == 0 || step_number == steps_) {
                {
                    const auto start = FastOpenMPSolverTimingStats::now();
                    fill_io_state();
                    writer_.write_snapshot(io_state_, time);
                    timing_.output += FastOpenMPSolverTimingStats::seconds_since(start);
                }

                {
                    const auto start = FastOpenMPSolverTimingStats::now();
                    run_sanity_checks(time, step_number);
                    timing_.sanity_checks += FastOpenMPSolverTimingStats::seconds_since(start);
                }
            }
        }

        timing_.print_summary();
    }

  private:
    void compute_rhs(const double *h, const double *hu, const double *hv, double *rhs_h,
                     double *rhs_hu, double *rhs_hv) {
        ++timing_.rhs_calls;

        {
            const auto start = FastOpenMPSolverTimingStats::now();
            fast_hll_muscl_bathy_omp::compute_x_fluxes(
                gv_, h, hu, hv, B_.data(), fxm_h_.data(), fxm_hu_.data(), fxm_hv_.data(),
                fxp_h_.data(), fxp_hu_.data(), fxp_hv_.data());
            timing_.x_fluxes += FastOpenMPSolverTimingStats::seconds_since(start);
        }

        {
            const auto start = FastOpenMPSolverTimingStats::now();
            fast_hll_muscl_bathy_omp::compute_y_fluxes(
                gv_, h, hu, hv, B_.data(), fym_h_.data(), fym_hu_.data(), fym_hv_.data(),
                fyp_h_.data(), fyp_hu_.data(), fyp_hv_.data());
            timing_.y_fluxes += FastOpenMPSolverTimingStats::seconds_since(start);
        }

        {
            const auto start = FastOpenMPSolverTimingStats::now();
            fast_hll_muscl_bathy_omp::apply_divergence(
                gv_, fxm_h_.data(), fxm_hu_.data(), fxm_hv_.data(), fxp_h_.data(), fxp_hu_.data(),
                fxp_hv_.data(), fym_h_.data(), fym_hu_.data(), fym_hv_.data(), fyp_h_.data(),
                fyp_hu_.data(), fyp_hv_.data(), rhs_h, rhs_hu, rhs_hv);
            timing_.finite_volume += FastOpenMPSolverTimingStats::seconds_since(start);
        }
    }

    void rk3_step(std::size_t step) {
        {
            const auto start = FastOpenMPSolverTimingStats::now();
            apply_reflecting_bc(h_.data(), hu_.data(), hv_.data());
            timing_.boundary_conditions += FastOpenMPSolverTimingStats::seconds_since(start);
        }

        compute_rhs(h_.data(), hu_.data(), hv_.data(), rhs_h_.data(), rhs_hu_.data(),
                    rhs_hv_.data());

        {
            const auto start = FastOpenMPSolverTimingStats::now();
#pragma omp parallel for schedule(static)
            for (std::ptrdiff_t k = 0; k < static_cast<std::ptrdiff_t>(n_total_); ++k) {
                h1_[k] = h_[k] + dt_ * rhs_h_[k];
                hu1_[k] = hu_[k] + dt_ * rhs_hu_[k];
                hv1_[k] = hv_[k] + dt_ * rhs_hv_[k];
            }
            timing_.time_integration += FastOpenMPSolverTimingStats::seconds_since(start);
        }

        {
            const auto start = FastOpenMPSolverTimingStats::now();
            fast_hll_muscl_bathy_omp::enforce_positivity(gv_, h1_.data(), hu1_.data(), hv1_.data());
            timing_.positivity_correction += FastOpenMPSolverTimingStats::seconds_since(start);
        }

        {
            const auto start = FastOpenMPSolverTimingStats::now();
            apply_reflecting_bc(h1_.data(), hu1_.data(), hv1_.data());
            timing_.boundary_conditions += FastOpenMPSolverTimingStats::seconds_since(start);
        }

        compute_rhs(h1_.data(), hu1_.data(), hv1_.data(), rhs_h_.data(), rhs_hu_.data(),
                    rhs_hv_.data());

        {
            const auto start = FastOpenMPSolverTimingStats::now();
#pragma omp parallel for schedule(static)
            for (std::ptrdiff_t k = 0; k < static_cast<std::ptrdiff_t>(n_total_); ++k) {
                h2_[k] = 0.75 * h_[k] + 0.25 * (h1_[k] + dt_ * rhs_h_[k]);
                hu2_[k] = 0.75 * hu_[k] + 0.25 * (hu1_[k] + dt_ * rhs_hu_[k]);
                hv2_[k] = 0.75 * hv_[k] + 0.25 * (hv1_[k] + dt_ * rhs_hv_[k]);
            }
            timing_.time_integration += FastOpenMPSolverTimingStats::seconds_since(start);
        }

        {
            const auto start = FastOpenMPSolverTimingStats::now();
            fast_hll_muscl_bathy_omp::enforce_positivity(gv_, h2_.data(), hu2_.data(), hv2_.data());
            timing_.positivity_correction += FastOpenMPSolverTimingStats::seconds_since(start);
        }

        {
            const auto start = FastOpenMPSolverTimingStats::now();
            apply_reflecting_bc(h2_.data(), hu2_.data(), hv2_.data());
            timing_.boundary_conditions += FastOpenMPSolverTimingStats::seconds_since(start);
        }

        compute_rhs(h2_.data(), hu2_.data(), hv2_.data(), rhs_h_.data(), rhs_hu_.data(),
                    rhs_hv_.data());

        {
            const auto start = FastOpenMPSolverTimingStats::now();
#pragma omp parallel for schedule(static)
            for (std::ptrdiff_t k = 0; k < static_cast<std::ptrdiff_t>(n_total_); ++k) {
                h_[k] = (1.0 / 3.0) * h_[k] + (2.0 / 3.0) * (h2_[k] + dt_ * rhs_h_[k]);

                hu_[k] = (1.0 / 3.0) * hu_[k] + (2.0 / 3.0) * (hu2_[k] + dt_ * rhs_hu_[k]);

                hv_[k] = (1.0 / 3.0) * hv_[k] + (2.0 / 3.0) * (hv2_[k] + dt_ * rhs_hv_[k]);
            }
            timing_.time_integration += FastOpenMPSolverTimingStats::seconds_since(start);
        }

        {
            const auto start = FastOpenMPSolverTimingStats::now();
            fast_hll_muscl_bathy_omp::enforce_positivity(gv_, h_.data(), hu_.data(), hv_.data());
            timing_.positivity_correction += FastOpenMPSolverTimingStats::seconds_since(start);
        }

        {
            const auto start = FastOpenMPSolverTimingStats::now();
            apply_reflecting_bc(h_.data(), hu_.data(), hv_.data());
            timing_.boundary_conditions += FastOpenMPSolverTimingStats::seconds_since(start);
        }

        (void)step;
    }

    void apply_reflecting_bc(double *h, double *hu, double *hv) const noexcept {
        const int nG = gv_.nG;
        const int Nx = gv_.Nx;
        const int Ny = gv_.Ny;
        const int s = gv_.stride;

#pragma omp parallel for schedule(static)
        for (int g = 0; g < nG; ++g) {
            const int iL = nG - 1 - g;
            const int iSrcL = nG + g;
            const int iR = nG + Nx + g;
            const int iSrcR = nG + Nx - 1 - g;

            for (int j = nG; j < nG + Ny; ++j) {
                const int kL = fast_hll_muscl_bathy_omp::idx(iL, j, s);
                const int kSrcL = fast_hll_muscl_bathy_omp::idx(iSrcL, j, s);

                h[kL] = h[kSrcL];
                hu[kL] = -hu[kSrcL];
                hv[kL] = hv[kSrcL];

                const int kR = fast_hll_muscl_bathy_omp::idx(iR, j, s);
                const int kSrcR = fast_hll_muscl_bathy_omp::idx(iSrcR, j, s);

                h[kR] = h[kSrcR];
                hu[kR] = -hu[kSrcR];
                hv[kR] = hv[kSrcR];
            }
        }

#pragma omp parallel for schedule(static)
        for (int g = 0; g < nG; ++g) {
            const int jB = nG - 1 - g;
            const int jSrcB = nG + g;
            const int jT = nG + Ny + g;
            const int jSrcT = nG + Ny - 1 - g;

            for (int i = 0; i < gv_.Nx_total; ++i) {
                const int kB = fast_hll_muscl_bathy_omp::idx(i, jB, s);
                const int kSrcB = fast_hll_muscl_bathy_omp::idx(i, jSrcB, s);

                h[kB] = h[kSrcB];
                hu[kB] = hu[kSrcB];
                hv[kB] = -hv[kSrcB];

                const int kT = fast_hll_muscl_bathy_omp::idx(i, jT, s);
                const int kSrcT = fast_hll_muscl_bathy_omp::idx(i, jSrcT, s);

                h[kT] = h[kSrcT];
                hu[kT] = hu[kSrcT];
                hv[kT] = -hv[kSrcT];
            }
        }
    }

    void apply_scalar_reflecting_bc(double *a) const noexcept {
        const int nG = gv_.nG;
        const int Nx = gv_.Nx;
        const int Ny = gv_.Ny;
        const int s = gv_.stride;

#pragma omp parallel for schedule(static)
        for (int g = 0; g < nG; ++g) {
            const int iL = nG - 1 - g;
            const int iSrcL = nG + g;
            const int iR = nG + Nx + g;
            const int iSrcR = nG + Nx - 1 - g;

            for (int j = nG; j < nG + Ny; ++j) {
                a[fast_hll_muscl_bathy_omp::idx(iL, j, s)] =
                    a[fast_hll_muscl_bathy_omp::idx(iSrcL, j, s)];

                a[fast_hll_muscl_bathy_omp::idx(iR, j, s)] =
                    a[fast_hll_muscl_bathy_omp::idx(iSrcR, j, s)];
            }
        }

#pragma omp parallel for schedule(static)
        for (int g = 0; g < nG; ++g) {
            const int jB = nG - 1 - g;
            const int jSrcB = nG + g;
            const int jT = nG + Ny + g;
            const int jSrcT = nG + Ny - 1 - g;

            for (int i = 0; i < gv_.Nx_total; ++i) {
                a[fast_hll_muscl_bathy_omp::idx(i, jB, s)] =
                    a[fast_hll_muscl_bathy_omp::idx(i, jSrcB, s)];

                a[fast_hll_muscl_bathy_omp::idx(i, jT, s)] =
                    a[fast_hll_muscl_bathy_omp::idx(i, jSrcT, s)];
            }
        }
    }

    static void copy_from_array2d(const Array2D &src, std::vector<double> &dst) {
        std::copy(src.data(), src.data() + src.size(), dst.begin());
    }

    static void copy_to_array2d(const std::vector<double> &src, Array2D &dst) {
        std::copy(src.begin(), src.end(), dst.data());
    }

    static void copy_from_state(const State &src, std::vector<double> &h, std::vector<double> &hu,
                                std::vector<double> &hv) {
        std::copy(src.h().data(), src.h().data() + src.h().size(), h.begin());
        std::copy(src.hu().data(), src.hu().data() + src.hu().size(), hu.begin());
        std::copy(src.hv().data(), src.hv().data() + src.hv().size(), hv.begin());
    }

    void fill_io_state() {
        double *dst_h = io_state_.h().data();
        double *dst_hu = io_state_.hu().data();
        double *dst_hv = io_state_.hv().data();

#pragma omp parallel for schedule(static)
        for (std::ptrdiff_t k = 0; k < static_cast<std::ptrdiff_t>(n_total_); ++k) {
            dst_h[k] = h_[k];
            dst_hu[k] = hu_[k];
            dst_hv[k] = hv_[k];
        }
    }

    void run_sanity_checks(double time, std::size_t step) {
        if (sanity_checks_.empty()) return;

        fill_io_state();

        for (auto &check : sanity_checks_) {
            check->evaluate(io_state_, grid_, time, step, sanity_writer_.get());
        }
    }

    static std::unique_ptr<SanityCheckNetCDFWriter> make_sanity_writer(const SimulationConfig &cfg,
                                                                       double dt) {
        if (!cfg.sanity_checks.debug) return nullptr;

        if (cfg.sanity_checks.output_path.empty()) {
            throw std::runtime_error(
                "sanity_checks.output_path must be set if debug mode is enabled");
        }

        return std::make_unique<SanityCheckNetCDFWriter>(
            cfg.sanity_checks.output_path.string(), cfg.time.time_unit, cfg.mesh.spatial_unit_h,
            cfg.time.save_every, dt, backend_name_from_cfg(cfg), "HLL", "MUSCL", "SSPRK3",
            "Reflecting Walls", bathymetry_name_from_cfg(cfg));
    }

  private:
    const SimulationConfig cfg_;

    Grid grid_;
    fast_hll_muscl_bathy_omp::GridView gv_;

    double dt_{};
    double end_time_{};
    std::size_t steps_{};
    std::size_t save_every_{};
    std::size_t n_total_{};

    std::vector<double> h_, hu_, hv_;
    std::vector<double> h1_, hu1_, hv1_;
    std::vector<double> h2_, hu2_, hv2_;
    std::vector<double> rhs_h_, rhs_hu_, rhs_hv_;
    std::vector<double> B_;

    std::vector<double> fxm_h_, fxm_hu_, fxm_hv_;
    std::vector<double> fxp_h_, fxp_hu_, fxp_hv_;
    std::vector<double> fym_h_, fym_hu_, fym_hv_;
    std::vector<double> fyp_h_, fyp_hu_, fyp_hv_;

    State io_state_;
    Array2D io_bathy_;
    NetCDFWriter writer_;
    std::vector<std::unique_ptr<SanityCheck>> sanity_checks_;
    std::unique_ptr<SanityCheckNetCDFWriter> sanity_writer_;

    FastOpenMPSolverTimingStats timing_;
};