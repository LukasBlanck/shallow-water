#pragma once

#ifndef USE_CUDA
#define USE_CUDA 0
#endif

#if USE_CUDA

#include "configs/config.hpp"
#include "configs/parse_config_helper.hpp"

#include "include/backend/cuda/cuda_kernels.cuh"
#include "include/backend/openmp/omp_kernels.hpp"
#include "include/bathymetry/apply_bathymetry.hpp"
#include "include/core/array2D.hpp"
#include "include/core/grid.hpp"
#include "include/core/state.hpp"
#include "include/initial_condition/apply_initial_condition.hpp"
#include "include/io/netCDF_writer.hpp"
#include "include/io/sanity_checks_netdcdf_writer.hpp"
#include "include/solver_assembly/time_utils.hpp"
#include "tests/sanity_checks/sanity_checks.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#define CUDA_SOLVER_CHECK(call)                                                                    \
    do {                                                                                           \
        const cudaError_t err__ = (call);                                                          \
        if (err__ != cudaSuccess) {                                                                \
            throw std::runtime_error(std::string("CUDA error in ") + __FILE__ + ":" +              \
                                     std::to_string(__LINE__) + " after " #call ": " +             \
                                     cudaGetErrorString(err__));                                   \
        }                                                                                          \
    } while (false)

struct FastCUDASolverTimingStats {
    using Clock = std::chrono::steady_clock;

    double cfl_check = 0.0;
    double host_to_device = 0.0;
    double device_to_host = 0.0;
    double boundary_conditions = 0.0;
    double fused_rhs = 0.0;
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
        return host_to_device + device_to_host + boundary_conditions + fused_rhs +
               time_integration + positivity_correction + output + sanity_checks;
    }

    void print_summary() const {
        const double total = measured_solver_time();

        std::cout << "\n========== Fused CUDA solver timing summary ==========\n";
        std::cout << "time steps           : " << time_steps << '\n';
        std::cout << "rhs calls            : " << rhs_calls << '\n';
        std::cout << "CFL check            : " << cfl_check << " s\n";
        std::cout << "measured solver time : " << total << " s\n";

        if (total <= 0.0) {
            std::cout << "No positive measured solver time.\n";
            std::cout << "====================================================\n";
            return;
        }

        const auto print_row = [total](const char *name, const double seconds) {
            const double percent = 100.0 * seconds / total;
            std::cout << std::left << std::setw(24) << name << std::right << std::setw(12)
                      << seconds << " s  " << std::setw(8) << percent << " %\n";
        };

        std::cout << std::fixed << std::setprecision(6);
        print_row("host to device", host_to_device);
        print_row("device to host", device_to_host);
        print_row("boundary conditions", boundary_conditions);
        print_row("fused rhs", fused_rhs);
        print_row("time integration", time_integration);
        print_row("positivity correction", positivity_correction);
        print_row("output", output);
        print_row("sanity checks", sanity_checks);
        std::cout << "====================================================\n";
    }
};

class FastHLLMUSCLBathyCUDASolver {
  public:
    explicit FastHLLMUSCLBathyCUDASolver(const SimulationConfig &cfg)
        : cfg_(cfg), grid_(cfg.mesh.Nx, cfg.mesh.Ny, cfg.mesh.Lx, cfg.mesh.Ly, cfg.mesh.nG),
          gv_{grid_.Nx(),       grid_.Ny(),       grid_.nG(), grid_.Nx_total(),
              grid_.Ny_total(), grid_.Ny_total(), grid_.dx(), grid_.dy()},
          cgv_{grid_.Nx(),       grid_.Ny(),       grid_.nG(), grid_.Nx_total(),
               grid_.Ny_total(), grid_.Ny_total(), grid_.dx(), grid_.dy()},
          dt_(cfg.time.end_time / static_cast<double>(cfg.time.time_steps)),
          end_time_(cfg.time.end_time), steps_(cfg.time.time_steps),
          save_every_(static_cast<std::size_t>(cfg.time.save_every)),
          n_total_(static_cast<std::size_t>(gv_.Nx_total) * static_cast<std::size_t>(gv_.Ny_total)),
          h_(n_total_), hu_(n_total_), hv_(n_total_), B_(n_total_), io_state_(grid_),
          io_bathy_(grid_.Nx_total(), grid_.Ny_total()),
          writer_(cfg.output.path.string(), grid_, cfg.mesh.spatial_unit_x, cfg.mesh.spatial_unit_y,
                  cfg.mesh.spatial_unit_h, cfg.time.time_unit, cfg.time.save_every),
          sanity_checks_(make_sanity_checks(cfg)), sanity_writer_(make_sanity_writer(cfg, dt_)) {

        if (grid_.nG() < 2) {
            throw std::runtime_error("FastHLLMUSCLBathyCUDASolver requires nG >= 2 for MUSCL.");
        }

        allocate_device();

        Array2D B_tmp(grid_.Nx_total(), grid_.Ny_total());
        apply_bathymetry(cfg_, grid_, B_tmp);
        apply_scalar_reflecting_bc(B_tmp.data());
        copy_from_array2d(B_tmp, B_);
        copy_to_array2d(B_, io_bathy_);
        writer_.write_bathymetry(io_bathy_);

        State U_tmp(grid_);
        apply_initial_condition(cfg_, grid_, U_tmp);
        copy_from_state(U_tmp, h_, hu_, hv_);
        apply_reflecting_bc_cpu(h_.data(), hu_.data(), hv_.data());

        for (auto &check : sanity_checks_) {
            fill_io_state_from_host();
            check->initialize(io_state_, grid_);
        }
    }

    ~FastHLLMUSCLBathyCUDASolver() { free_device_noexcept(); }

    FastHLLMUSCLBathyCUDASolver(const FastHLLMUSCLBathyCUDASolver &) = delete;
    FastHLLMUSCLBathyCUDASolver &operator=(const FastHLLMUSCLBathyCUDASolver &) = delete;

    void run() {
        std::cout << std::unitbuf;
        std::cerr << std::unitbuf;
        std::cout << "INFO: CUDA fused RHS solver active\n";

        double dt_stable = 0.0;
        {
            const auto start = FastCUDASolverTimingStats::now();
            dt_stable = fast_hll_muscl_bathy_omp::compute_stable_dt(gv_, h_.data(), hu_.data(),
                                                                    hv_.data(), cfg_.time.cfl);
            timing_.cfl_check += FastCUDASolverTimingStats::seconds_since(start);
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
            const auto start = FastCUDASolverTimingStats::now();
            fill_io_state_from_host();
            writer_.write_snapshot(io_state_, time, dt_, backend_name_from_cfg(cfg_), "HLL",
                                   "MUSCL", "SSPRK3", "Reflecting Walls",
                                   bathymetry_name_from_cfg(cfg_));
            timing_.output += FastCUDASolverTimingStats::seconds_since(start);
        }

        {
            const auto start = FastCUDASolverTimingStats::now();
            run_sanity_checks(time, 0);
            timing_.sanity_checks += FastCUDASolverTimingStats::seconds_since(start);
        }

        upload_initial_state_to_device();

        for (std::size_t step = 0; step < steps_; ++step) {
            rk3_step_cuda();
            time += dt_;

            const std::size_t step_number = step + 1;
            ++timing_.time_steps;

            if (estimate_eta && !eta_printed && step_number == eta_probe_step) {
                CUDA_SOLVER_CHECK(cudaDeviceSynchronize());
                const auto now = std::chrono::steady_clock::now();
                const double elapsed_sec = std::chrono::duration<double>(now - start_wall).count();
                const double eta_sec = estimate_eta_seconds(step_number, steps_, elapsed_sec);
                std::cout << "ETA = " << format_duration(eta_sec) << '\n';
                eta_printed = true;
            }

            if (step_number % save_every_ == 0 || step_number == steps_) {
                download_state_from_device();

                {
                    const auto start = FastCUDASolverTimingStats::now();
                    fill_io_state_from_host();
                    writer_.write_snapshot(io_state_, time);
                    timing_.output += FastCUDASolverTimingStats::seconds_since(start);
                }

                {
                    const auto start = FastCUDASolverTimingStats::now();
                    run_sanity_checks(time, step_number);
                    timing_.sanity_checks += FastCUDASolverTimingStats::seconds_since(start);
                }
            }
        }

        timing_.print_summary();
    }

  private:
    void allocate_device() {
        CUDA_SOLVER_CHECK(cudaMalloc(&d_h_, nbytes()));
        CUDA_SOLVER_CHECK(cudaMalloc(&d_hu_, nbytes()));
        CUDA_SOLVER_CHECK(cudaMalloc(&d_hv_, nbytes()));

        CUDA_SOLVER_CHECK(cudaMalloc(&d_h1_, nbytes()));
        CUDA_SOLVER_CHECK(cudaMalloc(&d_hu1_, nbytes()));
        CUDA_SOLVER_CHECK(cudaMalloc(&d_hv1_, nbytes()));

        CUDA_SOLVER_CHECK(cudaMalloc(&d_h2_, nbytes()));
        CUDA_SOLVER_CHECK(cudaMalloc(&d_hu2_, nbytes()));
        CUDA_SOLVER_CHECK(cudaMalloc(&d_hv2_, nbytes()));

        CUDA_SOLVER_CHECK(cudaMalloc(&d_rhs_h_, nbytes()));
        CUDA_SOLVER_CHECK(cudaMalloc(&d_rhs_hu_, nbytes()));
        CUDA_SOLVER_CHECK(cudaMalloc(&d_rhs_hv_, nbytes()));

        CUDA_SOLVER_CHECK(cudaMalloc(&d_B_, nbytes()));
        CUDA_SOLVER_CHECK(cudaMalloc(&d_nonfinite_, sizeof(int)));
    }

    void free_device_noexcept() noexcept {
        cudaFree(d_h_);
        cudaFree(d_hu_);
        cudaFree(d_hv_);

        cudaFree(d_h1_);
        cudaFree(d_hu1_);
        cudaFree(d_hv1_);

        cudaFree(d_h2_);
        cudaFree(d_hu2_);
        cudaFree(d_hv2_);

        cudaFree(d_rhs_h_);
        cudaFree(d_rhs_hu_);
        cudaFree(d_rhs_hv_);

        cudaFree(d_B_);
        cudaFree(d_nonfinite_);
    }

    std::size_t nbytes() const noexcept { return n_total_ * sizeof(double); }

    void upload_initial_state_to_device() {
        const auto start = FastCUDASolverTimingStats::now();

        CUDA_SOLVER_CHECK(cudaMemcpy(d_h_, h_.data(), nbytes(), cudaMemcpyHostToDevice));
        CUDA_SOLVER_CHECK(cudaMemcpy(d_hu_, hu_.data(), nbytes(), cudaMemcpyHostToDevice));
        CUDA_SOLVER_CHECK(cudaMemcpy(d_hv_, hv_.data(), nbytes(), cudaMemcpyHostToDevice));
        CUDA_SOLVER_CHECK(cudaMemcpy(d_B_, B_.data(), nbytes(), cudaMemcpyHostToDevice));

        timing_.host_to_device += FastCUDASolverTimingStats::seconds_since(start);
    }

    void download_state_from_device() {
        const auto start = FastCUDASolverTimingStats::now();

        CUDA_SOLVER_CHECK(cudaMemcpy(h_.data(), d_h_, nbytes(), cudaMemcpyDeviceToHost));
        CUDA_SOLVER_CHECK(cudaMemcpy(hu_.data(), d_hu_, nbytes(), cudaMemcpyDeviceToHost));
        CUDA_SOLVER_CHECK(cudaMemcpy(hv_.data(), d_hv_, nbytes(), cudaMemcpyDeviceToHost));

        timing_.device_to_host += FastCUDASolverTimingStats::seconds_since(start);
    }

    void launch_reflecting_bc(double *h, double *hu, double *hv) {
        const auto start = FastCUDASolverTimingStats::now();

        constexpr int block_1d = 256;
        const dim3 grid_x((gv_.Ny + block_1d - 1) / block_1d, gv_.nG);
        const dim3 grid_y((gv_.Nx_total + block_1d - 1) / block_1d, gv_.nG);

        fast_hll_muscl_bathy_cuda::apply_reflecting_x_kernel<<<grid_x, block_1d>>>(cgv_, h, hu, hv);
        CUDA_SOLVER_CHECK(cudaGetLastError());

        fast_hll_muscl_bathy_cuda::apply_reflecting_y_kernel<<<grid_y, block_1d>>>(cgv_, h, hu, hv);
        CUDA_SOLVER_CHECK(cudaGetLastError());

        CUDA_SOLVER_CHECK(cudaDeviceSynchronize());
        timing_.boundary_conditions += FastCUDASolverTimingStats::seconds_since(start);
    }

    void launch_fused_rhs(const double *h, const double *hu, const double *hv, double *rhs_h,
                          double *rhs_hu, double *rhs_hv) {
        ++timing_.rhs_calls;
        const auto start = FastCUDASolverTimingStats::now();

        const dim3 block = fast_hll_muscl_bathy_cuda::rhs_block();
        const dim3 grid = fast_hll_muscl_bathy_cuda::rhs_grid(cgv_);

        fast_hll_muscl_bathy_cuda::fused_rhs_kernel<<<grid, block>>>(cgv_, h, hu, hv, d_B_, rhs_h,
                                                                     rhs_hu, rhs_hv);
        CUDA_SOLVER_CHECK(cudaGetLastError());

        CUDA_SOLVER_CHECK(cudaDeviceSynchronize());
        timing_.fused_rhs += FastCUDASolverTimingStats::seconds_since(start);
    }

    void launch_rk_stage1() {
        const auto start = FastCUDASolverTimingStats::now();

        constexpr int block = 256;
        const dim3 grid = fast_hll_muscl_bathy_cuda::one_d_grid(n_total_, block);

        fast_hll_muscl_bathy_cuda::rk_stage1_kernel<<<grid, block>>>(
            n_total_, dt_, d_h_, d_hu_, d_hv_, d_rhs_h_, d_rhs_hu_, d_rhs_hv_, d_h1_, d_hu1_,
            d_hv1_);
        CUDA_SOLVER_CHECK(cudaGetLastError());

        CUDA_SOLVER_CHECK(cudaDeviceSynchronize());
        timing_.time_integration += FastCUDASolverTimingStats::seconds_since(start);
    }

    void launch_rk_stage2() {
        const auto start = FastCUDASolverTimingStats::now();

        constexpr int block = 256;
        const dim3 grid = fast_hll_muscl_bathy_cuda::one_d_grid(n_total_, block);

        fast_hll_muscl_bathy_cuda::rk_stage2_kernel<<<grid, block>>>(
            n_total_, dt_, d_h_, d_hu_, d_hv_, d_h1_, d_hu1_, d_hv1_, d_rhs_h_, d_rhs_hu_,
            d_rhs_hv_, d_h2_, d_hu2_, d_hv2_);
        CUDA_SOLVER_CHECK(cudaGetLastError());

        CUDA_SOLVER_CHECK(cudaDeviceSynchronize());
        timing_.time_integration += FastCUDASolverTimingStats::seconds_since(start);
    }

    void launch_rk_stage3() {
        const auto start = FastCUDASolverTimingStats::now();

        constexpr int block = 256;
        const dim3 grid = fast_hll_muscl_bathy_cuda::one_d_grid(n_total_, block);

        fast_hll_muscl_bathy_cuda::rk_stage3_kernel<<<grid, block>>>(
            n_total_, dt_, d_h_, d_hu_, d_hv_, d_h2_, d_hu2_, d_hv2_, d_rhs_h_, d_rhs_hu_,
            d_rhs_hv_);
        CUDA_SOLVER_CHECK(cudaGetLastError());

        CUDA_SOLVER_CHECK(cudaDeviceSynchronize());
        timing_.time_integration += FastCUDASolverTimingStats::seconds_since(start);
    }

    void launch_enforce_positivity(double *h, double *hu, double *hv) {
        const auto start = FastCUDASolverTimingStats::now();

        CUDA_SOLVER_CHECK(cudaMemset(d_nonfinite_, 0, sizeof(int)));

        const dim3 block = fast_hll_muscl_bathy_cuda::rhs_block();
        const dim3 grid = fast_hll_muscl_bathy_cuda::rhs_grid(cgv_);

        fast_hll_muscl_bathy_cuda::enforce_positivity_kernel<<<grid, block>>>(cgv_, h, hu, hv,
                                                                              d_nonfinite_);
        CUDA_SOLVER_CHECK(cudaGetLastError());

        int nonfinite = 0;
        CUDA_SOLVER_CHECK(
            cudaMemcpy(&nonfinite, d_nonfinite_, sizeof(int), cudaMemcpyDeviceToHost));
        CUDA_SOLVER_CHECK(cudaDeviceSynchronize());

        timing_.positivity_correction += FastCUDASolverTimingStats::seconds_since(start);

        if (nonfinite) {
            throw std::runtime_error("non-finite state detected in CUDA fast solver");
        }
    }

    void rk3_step_cuda() {
        launch_reflecting_bc(d_h_, d_hu_, d_hv_);
        launch_fused_rhs(d_h_, d_hu_, d_hv_, d_rhs_h_, d_rhs_hu_, d_rhs_hv_);
        launch_rk_stage1();
        launch_enforce_positivity(d_h1_, d_hu1_, d_hv1_);

        launch_reflecting_bc(d_h1_, d_hu1_, d_hv1_);
        launch_fused_rhs(d_h1_, d_hu1_, d_hv1_, d_rhs_h_, d_rhs_hu_, d_rhs_hv_);
        launch_rk_stage2();
        launch_enforce_positivity(d_h2_, d_hu2_, d_hv2_);

        launch_reflecting_bc(d_h2_, d_hu2_, d_hv2_);
        launch_fused_rhs(d_h2_, d_hu2_, d_hv2_, d_rhs_h_, d_rhs_hu_, d_rhs_hv_);
        launch_rk_stage3();
        launch_enforce_positivity(d_h_, d_hu_, d_hv_);

        launch_reflecting_bc(d_h_, d_hu_, d_hv_);
    }

    void apply_reflecting_bc_cpu(double *h, double *hu, double *hv) const noexcept {
        const int nG = gv_.nG;
        const int Nx = gv_.Nx;
        const int Ny = gv_.Ny;
        const int s = gv_.stride;

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

    void fill_io_state_from_host() {
        std::copy(h_.begin(), h_.end(), io_state_.h().data());
        std::copy(hu_.begin(), hu_.end(), io_state_.hu().data());
        std::copy(hv_.begin(), hv_.end(), io_state_.hv().data());
    }

    void run_sanity_checks(double time, std::size_t step) {
        if (sanity_checks_.empty()) return;

        fill_io_state_from_host();

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
    fast_hll_muscl_bathy_cuda::GridView cgv_;

    double dt_{};
    double end_time_{};
    std::size_t steps_{};
    std::size_t save_every_{};
    std::size_t n_total_{};

    std::vector<double> h_;
    std::vector<double> hu_;
    std::vector<double> hv_;
    std::vector<double> B_;

    double *d_h_{};
    double *d_hu_{};
    double *d_hv_{};

    double *d_h1_{};
    double *d_hu1_{};
    double *d_hv1_{};

    double *d_h2_{};
    double *d_hu2_{};
    double *d_hv2_{};

    double *d_rhs_h_{};
    double *d_rhs_hu_{};
    double *d_rhs_hv_{};

    double *d_B_{};
    int *d_nonfinite_{};

    State io_state_;
    Array2D io_bathy_;
    NetCDFWriter writer_;
    std::vector<std::unique_ptr<SanityCheck>> sanity_checks_;
    std::unique_ptr<SanityCheckNetCDFWriter> sanity_writer_;

    FastCUDASolverTimingStats timing_;
};

#undef CUDA_SOLVER_CHECK

#endif // USE_CUDA