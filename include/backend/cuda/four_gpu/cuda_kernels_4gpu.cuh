#pragma once

#include "include/constants.hpp"

#include <cstddef>

#ifndef USE_CUDA_4
#define USE_CUDA_4 0
#endif

#if USE_CUDA_4

#include <cuda_runtime.h>

#include <cmath>

namespace fast_hll_muscl_bathy_cuda_4 {

#if defined(__CUDACC__)
#define FAST_HD __host__ __device__
#define FAST_D __device__
#else
#define FAST_HD
#define FAST_D
#endif

#if defined(__GNUC__) || defined(__clang__)
#define FAST_RESTRICT __restrict__
#else
#define FAST_RESTRICT
#endif

struct GridView {
    int Nx{};
    int Ny{};
    int nG{};
    int Nx_total{};
    int Ny_total{};
    int stride{};
    double dx{};
    double dy{};
};

struct DeviceState {
    double *h{};
    double *hu{};
    double *hv{};
};

struct ConstDeviceState {
    const double *h{};
    const double *hu{};
    const double *hv{};
};

FAST_HD inline int idx(const int i, const int j, const int stride) noexcept {
    return i * stride + j;
}

FAST_D inline double positive_height(const double h) noexcept {
    return h > 0.0 ? h : 0.0;
}

FAST_D inline double safe_velocity(const double q, const double h) noexcept {
    return (h > constants::eps) ? (q / h) : 0.0;
}

FAST_D inline double minmod(const double a, const double b) noexcept {
    if (a * b <= 0.0) return 0.0;
    return copysign(fmin(fabs(a), fabs(b)), a);
}

FAST_D inline void hll_x(const double hL, const double huL, const double hvL,
                         const double hR, const double huR, const double hvR,
                         double &fh, double &fhu, double &fhv) noexcept {
    const double uL = safe_velocity(huL, hL);
    const double vL = safe_velocity(hvL, hL);
    const double uR = safe_velocity(huR, hR);
    const double vR = safe_velocity(hvR, hR);

    const double cL = sqrt(constants::g * hL);
    const double cR = sqrt(constants::g * hR);

    const double sL = fmin(uL - cL, uR - cR);
    const double sR = fmax(uL + cL, uR + cR);

    const double FL_h = huL;
    const double FL_hu = huL * uL + 0.5 * constants::g * hL * hL;
    const double FL_hv = huL * vL;

    const double FR_h = huR;
    const double FR_hu = huR * uR + 0.5 * constants::g * hR * hR;
    const double FR_hv = huR * vR;

    if (sL >= 0.0) {
        fh = FL_h;
        fhu = FL_hu;
        fhv = FL_hv;
    } else if (sR <= 0.0) {
        fh = FR_h;
        fhu = FR_hu;
        fhv = FR_hv;
    } else {
        const double inv = 1.0 / (sR - sL);
        fh = (sR * FL_h - sL * FR_h + sL * sR * (hR - hL)) * inv;
        fhu = (sR * FL_hu - sL * FR_hu + sL * sR * (huR - huL)) * inv;
        fhv = (sR * FL_hv - sL * FR_hv + sL * sR * (hvR - hvL)) * inv;
    }
}

FAST_D inline void hll_y(const double hB, const double huB, const double hvB,
                         const double hT, const double huT, const double hvT,
                         double &gh, double &ghu, double &ghv) noexcept {
    const double uB = safe_velocity(huB, hB);
    const double vB = safe_velocity(hvB, hB);
    const double uT = safe_velocity(huT, hT);
    const double vT = safe_velocity(hvT, hT);

    const double cB = sqrt(constants::g * hB);
    const double cT = sqrt(constants::g * hT);

    const double sB = fmin(vB - cB, vT - cT);
    const double sT = fmax(vB + cB, vT + cT);

    const double GB_h = hvB;
    const double GB_hu = hvB * uB;
    const double GB_hv = hvB * vB + 0.5 * constants::g * hB * hB;

    const double GT_h = hvT;
    const double GT_hu = hvT * uT;
    const double GT_hv = hvT * vT + 0.5 * constants::g * hT * hT;

    if (sB >= 0.0) {
        gh = GB_h;
        ghu = GB_hu;
        ghv = GB_hv;
    } else if (sT <= 0.0) {
        gh = GT_h;
        ghu = GT_hu;
        ghv = GT_hv;
    } else {
        const double inv = 1.0 / (sT - sB);
        gh = (sT * GB_h - sB * GT_h + sB * sT * (hT - hB)) * inv;
        ghu = (sT * GB_hu - sB * GT_hu + sB * sT * (huT - huB)) * inv;
        ghv = (sT * GB_hv - sB * GT_hv + sB * sT * (hvT - hvB)) * inv;
    }
}

FAST_D inline void compute_x_interface_flux(
    const GridView g,
    const double *FAST_RESTRICT h,
    const double *FAST_RESTRICT hu,
    const double *FAST_RESTRICT hv,
    const double *FAST_RESTRICT B,
    const int i,
    const int j,
    double &fm_h,
    double &fm_hu,
    double &fm_hv,
    double &fp_h,
    double &fp_hu,
    double &fp_hv) noexcept {

    const int s = g.stride;
    const int k = idx(i, j, s);
    const int kr = k + s;

    const double eta_im1 = h[k - s] + B[k - s];
    const double eta_i = h[k] + B[k];
    const double eta_ip1 = h[kr] + B[kr];
    const double eta_ip2 = h[kr + s] + B[kr + s];

    const double eta_minus = eta_i + 0.5 * minmod(eta_i - eta_im1, eta_ip1 - eta_i);
    const double eta_plus = eta_ip1 - 0.5 * minmod(eta_ip1 - eta_i, eta_ip2 - eta_ip1);

    const double hu_minus = hu[k] + 0.5 * minmod(hu[k] - hu[k - s], hu[kr] - hu[k]);
    const double hv_minus = hv[k] + 0.5 * minmod(hv[k] - hv[k - s], hv[kr] - hv[k]);
    const double hu_plus = hu[kr] - 0.5 * minmod(hu[kr] - hu[k], hu[kr + s] - hu[kr]);
    const double hv_plus = hv[kr] - 0.5 * minmod(hv[kr] - hv[k], hv[kr + s] - hv[kr]);

    const double b_minus = B[k];
    const double b_plus = B[kr];
    const double b_star = fmax(b_minus, b_plus);

    const double h_minus = positive_height(eta_minus - b_minus);
    const double h_plus = positive_height(eta_plus - b_plus);

    const double u_minus = safe_velocity(hu_minus, h_minus);
    const double v_minus = safe_velocity(hv_minus, h_minus);
    const double u_plus = safe_velocity(hu_plus, h_plus);
    const double v_plus = safe_velocity(hv_plus, h_plus);

    const double hs_minus = positive_height(eta_minus - b_star);
    const double hs_plus = positive_height(eta_plus - b_star);

    const double hL = hs_minus;
    const double huL = hs_minus * u_minus;
    const double hvL = hs_minus * v_minus;
    const double hR = hs_plus;
    const double huR = hs_plus * u_plus;
    const double hvR = hs_plus * v_plus;

    double fh{};
    double fhu{};
    double fhv{};

    hll_x(hL, huL, hvL, hR, huR, hvR, fh, fhu, fhv);

    const double corr_minus = 0.5 * constants::g * (h_minus * h_minus - hs_minus * hs_minus);
    const double corr_plus = 0.5 * constants::g * (h_plus * h_plus - hs_plus * hs_plus);

    fm_h = fh;
    fm_hu = fhu + corr_minus;
    fm_hv = fhv;

    fp_h = fh;
    fp_hu = fhu + corr_plus;
    fp_hv = fhv;
}

FAST_D inline void compute_y_interface_flux(
    const GridView g,
    const double *FAST_RESTRICT h,
    const double *FAST_RESTRICT hu,
    const double *FAST_RESTRICT hv,
    const double *FAST_RESTRICT B,
    const int i,
    const int j,
    double &fm_h,
    double &fm_hu,
    double &fm_hv,
    double &fp_h,
    double &fp_hu,
    double &fp_hv) noexcept {

    const int s = g.stride;
    const int k = idx(i, j, s);
    const int kt = k + 1;

    const double eta_jm1 = h[k - 1] + B[k - 1];
    const double eta_j = h[k] + B[k];
    const double eta_jp1 = h[kt] + B[kt];
    const double eta_jp2 = h[kt + 1] + B[kt + 1];

    const double eta_minus = eta_j + 0.5 * minmod(eta_j - eta_jm1, eta_jp1 - eta_j);
    const double eta_plus = eta_jp1 - 0.5 * minmod(eta_jp1 - eta_j, eta_jp2 - eta_jp1);

    const double hu_minus = hu[k] + 0.5 * minmod(hu[k] - hu[k - 1], hu[kt] - hu[k]);
    const double hv_minus = hv[k] + 0.5 * minmod(hv[k] - hv[k - 1], hv[kt] - hv[k]);
    const double hu_plus = hu[kt] - 0.5 * minmod(hu[kt] - hu[k], hu[kt + 1] - hu[kt]);
    const double hv_plus = hv[kt] - 0.5 * minmod(hv[kt] - hv[k], hv[kt + 1] - hv[kt]);

    const double b_minus = B[k];
    const double b_plus = B[kt];
    const double b_star = fmax(b_minus, b_plus);

    const double h_minus = positive_height(eta_minus - b_minus);
    const double h_plus = positive_height(eta_plus - b_plus);

    const double u_minus = safe_velocity(hu_minus, h_minus);
    const double v_minus = safe_velocity(hv_minus, h_minus);
    const double u_plus = safe_velocity(hu_plus, h_plus);
    const double v_plus = safe_velocity(hv_plus, h_plus);

    const double hs_minus = positive_height(eta_minus - b_star);
    const double hs_plus = positive_height(eta_plus - b_star);

    const double hB = hs_minus;
    const double huB = hs_minus * u_minus;
    const double hvB = hs_minus * v_minus;
    const double hT = hs_plus;
    const double huT = hs_plus * u_plus;
    const double hvT = hs_plus * v_plus;

    double gh{};
    double ghu{};
    double ghv{};

    hll_y(hB, huB, hvB, hT, huT, hvT, gh, ghu, ghv);

    const double corr_minus = 0.5 * constants::g * (h_minus * h_minus - hs_minus * hs_minus);
    const double corr_plus = 0.5 * constants::g * (h_plus * h_plus - hs_plus * hs_plus);

    fm_h = gh;
    fm_hu = ghu;
    fm_hv = ghv + corr_minus;

    fp_h = gh;
    fp_hu = ghu;
    fp_hv = ghv + corr_plus;
}

__global__ inline void apply_reflecting_x_kernel(
    const GridView g,
    double *FAST_RESTRICT h,
    double *FAST_RESTRICT hu,
    double *FAST_RESTRICT hv) {

    const int jj = blockIdx.x * blockDim.x + threadIdx.x;
    const int gg = blockIdx.y;

    if (jj >= g.Ny || gg >= g.nG) return;

    const int j = g.nG + jj;

    const int iL = g.nG - 1 - gg;
    const int iSrcL = g.nG + gg;

    const int iR = g.nG + g.Nx + gg;
    const int iSrcR = g.nG + g.Nx - 1 - gg;

    const int kL = idx(iL, j, g.stride);
    const int kSrcL = idx(iSrcL, j, g.stride);

    const int kR = idx(iR, j, g.stride);
    const int kSrcR = idx(iSrcR, j, g.stride);

    h[kL] = h[kSrcL];
    hu[kL] = -hu[kSrcL];
    hv[kL] = hv[kSrcL];

    h[kR] = h[kSrcR];
    hu[kR] = -hu[kSrcR];
    hv[kR] = hv[kSrcR];
}


__global__ inline void apply_physical_reflecting_x_kernel(
    const GridView g,
    double *FAST_RESTRICT h,
    double *FAST_RESTRICT hu,
    double *FAST_RESTRICT hv,
    const int reflect_left,
    const int reflect_right) {

    const int jj = blockIdx.x * blockDim.x + threadIdx.x;
    const int gg = blockIdx.y;

    if (jj >= g.Ny || gg >= g.nG) return;

    const int j = g.nG + jj;

    if (reflect_left) {
        const int iL = g.nG - 1 - gg;
        const int iSrcL = g.nG + gg;

        const int kL = idx(iL, j, g.stride);
        const int kSrcL = idx(iSrcL, j, g.stride);

        h[kL] = h[kSrcL];
        hu[kL] = -hu[kSrcL];
        hv[kL] = hv[kSrcL];
    }

    if (reflect_right) {
        const int iR = g.nG + g.Nx + gg;
        const int iSrcR = g.nG + g.Nx - 1 - gg;

        const int kR = idx(iR, j, g.stride);
        const int kSrcR = idx(iSrcR, j, g.stride);

        h[kR] = h[kSrcR];
        hu[kR] = -hu[kSrcR];
        hv[kR] = hv[kSrcR];
    }
}

__global__ inline void apply_reflecting_y_kernel(
    const GridView g,
    double *FAST_RESTRICT h,
    double *FAST_RESTRICT hu,
    double *FAST_RESTRICT hv) {

    const int ii = blockIdx.x * blockDim.x + threadIdx.x;
    const int gg = blockIdx.y;

    if (ii >= g.Nx_total || gg >= g.nG) return;

    const int i = ii;

    const int jB = g.nG - 1 - gg;
    const int jSrcB = g.nG + gg;

    const int jT = g.nG + g.Ny + gg;
    const int jSrcT = g.nG + g.Ny - 1 - gg;

    const int kB = idx(i, jB, g.stride);
    const int kSrcB = idx(i, jSrcB, g.stride);

    const int kT = idx(i, jT, g.stride);
    const int kSrcT = idx(i, jSrcT, g.stride);

    h[kB] = h[kSrcB];
    hu[kB] = hu[kSrcB];
    hv[kB] = -hv[kSrcB];

    h[kT] = h[kSrcT];
    hu[kT] = hu[kSrcT];
    hv[kT] = -hv[kSrcT];
}

__global__ inline void fused_rhs_kernel(
    const GridView g,
    const double *FAST_RESTRICT h,
    const double *FAST_RESTRICT hu,
    const double *FAST_RESTRICT hv,
    const double *FAST_RESTRICT B,
    double *FAST_RESTRICT rhs_h,
    double *FAST_RESTRICT rhs_hu,
    double *FAST_RESTRICT rhs_hv) {

    const int jj = blockIdx.x * blockDim.x + threadIdx.x;
    const int ii = blockIdx.y * blockDim.y + threadIdx.y;

    if (ii >= g.Nx || jj >= g.Ny) return;

    const int i = g.nG + ii;
    const int j = g.nG + jj;
    const int k = idx(i, j, g.stride);

    double xr_m_h{}, xr_m_hu{}, xr_m_hv{}, xr_p_h{}, xr_p_hu{}, xr_p_hv{};
    double xl_m_h{}, xl_m_hu{}, xl_m_hv{}, xl_p_h{}, xl_p_hu{}, xl_p_hv{};
    double yt_m_h{}, yt_m_hu{}, yt_m_hv{}, yt_p_h{}, yt_p_hu{}, yt_p_hv{};
    double yb_m_h{}, yb_m_hu{}, yb_m_hv{}, yb_p_h{}, yb_p_hu{}, yb_p_hv{};

    compute_x_interface_flux(g, h, hu, hv, B, i, j,
                             xr_m_h, xr_m_hu, xr_m_hv,
                             xr_p_h, xr_p_hu, xr_p_hv);

    compute_x_interface_flux(g, h, hu, hv, B, i - 1, j,
                             xl_m_h, xl_m_hu, xl_m_hv,
                             xl_p_h, xl_p_hu, xl_p_hv);

    compute_y_interface_flux(g, h, hu, hv, B, i, j,
                             yt_m_h, yt_m_hu, yt_m_hv,
                             yt_p_h, yt_p_hu, yt_p_hv);

    compute_y_interface_flux(g, h, hu, hv, B, i, j - 1,
                             yb_m_h, yb_m_hu, yb_m_hv,
                             yb_p_h, yb_p_hu, yb_p_hv);

    const double inv_dx = 1.0 / g.dx;
    const double inv_dy = 1.0 / g.dy;

    rhs_h[k] = -inv_dx * (xr_m_h - xl_p_h) - inv_dy * (yt_m_h - yb_p_h);
    rhs_hu[k] = -inv_dx * (xr_m_hu - xl_p_hu) - inv_dy * (yt_m_hu - yb_p_hu);
    rhs_hv[k] = -inv_dx * (xr_m_hv - xl_p_hv) - inv_dy * (yt_m_hv - yb_p_hv);
}

__global__ inline void rk_stage1_kernel(
    const std::size_t n,
    const double dt,
    const double *FAST_RESTRICT h,
    const double *FAST_RESTRICT hu,
    const double *FAST_RESTRICT hv,
    const double *FAST_RESTRICT rhs_h,
    const double *FAST_RESTRICT rhs_hu,
    const double *FAST_RESTRICT rhs_hv,
    double *FAST_RESTRICT h1,
    double *FAST_RESTRICT hu1,
    double *FAST_RESTRICT hv1) {

    const std::size_t k =
        static_cast<std::size_t>(blockIdx.x) * static_cast<std::size_t>(blockDim.x) +
        static_cast<std::size_t>(threadIdx.x);

    if (k >= n) return;

    h1[k] = h[k] + dt * rhs_h[k];
    hu1[k] = hu[k] + dt * rhs_hu[k];
    hv1[k] = hv[k] + dt * rhs_hv[k];
}

__global__ inline void rk_stage2_kernel(
    const std::size_t n,
    const double dt,
    const double *FAST_RESTRICT h0,
    const double *FAST_RESTRICT hu0,
    const double *FAST_RESTRICT hv0,
    const double *FAST_RESTRICT h1,
    const double *FAST_RESTRICT hu1,
    const double *FAST_RESTRICT hv1,
    const double *FAST_RESTRICT rhs_h,
    const double *FAST_RESTRICT rhs_hu,
    const double *FAST_RESTRICT rhs_hv,
    double *FAST_RESTRICT h2,
    double *FAST_RESTRICT hu2,
    double *FAST_RESTRICT hv2) {

    const std::size_t k =
        static_cast<std::size_t>(blockIdx.x) * static_cast<std::size_t>(blockDim.x) +
        static_cast<std::size_t>(threadIdx.x);

    if (k >= n) return;

    h2[k] = 0.75 * h0[k] + 0.25 * (h1[k] + dt * rhs_h[k]);
    hu2[k] = 0.75 * hu0[k] + 0.25 * (hu1[k] + dt * rhs_hu[k]);
    hv2[k] = 0.75 * hv0[k] + 0.25 * (hv1[k] + dt * rhs_hv[k]);
}

__global__ inline void rk_stage3_kernel(
    const std::size_t n,
    const double dt,
    double *FAST_RESTRICT h0,
    double *FAST_RESTRICT hu0,
    double *FAST_RESTRICT hv0,
    const double *FAST_RESTRICT h2,
    const double *FAST_RESTRICT hu2,
    const double *FAST_RESTRICT hv2,
    const double *FAST_RESTRICT rhs_h,
    const double *FAST_RESTRICT rhs_hu,
    const double *FAST_RESTRICT rhs_hv) {

    const std::size_t k =
        static_cast<std::size_t>(blockIdx.x) * static_cast<std::size_t>(blockDim.x) +
        static_cast<std::size_t>(threadIdx.x);

    if (k >= n) return;

    h0[k] = (1.0 / 3.0) * h0[k] + (2.0 / 3.0) * (h2[k] + dt * rhs_h[k]);
    hu0[k] = (1.0 / 3.0) * hu0[k] + (2.0 / 3.0) * (hu2[k] + dt * rhs_hu[k]);
    hv0[k] = (1.0 / 3.0) * hv0[k] + (2.0 / 3.0) * (hv2[k] + dt * rhs_hv[k]);
}

__global__ inline void enforce_positivity_kernel(
    const GridView g,
    double *FAST_RESTRICT h,
    double *FAST_RESTRICT hu,
    double *FAST_RESTRICT hv,
    int *FAST_RESTRICT nonfinite,
    const double h_floor = 1.0e-8) {

    const int jj = blockIdx.x * blockDim.x + threadIdx.x;
    const int ii = blockIdx.y * blockDim.y + threadIdx.y;

    if (ii >= g.Nx || jj >= g.Ny) return;

    const int i = g.nG + ii;
    const int j = g.nG + jj;
    const int k = idx(i, j, g.stride);

    if (h[k] < h_floor) {
        h[k] = h_floor;
        hu[k] = 0.0;
        hv[k] = 0.0;
    }

    if (!isfinite(h[k]) || !isfinite(hu[k]) || !isfinite(hv[k])) {
        atomicExch(nonfinite, 1);
    }
}

inline dim3 rhs_block() noexcept {
    return dim3(16, 16);
}

inline dim3 rhs_grid(const GridView &g) noexcept {
    const dim3 block = rhs_block();
    return dim3(
        static_cast<unsigned int>((g.Ny + block.x - 1) / block.x),
        static_cast<unsigned int>((g.Nx + block.y - 1) / block.y));
}

inline dim3 one_d_grid(const std::size_t n, const int block = 256) noexcept {
    return dim3(
        static_cast<unsigned int>(
            (n + static_cast<std::size_t>(block) - 1) /
            static_cast<std::size_t>(block)));
}

#undef FAST_HD
#undef FAST_D
#undef FAST_RESTRICT

} // namespace fast_hll_muscl_bathy_cuda_4

#else  // USE_CUDA_4

namespace fast_hll_muscl_bathy_cuda_4 {

// Intentionally empty in non-CUDA builds.
//
// This makes the header safe to include when USE_CUDA_4=0.

} // namespace fast_hll_muscl_bathy_cuda_4

#endif // USE_CUDA_4