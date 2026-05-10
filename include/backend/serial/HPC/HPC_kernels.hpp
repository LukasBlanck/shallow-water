#pragma once

#include "include/constants.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace fast_hll_muscl_bathy {

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

inline int idx(const int i, const int j, const int stride) noexcept {
    return i * stride + j;
}

inline double positive_height(const double h) noexcept {
    return h > 0.0 ? h : 0.0;
}

inline double safe_velocity(const double q, const double h) noexcept {
    return (h > constants::eps) ? (q / h) : 0.0;
}

inline double minmod(const double a, const double b) noexcept {
    if (a * b <= 0.0) return 0.0;
    return (a > 0.0 ? 1.0 : -1.0) * std::min(std::abs(a), std::abs(b));
}

inline void hll_x(const double hL, const double huL, const double hvL, const double hR,
                  const double huR, const double hvR, double &fh, double &fhu,
                  double &fhv) noexcept {
    const double uL = safe_velocity(huL, hL);
    const double vL = safe_velocity(hvL, hL);
    const double uR = safe_velocity(huR, hR);
    const double vR = safe_velocity(hvR, hR);

    const double cL = std::sqrt(constants::g * hL);
    const double cR = std::sqrt(constants::g * hR);

    const double sL = std::min(uL - cL, uR - cR);
    const double sR = std::max(uL + cL, uR + cR);

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

inline void hll_y(const double hB, const double huB, const double hvB, const double hT,
                  const double huT, const double hvT, double &gh, double &ghu,
                  double &ghv) noexcept {
    const double uB = safe_velocity(huB, hB);
    const double vB = safe_velocity(hvB, hB);
    const double uT = safe_velocity(huT, hT);
    const double vT = safe_velocity(hvT, hT);

    const double cB = std::sqrt(constants::g * hB);
    const double cT = std::sqrt(constants::g * hT);

    const double sB = std::min(vB - cB, vT - cT);
    const double sT = std::max(vB + cB, vT + cT);

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

inline void compute_x_fluxes(const GridView &g, const double *FAST_RESTRICT h,
                             const double *FAST_RESTRICT hu, const double *FAST_RESTRICT hv,
                             const double *FAST_RESTRICT B, double *FAST_RESTRICT fxm_h,
                             double *FAST_RESTRICT fxm_hu, double *FAST_RESTRICT fxm_hv,
                             double *FAST_RESTRICT fxp_h, double *FAST_RESTRICT fxp_hu,
                             double *FAST_RESTRICT fxp_hv) {
    const int s = g.stride;
    const int nG = g.nG;

    for (int i = nG - 1; i < nG + g.Nx; ++i) {
        for (int j = nG; j < nG + g.Ny; ++j) {
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
            const double b_star = std::max(b_minus, b_plus);

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

            double fh, fhu, fhv;
            hll_x(hL, huL, hvL, hR, huR, hvR, fh, fhu, fhv);

            const double corr_minus =
                0.5 * constants::g * (h_minus * h_minus - hs_minus * hs_minus);
            const double corr_plus = 0.5 * constants::g * (h_plus * h_plus - hs_plus * hs_plus);

            fxm_h[k] = fh;
            fxm_hu[k] = fhu + corr_minus;
            fxm_hv[k] = fhv;

            fxp_h[k] = fh;
            fxp_hu[k] = fhu + corr_plus;
            fxp_hv[k] = fhv;
        }
    }
}

inline void compute_y_fluxes(const GridView &g, const double *FAST_RESTRICT h,
                             const double *FAST_RESTRICT hu, const double *FAST_RESTRICT hv,
                             const double *FAST_RESTRICT B, double *FAST_RESTRICT fym_h,
                             double *FAST_RESTRICT fym_hu, double *FAST_RESTRICT fym_hv,
                             double *FAST_RESTRICT fyp_h, double *FAST_RESTRICT fyp_hu,
                             double *FAST_RESTRICT fyp_hv) {
    const int s = g.stride;
    const int nG = g.nG;

    for (int i = nG; i < nG + g.Nx; ++i) {
        for (int j = nG - 1; j < nG + g.Ny; ++j) {
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
            const double b_star = std::max(b_minus, b_plus);

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

            double gh, ghu, ghv;
            hll_y(hB, huB, hvB, hT, huT, hvT, gh, ghu, ghv);

            const double corr_minus =
                0.5 * constants::g * (h_minus * h_minus - hs_minus * hs_minus);
            const double corr_plus = 0.5 * constants::g * (h_plus * h_plus - hs_plus * hs_plus);

            fym_h[k] = gh;
            fym_hu[k] = ghu;
            fym_hv[k] = ghv + corr_minus;

            fyp_h[k] = gh;
            fyp_hu[k] = ghu;
            fyp_hv[k] = ghv + corr_plus;
        }
    }
}

inline void apply_divergence(const GridView &g, const double *FAST_RESTRICT fxm_h,
                             const double *FAST_RESTRICT fxm_hu, const double *FAST_RESTRICT fxm_hv,
                             const double *FAST_RESTRICT fxp_h, const double *FAST_RESTRICT fxp_hu,
                             const double *FAST_RESTRICT fxp_hv, const double *FAST_RESTRICT fym_h,
                             const double *FAST_RESTRICT fym_hu, const double *FAST_RESTRICT fym_hv,
                             const double *FAST_RESTRICT fyp_h, const double *FAST_RESTRICT fyp_hu,
                             const double *FAST_RESTRICT fyp_hv, double *FAST_RESTRICT rhs_h,
                             double *FAST_RESTRICT rhs_hu, double *FAST_RESTRICT rhs_hv) noexcept {
    const int s = g.stride;
    const double inv_dx = 1.0 / g.dx;
    const double inv_dy = 1.0 / g.dy;

    for (int i = g.nG; i < g.nG + g.Nx; ++i) {
        for (int j = g.nG; j < g.nG + g.Ny; ++j) {
            const int k = idx(i, j, s);
            const int kmx = k - s;
            const int kmy = k - 1;

            rhs_h[k] = -inv_dx * (fxm_h[k] - fxp_h[kmx]) - inv_dy * (fym_h[k] - fyp_h[kmy]);
            rhs_hu[k] = -inv_dx * (fxm_hu[k] - fxp_hu[kmx]) - inv_dy * (fym_hu[k] - fyp_hu[kmy]);
            rhs_hv[k] = -inv_dx * (fxm_hv[k] - fxp_hv[kmx]) - inv_dy * (fym_hv[k] - fyp_hv[kmy]);
        }
    }
}

inline double compute_stable_dt(const GridView &g, const double *FAST_RESTRICT h,
                                const double *FAST_RESTRICT hu, const double *FAST_RESTRICT hv,
                                const double cfl) noexcept {
    double max_speed = 0.0;
    for (int i = g.nG; i < g.nG + g.Nx; ++i) {
        for (int j = g.nG; j < g.nG + g.Ny; ++j) {
            const int k = idx(i, j, g.stride);
            const double hk = h[k];
            const double u = safe_velocity(hu[k], hk);
            const double v = safe_velocity(hv[k], hk);
            const double c = std::sqrt(constants::g * std::max(0.0, hk));
            max_speed = std::max(max_speed, std::abs(u) + c);
            max_speed = std::max(max_speed, std::abs(v) + c);
        }
    }
    if (max_speed <= constants::eps) return 1.0e100;
    return cfl * std::min(g.dx, g.dy) / max_speed;
}

inline void enforce_positivity(const GridView &g, double *FAST_RESTRICT h, double *FAST_RESTRICT hu,
                               double *FAST_RESTRICT hv, const double h_floor = 1.0e-8) {
    for (int i = g.nG; i < g.nG + g.Nx; ++i) {
        for (int j = g.nG; j < g.nG + g.Ny; ++j) {
            const int k = idx(i, j, g.stride);
            if (h[k] < h_floor) {
                h[k] = h_floor;
                hu[k] = 0.0;
                hv[k] = 0.0;
            }
            if (!std::isfinite(h[k]) || !std::isfinite(hu[k]) || !std::isfinite(hv[k])) {
                throw std::runtime_error("non-finite state detected in fast solver");
            }
        }
    }
}

#undef FAST_RESTRICT

} // namespace fast_hll_muscl_bathy
