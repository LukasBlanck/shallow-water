#pragma once

#include "include/constants.hpp"
#include "include/core/state.hpp"
#include <iomanip>

#include <sstream>

inline double compute_stable_dt(const State &U, const Grid &grid, double cfl) {
    double dt_limit = std::numeric_limits<double>::max();

    for (int i = grid.nG(); i < grid.nG() + grid.Nx(); ++i) {
        for (int j = grid.nG(); j < grid.nG() + grid.Ny(); ++j) {
            const double h = U.h()(i, j);
            const double hu = U.hu()(i, j);
            const double hv = U.hv()(i, j);

            if (h <= 0.0) {
                throw std::runtime_error(
                    "Could not perform stability check for dt because h <= 0!");
            }

            const double u = hu / h;
            const double v = hv / h;

            const double dt_cell = cfl / ((std::abs(u) + std::sqrt(constants::g * h)) / grid.dx() +
                                          (std::abs(v) + std::sqrt(constants::g * h)) / grid.dy());

            dt_limit = std::min(dt_limit, dt_cell);
        }
    }

    return dt_limit;
}

inline double estimate_eta_seconds(std::size_t completed_steps, std::size_t total_steps,
                                   double elapsed_seconds) {
    if (completed_steps == 0 || completed_steps >= total_steps) {
        return 0.0;
    }

    const double sec_per_step = elapsed_seconds / static_cast<double>(completed_steps);

    return sec_per_step * static_cast<double>(total_steps - completed_steps);
}

inline std::string format_duration(double seconds) {
    if (seconds < 0.0) {
        seconds = 0.0;
    }

    const auto total_seconds = static_cast<long long>(std::llround(seconds));
    const long long hours = total_seconds / 3600;
    const long long minutes = (total_seconds % 3600) / 60;
    const long long secs = total_seconds % 60;

    std::ostringstream oss;
    if (hours > 0) {
        oss << hours << "h " << std::setw(2) << std::setfill('0') << minutes << "m " << std::setw(2)
            << std::setfill('0') << secs << "s";
    } else if (minutes > 0) {
        oss << minutes << "m " << std::setw(2) << std::setfill('0') << secs << "s";
    } else {
        oss << secs << "s";
    }
    return oss.str();
}