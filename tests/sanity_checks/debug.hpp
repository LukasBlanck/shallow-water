#pragma once

#include "tests/sanity_checks/sanity_checks.hpp"
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

class Debug : public SanityCheck {
  public: // write for positivity the min of all h per time step
    const char *name() const override { return "Debug"; };
    void initialize(const State &U0, const Grid &grid) override {
        initial_mass_ = compute_mass_and_h_min(U0, grid).first;
        initialized_ = true;
    };
    void evaluate(const State &U, const Grid &grid, double time, std::size_t step,
                  SanityCheckNetCDFWriter *writer) override {
        if (!initialized_) {
            initialize(U, grid);
        }
        if (writer == nullptr) {
            throw std::runtime_error("Debug sanity check requires a sanity writer");
        }

        const std::pair<double, double> mass_and_hmin = compute_mass_and_h_min(U, grid);
        const double mass = mass_and_hmin.first;
        const double h_min = mass_and_hmin.second;

        const double abs_err = std::abs(mass - initial_mass_);
        const double rel_err_mass =
            (std::abs(initial_mass_) > 0.0) ? abs_err / std::abs(initial_mass_) : abs_err;

        writer->write(step, time, rel_err_mass, h_min);
    };

  private:
    double initial_mass_;
    bool initialized_;

    std::pair<double, double> compute_mass_and_h_min(const State &U, const Grid &grid) const {
        const int nG = grid.nG();
        const double dxdy = grid.dx() * grid.dy();
        double mass = 0.0;

        double h_limit = std::numeric_limits<double>::max();

        for (int i = nG; i < nG + grid.Nx(); ++i) {
            for (int j = nG; j < nG + grid.Ny(); ++j) {
                mass += U.h()(i, j) * dxdy;
                if (U.h()(i, j) < h_limit) {
                    h_limit = U.h()(i, j);
                }
            }
        }

        return {mass, h_limit};
    }
};