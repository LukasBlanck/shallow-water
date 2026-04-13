#pragma once

#include "include/core/grid.hpp"
#include "include/core/state.hpp"
#include "include/io/sanity_checks_netdcdf_writer.hpp"
#include "tests/sanity_checks/sanity_checks.hpp"

#include <cmath>
#include <cstddef>

class MassConservation : public SanityCheck {
  public:
    const char *name() const override { return "MassConservation"; }

    void initialize(const State &U0, const Grid &grid) override {
        initial_mass_ = compute_mass(U0, grid);
        initialized_ = true;
    }

    void evaluate(const State &U, const Grid &grid, double time, std::size_t step,
                  SanityCheckNetCDFWriter &writer) override {
        if (!initialized_) {
            initialize(U, grid);
        }

        const double mass = compute_mass(U, grid);
        const double abs_err = std::abs(mass - initial_mass_);
        const double rel_err =
            (std::abs(initial_mass_) > 0.0) ? abs_err / std::abs(initial_mass_) : abs_err;

        writer.write_mass_conservation(step, time, mass, abs_err, rel_err);
    }

  private:
    double compute_mass(const State &U, const Grid &grid) const {
        const int nG = grid.nG();
        const double dxdy = grid.dx() * grid.dy();
        double mass = 0.0;

        for (int i = nG; i < nG + grid.Nx(); ++i) {
            for (int j = nG; j < nG + grid.Ny(); ++j) {
                mass += U.h()(i, j) * dxdy;
            }
        }
        return mass;
    }

    double initial_mass_{0.0};
    bool initialized_{false};
};