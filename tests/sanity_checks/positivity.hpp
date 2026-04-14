#pragma once

#include "include/core/grid.hpp"
#include "include/core/state.hpp"
#include "include/io/sanity_checks_netdcdf_writer.hpp"
#include "tests/sanity_checks/sanity_checks.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>

class Positivity : public SanityCheck {
  public:
    const char *name() const override { return "Positivity"; }

    void initialize(const State &U0, const Grid &grid) override {
        check_nonnegative(U0, grid, 0.0, 0);
    }

    void evaluate(const State &U, const Grid &grid, double time, std::size_t step,
                  SanityCheckNetCDFWriter *writer) override {
        check_nonnegative(U, grid, time, step);
    }

  private:
    void check_nonnegative(const State &U, const Grid &grid, double time, std::size_t step) const {
        const int nG = grid.nG();

        for (int i = nG; i < nG + grid.Nx(); ++i) {
            for (int j = nG; j < nG + grid.Ny(); ++j) {
                const double h = U.h()(i, j);

                if (h < 0.0) {
                    throw std::runtime_error(
                        "Positivity check failed: negative water height h = " + std::to_string(h) +
                        " at cell (" + std::to_string(i - nG) + ", " + std::to_string(j - nG) +
                        "), step = " + std::to_string(step) + ", time = " + std::to_string(time));
                }
            }
        }
    }
};