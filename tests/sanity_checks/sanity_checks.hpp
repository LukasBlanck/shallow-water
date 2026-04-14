#pragma once

#include "configs/config.hpp"
#include "include/core/grid.hpp"
#include "include/core/state.hpp"
#include "include/io/sanity_checks_netdcdf_writer.hpp"

#include <cstddef>
#include <memory>
#include <vector>

class SanityCheck {
  public:
    virtual ~SanityCheck() = default;

    virtual const char *name() const = 0;
    virtual void initialize(const State &U0, const Grid &grid) = 0;
    virtual void evaluate(const State &U, const Grid &grid, double time, std::size_t step,
                          SanityCheckNetCDFWriter *writer) = 0;
};

std::vector<std::unique_ptr<SanityCheck>> make_sanity_checks(const SimulationConfig &cfg);