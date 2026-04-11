#pragma once

#include "include/core/grid.hpp"
#include "include/core/state.hpp"

class InitialCondition {
  public:
    virtual ~InitialCondition() = default;
    virtual void apply(const Grid &grid, State &U) const = 0;
};