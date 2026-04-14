#pragma once

#include "configs/config.hpp"
#include "include/core/grid.hpp"
#include "include/core/state.hpp"
#include "include/initial_condition/dam_break.hpp"
#include "include/initial_condition/dam_break_radial.hpp"
#include "include/initial_condition/gauss_initial.hpp"
#include "include/initial_condition/still_water.hpp"

#include <stdexcept>

inline void apply_initial_condition(const SimulationConfig &cfg, const Grid &grid, State &U) {
    switch (cfg.initial_condition.type) {
    case InitialConditionType::GaussInitial: {
        GaussInitial ic(cfg.initial_condition.peak_height, cfg.initial_condition.sigma_x,
                        cfg.initial_condition.sigma_y, cfg.initial_condition.x0,
                        cfg.initial_condition.y0, cfg.initial_condition.h0);
        ic.apply(grid, U);
        return;
    }
    case InitialConditionType::StillWater: {
        StillWater ic(cfg.initial_condition.h0);
        ic.apply(grid, U);
        return;
    }
    case InitialConditionType::DamBreak: {
        DamBreak ic(cfg.initial_condition.dam_height, cfg.initial_condition.h0,
                    cfg.initial_condition.dam_x);
        ic.apply(grid, U);
        return;
    }
    case InitialConditionType::DamBreakRadial: {
        DamBreakRadial ic(cfg.initial_condition.dam_height, cfg.initial_condition.h0,
                          cfg.initial_condition.dam_x0, cfg.initial_condition.dam_y0,
                          cfg.initial_condition.dam_radius);
        ic.apply(grid, U);
        return;
    }
    }

    throw std::runtime_error("Unsupported initial condition");
}