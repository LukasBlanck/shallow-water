#pragma once
#include "include/core/grid.hpp"
#include "include/core/state.hpp"
#include "include/initial_condition/initial_condition.hpp"

class still_water : public InitialCondition {
  public:
    explicit still_water(double h0) : h0_(h0) {};

    void apply(const Grid &grid, State &U) const override {
        for (int i = 0; i < grid.Nx_total(); i++) {
            for (int j = 0; j < grid.Ny_total(); j++) {
                U.h()(i, j) = h0_;
                U.hu()(i, j) = 0.0;
                U.hv()(i, j) = 0.0;
            }
        }
    }

  private:
    double h0_;
};