#pragma once
#include "include/initial_condition/initial_condition.hpp"

class StillWater : public InitialCondition {
  public:
    explicit StillWater(double h0) : h0_(h0) {};

    void apply(const Grid &grid, State &U) const override {
        const int nG = grid.nG();
        for (int i = 0; i < grid.Nx(); i++) {
            for (int j = 0; j < grid.Ny(); j++) {
                U.h()(i + nG, j + nG) = h0_;
                U.hu()(i + nG, j + nG) = 0.0;
                U.hv()(i + nG, j + nG) = 0.0;
            }
        }
    }

  private:
    double h0_;
};