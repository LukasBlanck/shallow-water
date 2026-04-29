#pragma once

#include "include/boundary/boundary_condition.hpp"

class ReflectingWalls : public BoundaryCondition {
  public:
    ReflectingWalls(Grid &grid) : BoundaryCondition(grid) {}

    void apply_BC(State &U) override {
        // fill top and bottom
        for (int g = 0; g < nG_; ++g) {
            int jg_bot = nG_ - 1 - g; // ghost row at bottom
            int ji_bot = nG_ + g;     // mirrored interior row

            int jg_top = Ny_total_ - nG_ + g;     // ghost row at top
            int ji_top = Ny_total_ - nG_ - 1 - g; // mirrored interior row

            for (int i = 0; i < Nx_total_; ++i) {
                // bottom
                U.h()(i, jg_bot) = U.h()(i, ji_bot);
                U.hu()(i, jg_bot) = U.hu()(i, ji_bot);
                U.hv()(i, jg_bot) = -U.hv()(i, ji_bot);

                // top
                U.h()(i, jg_top) = U.h()(i, ji_top);
                U.hu()(i, jg_top) = U.hu()(i, ji_top);
                U.hv()(i, jg_top) = -U.hv()(i, ji_top);
            }
        }

        // fill left and right
        for (int g = 0; g < nG_; ++g) {
            int ig_left = nG_ - 1 - g;
            int ii_left = nG_ + g;

            int ig_right = Nx_total_ - nG_ + g;
            int ii_right = Nx_total_ - nG_ - 1 - g;

            for (int j = 0; j < Ny_total_; ++j) {
                // left wall
                U.h()(ig_left, j) = U.h()(ii_left, j);
                U.hu()(ig_left, j) = -U.hu()(ii_left, j);
                U.hv()(ig_left, j) = U.hv()(ii_left, j);

                // right wall
                U.h()(ig_right, j) = U.h()(ii_right, j);
                U.hu()(ig_right, j) = -U.hu()(ii_right, j);
                U.hv()(ig_right, j) = U.hv()(ii_right, j);
            }
        }
    }
    void apply_BC(Array2D &B) override {

        // fill top and bottom

        for (int g = 0; g < nG_; ++g) {

            int jg_bot = nG_ - 1 - g;

            int ji_bot = nG_ + g;

            int jg_top = Ny_total_ - nG_ + g;

            int ji_top = Ny_total_ - nG_ - 1 - g;

            for (int i = 0; i < Nx_total_; ++i) {

                B(i, jg_bot) = B(i, ji_bot);

                B(i, jg_top) = B(i, ji_top);
            }
        }

        // fill left and right

        for (int g = 0; g < nG_; ++g) {

            int ig_left = nG_ - 1 - g;

            int ii_left = nG_ + g;

            int ig_right = Nx_total_ - nG_ + g;

            int ii_right = Nx_total_ - nG_ - 1 - g;

            for (int j = 0; j < Ny_total_; ++j) {

                B(ig_left, j) = B(ii_left, j);

                B(ig_right, j) = B(ii_right, j);
            }
        }
    }
};