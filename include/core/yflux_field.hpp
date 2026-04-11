#pragma once

#include "include/core/YFlux.hpp"

class YFluxField {
  public:
    YFluxField(const Grid &grid) : h_(grid), hu_(grid), hv_(grid) {}

    YFlux &h() { return h_; }
    YFlux &hu() { return hu_; }
    YFlux &hv() { return hv_; }

    const YFlux &h() const { return h_; }
    const YFlux &hu() const { return hu_; }
    const YFlux &hv() const { return hv_; }

  private:
    YFlux h_;
    YFlux hu_;
    YFlux hv_;
};