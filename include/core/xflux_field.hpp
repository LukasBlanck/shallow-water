#pragma once

#include "include/core/xflux.hpp"

class XFluxField {
  public:
    XFluxField(const Grid &grid) : h_(grid), hu_(grid), hv_(grid) {}

    XFlux &h() { return h_; }
    XFlux &hu() { return hu_; }
    XFlux &hv() { return hv_; }

    const XFlux &h() const { return h_; }
    const XFlux &hu() const { return hu_; }
    const XFlux &hv() const { return hv_; }

  private:
    XFlux h_;
    XFlux hu_;
    XFlux hv_;
};