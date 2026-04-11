#pragma once

class Grid {
  public:
    Grid(int Nx, int Ny, double Lx, double Ly, int nG)
        : Nx_(Nx), Ny_(Ny), Lx_(Lx), Ly_(Ly), nG_(nG) {};

    int Nx() const noexcept { return Nx_; }
    int Ny() const noexcept { return Ny_; }

    int Nx_total() const noexcept { return Nx_ + 2 * nG_; }
    int Ny_total() const noexcept { return Ny_ + 2 * nG_; }

    int nG() const noexcept { return nG_; }

    double Lx() const noexcept { return Lx_; }
    double Ly() const noexcept { return Ly_; }

    double dx() const noexcept { return Lx_ / Nx_; }  // implicit cast from double Lx_
    double dy() const noexcept { return Ly_ / Ny_; }

    double x_center(int i) const noexcept { return (i + 0.5) * dx(); }
    double y_center(int j) const noexcept { return (j + 0.5) * dy(); }

  private:
    int Nx_; // number of cells
    int Ny_;
    double Lx_; // length of domain
    double Ly_;
    int nG_; // number of ghost cells
};