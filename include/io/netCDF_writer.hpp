#pragma once

#include "include/core/grid.hpp"
#include "include/core/state.hpp"

#include <cstddef>
#include <limits>
#include <netcdf>
#include <string>
#include <vector>

class NetCDFWriter {
  public:
    NetCDFWriter(const std::string &path, const Grid &grid, const std::string &spatial_unit_x,
                 const std::string &spatial_unit_y, const std::string &spatial_unit_h,
                 const std::string &time_unit, int save_every);

    void write_snapshot(const State &U, double time,
                        double dt = std::numeric_limits<double>::quiet_NaN(),
                        const std::string &riemann_solver = "",
                        const std::string &reconstruction = "",
                        const std::string &time_integrator = "");

  private:
    void define_file_structure();
    void write_coordinates();
    void pack_interior_into(const Array2D &A, std::vector<double> &buf) const;

    const Grid &grid_;
    netCDF::NcFile file_;

    netCDF::NcDim time_dim_, x_dim_, y_dim_;
    netCDF::NcVar time_var_, x_var_, y_var_;
    netCDF::NcVar h_var_, hu_var_, hv_var_;

    std::size_t nx_, ny_;
    std::size_t next_record_ = 0;

    std::vector<double> h_buf_, hu_buf_, hv_buf_;

    bool metadata_written_ = false;

    std::string spatial_unit_x_;
    std::string spatial_unit_y_;
    std::string spatial_unit_h_;
    std::string time_unit_;
    int save_every_{1};
};