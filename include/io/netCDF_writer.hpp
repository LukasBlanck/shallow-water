#pragma once
#include <netcdf>
#include <string>
#include <vector>

class Grid;
class State;
class Array2D;

class NetCDFWriter {
  public:
    NetCDFWriter(const std::string &path, const Grid &grid);

    void write_snapshot(const State &U, double time);

  private:
    void define_file_structure();
    void write_coordinates();
    void pack_interior_into(const Array2D &A, std::vector<double> &buf) const;

  private:
    const Grid &grid_;

    netCDF::NcFile file_;

    netCDF::NcDim time_dim_;
    netCDF::NcDim x_dim_;
    netCDF::NcDim y_dim_;

    netCDF::NcVar time_var_;
    netCDF::NcVar x_var_;
    netCDF::NcVar y_var_;
    netCDF::NcVar h_var_;
    netCDF::NcVar hu_var_;
    netCDF::NcVar hv_var_;

    std::size_t next_record_ = 0;

    std::size_t nx_ = 0;
    std::size_t ny_ = 0;

    std::vector<double> h_buf_;
    std::vector<double> hu_buf_;
    std::vector<double> hv_buf_;
};