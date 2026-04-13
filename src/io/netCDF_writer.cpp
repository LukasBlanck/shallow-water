#include "include/io/netCDF_writer.hpp"
#include "include/core/array2D.hpp"
#include "include/core/grid.hpp"
#include "include/core/state.hpp"

#include <cmath>
#include <filesystem>
#include <string>

NetCDFWriter::NetCDFWriter(const std::string &path, const Grid &grid,
                           const std::string &spatial_unit_x, const std::string &spatial_unit_y,
                           const std::string &spatial_unit_h, const std::string &time_unit,
                           int save_every)
    : grid_(grid), file_([&]() {
          std::filesystem::path p(path);
          if (p.has_parent_path()) {
              std::filesystem::create_directories(p.parent_path());
          }
          return netCDF::NcFile(path, netCDF::NcFile::replace);
      }()),
      nx_(static_cast<std::size_t>(grid_.Nx())), ny_(static_cast<std::size_t>(grid_.Ny())),
      h_buf_(nx_ * ny_), hu_buf_(nx_ * ny_), hv_buf_(nx_ * ny_), spatial_unit_x_(spatial_unit_x),
      spatial_unit_y_(spatial_unit_y), spatial_unit_h_(spatial_unit_h), time_unit_(time_unit),
      save_every_(save_every) {
    define_file_structure();
    write_coordinates();
}

void NetCDFWriter::define_file_structure() {
    time_dim_ = file_.addDim("time");
    x_dim_ = file_.addDim("x", grid_.Nx());
    y_dim_ = file_.addDim("y", grid_.Ny());

    time_var_ = file_.addVar("time", netCDF::ncDouble, time_dim_);
    x_var_ = file_.addVar("x", netCDF::ncDouble, x_dim_);
    y_var_ = file_.addVar("y", netCDF::ncDouble, y_dim_);

    h_var_ = file_.addVar("h", netCDF::ncDouble, {time_dim_, x_dim_, y_dim_});
    hu_var_ = file_.addVar("hu", netCDF::ncDouble, {time_dim_, x_dim_, y_dim_});
    hv_var_ = file_.addVar("hv", netCDF::ncDouble, {time_dim_, x_dim_, y_dim_});

    if (!spatial_unit_x_.empty()) {
        x_var_.putAtt("units", spatial_unit_x_);
    }
    if (!spatial_unit_y_.empty()) {
        y_var_.putAtt("units", spatial_unit_y_);
    }
    if (!spatial_unit_h_.empty()) {
        h_var_.putAtt("units", spatial_unit_h_);
    }
    if (!time_unit_.empty()) {
        time_var_.putAtt("units", time_unit_);
    }

    file_.putAtt("save_every", netCDF::ncInt, save_every_);
}

void NetCDFWriter::write_coordinates() {
    std::vector<double> x(nx_);
    std::vector<double> y(ny_);

    for (int i = 0; i < grid_.Nx(); ++i) {
        x[static_cast<std::size_t>(i)] = grid_.x_center(i);
    }
    for (int j = 0; j < grid_.Ny(); ++j) {
        y[static_cast<std::size_t>(j)] = grid_.y_center(j);
    }

    x_var_.putVar(x.data());
    y_var_.putVar(y.data());
}

void NetCDFWriter::pack_interior_into(const Array2D &A, std::vector<double> &buf) const {
    const int ng = grid_.nG();

    for (int i = 0; i < grid_.Nx(); ++i) {
        const std::size_t ii = static_cast<std::size_t>(i);
        for (int j = 0; j < grid_.Ny(); ++j) {
            const std::size_t jj = static_cast<std::size_t>(j);
            buf[ii * ny_ + jj] = A(i + ng, j + ng);
        }
    }
}

void NetCDFWriter::write_snapshot(const State &U, double time, double dt,
                                  const std::string &riemann_solver,
                                  const std::string &reconstruction,
                                  const std::string &time_integrator) {
    if (!metadata_written_ && (!std::isnan(dt) || !riemann_solver.empty() ||
                               !reconstruction.empty() || !time_integrator.empty())) {

        if (!std::isnan(dt)) {
            file_.putAtt("dt", netCDF::ncDouble, dt);
        }
        if (!riemann_solver.empty()) {
            file_.putAtt("riemann_solver", riemann_solver);
        }
        if (!reconstruction.empty()) {
            file_.putAtt("reconstruction", reconstruction);
        }
        if (!time_integrator.empty()) {
            file_.putAtt("time_integrator", time_integrator);
        }

        metadata_written_ = true;
    }

    time_var_.putVar(std::vector<size_t>{next_record_}, time);

    pack_interior_into(U.h(), h_buf_);
    pack_interior_into(U.hu(), hu_buf_);
    pack_interior_into(U.hv(), hv_buf_);

    const std::vector<size_t> start{next_record_, 0, 0};
    const std::vector<size_t> count{1, nx_, ny_};

    h_var_.putVar(start, count, h_buf_.data());
    hu_var_.putVar(start, count, hu_buf_.data());
    hv_var_.putVar(start, count, hv_buf_.data());

    ++next_record_;
}