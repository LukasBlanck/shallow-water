#pragma once

#include "netcdf"

class SanityCheckNetCDFWriter {
  public:
    SanityCheckNetCDFWriter(const std::string &path, const std::string &time_unit,
                            const std::string &h_unit, int save_every, double dt,
                            const std::string &riemann_solver, const std::string &reconstruction,
                            const std::string &time_integrator,
                            const std::string &boundary_condition, const std::string &bathymetry);

    void write(std::size_t step, double time, double rel_err, double h_min);

  private:
    void define_file_structure();

    netCDF::NcFile file_;

    std::string time_unit_;
    std::string h_unit_;
    int save_every_;
    double dt_;
    std::string riemann_solver_;
    std::string reconstruction_;
    std::string time_integrator_;
    std::string boundary_condition_;
    std::string bathymetry_;

    netCDF::NcDim record_dim_;

    netCDF::NcVar step_var_;
    netCDF::NcVar time_var_;
    netCDF::NcVar mass_rel_err_var_;
    netCDF::NcVar h_min_var_;

    std::size_t next_record_{0};
};