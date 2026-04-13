#pragma once

#include <cstddef>
#include <netcdf>
#include <string>

class SanityCheckNetCDFWriter {
  public:
    SanityCheckNetCDFWriter(const std::string &path,
                            const std::string &time_unit,
                            int save_every,
                            double dt,
                            const std::string &riemann_solver,
                            const std::string &reconstruction,
                            const std::string &time_integrator);

    void write_mass_conservation(std::size_t step, double time, double mass,
                                 double abs_err, double rel_err);

  private:
    void define_file_structure();

    netCDF::NcFile file_;

    netCDF::NcDim record_dim_;
    netCDF::NcVar step_var_;
    netCDF::NcVar time_var_;
    netCDF::NcVar mass_var_;
    netCDF::NcVar mass_abs_err_var_;
    netCDF::NcVar mass_rel_err_var_;

    std::size_t next_record_{0};

    std::string time_unit_;
    int save_every_{1};
    double dt_{};
    std::string riemann_solver_;
    std::string reconstruction_;
    std::string time_integrator_;
};