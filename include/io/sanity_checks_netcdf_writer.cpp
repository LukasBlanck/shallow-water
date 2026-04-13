#include "include/io/sanity_checks_netdcdf_writer.hpp"

#include <filesystem>
#include <vector>

SanityCheckNetCDFWriter::SanityCheckNetCDFWriter(const std::string &path,
                                                 const std::string &time_unit, int save_every,
                                                 double dt, const std::string &riemann_solver,
                                                 const std::string &reconstruction,
                                                 const std::string &time_integrator)
    : file_([&]() {
          std::filesystem::path p(path);
          if (p.has_parent_path()) {
              std::filesystem::create_directories(p.parent_path());
          }
          return netCDF::NcFile(path, netCDF::NcFile::replace);
      }()),
      time_unit_(time_unit), save_every_(save_every), dt_(dt), riemann_solver_(riemann_solver),
      reconstruction_(reconstruction), time_integrator_(time_integrator) {
    define_file_structure();
}

void SanityCheckNetCDFWriter::define_file_structure() {
    record_dim_ = file_.addDim("record");

    step_var_ = file_.addVar("step", netCDF::ncUint64, record_dim_);
    time_var_ = file_.addVar("time", netCDF::ncDouble, record_dim_);
    mass_var_ = file_.addVar("mass", netCDF::ncDouble, record_dim_);
    mass_abs_err_var_ = file_.addVar("mass_abs_err", netCDF::ncDouble, record_dim_);
    mass_rel_err_var_ = file_.addVar("mass_rel_err", netCDF::ncDouble, record_dim_);

    if (!time_unit_.empty()) {
        time_var_.putAtt("units", time_unit_);
    }

    mass_var_.putAtt("long_name", "total_mass");
    mass_abs_err_var_.putAtt("long_name", "absolute_mass_error");
    mass_rel_err_var_.putAtt("long_name", "relative_mass_error");
    mass_rel_err_var_.putAtt("units", "1");

    file_.putAtt("save_every", netCDF::ncInt, save_every_);
    file_.putAtt("dt", netCDF::ncDouble, dt_);
    file_.putAtt("riemann_solver", riemann_solver_);
    file_.putAtt("reconstruction", reconstruction_);
    file_.putAtt("time_integrator", time_integrator_);

    file_.putAtt("mass_conservation_enabled", netCDF::ncInt, 1);
    // file_.putAtt("convergence_enabled", netCDF::ncInt, 0);
}

void SanityCheckNetCDFWriter::write_mass_conservation(std::size_t step, double time, double mass,
                                                      double abs_err, double rel_err) {
    const std::vector<std::size_t> start{next_record_};

    const auto step_u64 = static_cast<unsigned long long>(step);

    step_var_.putVar(start, &step_u64);
    time_var_.putVar(start, &time);
    mass_var_.putVar(start, &mass);
    mass_abs_err_var_.putVar(start, &abs_err);
    mass_rel_err_var_.putVar(start, &rel_err);

    ++next_record_;
}