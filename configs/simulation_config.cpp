#include "configs/simulation_config.hpp"
#include "configs/config.hpp"

#include <cstdint>
#include <extern/toml.hpp>
#include <stdexcept>
#include <string>

namespace {

BackendType parse_backend_type(const std::string &s) {
    if (s == "serial") return BackendType::Serial;
    throw std::runtime_error("Unknown backend type: " + s);
}

BoundaryType parse_boundary_type(const std::string &s) {
    if (s == "ReflectingWalls") return BoundaryType::ReflectingWalls;
    throw std::runtime_error("Unknown boundary type: " + s);
}

ReconstructionType parse_reconstruction_type(const std::string &s) {
    if (s == "PiecewiseConst") return ReconstructionType::PiecewiseConst;
    throw std::runtime_error("Unknown reconstruction type: " + s);
}

RiemannType parse_riemann_type(const std::string &s) {
    if (s == "Rusanov") return RiemannType::Rusanov;
    throw std::runtime_error("Unknown riemann solver: " + s);
}

TimeIntegratorType parse_time_integrator_type(const std::string &s) {
    if (s == "SSPRK3") return TimeIntegratorType::SSPRK3;
    throw std::runtime_error("Unknown time integrator: " + s);
}

InitialConditionType parse_initial_condition_type(const std::string &s) {
    if (s == "GaussInitial")
        return InitialConditionType::GaussInitial;
    else if (s == "StillWater")
        return InitialConditionType::StillWater;
    throw std::runtime_error("Unknown initial condition: " + s);
}

BathymetryType parse_bathymetry_type(const std::string &s) {
    if (s == "Flat") return BathymetryType::Flat;
    throw std::runtime_error("Unknown bathymetry type: " + s);
}

template <typename T>
T require_value(const toml::node_view<const toml::node> &node, const std::string &key_path) {
    if (auto v = node.template value<T>()) {
        return *v;
    }
    throw std::runtime_error("Missing or invalid TOML value for key '" + key_path + "'");
}

} // namespace

SimulationConfig load_config(const std::filesystem::path &path) {
    toml::table tbl;

    try {
        tbl = toml::parse_file(path.string());
    } catch (const toml::parse_error &err) {
        throw std::runtime_error("Failed to parse TOML file '" + path.string() +
                                 "': " + std::string(err.description()));
    }

    SimulationConfig cfg{};

    const toml::table *mesh = tbl["mesh"].as_table();
    if (!mesh) throw std::runtime_error("Missing or invalid [mesh] table");

    cfg.mesh.Lx = require_value<double>((*mesh)["Lx"], "mesh.Lx");
    cfg.mesh.Ly = require_value<double>((*mesh)["Ly"], "mesh.Ly");
    cfg.mesh.nG = static_cast<int>(require_value<std::int64_t>((*mesh)["nG"], "mesh.nG"));
    cfg.mesh.Nx = static_cast<int>(require_value<std::int64_t>((*mesh)["Nx"], "mesh.Nx"));
    cfg.mesh.Ny = static_cast<int>(require_value<std::int64_t>((*mesh)["Ny"], "mesh.Ny"));
    cfg.mesh.spatial_unit_x =
        require_value<std::string>((*mesh)["spatial_unit_x"], "mesh.spatial_unit_x");
    cfg.mesh.spatial_unit_y =
        require_value<std::string>((*mesh)["spatial_unit_y"], "mesh.spatial_unit_y");
    cfg.mesh.spatial_unit_h =
        require_value<std::string>((*mesh)["spatial_unit_h"], "mesh.spatial_unit_h");

    const toml::table *bathy = tbl["bathymetry"].as_table();
    if (!bathy) throw std::runtime_error("Missing or invalid [bathymetry] table");

    cfg.bathymetry.type =
        parse_bathymetry_type(require_value<std::string>((*bathy)["type"], "bathymetry.type"));

    const toml::table *boundary = tbl["boundary"].as_table();
    if (!boundary) throw std::runtime_error("Missing or invalid [boundary] table");

    cfg.boundary.type =
        parse_boundary_type(require_value<std::string>((*boundary)["type"], "boundary.type"));

    const toml::table *ic = tbl["initial_condition"].as_table();
    if (!ic) throw std::runtime_error("Missing or invalid [initial_condition] table");

    cfg.initial_condition.type = parse_initial_condition_type(
        require_value<std::string>((*ic)["type"], "initial_condition.type"));

    cfg.initial_condition.peak_height =
        require_value<double>((*ic)["peak_height"], "initial_condition.peak_height");
    cfg.initial_condition.x0 = require_value<double>((*ic)["x0"], "initial_condition.x0");
    cfg.initial_condition.y0 = require_value<double>((*ic)["y0"], "initial_condition.y0");
    cfg.initial_condition.sigma_x =
        require_value<double>((*ic)["sigma_x"], "initial_condition.sigma_x");
    cfg.initial_condition.sigma_y =
        require_value<double>((*ic)["sigma_y"], "initial_condition.sigma_y");
    cfg.initial_condition.h0 = require_value<double>((*ic)["h0"], "initial_condition.h0");

    const toml::table *solver = tbl["solver"].as_table();
    if (!solver) throw std::runtime_error("Missing or invalid [solver] table");

    cfg.solver.reconstruction = parse_reconstruction_type(
        require_value<std::string>((*solver)["reconstruction"], "solver.reconstruction"));
    cfg.solver.riemann =
        parse_riemann_type(require_value<std::string>((*solver)["riemann"], "solver.riemann"));
    cfg.solver.time =
        parse_time_integrator_type(require_value<std::string>((*solver)["time"], "solver.time"));
    cfg.solver.limiter = require_value<std::string>((*solver)["limiter"], "solver.limiter");
    cfg.solver.positivity_preserving =
        require_value<bool>((*solver)["positivity_preserving"], "solver.positivity_preserving");

    const toml::table *time = tbl["time"].as_table();
    if (!time) throw std::runtime_error("Missing or invalid [time] table");

    cfg.time.end_time = require_value<double>((*time)["end_time"], "time.end_time");
    cfg.time.time_steps = static_cast<std::size_t>(
        require_value<std::int64_t>((*time)["time_steps"], "time.time_steps"));
    cfg.time.cfl = require_value<double>((*time)["cfl"], "time.cfl");
    cfg.time.save_every =
        static_cast<int>(require_value<std::int64_t>((*time)["save_every"], "time.save_every"));
    cfg.time.time_unit = require_value<std::string>((*time)["time_unit"], "time.time_unit");

    const toml::table *backend = tbl["backend"].as_table();
    if (!backend) throw std::runtime_error("Missing or invalid [backend] table");

    cfg.backend.type =
        parse_backend_type(require_value<std::string>((*backend)["type"], "backend.type"));
    cfg.backend.threads =
        static_cast<int>(require_value<std::int64_t>((*backend)["threads"], "backend.threads"));

    const toml::table *output = tbl["output"].as_table();
    if (!output) throw std::runtime_error("Missing or invalid [output] table");

    cfg.output.path = require_value<std::string>((*output)["path"], "output.path");

    const toml::table *sanity = tbl["sanity_checks"].as_table();
    if (sanity) {
        cfg.sanity_checks.mass_conservation =
            (*sanity)["mass_conservation"].value<bool>().value_or(false);

        cfg.sanity_checks.positivity = (*sanity)["positivity"].value<bool>().value_or(false);

        cfg.sanity_checks.debug = (*sanity)["debug"].value<bool>().value_or(false);

        cfg.sanity_checks.output_path =
            require_value<std::string>((*sanity)["output_path"], "output_path");
    }

    return cfg;
}