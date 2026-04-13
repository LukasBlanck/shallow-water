import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

FILE_PATH = "output/sanity_checks.nc"
SHOW_INFO_BOX = True
USE_LOG_SCALE_FOR_ERRORS = True


def get_bool_attr(ds, name, default=False):
    value = ds.attrs.get(name, int(default))
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return value != 0
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return default


ds = xr.open_dataset(FILE_PATH, engine="netcdf4")

# ---------------- Read shared metadata ----------------
time_unit = ds["time"].attrs.get("units", "not specified") if "time" in ds else "not specified"
dt_used = ds.attrs.get("dt", None)
riemann_solver = ds.attrs.get("riemann_solver", "not specified")
reconstruction = ds.attrs.get("reconstruction", "not specified")
time_integrator = ds.attrs.get("time_integrator", "not specified")
save_every = ds.attrs.get("save_every", "not specified")

mass_conservation_enabled = get_bool_attr(ds, "mass_conservation_enabled", False)
convergence_enabled = get_bool_attr(ds, "convergence_enabled", False)

# ---------------- Decide what to plot ----------------
available_plots = []

if mass_conservation_enabled and "mass_rel_err" in ds and "time" in ds:
    available_plots.append("mass_conservation")

if convergence_enabled:
    # adapt these variable names to your actual convergence writer once implemented
    if "convergence_err" in ds and "time" in ds:
        available_plots.append("convergence")

if not available_plots:
    raise RuntimeError(
        "No supported sanity-check data found in output/sanity_checks.nc"
    )

nplots = len(available_plots)
fig, axes = plt.subplots(nplots, 1, figsize=(9, 5 * nplots), squeeze=False)
axes = axes.ravel()

for ax, plot_kind in zip(axes, available_plots):
    if plot_kind == "mass_conservation":
        t = ds["time"].values
        rel_err = ds["mass_rel_err"].values

        if USE_LOG_SCALE_FOR_ERRORS:
            ax.semilogy(t, rel_err, marker="o", linewidth=1.5, markersize=4)
        else:
            ax.plot(t, rel_err, marker="o", linewidth=1.5, markersize=4)

        ax.set_xlabel(f"time [{time_unit}]".strip())
        ax.set_ylabel("relative mass error [1]")
        ax.set_title("Mass conservation")
        ax.grid(True, which="both", alpha=0.3)

        if SHOW_INFO_BOX:
            info_text = "\n".join([
                f"dt = {dt_used:.6g} {time_unit}" if dt_used is not None else "dt = n/a",
                f"Riemann: {riemann_solver}",
                f"Recon: {reconstruction}",
                f"Time int.: {time_integrator}",
            ])

            ax.text(
                0.02, 0.98, info_text,
                transform=ax.transAxes,
                va="top",
                ha="left",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.85)
            )

    elif plot_kind == "convergence":
        t = ds["time"].values
        conv_err = ds["convergence_err"].values
        conv_unit = ds["convergence_err"].attrs.get("units", "1")

        if USE_LOG_SCALE_FOR_ERRORS:
            ax.semilogy(t, conv_err, marker="o", linewidth=1.5, markersize=4)
        else:
            ax.plot(t, conv_err, marker="o", linewidth=1.5, markersize=4)

        ax.set_xlabel(f"time [{time_unit}]".strip())
        ax.set_ylabel(f"convergence error [{conv_unit}]".strip())
        ax.set_title("Convergence")
        ax.grid(True, which="both", alpha=0.3)

        if SHOW_INFO_BOX:
            info_text = "\n".join([
                f"dt = {dt_used:.6g} {time_unit}" if dt_used is not None else "dt = n/a",
                f"Riemann: {riemann_solver}",
                f"Recon: {reconstruction}",
                f"Time int.: {time_integrator}",
            ])

            ax.text(
                0.02, 0.98, info_text,
                transform=ax.transAxes,
                va="top",
                ha="left",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.85)
            )

plt.tight_layout()
plt.show()