import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

FILE_PATH = "output/sanity_checks.nc"

SHOW_INFO_BOX = True
USE_LOG_SCALE_FOR_ERRORS = True

SAVE_PLOTS = False
OUTPUT_PATH = "output/sanity_checks_plots.png"

# ---------------- CLI printing options ----------------
PRINT_DATASET_SUMMARY_TO_CLI = False
PRINT_SERIES_TO_CLI = True
CLI_SERIES_TO_PRINT = "h_min"   # choose: "h_min" or "mass_rel_err"


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

# ---------------- Optional dataset summary ----------------
if PRINT_DATASET_SUMMARY_TO_CLI:
    print("\n===== NetCDF file contents =====")
    print(ds)
    print("================================\n")

# ---------------- Read shared metadata ----------------
time_unit = ds["time"].attrs.get("units", "not specified") if "time" in ds else "not specified"
h_unit = ds["h_min"].attrs.get("units", "not specified") if "h_min" in ds else "not specified"

dt_used = ds.attrs.get("dt", None)
backend = ds.attrs.get("backend", "not specified")
riemann_solver = ds.attrs.get("riemann_solver", "not specified")
reconstruction = ds.attrs.get("reconstruction", "not specified")
time_integrator = ds.attrs.get("time_integrator", "not specified")
boundary_condition = ds.attrs.get("boundary_condition", "not specified")
bathymetry = ds.attrs.get("bathymetry", "not specified")
save_every = ds.attrs.get("save_every", "not specified")

mass_conservation_enabled = get_bool_attr(ds, "mass_conservation_enabled", False)
positivity_enabled = get_bool_attr(ds, "positivity_enabled", False)

# ---------------- Print metadata to CLI ----------------
print("Simulation metadata:")
print(f"  dt              = {dt_used}" if dt_used is not None else "  dt              = n/a")
print(f"  Backend         = {backend}")
print(f"  Riemann solver  = {riemann_solver}")
print(f"  reconstruction  = {reconstruction}")
print(f"  time integrator = {time_integrator}")
print(f"  BC              = {boundary_condition}")
print(f"  Bathymetry      = {bathymetry}")
print(f"  save_every      = {save_every}")
print()

# ---------------- Optional explicit series output ----------------
if PRINT_SERIES_TO_CLI:
    allowed_series = {"h_min", "mass_rel_err"}
    if CLI_SERIES_TO_PRINT not in allowed_series:
        raise ValueError(
            f"CLI_SERIES_TO_PRINT must be one of {allowed_series}, got '{CLI_SERIES_TO_PRINT}'"
        )

    if CLI_SERIES_TO_PRINT not in ds:
        raise RuntimeError(f"Variable '{CLI_SERIES_TO_PRINT}' not found in {FILE_PATH}")

    t = ds["time"].values if "time" in ds else None
    values = ds[CLI_SERIES_TO_PRINT].values
    steps = ds["step"].values if "step" in ds else None
    unit = ds[CLI_SERIES_TO_PRINT].attrs.get("units", "not specified")

    print(f"\n===== Explicit CLI output: {CLI_SERIES_TO_PRINT} =====")
    print(f"unit = {unit}")

    if t is not None and steps is not None:
        print(f"{'idx':>5} {'step':>10} {'time':>15} {'value':>20}")
        print("-" * 56)
        for i, (s, ti, val) in enumerate(zip(steps, t, values)):
            print(f"{i:5d} {int(s):10d} {ti:15.8f} {val:20.12e}")
    elif t is not None:
        print(f"{'idx':>5} {'time':>15} {'value':>20}")
        print("-" * 44)
        for i, (ti, val) in enumerate(zip(t, values)):
            print(f"{i:5d} {ti:15.8f} {val:20.12e}")
    else:
        print(f"{'idx':>5} {'value':>20}")
        print("-" * 28)
        for i, val in enumerate(values):
            print(f"{i:5d} {val:20.12e}")

    print("=" * 56 + "\n")

# ---------------- Decide what to plot ----------------
available_plots = []

if mass_conservation_enabled and "mass_rel_err" in ds and "time" in ds:
    available_plots.append("mass_conservation")

if positivity_enabled and "h_min" in ds and "time" in ds:
    available_plots.append("positivity")

if not available_plots:
    raise RuntimeError(f"No supported sanity-check data found in {FILE_PATH}")

# ---------------- Build info box text ----------------
info_lines = [
    f"dt = {dt_used:.6g} {time_unit}" if dt_used is not None else "dt = n/a",
    f"save_every = {save_every}",
    f"Backend: {backend}",
    f"Riemann: {riemann_solver}",
    f"Recon: {reconstruction}",
    f"Time int.: {time_integrator}",
    f"BC: {boundary_condition}",
    f"Bathymetry: {bathymetry}",
]
info_text = "\n".join(info_lines)

# ---------------- Create figure with dedicated info subplot ----------------
nplots = len(available_plots)
fig = plt.figure(figsize=(9, 2.4 + 5 * nplots))
gs = fig.add_gridspec(nplots + 1, 1, height_ratios=[1.6] + [5] * nplots)

info_ax = fig.add_subplot(gs[0])
info_ax.axis("off")

if SHOW_INFO_BOX:
    info_ax.text(
        0.01,
        0.95,
        info_text,
        ha="left",
        va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )

axes = [fig.add_subplot(gs[i + 1]) for i in range(nplots)]

for ax, plot_kind in zip(axes, available_plots):
    if plot_kind == "mass_conservation":
        t = ds["time"].values
        rel_err = ds["mass_rel_err"].values

        if USE_LOG_SCALE_FOR_ERRORS and np.any(rel_err > 0):
            ax.semilogy(t, rel_err, marker="o", linewidth=1.5, markersize=4)
        else:
            ax.plot(t, rel_err, marker="o", linewidth=1.5, markersize=4)

        ax.set_xlabel(f"time [{time_unit}]".strip())
        ax.set_ylabel("relative mass error")
        ax.set_title("Mass conservation")
        ax.grid(True, which="both", alpha=0.3)

    elif plot_kind == "positivity":
        t = ds["time"].values
        h_min = ds["h_min"].values

        ax.plot(t, h_min, marker="o", linewidth=1.5, markersize=4)

        ax.set_xlabel(f"time [{time_unit}]".strip())
        ax.set_ylabel(f"h_min [{h_unit}]".strip())
        ax.set_title("Minimum water height")
        ax.grid(True, which="both", alpha=0.3)

plt.tight_layout()

if SAVE_PLOTS:
    output_dir = os.path.dirname(OUTPUT_PATH)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    fig.savefig(OUTPUT_PATH, dpi=300, bbox_inches="tight")
    print(f"Saved plot to: {OUTPUT_PATH}")

plt.show()