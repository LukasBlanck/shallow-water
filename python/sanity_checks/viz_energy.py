import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# ============================================================
# User settings
# ============================================================

SANITY_FILE_PATH = "output/sanity_checks.nc"

# IMPORTANT:
# This must be the file written by NetCDFWriter, i.e. the one containing
# h(time,x,y), hu(time,x,y), hv(time,x,y), B(x,y).
FULL_STATE_FILE_PATH = "output/simulation.nc"

SAVE_PLOTS = False
OUTPUT_DIR = "output/diffusivity_plots"

G = 9.81
H_DRY = 1e-12

USE_LOG_SCALE_FOR_ERRORS = True
SHOW_INFO_BOX = True

PRINT_DATASET_SUMMARY_TO_CLI = False
PRINT_DIFFUSIVITY_SUMMARY_TO_CLI = True

SNAPSHOT_INDICES = [10, 20, 50]  # None means first, middle, last


# ============================================================
# Helper functions
# ============================================================

def get_bool_attr(ds, name, default=False):
    value = ds.attrs.get(name, int(default))
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return value != 0
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return default


def infer_spacing_from_coords(ds):
    dx = ds.attrs.get("dx", None)
    dy = ds.attrs.get("dy", None)

    if dx is None and "x" in ds.coords and ds["x"].size > 1:
        dx = float(np.mean(np.diff(ds["x"].values)))

    if dy is None and "y" in ds.coords and ds["y"].size > 1:
        dy = float(np.mean(np.diff(ds["y"].values)))

    if dx is None:
        dx = 1.0
    if dy is None:
        dy = 1.0

    return float(dx), float(dy)


def relative_change(series, eps=1e-300):
    series = np.asarray(series, dtype=float)
    denom = max(abs(series[0]), eps)
    return (series - series[0]) / denom


def safe_velocity(q, h, h_dry=H_DRY):
    out = np.zeros_like(q, dtype=float)
    mask = h > h_dry
    out[mask] = q[mask] / h[mask]
    return out


def total_variation_2d_series(h):
    tv = np.empty(h.shape[0], dtype=float)

    for n in range(h.shape[0]):
        snap = np.asarray(h[n])
        tv_x = np.abs(np.diff(snap, axis=0)).sum()
        tv_y = np.abs(np.diff(snap, axis=1)).sum()
        tv[n] = tv_x + tv_y

    return tv


def get_snapshot_indices(nt):
    if SNAPSHOT_INDICES is not None:
        return [idx for idx in SNAPSHOT_INDICES if 0 <= idx < nt]

    if nt == 1:
        return [0]

    return sorted(set([0, nt // 2, nt - 1]))


def save_or_show(fig, filename=None):
    plt.tight_layout()

    if SAVE_PLOTS and filename is not None:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        path = os.path.join(OUTPUT_DIR, filename)
        fig.savefig(path, dpi=300, bbox_inches="tight")
        print(f"Saved plot to: {path}")

    plt.show()


def require_vars(ds, names, file_label):
    missing = [name for name in names if name not in ds]
    if missing:
        raise RuntimeError(
            f"{file_label} is missing variables {missing}. "
            f"Available variables are {list(ds.data_vars)}."
        )


# ============================================================
# Load sanity file
# ============================================================

sanity_ds = xr.open_dataset(SANITY_FILE_PATH, engine="netcdf4")

if PRINT_DATASET_SUMMARY_TO_CLI:
    print("\n===== Sanity NetCDF contents =====")
    print(sanity_ds)
    print("==================================\n")

require_vars(sanity_ds, ["time", "mass_rel_err", "h_min"], SANITY_FILE_PATH)

t_sanity = sanity_ds["time"].values
mass_rel_err = sanity_ds["mass_rel_err"].values
h_min_series = sanity_ds["h_min"].values

time_unit = sanity_ds["time"].attrs.get("units", "not specified")

dt_used = sanity_ds.attrs.get("dt", None)
backend = sanity_ds.attrs.get("backend", "not specified")
riemann_solver = sanity_ds.attrs.get("riemann_solver", "not specified")
reconstruction = sanity_ds.attrs.get("reconstruction", "not specified")
time_integrator = sanity_ds.attrs.get("time_integrator", "not specified")
boundary_condition = sanity_ds.attrs.get("boundary_condition", "not specified")
bathymetry = sanity_ds.attrs.get("bathymetry", "not specified")
save_every = sanity_ds.attrs.get("save_every", "not specified")

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


# ============================================================
# Plot sanity checks first
# ============================================================

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

fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

ax = axes[0]
if USE_LOG_SCALE_FOR_ERRORS and np.any(np.abs(mass_rel_err) > 0):
    ax.semilogy(t_sanity, np.abs(mass_rel_err), marker="o", linewidth=1.5, markersize=4)
    ax.set_ylabel("|relative mass error|")
else:
    ax.plot(t_sanity, mass_rel_err, marker="o", linewidth=1.5, markersize=4)
    ax.set_ylabel("relative mass error")

ax.set_title("Mass conservation")
ax.grid(True, which="both", alpha=0.3)

if SHOW_INFO_BOX:
    ax.text(
        0.01,
        0.97,
        info_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )

ax = axes[1]
ax.plot(t_sanity, h_min_series, marker="o", linewidth=1.5, markersize=4)
ax.set_xlabel(f"time [{time_unit}]".strip())
ax.set_ylabel("h_min")
ax.set_title("Minimum water height")
ax.grid(True, alpha=0.3)

save_or_show(fig, "01_sanity_checks.png")


# ============================================================
# Load full state file
# ============================================================

if not os.path.exists(FULL_STATE_FILE_PATH):
    raise RuntimeError(
        f"Full state file not found: {FULL_STATE_FILE_PATH}\n"
        "The diffusivity diagnostics need the NetCDF file written by NetCDFWriter, "
        "not only sanity_checks.nc."
    )

state_ds = xr.open_dataset(FULL_STATE_FILE_PATH, engine="netcdf4")

if PRINT_DATASET_SUMMARY_TO_CLI:
    print("\n===== Full state NetCDF contents =====")
    print(state_ds)
    print("======================================\n")

require_vars(state_ds, ["time", "h", "hu", "hv"], FULL_STATE_FILE_PATH)

has_B = "B" in state_ds

h = state_ds["h"].values          # shape: time, x, y
hu = state_ds["hu"].values        # shape: time, x, y
hv = state_ds["hv"].values        # shape: time, x, y
t = state_ds["time"].values

if has_B:
    B_2d = state_ds["B"].values   # shape: x, y
    B = np.broadcast_to(B_2d, h.shape)
else:
    B = None

dx, dy = infer_spacing_from_coords(state_ds)
cell_area = dx * dy

nt = h.shape[0]
spatial_axes = (1, 2)


# ============================================================
# Compute diffusivity diagnostics
# ============================================================

mass = np.sum(h, axis=spatial_axes) * cell_area
mass_rel_change_from_state = relative_change(mass)

h_min = np.min(h, axis=spatial_axes)
h_max = np.max(h, axis=spatial_axes)
h_mean = np.mean(h, axis=spatial_axes)
h_std = np.std(h, axis=spatial_axes)

h_amplitude = h_max - h_min
h_amplitude_rel_change = relative_change(h_amplitude)

tv_h = total_variation_2d_series(h)
tv_h_rel_change = relative_change(tv_h)

u = safe_velocity(hu, h)
v = safe_velocity(hv, h)
speed = np.sqrt(u**2 + v**2)

momentum_abs_x = np.sum(np.abs(hu), axis=spatial_axes) * cell_area
momentum_abs_y = np.sum(np.abs(hv), axis=spatial_axes) * cell_area
momentum_abs_total = np.sum(np.sqrt(hu**2 + hv**2), axis=spatial_axes) * cell_area

momentum_x_total = np.sum(hu, axis=spatial_axes) * cell_area
momentum_y_total = np.sum(hv, axis=spatial_axes) * cell_area

speed_max = np.max(speed, axis=spatial_axes)
speed_mean = np.mean(speed, axis=spatial_axes)

kinetic_density = np.zeros_like(h)
mask = h > H_DRY
kinetic_density[mask] = 0.5 * (hu[mask]**2 + hv[mask]**2) / h[mask]
kinetic_energy = np.sum(kinetic_density, axis=spatial_axes) * cell_area

potential_energy = np.sum(0.5 * G * h**2, axis=spatial_axes) * cell_area

if B is not None:
    bathymetry_energy = np.sum(G * h * B, axis=spatial_axes) * cell_area
else:
    bathymetry_energy = None

total_energy = potential_energy + kinetic_energy
if bathymetry_energy is not None:
    total_energy = total_energy + bathymetry_energy

energy_rel_change = relative_change(total_energy)

if PRINT_DIFFUSIVITY_SUMMARY_TO_CLI:
    print("Diffusivity / damping diagnostics:")
    print(f"  Full state file             = {FULL_STATE_FILE_PATH}")
    print(f"  h shape                     = {h.shape}")
    print(f"  dx, dy                      = {dx}, {dy}")
    print()
    print(f"  final sanity mass rel. err  = {mass_rel_err[-1]:.12e}")
    print(f"  final state mass rel. chg   = {mass_rel_change_from_state[-1]:.12e}")
    print()
    print(f"  initial h amplitude         = {h_amplitude[0]:.12e}")
    print(f"  final h amplitude           = {h_amplitude[-1]:.12e}")
    print(f"  amplitude rel. change       = {h_amplitude_rel_change[-1]:.12e}")
    print()
    print(f"  initial TV(h) proxy         = {tv_h[0]:.12e}")
    print(f"  final TV(h) proxy           = {tv_h[-1]:.12e}")
    print(f"  TV(h) rel. change           = {tv_h_rel_change[-1]:.12e}")
    print()
    print(f"  initial total energy        = {total_energy[0]:.12e}")
    print(f"  final total energy          = {total_energy[-1]:.12e}")
    print(f"  energy rel. change          = {energy_rel_change[-1]:.12e}")
    print()
    print(f"  initial sum |momentum|      = {momentum_abs_total[0]:.12e}")
    print(f"  final sum |momentum|        = {momentum_abs_total[-1]:.12e}")
    print(f"  rel. change sum |momentum|  = {relative_change(momentum_abs_total)[-1]:.12e}")
    print()
    print(f"  initial kinetic energy      = {kinetic_energy[0]:.12e}")
    print(f"  final kinetic energy        = {kinetic_energy[-1]:.12e}")
    print(f"  rel. change kinetic energy  = {relative_change(kinetic_energy)[-1]:.12e}")
    print()


# ============================================================
# Plot height diffusion indicators
# ============================================================

fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)

axes[0].plot(t, h_amplitude, marker="o", linewidth=1.5, markersize=4)
axes[0].set_ylabel("max(h) - min(h)")
axes[0].set_title("Height amplitude")
axes[0].grid(True, alpha=0.3)

axes[1].plot(t, h_std, marker="o", linewidth=1.5, markersize=4)
axes[1].set_ylabel("std(h)")
axes[1].set_title("Spatial standard deviation of h")
axes[1].grid(True, alpha=0.3)

axes[2].plot(t, tv_h, marker="o", linewidth=1.5, markersize=4)
axes[2].set_xlabel(f"time [{time_unit}]".strip())
axes[2].set_ylabel("TV(h) proxy")
axes[2].set_title("Total variation proxy of h")
axes[2].grid(True, alpha=0.3)

save_or_show(fig, "02_height_diffusion_indicators.png")


# ============================================================
# Plot energy diagnostics
# ============================================================

fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)

axes[0].plot(t, potential_energy, marker="o", linewidth=1.5, markersize=4, label="potential")
axes[0].plot(t, kinetic_energy, marker="o", linewidth=1.5, markersize=4, label="kinetic")

if bathymetry_energy is not None:
    axes[0].plot(t, bathymetry_energy, marker="o", linewidth=1.5, markersize=4, label="bathymetry")

axes[0].set_ylabel("energy")
axes[0].set_title("Energy components")
axes[0].grid(True, alpha=0.3)
axes[0].legend()

axes[1].plot(t, total_energy, marker="o", linewidth=1.5, markersize=4)
axes[1].set_ylabel("total energy")
axes[1].set_title("Total mechanical energy proxy")
axes[1].grid(True, alpha=0.3)

axes[2].plot(t, energy_rel_change, marker="o", linewidth=1.5, markersize=4)
axes[2].set_xlabel(f"time [{time_unit}]".strip())
axes[2].set_ylabel("relative change")
axes[2].set_title("Relative total energy change")
axes[2].grid(True, alpha=0.3)

save_or_show(fig, "03_energy_diagnostics.png")


# ============================================================
# Plot momentum diagnostics
# ============================================================

fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)

axes[0].plot(t, momentum_abs_x, marker="o", linewidth=1.5, markersize=4, label="sum |hu|")
axes[0].plot(t, momentum_abs_y, marker="o", linewidth=1.5, markersize=4, label="sum |hv|")
axes[0].plot(t, momentum_abs_total, marker="o", linewidth=1.5, markersize=4, label="sum sqrt(hu²+hv²)")
axes[0].set_ylabel("absolute momentum")
axes[0].set_title("Absolute momentum content")
axes[0].grid(True, alpha=0.3)
axes[0].legend()

axes[1].plot(t, momentum_x_total, marker="o", linewidth=1.5, markersize=4, label="sum hu")
axes[1].plot(t, momentum_y_total, marker="o", linewidth=1.5, markersize=4, label="sum hv")
axes[1].set_ylabel("net momentum")
axes[1].set_title("Net momentum")
axes[1].grid(True, alpha=0.3)
axes[1].legend()

axes[2].plot(t, speed_max, marker="o", linewidth=1.5, markersize=4, label="max speed")
axes[2].plot(t, speed_mean, marker="o", linewidth=1.5, markersize=4, label="mean speed")
axes[2].set_xlabel(f"time [{time_unit}]".strip())
axes[2].set_ylabel("speed")
axes[2].set_title("Velocity magnitude")
axes[2].grid(True, alpha=0.3)
axes[2].legend()

save_or_show(fig, "04_momentum_diagnostics.png")


# ============================================================
# Plot spatial snapshots
# ============================================================

snapshot_indices = get_snapshot_indices(nt)

for idx in snapshot_indices:
    time_label = f"t = {t[idx]:.6g}"

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    im = axes[0, 0].imshow(h[idx].T, origin="lower", aspect="auto")
    axes[0, 0].set_title(f"h, {time_label}")
    plt.colorbar(im, ax=axes[0, 0])

    im = axes[0, 1].imshow(hu[idx].T, origin="lower", aspect="auto")
    axes[0, 1].set_title(f"hu, {time_label}")
    plt.colorbar(im, ax=axes[0, 1])

    im = axes[1, 0].imshow(hv[idx].T, origin="lower", aspect="auto")
    axes[1, 0].set_title(f"hv, {time_label}")
    plt.colorbar(im, ax=axes[1, 0])

    im = axes[1, 1].imshow(speed[idx].T, origin="lower", aspect="auto")
    axes[1, 1].set_title(f"speed, {time_label}")
    plt.colorbar(im, ax=axes[1, 1])

    save_or_show(fig, f"05_snapshot_{idx:04d}.png")


# ============================================================
# Interpretation guide
# ============================================================

print("Interpretation guide:")
print("  Small mass error:")
print("    Good conservation, but not proof of low diffusion.")
print()
print("  Decreasing max(h)-min(h), std(h), or TV(h):")
print("    Indicates smoothing of the height field.")
print()
print("  Decreasing kinetic energy or sum |momentum|:")
print("    Indicates damping of the velocity / momentum field.")
print()
print("  Mass conserved but h/hu/hv snapshots visibly spreading:")
print("    Conservative but numerically diffusive behavior.")