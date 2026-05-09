#!/bin/bash
import os
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import cm
from matplotlib.colors import Normalize
import numpy as np

FILE_PATH = "output/simulation.nc"
VAR_NAME = "eta"  # eta = h + B

SAVE_VIDEO = True
VIDEO_PATH = "output/simulation_3d.mp4"   # or "output/simulation_3d.gif"

ADAPTIVE_FPS = True
FPS = 150   # used only if ADAPTIVE_FPS = False

CMAP = "viridis"

# Show bathymetry B as a dark bottom surface
SHOW_B_BOTTOM = True

# Performance knob:
# 1 = use every cell
# 2 = every second cell
# 4 = every fourth cell
SPATIAL_STRIDE = 1

# Fixed z/color range from initial frame only
BUFFER_FRAC = 0.08
MIN_BUFFER = 1e-12

# Camera
ELEV = 30
AZIM = -135

# Optional info box in plot
SHOW_INFO_BOX = True


def time_unit_to_seconds_factor(unit):
    u = str(unit).strip().lower()

    mapping = {
        "s": 1.0,
        "sec": 1.0,
        "secs": 1.0,
        "second": 1.0,
        "seconds": 1.0,

        "ms": 1e-3,
        "millisecond": 1e-3,
        "milliseconds": 1e-3,

        "us": 1e-6,
        "µs": 1e-6,
        "microsecond": 1e-6,
        "microseconds": 1e-6,

        "min": 60.0,
        "mins": 60.0,
        "minute": 60.0,
        "minutes": 60.0,

        "h": 3600.0,
        "hr": 3600.0,
        "hrs": 3600.0,
        "hour": 3600.0,
        "hours": 3600.0,
    }

    if u not in mapping:
        raise ValueError(f"Unsupported time unit for adaptive FPS: {unit!r}")

    return mapping[u]


# ---------------- Load dataset ----------------
ds = xr.open_dataset(FILE_PATH, engine="netcdf4")

# ---------------- Metadata: read once, print once ----------------
dt_used = ds.attrs.get("dt", None)
backend = ds.attrs.get("backend", "n/a")
riemann_solver = ds.attrs.get("riemann_solver", "n/a")
reconstruction = ds.attrs.get("reconstruction", "n/a")
time_integrator = ds.attrs.get("time_integrator", "n/a")
boundary_condition = ds.attrs.get("boundary_condition", "n/a")
bathymetry = ds.attrs.get("bathymetry", "n/a")

print("Simulation metadata:")
print(f"  dt              = {dt_used}" if dt_used is not None else "  dt              = n/a")
print(f"  Backend         = {backend}")
print(f"  Riemann solver  = {riemann_solver}")
print(f"  reconstruction  = {reconstruction}")
print(f"  time integrator = {time_integrator}")
print(f"  BC              = {boundary_condition}")
print(f"  Bathymetry      = {bathymetry}")

time_unit = ds["time"].attrs.get("units", "s")
save_every = ds.attrs.get("save_every", "n/a")

info_text = "\n".join([
    f"dt = {dt_used:.6g} {time_unit}" if dt_used is not None else "dt = n/a",
    f"save_every = {save_every}",
    f"Backend: {backend}",
    f"Riemann: {riemann_solver}",
    f"Recon: {reconstruction}",
    f"Time int.: {time_integrator}",
    f"BC: {boundary_condition}",
    f"Bathymetry: {bathymetry}",
])

x = ds["x"].values
y = ds["y"].values
t = ds["time"].values

if "h" not in ds:
    raise RuntimeError(f"Variable 'h' not found in {FILE_PATH}")

if "B" not in ds:
    raise RuntimeError(f"Variable 'B' not found in {FILE_PATH}")

# Stored h as (time, x, y), B as (x, y).
# Transpose to plotting layout: h -> (time, y, x), B -> (y, x).
# Efficient: spatially stride before loading, then broadcast B over time.
h_plot = (
    ds["h"]
    .isel(x=slice(None, None, SPATIAL_STRIDE), y=slice(None, None, SPATIAL_STRIDE))
    .transpose("time", "y", "x")
    .load()
    .values
)

B_plot = (
    ds["B"]
    .isel(x=slice(None, None, SPATIAL_STRIDE), y=slice(None, None, SPATIAL_STRIDE))
    .transpose("y", "x")
    .load()
    .values
)

data_plot = h_plot + B_plot[None, :, :]  # eta = h + B

if data_plot.shape[0] == 0:
    raise RuntimeError("No time frames found in dataset.")

if ADAPTIVE_FPS:
    if len(t) < 2:
        raise RuntimeError("Adaptive FPS needs at least two saved frames.")

    dt_frame = float(t[1] - t[0])
    if dt_frame <= 0.0:
        raise RuntimeError("Non-positive saved time difference found in dataset.")

    seconds_per_time_unit = time_unit_to_seconds_factor(time_unit)
    dt_frame_seconds = dt_frame * seconds_per_time_unit
    FPS_USED = 1.0 / dt_frame_seconds
else:
    FPS_USED = float(FPS)

print(f"  save_every      = {save_every}")
print(f"  time unit       = {time_unit}")
print(f"  fps used        = {FPS_USED:.6g}")

# Optional spatial downsampling for speed
x_plot = x[::SPATIAL_STRIDE]
y_plot = y[::SPATIAL_STRIDE]

X, Y = np.meshgrid(x_plot, y_plot)

# ---------------- Fixed range from initial frame ----------------
z0 = data_plot[0]

z0_min = float(np.nanmin(z0))
z0_max = float(np.nanmax(z0))

if SHOW_B_BOTTOM:
    z0_min = min(z0_min, float(np.nanmin(B_plot)))
    z0_max = max(z0_max, float(np.nanmax(B_plot)))

z0_span = z0_max - z0_min
buffer = max(BUFFER_FRAC * z0_span, MIN_BUFFER)

vmin = z0_min - buffer
vmax = z0_max + buffer

norm = Normalize(vmin=vmin, vmax=vmax)

# ---------------- Figure / axes ----------------
fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(111, projection="3d")

x_unit = ds["x"].attrs.get("units", "")
y_unit = ds["y"].attrs.get("units", "")
eta_unit = ds["h"].attrs.get("units", "")

ax.set_xlabel(f"x [{x_unit}]".strip())
ax.set_ylabel(f"y [{y_unit}]".strip())
ax.set_zlabel(rf"$\eta$ [{eta_unit}]".strip())

ax.set_xlim(float(x_plot[0]), float(x_plot[-1]))
ax.set_ylim(float(y_plot[0]), float(y_plot[-1]))
ax.set_zlim(vmin, vmax)
ax.view_init(elev=ELEV, azim=AZIM)

title = ax.set_title(
    rf"$\eta = h + B$ at t={float(t[0]):.6f} {time_unit}"
)

if SHOW_INFO_BOX:
    ax.text2D(
        0.02, 0.98, info_text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85)
    )

# Static bathymetry bottom surface
if SHOW_B_BOTTOM:
    ax.plot_surface(
        X,
        Y,
        B_plot,
        color="black",
        alpha=0.45,
        linewidth=0,
        antialiased=False,
        shade=True
    )

# Initial eta surface
surf = ax.plot_surface(
    X,
    Y,
    data_plot[0],
    cmap=CMAP,
    norm=norm,
    linewidth=0,
    antialiased=False,
    shade=False
)

# Independent colorbar with fixed normalization
mappable = cm.ScalarMappable(norm=norm, cmap=CMAP)
mappable.set_array([])
cbar = fig.colorbar(mappable, ax=ax, shrink=0.75, pad=0.12)
cbar.set_label(rf"$\eta$ [{eta_unit}]".strip())


def update(frame_idx):
    global surf
    surf.remove()

    surf = ax.plot_surface(
        X,
        Y,
        data_plot[frame_idx],
        cmap=CMAP,
        norm=norm,
        linewidth=0,
        antialiased=False,
        shade=False
    )

    title.set_text(
        rf"$\eta = h + B$ at t={float(t[frame_idx]):.6f} {time_unit}"
    )
    return surf, title


ani = FuncAnimation(
    fig,
    update,
    frames=data_plot.shape[0],
    interval=1000 / FPS_USED,
    blit=False,   # 3D matplotlib does not benefit from blit
    repeat=True
)

plt.tight_layout()

if SAVE_VIDEO:
    if VIDEO_PATH.endswith(".gif"):
        ani.save(VIDEO_PATH, writer="pillow", fps=FPS_USED)
    else:
        ani.save(VIDEO_PATH, writer="ffmpeg", fps=FPS_USED)

plt.show()