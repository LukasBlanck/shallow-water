import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import cm
from matplotlib.colors import Normalize
import numpy as np

FILE_PATH = "output/simulation.nc"
VAR_NAME = "h"

SAVE_VIDEO = True
VIDEO_PATH = "output/simulation_3d.mp4"   # or "output/simulation_3d.gif"
FPS = 150

CMAP = "viridis"

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


# ---------------- Load dataset ----------------
ds = xr.open_dataset(FILE_PATH, engine="netcdf4")

# ---------------- Metadata: read once, print once ----------------
dt_used = ds.attrs.get("dt", None)
riemann_solver = ds.attrs.get("riemann_solver", "n/a")
reconstruction = ds.attrs.get("reconstruction", "n/a")
time_integrator = ds.attrs.get("time_integrator", "n/a")

print("Simulation metadata:")
print(f"  dt              = {dt_used}" if dt_used is not None else "  dt              = n/a")
print(f"  Riemann solver  = {riemann_solver}")
print(f"  reconstruction  = {reconstruction}")
print(f"  time integrator = {time_integrator}")

info_text = "\n".join([
    f"dt = {dt_used}" if dt_used is not None else "dt = n/a",
    f"Riemann: {riemann_solver}",
    f"Recon: {reconstruction}",
    f"Time int.: {time_integrator}",
])

x = ds["x"].values
y = ds["y"].values
t = ds["time"].values

# Stored as (time, x, y) -> transpose to (time, y, x)
data = ds[VAR_NAME].transpose("time", "y", "x").load().values

if data.shape[0] == 0:
    raise RuntimeError("No time frames found in dataset.")

# Optional spatial downsampling for speed
x_plot = x[::SPATIAL_STRIDE]
y_plot = y[::SPATIAL_STRIDE]
data_plot = data[:, ::SPATIAL_STRIDE, ::SPATIAL_STRIDE]

X, Y = np.meshgrid(x_plot, y_plot)

# ---------------- Fixed range from initial frame ----------------
z0 = data_plot[0]
z0_min = float(np.nanmin(z0))
z0_max = float(np.nanmax(z0))
z0_span = z0_max - z0_min
buffer = max(BUFFER_FRAC * z0_span, MIN_BUFFER)

vmin = z0_min - buffer
vmax = z0_max + buffer

norm = Normalize(vmin=vmin, vmax=vmax)

# ---------------- Figure / axes ----------------
fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(111, projection="3d")

ax.set_xlabel(f"x [{ds['x'].attrs.get('units', '')}]".strip())
ax.set_ylabel(f"y [{ds['y'].attrs.get('units', '')}]".strip())
ax.set_zlabel(f"{VAR_NAME} [{ds[VAR_NAME].attrs.get('units', '')}]".strip())

ax.set_xlim(float(x_plot[0]), float(x_plot[-1]))
ax.set_ylim(float(y_plot[0]), float(y_plot[-1]))
ax.set_zlim(vmin, vmax)
ax.view_init(elev=ELEV, azim=AZIM)

title = ax.set_title(f"{VAR_NAME} at t={float(t[0]):.6f}")

if SHOW_INFO_BOX:
    ax.text2D(
        0.02, 0.98, info_text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85)
    )

# Initial surface
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
cbar.set_label(f"{VAR_NAME} [{ds[VAR_NAME].attrs.get('units', '')}]".strip())


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

    title.set_text(f"{VAR_NAME} at t={float(t[frame_idx]):.6f}")
    return surf, title


ani = FuncAnimation(
    fig,
    update,
    frames=data_plot.shape[0],
    interval=1000 / FPS,
    blit=False,   # 3D matplotlib does not benefit from blit
    repeat=True
)

plt.tight_layout()

if SAVE_VIDEO:
    if VIDEO_PATH.endswith(".gif"):
        ani.save(VIDEO_PATH, writer="pillow", fps=FPS)
    else:
        ani.save(VIDEO_PATH, writer="ffmpeg", fps=FPS)

plt.show()