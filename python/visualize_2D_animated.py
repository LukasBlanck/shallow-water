import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

FILE_PATH = "output/simulation.nc"
VAR_NAME = "h"

SAVE_VIDEO = True
VIDEO_PATH = "output/simulation.mp4"   # or "output/simulation.gif"
FPS = 40

SHOW_GRID = True
CMAP = "viridis"

# fixed range from initial frame only
BUFFER_FRAC = 0.08   # 8% buffer around initial min/max
MIN_BUFFER = 1e-12   # fallback if initial field is constant


def centers_to_edges(c):
    c = np.asarray(c, dtype=float)
    e = np.empty(c.size + 1, dtype=float)

    if c.size == 1:
        e[0] = c[0] - 0.5
        e[1] = c[0] + 0.5
        return e

    e[1:-1] = 0.5 * (c[:-1] + c[1:])
    e[0] = c[0] - 0.5 * (c[1] - c[0])
    e[-1] = c[-1] + 0.5 * (c[-1] - c[-2])
    return e


# ---------------- Load dataset ----------------
ds = xr.open_dataset(FILE_PATH, engine="netcdf4")

x = ds["x"].values
y = ds["y"].values
t = ds["time"].values

# Stored as (time, x, y) in file -> transpose to (time, y, x) for plotting
data = ds[VAR_NAME].transpose("time", "y", "x").load().values

n_frames = data.shape[0]
if n_frames == 0:
    raise RuntimeError("No time frames found in dataset.")

# ---------------- Fixed color range from initial frame ----------------
z0 = data[0]
z0_min = float(np.nanmin(z0))
z0_max = float(np.nanmax(z0))
z0_span = z0_max - z0_min
buffer = max(BUFFER_FRAC * z0_span, MIN_BUFFER)

vmin = z0_min - buffer
vmax = z0_max + buffer

# ---------------- Plot geometry ----------------
x_edges = centers_to_edges(x)
y_edges = centers_to_edges(y)

extent = [x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]]

fig, ax = plt.subplots(figsize=(7, 5))

im = ax.imshow(
    data[0],
    origin="lower",
    extent=extent,
    cmap=CMAP,
    vmin=vmin,
    vmax=vmax,
    interpolation="nearest",
    aspect="equal",
    animated=True
)

cbar = fig.colorbar(
    im,
    ax=ax,
    label=f"{VAR_NAME} [{ds[VAR_NAME].attrs.get('units', '')}]".strip()
)

ax.set_xlabel(f"x [{ds['x'].attrs.get('units', '')}]".strip())
ax.set_ylabel(f"y [{ds['y'].attrs.get('units', '')}]".strip())
title = ax.set_title(f"{VAR_NAME} at t={float(t[0]):.6f}")

if SHOW_GRID:
    ax.set_xticks(x_edges)
    ax.set_yticks(y_edges)
    ax.grid(True, color="white", linewidth=0.5)
    ax.set_axisbelow(False)

# lock axes once
ax.set_xlim(x_edges[0], x_edges[-1])
ax.set_ylim(y_edges[0], y_edges[-1])


def update(frame_idx):
    im.set_data(data[frame_idx])
    title.set_text(f"{VAR_NAME} at t={float(t[frame_idx]):.6f}")
    return im, title


ani = FuncAnimation(
    fig,
    update,
    frames=n_frames,
    interval=1000 / FPS,
    blit=True,
    repeat=True
)

plt.tight_layout()

if SAVE_VIDEO:
    if VIDEO_PATH.endswith(".gif"):
        ani.save(VIDEO_PATH, writer="pillow", fps=FPS)
    else:
        ani.save(VIDEO_PATH, writer="ffmpeg", fps=FPS)

plt.show()