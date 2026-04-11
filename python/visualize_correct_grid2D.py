import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

VERBOSE = True
VAR_NAME = "h"
TIME_INDEX = 0
FILE_PATH = "output/simulation.nc"

ds = xr.open_dataset(FILE_PATH, engine="netcdf4", chunks={"time": 1})

frame = ds[VAR_NAME].isel(time=TIME_INDEX).transpose("y", "x")
Z = frame.values

x = ds["x"].values
y = ds["y"].values

def centers_to_edges(c):
    c = np.asarray(c)
    e = np.empty(c.size + 1)

    if c.size == 1:
        e[0] = c[0] - 0.5
        e[1] = c[0] + 0.5
        return e

    e[1:-1] = 0.5 * (c[:-1] + c[1:])
    e[0] = c[0] - 0.5 * (c[1] - c[0])
    e[-1] = c[-1] + 0.5 * (c[-1] - c[-2])
    return e

x_edges = centers_to_edges(x)
y_edges = centers_to_edges(y)

Xe, Ye = np.meshgrid(x_edges, y_edges)

fig, ax = plt.subplots()

pcm = ax.pcolormesh(Xe, Ye, Z, shading="flat", cmap="viridis")

fig.colorbar(
    pcm,
    ax=ax,
    label=f"{VAR_NAME} [{ds[VAR_NAME].attrs.get('units', '')}]".strip()
)

ax.set_xlabel(f"x [{ds['x'].attrs.get('units', '')}]".strip())
ax.set_ylabel(f"y [{ds['y'].attrs.get('units', '')}]".strip())
ax.set_title(f"{VAR_NAME} at t={float(ds['time'].isel(time=TIME_INDEX)):.3f}")
ax.set_aspect("equal")

# grid lines at cell boundaries
ax.set_xticks(x_edges)
ax.set_yticks(y_edges)
ax.grid(True, color="white", linewidth=0.5)
ax.set_axisbelow(False)  # draw grid on top

plt.show()