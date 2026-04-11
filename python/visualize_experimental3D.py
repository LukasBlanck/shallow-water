import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

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

x_step = np.repeat(x_edges, 2)[1:-1]
y_step = np.repeat(y_edges, 2)[1:-1]
Z_step = np.repeat(np.repeat(Z, 2, axis=1), 2, axis=0)

X_step, Y_step = np.meshgrid(x_step, y_step)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

surf = ax.plot_surface(X_step, Y_step, Z_step, cmap="viridis", linewidth=0)

fig.colorbar(
    surf,
    ax=ax,
    shrink=0.6,
    label=f"{VAR_NAME} [{ds[VAR_NAME].attrs.get('units', '')}]".strip()
)

ax.set_xlabel(f"x [{ds['x'].attrs.get('units', '')}]".strip())
ax.set_ylabel(f"y [{ds['y'].attrs.get('units', '')}]".strip())
ax.set_zlabel(f"{VAR_NAME} [{ds[VAR_NAME].attrs.get('units', '')}]".strip())
ax.set_title(f"{VAR_NAME} at t={float(ds['time'].isel(time=TIME_INDEX)):.3f}")

plt.show()