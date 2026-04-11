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

x = ds["x"].values   # cell centers
y = ds["y"].values   # cell centers

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

dx_cells = np.diff(x_edges)
dy_cells = np.diff(y_edges)

X0, Y0 = np.meshgrid(x_edges[:-1], y_edges[:-1], indexing="xy")
DX, DY = np.meshgrid(dx_cells, dy_cells, indexing="xy")

xpos = X0.ravel()
ypos = Y0.ravel()
zpos = np.zeros_like(xpos)

dx_plot = DX.ravel()
dy_plot = DY.ravel()
dz_plot = Z.ravel()

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

norm = plt.Normalize(vmin=np.nanmin(Z), vmax=np.nanmax(Z))
colors = plt.cm.viridis(norm(dz_plot))

ax.bar3d(
    xpos, ypos, zpos,
    dx_plot, dy_plot, dz_plot,
    color=colors,
    shade=True,
    zsort="average"
)

mappable = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
mappable.set_array([])

fig.colorbar(
    mappable,
    ax=ax,
    shrink=0.6,
    label=f"{VAR_NAME} [{ds[VAR_NAME].attrs.get('units', '')}]".strip()
)

ax.set_xlabel(f"x [{ds['x'].attrs.get('units', '')}]".strip())
ax.set_ylabel(f"y [{ds['y'].attrs.get('units', '')}]".strip())
ax.set_zlabel(f"{VAR_NAME} [{ds[VAR_NAME].attrs.get('units', '')}]".strip())
ax.set_title(f"{VAR_NAME} at t={float(ds['time'].isel(time=TIME_INDEX)):.3f}")

plt.show()