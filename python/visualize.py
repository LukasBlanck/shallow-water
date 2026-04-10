import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

VERBOSE = True
VAR_NAME = "h"
TIME_INDEX = 0
FILE_PATH = "output/simulation.nc"

ds = xr.open_dataset(FILE_PATH, engine="netcdf4", chunks={"time": 1})

# pick variable and timestep
frame = ds[VAR_NAME].isel(time=TIME_INDEX).transpose("y", "x")

# coordinates
x = ds["x"].values
y = ds["y"].values
X, Y = np.meshgrid(x, y)

Z = frame.values

if VERBOSE:
    print("=" * 60)
    print("DATASET INFO")
    print("=" * 60)
    print(f"File: {FILE_PATH}")
    print(f"Variables: {list(ds.data_vars)}")
    print(f"Coordinates: {list(ds.coords)}")
    print(f"Dimensions: {dict(ds.sizes)}")
    print()

    nx = ds.sizes["x"]
    ny = ds.sizes["y"]
    nt = ds.sizes.get("time", 1)

    dx = float(x[1] - x[0]) if len(x) > 1 else 0.0
    dy = float(y[1] - y[0]) if len(y) > 1 else 0.0

    print("=" * 60)
    print("MESH INFO")
    print("=" * 60)
    print(f"Number of x cells: {nx}")
    print(f"Number of y cells: {ny}")
    print(f"Total number of cells: {nx * ny}")
    print(f"Number of time steps: {nt}")
    print(f"x range: [{float(x.min()):.6g}, {float(x.max()):.6g}]")
    print(f"y range: [{float(y.min()):.6g}, {float(y.max()):.6g}]")
    print(f"dx: {dx:.6g}")
    print(f"dy: {dy:.6g}")
    print()

    print("=" * 60)
    print("UNITS / METADATA")
    print("=" * 60)
    print(f"{VAR_NAME} units: {ds[VAR_NAME].attrs.get('units', 'not specified')}")
    print(f"x units: {ds['x'].attrs.get('units', 'not specified')}")
    print(f"y units: {ds['y'].attrs.get('units', 'not specified')}")
    print(f"time units: {ds['time'].attrs.get('units', 'not specified')}")
    print(f"title: {ds.attrs.get('title', 'not specified')}")
    print(f"model: {ds.attrs.get('model', 'not specified')}")
    print()

    print("=" * 60)
    print("SELECTED SNAPSHOT")
    print("=" * 60)
    print(f"Variable: {VAR_NAME}")
    print(f"Time index: {TIME_INDEX}")
    print(f"Physical time: {float(ds['time'].isel(time=TIME_INDEX)):.6g}")
    print(f"Shape of plotted matrix (y, x): {Z.shape}")
    print(f"Min value: {np.nanmin(Z):.6g}")
    print(f"Max value: {np.nanmax(Z):.6g}")
    print(f"Mean value: {np.nanmean(Z):.6g}")
    print()

    print("=" * 60)
    print(f"{VAR_NAME} VALUES AS MATRIX")
    print("=" * 60)
    with np.printoptions(precision=4, suppress=True, linewidth=140, threshold=400):
        print(Z)
    print()

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

surf = ax.plot_surface(X, Y, Z, cmap="viridis")

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