import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

# Plots initial state correctly with important metadat output

VERBOSE = True
SHOW_INFO_BOX = True

VAR_NAME = "h"
TIME_INDEX = 0
FILE_PATH = "output/simulation.nc"


def centers_to_edges(c):
    c = np.asarray(c)
    if c.size == 1:
        # fallback if only one cell exists
        return np.array([c[0] - 0.5, c[0] + 0.5])

    edges = np.empty(c.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (c[:-1] + c[1:])
    edges[0] = c[0] - 0.5 * (c[1] - c[0])
    edges[-1] = c[-1] + 0.5 * (c[-1] - c[-2])
    return edges


ds = xr.open_dataset(FILE_PATH, engine="netcdf4", chunks={"time": 1})

# ---------------- Read metadata ----------------
dt_used = ds.attrs.get("dt", None)
riemann_solver = ds.attrs.get("riemann_solver", "not specified")
reconstruction = ds.attrs.get("reconstruction", "not specified")
time_integrator = ds.attrs.get("time_integrator", "not specified")
title_attr = ds.attrs.get("title", "not specified")
model_attr = ds.attrs.get("model", "not specified")

# pick variable and timestep
frame = ds[VAR_NAME].isel(time=TIME_INDEX).transpose("y", "x")

# coordinates = cell centers from NetCDF
x = ds["x"].values
y = ds["y"].values

# cell values
Z = frame.values   # shape (ny, nx)

# ---------------- Build correct cell grid ----------------
x_edges = centers_to_edges(x)
y_edges = centers_to_edges(y)

x_edges[np.abs(x_edges) < 1e-14] = 0.0
y_edges[np.abs(y_edges) < 1e-14] = 0.0

dx = float(x[1] - x[0]) if len(x) > 1 else 0.0
dy = float(y[1] - y[0]) if len(y) > 1 else 0.0

# "Stepped" coordinates for piecewise-constant cell visualization
x_plot = np.repeat(x_edges, 2)[1:-1]   # length 2*nx
y_plot = np.repeat(y_edges, 2)[1:-1]   # length 2*ny

X_plot, Y_plot = np.meshgrid(x_plot, y_plot)
Z_plot = np.repeat(np.repeat(Z, 2, axis=0), 2, axis=1)

# ---------------- Useful derived quantities ----------------
nx = ds.sizes["x"]
ny = ds.sizes["y"]
nt = ds.sizes.get("time", 1)
physical_time = float(ds["time"].isel(time=TIME_INDEX))

if VERBOSE:
    print("=" * 60)
    print("DATASET INFO")
    print("=" * 60)
    print(f"File: {FILE_PATH}")
    print(f"Variables: {list(ds.data_vars)}")
    print(f"Coordinates: {list(ds.coords)}")
    print(f"Dimensions: {dict(ds.sizes)}")
    print()

    print("=" * 60)
    print("MESH INFO")
    print("=" * 60)
    print(f"Number of x cells: {nx}")
    print(f"Number of y cells: {ny}")
    print(f"Total number of cells: {nx * ny}")
    print(f"Number of time steps: {nt}")
    print(f"x-center range: [{float(x.min()):.6g}, {float(x.max()):.6g}]")
    print(f"y-center range: [{float(y.min()):.6g}, {float(y.max()):.6g}]")
    print(f"x-edge range:   [{float(x_edges.min()):.6g}, {float(x_edges.max()):.6g}]")
    print(f"y-edge range:   [{float(y_edges.min()):.6g}, {float(y_edges.max()):.6g}]")
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
    print(f"save_every = {ds.attrs.get("save_every", "not specified")}")
    print()
    print("Simulation setup:")
    print(f"dt: {dt_used if dt_used is not None else 'not specified'}")
    print(f"Riemann solver: {riemann_solver}")
    print(f"Reconstruction: {reconstruction}")
    print(f"Time integrator: {time_integrator}")
    print()

    print("=" * 60)
    print("SELECTED SNAPSHOT")
    print("=" * 60)
    print(f"Variable: {VAR_NAME}")
    print(f"Time index: {TIME_INDEX}")
    print(f"Physical time: {physical_time:.6g}")
    print(f"Original cell matrix shape (y, x): {Z.shape}")
    print(f"Expanded plot matrix shape: {Z_plot.shape}")
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

fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(111, projection="3d")

surf = ax.plot_surface(
    X_plot,
    Y_plot,
    Z_plot,
    cmap="viridis",
    linewidth=0.15,
    edgecolor="k",
    antialiased=False,
    shade=False
)

fig.colorbar(
    surf,
    ax=ax,
    shrink=0.6,
    pad=0.12,
    label=f"{VAR_NAME} [{ds[VAR_NAME].attrs.get('units', '')}]".strip()
)

ax.set_xlabel(f"x [{ds['x'].attrs.get('units', '')}]".strip())
ax.set_ylabel(f"y [{ds['y'].attrs.get('units', '')}]".strip())
ax.set_zlabel(f"{VAR_NAME} [{ds[VAR_NAME].attrs.get('units', '')}]".strip())
ax.set_title(f"{VAR_NAME} at t={physical_time:.3f} {ds['time'].attrs.get('units', 'not specified')}")

# use actual cell edges as plot bounds
ax.set_xlim(float(x_edges[0]), float(x_edges[-1]))
ax.set_ylim(float(y_edges[0]), float(y_edges[-1]))

if SHOW_INFO_BOX:
    info_text = "\n".join([
        f"t = {physical_time:.6g} {ds['time'].attrs.get('units', 'not specified')}",
        f"dt = {dt_used:.6g} {ds['time'].attrs.get('units', 'not specified')}" if dt_used is not None else "dt ='n/a'",
        f"Riemann: {riemann_solver}",
        f"Recon: {reconstruction}",
        f"Time int.: {time_integrator}",
    ])

    ax.text2D(
        0.02, 0.98, info_text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85)
    )

plt.tight_layout()
plt.show()