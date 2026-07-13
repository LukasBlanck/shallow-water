#!/usr/bin/env python3

from pathlib import Path
import argparse

import pandas as pd
import matplotlib.pyplot as plt


COMPONENTS = [
    "cfl_s",
    "host_to_device_s",
    "device_to_host_s",
    "boundary_conditions_s",
    "halo_exchange_s",
    "fused_rhs_s",
    "time_integration_s",
    "positivity_correction_s",
    "output_s",
    "sanity_checks_s",
]


PRETTY_NAMES = {
    "cfl_s": "CFL check",
    "host_to_device_s": "Host to device",
    "device_to_host_s": "Device to host",
    "boundary_conditions_s": "Boundary conditions",
    "halo_exchange_s": "Halo exchange",
    "fused_rhs_s": "Fused RHS",
    "time_integration_s": "Time integration",
    "positivity_correction_s": "Positivity correction",
    "output_s": "Output",
    "sanity_checks_s": "Sanity checks",
}


# Bar width for component plots.
# Since the x-axis uses actual Nx values, this is also measured in Nx units.
# Try 20 for thinner bars, 35 for medium bars, 60 for thicker bars.
BAR_WIDTH = 50


def read_results(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    numeric_cols = [
        "gpus",
        "Nx",
        "Ny",
        "cells",
        "repeat",
        "time_steps",
        "rhs_calls",
        "solver_time_s",
        "success",
        *COMPONENTS,
    ]

    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df[df["success"] == 1].copy()

    if df.empty:
        raise RuntimeError("No successful benchmark rows found in the CSV file.")

    return df


def aggregate(df: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["backend", "gpus", "Nx", "Ny", "cells"]

    value_cols = [
        "solver_time_s",
        "time_steps",
        "rhs_calls",
        *COMPONENTS,
    ]

    grouped_mean = (
        df.groupby(group_cols, as_index=False)[value_cols]
        .mean(numeric_only=True)
    )

    grouped_std = (
        df.groupby(group_cols, as_index=False)[value_cols]
        .std(numeric_only=True)
        .rename(columns={col: f"{col}_std" for col in value_cols})
    )

    grouped_count = (
        df.groupby(group_cols, as_index=False)["solver_time_s"]
        .count()
        .rename(columns={"solver_time_s": "n_successful_repeats"})
    )

    grouped = grouped_mean.merge(grouped_std, on=group_cols, how="left")
    grouped = grouped.merge(grouped_count, on=group_cols, how="left")

    # Standard error of the mean: std / sqrt(n)
    for col in value_cols:
        std_col = f"{col}_std"
        sem_col = f"{col}_sem"

        if std_col in grouped.columns:
            grouped[sem_col] = grouped[std_col] / grouped["n_successful_repeats"] ** 0.5

    grouped = grouped.sort_values(["Nx", "gpus"])

    return grouped


def plot_solver_time(df: pd.DataFrame, outdir: Path) -> None:
    plt.figure()

    for backend, part in df.groupby("backend"):
        part = part.sort_values("cells")

        yerr = None
        if "solver_time_s_std" in part.columns:
            yerr = part["solver_time_s_std"]

        plt.errorbar(
            part["cells"],
            part["solver_time_s"],
            yerr=yerr,
            marker="o",
            capsize=4,
            label=backend,
        )

    plt.xlabel("Grid cells Nx × Ny")
    plt.ylabel("Measured solver time [s]")
    plt.title("Scaling of solver runtime")
    plt.xscale("log")
    plt.yscale("log")
    plt.grid(True, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "solver_time_vs_cells_with_error.png", dpi=200)
    plt.close()


def plot_solver_time_vs_n(df: pd.DataFrame, outdir: Path) -> None:
    plt.figure()

    for backend, part in df.groupby("backend"):
        part = part.sort_values("Nx")

        yerr = None
        if "solver_time_s_std" in part.columns:
            yerr = part["solver_time_s_std"]

        plt.errorbar(
            part["Nx"],
            part["solver_time_s"],
            yerr=yerr,
            marker="o",
            capsize=4,
            label=backend,
        )

    plt.xlabel("N for quadratic grid N × N")
    plt.ylabel("Measured solver time [s]")
    plt.title("Solver runtime for quadratic grids, mean ± std")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "solver_time_vs_N_with_error.png", dpi=200)
    plt.close()


def plot_speedup(df: pd.DataFrame, outdir: Path) -> None:
    one = df[df["backend"] == "CUDA"][
        ["Nx", "solver_time_s", "solver_time_s_std", "n_successful_repeats"]
    ].rename(
        columns={
            "solver_time_s": "time_1gpu",
            "solver_time_s_std": "std_1gpu",
            "n_successful_repeats": "n_1gpu",
        }
    )

    four = df[df["backend"] == "CUDA_4"][
        ["Nx", "solver_time_s", "solver_time_s_std", "n_successful_repeats"]
    ].rename(
        columns={
            "solver_time_s": "time_4gpu",
            "solver_time_s_std": "std_4gpu",
            "n_successful_repeats": "n_4gpu",
        }
    )

    merged = one.merge(four, on="Nx", how="inner")

    if merged.empty:
        print("Skipping speedup plot: need matching CUDA and CUDA_4 rows.")
        return

    merged["speedup_4gpu_vs_1gpu"] = merged["time_1gpu"] / merged["time_4gpu"]
    merged["parallel_efficiency"] = merged["speedup_4gpu_vs_1gpu"] / 4.0

    # Error propagation for speedup S = T1 / T4:
    # sigma_S = S * sqrt((sigma_T1/T1)^2 + (sigma_T4/T4)^2)
    merged["speedup_std"] = merged["speedup_4gpu_vs_1gpu"] * (
        (merged["std_1gpu"] / merged["time_1gpu"]) ** 2
        + (merged["std_4gpu"] / merged["time_4gpu"]) ** 2
    ) ** 0.5

    merged["efficiency_std"] = merged["speedup_std"] / 4.0

    plt.figure()
    plt.errorbar(
        merged["Nx"],
        merged["speedup_4gpu_vs_1gpu"],
        yerr=merged["speedup_std"],
        marker="o",
        capsize=4,
        label="Measured speedup",
    )
    plt.axhline(1.0, linestyle="--", label="No speedup")
    plt.axhline(4.0, linestyle="--", label="Ideal 4× speedup")
    plt.xlabel("N for quadratic grid N × N")
    plt.ylabel("Speedup: CUDA time / CUDA_4 time")
    plt.title("Four-GPU speedup relative to one GPU, mean ± propagated std")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "speedup_4gpu_vs_1gpu_with_error.png", dpi=200)
    plt.close()

    plt.figure()
    plt.errorbar(
        merged["Nx"],
        merged["parallel_efficiency"],
        yerr=merged["efficiency_std"],
        marker="o",
        capsize=4,
        label="Parallel efficiency",
    )
    plt.axhline(1.0, linestyle="--", label="Ideal efficiency")
    plt.xlabel("N for quadratic grid N × N")
    plt.ylabel("Efficiency = speedup / 4")
    plt.title("Four-GPU parallel efficiency, mean ± propagated std")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "parallel_efficiency_with_error.png", dpi=200)
    plt.close()

    merged.to_csv(outdir / "speedup_table_with_error.csv", index=False)


def plot_components_absolute(df: pd.DataFrame, outdir: Path) -> None:
    for backend, part in df.groupby("backend"):
        part = part.sort_values("Nx")

        available = [
            col for col in COMPONENTS
            if col in part.columns and part[col].notna().any()
        ]

        if not available:
            continue

        plt.figure(figsize=(11, 5))

        bottom = None
        for col in available:
            values = part[col].fillna(0.0)
            if bottom is None:
                plt.bar(
                    part["Nx"],
                    values,
                    width=BAR_WIDTH,
                    label=PRETTY_NAMES.get(col, col),
                )
                bottom = values.copy()
            else:
                plt.bar(
                    part["Nx"],
                    values,
                    width=BAR_WIDTH,
                    bottom=bottom,
                    label=PRETTY_NAMES.get(col, col),
                )
                bottom = bottom + values

        plt.xlabel("N for quadratic grid N × N")
        plt.ylabel("Component time [s]")
        plt.title(f"Timing components: {backend}")
        plt.legend(
            fontsize=11,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            borderaxespad=0.0,
        )
        plt.tight_layout(rect=[0, 0, 0.99, 1])
        plt.savefig(outdir / f"components_absolute_{backend}.png", dpi=200)
        plt.close()


def plot_components_percent(df: pd.DataFrame, outdir: Path) -> None:
    for backend, part in df.groupby("backend"):
        part = part.sort_values("Nx").copy()

        available = [
            col for col in COMPONENTS
            if col in part.columns and part[col].notna().any()
        ]

        if not available:
            continue

        for col in available:
            part[col] = part[col].fillna(0.0)

        total = part[available].sum(axis=1)
        total = total.replace(0.0, pd.NA)

        percent = part[available].div(total, axis=0) * 100.0

        plt.figure(figsize=(11, 5))

        bottom = None
        for col in available:
            values = percent[col].fillna(0.0)
            if bottom is None:
                plt.bar(
                    part["Nx"],
                    values,
                    width=BAR_WIDTH,
                    label=PRETTY_NAMES.get(col, col),
                )
                bottom = values.copy()
            else:
                plt.bar(
                    part["Nx"],
                    values,
                    width=BAR_WIDTH,
                    bottom=bottom,
                    label=PRETTY_NAMES.get(col, col),
                )
                bottom = bottom + values

        plt.xlabel("N for quadratic grid N × N")
        plt.ylabel("Share of measured components [%]")
        plt.title(f"Relative timing components: {backend}")
        plt.legend(
            fontsize=11,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            borderaxespad=0.0,
        )
        plt.tight_layout(rect=[0, 0, 0.99, 1])
        plt.savefig(outdir / f"components_percent_{backend}.png", dpi=200)
        plt.close()


def print_summary(df: pd.DataFrame) -> None:
    cols = [
        "backend",
        "gpus",
        "Nx",
        "cells",
        "n_successful_repeats",
        "solver_time_s",
        "solver_time_s_std",
        "solver_time_s_sem",
        "boundary_conditions_s",
        "halo_exchange_s",
        "fused_rhs_s",
        "time_integration_s",
        "positivity_correction_s",
    ]

    cols = [col for col in cols if col in df.columns]

    print()
    print("Averaged benchmark data:")
    print(df[cols].to_string(index=False))

    one = df[df["backend"] == "CUDA"][
        ["Nx", "solver_time_s", "solver_time_s_std"]
    ].rename(
        columns={
            "solver_time_s": "time_1gpu",
            "solver_time_s_std": "std_1gpu",
        }
    )

    four = df[df["backend"] == "CUDA_4"][
        ["Nx", "solver_time_s", "solver_time_s_std"]
    ].rename(
        columns={
            "solver_time_s": "time_4gpu",
            "solver_time_s_std": "std_4gpu",
        }
    )

    merged = one.merge(four, on="Nx", how="inner")

    if not merged.empty:
        merged["speedup"] = merged["time_1gpu"] / merged["time_4gpu"]
        merged["speedup_std"] = merged["speedup"] * (
            (merged["std_1gpu"] / merged["time_1gpu"]) ** 2
            + (merged["std_4gpu"] / merged["time_4gpu"]) ** 2
        ) ** 0.5
        merged["efficiency"] = merged["speedup"] / 4.0
        merged["efficiency_std"] = merged["speedup_std"] / 4.0

        print()
        print("Speedup table with propagated error:")
        print(merged.to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "csv",
        nargs="?",
        default="benchmarks/scaling/scaling_results.csv",
        help="Path to scaling_results.csv",
    )
    parser.add_argument(
        "--outdir",
        default="benchmarks/scaling/plots",
        help="Output directory for plots",
    )
    args = parser.parse_args()

    csv_path = Path(args.csv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df_raw = read_results(csv_path)
    df = aggregate(df_raw)

    print_summary(df)

    plot_solver_time(df, outdir)
    plot_solver_time_vs_n(df, outdir)
    plot_speedup(df, outdir)
    plot_components_absolute(df, outdir)
    plot_components_percent(df, outdir)

    print()
    print(f"Plots written to: {outdir}")


if __name__ == "__main__":
    main()