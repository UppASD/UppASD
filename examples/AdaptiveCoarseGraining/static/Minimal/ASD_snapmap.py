#!/usr/bin/env python3
"""Plot a 2D minimap of the z-component of UppASD moments."""

import argparse
import glob
import os

import matplotlib.pyplot as plt
import numpy as np


def find_single_file(pattern, override=None):
    """Return the first file matching a glob pattern or an explicit override."""
    if override:
        return override

    matches = sorted(glob.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No file matches '{pattern}'")
    return matches[0]


def load_coordinates(path):
    """Load x, y, z coordinates from a coord output file."""
    coord = np.loadtxt(path, usecols=(1, 2, 3), comments="#")
    return np.atleast_2d(coord)


def load_moments(path, ensemble=1):
    """Load moment vectors for the selected ensemble at the final iteration."""
    restart = np.loadtxt(path, comments="#")
    restart = np.atleast_2d(restart)

    if restart.shape[1] < 7:
        raise ValueError(
            "Restart file must contain at least 7 columns: iter, ensemble, atom, |M|, Mx, My, Mz"
        )

    ensemble_mask = restart[:, 1].astype(int) == ensemble
    if not np.any(ensemble_mask):
        raise ValueError(f"No rows found for ensemble {ensemble} in {path}")

    selected = restart[ensemble_mask]
    last_iteration = np.max(selected[:, 0])
    selected = selected[selected[:, 0] == last_iteration]
    atom_index = selected[:, 2].astype(int)
    order = np.argsort(atom_index)
    return selected[order, 4:7]


def choose_plot_plane(coord):
    """Pick the two coordinate axes with the largest spatial extent."""
    spans = np.ptp(coord, axis=0)
    axes = np.argsort(spans)[-2:]
    axes.sort()
    labels = ["x", "y", "z"]
    return axes, labels[axes[0]], labels[axes[1]]


def marker_size(x, y):
    """Choose a square marker size based on the point spacing."""
    unique_x = np.unique(np.round(x, 10))
    unique_y = np.unique(np.round(y, 10))
    spacings = []

    if unique_x.size > 1:
        dx = np.diff(unique_x)
        spacings.extend(dx[dx > 0])
    if unique_y.size > 1:
        dy = np.diff(unique_y)
        spacings.extend(dy[dy > 0])

    if not spacings:
        return 120.0

    min_spacing = float(np.min(spacings))
    return max(90.0, min(450.0, 6000.0 * min_spacing))


def build_title(jobname, plane_x, plane_y, ensemble):
    """Build a descriptive plot title for the selected output set."""
    title = f"$m_z$ map ({plane_x}{plane_y}-plane)"
    if jobname:
        title += f"\n{jobname}"
    if ensemble != 1:
        title += f", ensemble {ensemble}"
    return title


def main():
    """Parse arguments, load UppASD output, and save the minimap image."""
    parser = argparse.ArgumentParser(
        description="Plot the z-component of moments from UppASD coord/restart output files."
    )
    parser.add_argument("--coord", help="Path to coord.*.out file")
    parser.add_argument("--restart", help="Path to restart.*.out file")
    parser.add_argument("--ensemble", type=int, default=1, help="Ensemble index to plot")
    parser.add_argument(
        "--output",
        help="Output PNG path (default: mzmap.<jobname>.png next to the restart file)",
    )
    args = parser.parse_args()

    coordfile = find_single_file("coord.*.out", args.coord)
    restartfile = find_single_file("restart.*.out", args.restart)

    coord = load_coordinates(coordfile)
    mom = load_moments(restartfile, ensemble=args.ensemble)

    if coord.shape[0] != mom.shape[0]:
        raise ValueError(
            f"Coordinate count ({coord.shape[0]}) does not match moment count ({mom.shape[0]})"
        )

    axes, plane_x, plane_y = choose_plot_plane(coord)
    x = coord[:, axes[0]]
    y = coord[:, axes[1]]
    mz = mom[:, 2]

    basename = os.path.basename(restartfile)
    parts = basename.split(".")
    jobname = parts[1] if len(parts) > 2 else os.path.splitext(basename)[0]
    output = args.output or f"mzmap.{jobname}.png"

    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.labelsize": 13,
            "axes.titlesize": 15,
            "figure.dpi": 140,
        }
    )

    fig, ax = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
    scatter = ax.scatter(
        x,
        y,
        c=mz,
        s=marker_size(x, y),
        cmap="coolwarm",
        vmin=-1.0,
        vmax=1.0,
        marker="s",
        edgecolors="black",
        linewidths=0.15,
    )

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(fr"$r_{plane_x}$")
    ax.set_ylabel(fr"$r_{plane_y}$")
    ax.set_title(build_title(jobname, plane_x, plane_y, args.ensemble))
    ax.set_facecolor("#f4f1ea")
    ax.grid(color="white", linewidth=0.6)

    cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label(r"$m_z$")

    fig.savefig(output, bbox_inches="tight")
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()

