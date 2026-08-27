"""GPU-crossover campaign plots (WP-07 section E, plots 1-2 and 5-9).

Plots 3 ("OpenMP scaling") and 4 ("OpenMP efficiency") are already produced
by `analysis.omp_scaling_report.plot_omp_scaling` (WP-05) -- this module does
not reproduce them, it only adds the plots that need
`analysis.campaign_summary.build_case_summary`'s GPU-crossover output:

1. steady step time versus `N_atom` (`plot_step_time_vs_natom`);
2. steady step time versus `N_directed` for neighbour workloads
   (`plot_step_time_vs_ndirected`);
5. GPU/BESTCPU speedup versus `N_atom` (`plot_speedup_vs_natom`);
6. GPU/BESTCPU speedup versus `N_directed` (`plot_speedup_vs_ndirected`);
7. SINGLE versus DOUBLE GPU performance (`plot_single_vs_double`);
8. crossover summary (`plot_crossover_summary`);
9. FFT/dipole time versus grid size (`plot_fft_vs_grid_size`).

Every function takes one case/variant's `build_case_summary` output and an
output directory, and returns the path it wrote (or `None` when the summary
has nothing relevant to plot -- e.g. `plot_fft_vs_grid_size` on a pure
neighbour-list case). `plot_campaign_report` runs all of them and returns the
non-`None` paths, mirroring `omp_scaling_report.plot_omp_scaling`'s
"list of written paths" contract.

Median values are plotted with MAD error bars wherever the summary carries
one (blueprint section F: "Carry MAD or an equivalent dispersion indication
into plots/tables"); axes are logarithmic by default since problem size spans
orders of magnitude.
"""

from __future__ import annotations

import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (backend must be selected first)

from analysis import crossover as crossover_mod

_PRECISION_COLORS = {"SINGLE": "tab:orange", "DOUBLE": "tab:blue"}
_THRESHOLD_STYLES = {1.0: ":", 1.25: "--", 2.0: "-.", 5.0: (0, (1, 1))}


def _all_rows(summary):
    rows = []
    for block in summary["by_precision"].values():
        rows.extend(block["rows"])
    return rows


def _rows_by_precision(summary, precision):
    block = summary["by_precision"].get(precision)
    return block["rows"] if block else []


def _identity_slug(summary):
    identity = summary["identity"]
    return f"{identity['case_id']}__{identity['variant_id']}"


def plot_step_time_vs_natom(summary, output_dir):
    """Plot 1: steady-state step time (GPU, CPU-BEST, CPU-1T) versus `N_atom`."""
    rows = sorted(_all_rows(summary), key=lambda row: row["natom"])
    if not rows:
        return None
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    figure, axis = plt.subplots(figsize=(6, 4.5))
    for precision in sorted({row["requested_precision"] for row in rows}):
        precision_rows = [row for row in rows if row["requested_precision"] == precision]
        axis.errorbar(
            [row["natom"] for row in precision_rows],
            [row["t_gpu_seconds"] for row in precision_rows],
            yerr=[row.get("t_gpu_mad") or 0 for row in precision_rows],
            marker="o", label=f"GPU {precision}", color=_PRECISION_COLORS.get(precision),
        )
    reference_rows = _rows_by_precision(summary, sorted(summary["by_precision"])[0])
    axis.errorbar(
        [row["natom"] for row in reference_rows],
        [row["t_cpu_best_seconds"] for row in reference_rows],
        yerr=[row.get("t_cpu_best_mad") or 0 for row in reference_rows],
        marker="s", label="CPU-BEST", color="black",
    )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("N_atom")
    axis.set_ylabel("steady step time (s)")
    axis.set_title(f"{_identity_slug(summary)}: step time vs N_atom")
    axis.legend()
    figure.tight_layout()
    path = output_dir / f"step_time_vs_natom__{_identity_slug(summary)}.png"
    figure.savefig(path)
    plt.close(figure)
    return path


def plot_step_time_vs_ndirected(summary, output_dir):
    """Plot 2: steady-state step time versus `N_directed`, neighbour workloads only."""
    rows = sorted(
        (row for row in _all_rows(summary) if row.get("directed_interactions") is not None),
        key=lambda row: row["directed_interactions"],
    )
    if not rows:
        return None
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    figure, axis = plt.subplots(figsize=(6, 4.5))
    for precision in sorted({row["requested_precision"] for row in rows}):
        precision_rows = [row for row in rows if row["requested_precision"] == precision]
        axis.errorbar(
            [row["directed_interactions"] for row in precision_rows],
            [row["t_gpu_seconds"] for row in precision_rows],
            yerr=[row.get("t_gpu_mad") or 0 for row in precision_rows],
            marker="o", label=f"GPU {precision}", color=_PRECISION_COLORS.get(precision),
        )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("N_directed (directed interactions)")
    axis.set_ylabel("steady step time (s)")
    axis.set_title(f"{_identity_slug(summary)}: step time vs N_directed")
    axis.legend()
    figure.tight_layout()
    path = output_dir / f"step_time_vs_ndirected__{_identity_slug(summary)}.png"
    figure.savefig(path)
    plt.close(figure)
    return path


def _plot_speedup_vs(summary, output_dir, *, x_field, x_label, filename_stub):
    rows = sorted(
        (row for row in _all_rows(summary) if row.get(x_field) is not None),
        key=lambda row: row[x_field],
    )
    if not rows:
        return None
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    figure, axis = plt.subplots(figsize=(6, 4.5))
    for precision in sorted({row["requested_precision"] for row in rows}):
        precision_rows = [row for row in rows if row["requested_precision"] == precision]
        axis.plot(
            [row[x_field] for row in precision_rows],
            [row["s_gpu_bestcpu"] for row in precision_rows],
            marker="o", label=f"S_GPU/BESTCPU {precision}", color=_PRECISION_COLORS.get(precision),
        )
    for threshold, style in _THRESHOLD_STYLES.items():
        axis.axhline(threshold, linestyle=style, color="gray", linewidth=1, label=f"{threshold}x")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel(x_label)
    axis.set_ylabel("S_GPU/BESTCPU")
    axis.set_title(f"{_identity_slug(summary)}: GPU/BESTCPU speedup vs {x_label}")
    axis.legend(fontsize="small")
    figure.tight_layout()
    path = output_dir / f"{filename_stub}__{_identity_slug(summary)}.png"
    figure.savefig(path)
    plt.close(figure)
    return path


def plot_speedup_vs_natom(summary, output_dir):
    """Plot 5: GPU/BESTCPU speedup versus `N_atom`."""
    return _plot_speedup_vs(
        summary, output_dir, x_field="natom", x_label="N_atom", filename_stub="speedup_vs_natom"
    )


def plot_speedup_vs_ndirected(summary, output_dir):
    """Plot 6: GPU/BESTCPU speedup versus `N_directed`."""
    return _plot_speedup_vs(
        summary, output_dir, x_field="directed_interactions", x_label="N_directed",
        filename_stub="speedup_vs_ndirected",
    )


def plot_single_vs_double(summary, output_dir):
    """Plot 7: SINGLE versus DOUBLE GPU step time versus `N_atom`."""
    single_rows = sorted(_rows_by_precision(summary, "SINGLE"), key=lambda row: row["natom"])
    double_rows = sorted(_rows_by_precision(summary, "DOUBLE"), key=lambda row: row["natom"])
    if not single_rows or not double_rows:
        return None
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    figure, axis = plt.subplots(figsize=(6, 4.5))
    axis.errorbar(
        [row["natom"] for row in single_rows], [row["t_gpu_seconds"] for row in single_rows],
        yerr=[row.get("t_gpu_mad") or 0 for row in single_rows],
        marker="o", label="GPU SINGLE", color=_PRECISION_COLORS["SINGLE"],
    )
    axis.errorbar(
        [row["natom"] for row in double_rows], [row["t_gpu_seconds"] for row in double_rows],
        yerr=[row.get("t_gpu_mad") or 0 for row in double_rows],
        marker="o", label="GPU DOUBLE", color=_PRECISION_COLORS["DOUBLE"],
    )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("N_atom")
    axis.set_ylabel("steady step time (s)")
    axis.set_title(f"{_identity_slug(summary)}: SINGLE vs DOUBLE GPU performance")
    axis.legend()
    figure.tight_layout()
    path = output_dir / f"single_vs_double__{_identity_slug(summary)}.png"
    figure.savefig(path)
    plt.close(figure)
    return path


def plot_crossover_summary(summary, output_dir):
    """Plot 8: crossover-threshold summary, one row of markers per precision.

    Interpolated crossings are drawn as filled circles at their (rounded)
    `crossover_natom`; a threshold classified `below_tested_range` or
    `above_tested_range` is drawn as an open triangle at the tested-range edge
    it falls outside of, so the plot never implies a measured point where none
    exists.
    """
    precisions = sorted(summary["by_precision"])
    if not precisions:
        return None
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    figure, axis = plt.subplots(figsize=(6, 2 + len(precisions)))
    for row_index, precision in enumerate(precisions):
        block = summary["by_precision"][precision]
        for result in block["crossover"].values():
            threshold = result["threshold"]
            if result["status"] == crossover_mod.WITHIN_TESTED_RANGE:
                marker = "o" if result["interpolated"] else "D"
                axis.scatter(result["crossover_natom"], row_index, marker=marker, color="black", zorder=3)
                axis.annotate(
                    f"{threshold}x", (result["crossover_natom"], row_index),
                    textcoords="offset points", xytext=(0, 8), ha="center", fontsize="small",
                )
            else:
                edge = result["tested_range"][0] if result["status"] == crossover_mod.BELOW_TESTED_RANGE \
                    else result["tested_range"][1]
                axis.scatter(edge, row_index, marker="^", facecolors="none", edgecolors="gray", zorder=3)
                axis.annotate(
                    f"{threshold}x {result['status']}", (edge, row_index),
                    textcoords="offset points", xytext=(0, 8), ha="center", fontsize="x-small", color="gray",
                )
    axis.set_xscale("log")
    axis.set_yticks(range(len(precisions)))
    axis.set_yticklabels(precisions)
    axis.set_xlabel("N_atom")
    axis.set_title(f"{_identity_slug(summary)}: crossover summary")
    figure.tight_layout()
    path = output_dir / f"crossover_summary__{_identity_slug(summary)}.png"
    figure.savefig(path)
    plt.close(figure)
    return path


def plot_fft_vs_grid_size(summary, output_dir):
    """Plot 9: FFT/dipole step time versus `fft_grid_points`."""
    rows = sorted(
        (row for row in _all_rows(summary) if row.get("fft_grid_points") is not None),
        key=lambda row: row["fft_grid_points"],
    )
    if not rows:
        return None
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    figure, axis = plt.subplots(figsize=(6, 4.5))
    for precision in sorted({row["requested_precision"] for row in rows}):
        precision_rows = [row for row in rows if row["requested_precision"] == precision]
        axis.errorbar(
            [row["fft_grid_points"] for row in precision_rows],
            [row["t_gpu_seconds"] for row in precision_rows],
            yerr=[row.get("t_gpu_mad") or 0 for row in precision_rows],
            marker="o", label=f"GPU {precision}", color=_PRECISION_COLORS.get(precision),
        )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("FFT grid points")
    axis.set_ylabel("steady step time (s)")
    axis.set_title(f"{_identity_slug(summary)}: FFT/dipole time vs grid size")
    axis.legend()
    figure.tight_layout()
    path = output_dir / f"fft_time_vs_grid__{_identity_slug(summary)}.png"
    figure.savefig(path)
    plt.close(figure)
    return path


_PLOT_BUILDERS = (
    plot_step_time_vs_natom,
    plot_step_time_vs_ndirected,
    plot_speedup_vs_natom,
    plot_speedup_vs_ndirected,
    plot_single_vs_double,
    plot_crossover_summary,
    plot_fft_vs_grid_size,
)


def plot_campaign_report(summary, output_dir):
    """Run every plot builder in this module; return the non-`None` paths written."""
    written = [builder(summary, output_dir) for builder in _PLOT_BUILDERS]
    return [path for path in written if path is not None]
