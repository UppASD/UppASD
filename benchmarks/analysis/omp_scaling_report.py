"""OpenMP scaling reports: machine-readable summaries and plots (WP-05
section H).

Consumes `aggregate` records through `harness.omp_scaling` and produces
understanding; this module never re-runs a measurement and never produces a
new result-record kind of its own (see `benchmarks/analysis/README.md`'s
directory responsibility: "Consumes records; never produces them").

Callers supply one "size family" per problem size: the `aggregate` records of
that size's OpenMP scaling sweep (varying only `omp_threads`) plus the
`natom` workload size, since `aggregate` records themselves do not carry
`natom` (only raw sample records do) -- the caller reads it off any of that
size's own raw sample records or the case's size ladder.
"""

from __future__ import annotations

import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (backend must be selected first)

from harness import omp_scaling


class OmpScalingReportError(ValueError):
    """A scaling report could not be built from the given size families."""


def build_scaling_summary(size_families, *, require_environment_valid=False):
    """Build the machine-readable OpenMP-scaling summary (blueprint section H).

    `size_families` is an iterable of `{"size_id": ..., "natom": ...,
    "aggregates": [...]}`, one entry per problem size, where `aggregates` is
    that size's OpenMP-scaling family (blueprint section E/F, exactly what
    `omp_scaling.compute_omp_speedup_table` expects). Returns a plain dict:
    per-size speedup/efficiency rows and CPU-1T/CPU-BEST, plus one
    `p_best_vs_natom` table across every size, sorted by `natom` -- never
    assuming size ladders are supplied in size order.
    """
    size_families = list(size_families)
    if not size_families:
        raise OmpScalingReportError("no size families given")

    per_size = []
    p_best_vs_natom = []
    for family in size_families:
        aggregates = family["aggregates"]
        rows = omp_scaling.compute_omp_speedup_table(
            aggregates, require_environment_valid=require_environment_valid
        )
        cpu_1t = omp_scaling.select_cpu_1t(aggregates, require_environment_valid=require_environment_valid)
        cpu_best = omp_scaling.determine_p_best(aggregates, require_environment_valid=require_environment_valid)
        per_size.append({
            "size_id": family["size_id"],
            "natom": family["natom"],
            "rows": rows,
            "cpu_1t_aggregate_id": cpu_1t["aggregate_id"],
            "cpu_best_aggregate_id": cpu_best["aggregate_id"],
            "p_best": cpu_best["omp_threads"],
            "t_best": cpu_best["median"],
        })
        p_best_vs_natom.append({
            "size_id": family["size_id"], "natom": family["natom"], "p_best": cpu_best["omp_threads"],
        })

    p_best_vs_natom.sort(key=lambda row: row["natom"])
    return {"sizes": per_size, "p_best_vs_natom": p_best_vs_natom}


def write_scaling_summary(summary, path):
    """Write `build_scaling_summary`'s output as indented, sorted-key JSON."""
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return path


def plot_omp_scaling(summary, output_dir):
    """Render the blueprint section H plots as PNGs under `output_dir`:
    one step-time/speedup/efficiency-vs-threads figure per size, plus one
    p_best-vs-N figure across every size. Returns the list of written paths.
    """
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    written = []

    for size in summary["sizes"]:
        threads = [row["omp_threads"] for row in size["rows"]]
        t_p = [row["t_p"] for row in size["rows"]]
        speedup = [row["speedup"] for row in size["rows"]]
        efficiency = [row["efficiency"] for row in size["rows"]]

        figure, axes = plt.subplots(1, 3, figsize=(15, 4))
        axes[0].plot(threads, t_p, marker="o")
        axes[0].set_xlabel("OMP threads")
        axes[0].set_ylabel("t_p (s)")
        axes[0].set_title(f"{size['size_id']}: step time vs threads")

        axes[1].plot(threads, speedup, marker="o", label="measured")
        axes[1].plot(threads, threads, linestyle="--", color="gray", label="ideal (S=p)")
        axes[1].set_xlabel("OMP threads")
        axes[1].set_ylabel("S_OMP(p)")
        axes[1].set_title("OpenMP speedup")
        axes[1].legend()

        axes[2].plot(threads, efficiency, marker="o")
        axes[2].axhline(1.0, linestyle="--", color="gray")
        axes[2].set_xlabel("OMP threads")
        axes[2].set_ylabel("E_OMP(p)")
        axes[2].set_title("OpenMP efficiency")

        figure.tight_layout()
        path = output_dir / f"omp_scaling__{size['size_id']}.png"
        figure.savefig(path)
        plt.close(figure)
        written.append(path)

    natoms = [row["natom"] for row in summary["p_best_vs_natom"]]
    p_bests = [row["p_best"] for row in summary["p_best_vs_natom"]]
    figure, axis = plt.subplots(figsize=(5, 4))
    axis.plot(natoms, p_bests, marker="o")
    axis.set_xlabel("N_atom")
    axis.set_ylabel("p_best")
    axis.set_title("CPU-BEST thread count vs problem size")
    figure.tight_layout()
    path = output_dir / "p_best_vs_natom.png"
    figure.savefig(path)
    plt.close(figure)
    written.append(path)

    return written
