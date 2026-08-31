"""Pair-interaction vs. GPU-convolution comparison plots.

Consumes the aggregate records `harness/convolution_driver.py` writes: for
every GPU (case, variant, size, precision) cell whose case passes
`harness.gpu_convolution.check_convolution_gate`, that driver writes two
`steady_step_seconds` aggregates sharing every `harness/aggregate.py`
`_CELL_FIELDS` value -- the only difference is a `__pair__`/`__conv__`
substring in `aggregate_id` itself, since the schema's aggregate cell
fields have no axis for "which Hamiltonian evaluator" (unlike
`requested_precision`, which is a real field). This module's loaders
resolve that split from `aggregate_id`, not from a schema field.

Three figures per case/variant:

1. `plot_step_time_comparison` -- grouped bars, CPU-BEST plus GPU
   SINGLE/DOUBLE x pair/conv, one group per size (log y-axis: convolution
   can be an order of magnitude away from pair, see conv_bench.txt's own
   dhcp_nd_long_range point, ~13x Hamiltonian-only).
2. `plot_convolution_speedup` -- conv-vs-pair speedup (`t_pair/t_conv`)
   and conv-vs-CPU-BEST speedup, one bar pair per (size, precision).
3. `plot_omp_scaling_with_conv_baselines` -- the CPU 1-8 thread sweep's
   step-time panel (as `analysis/wp10_plots.py`'s
   `plot_omp_scaling_with_gpu_baseline` already draws for WP-10) with all
   four GPU baselines overlaid: pair SINGLE/DOUBLE and conv SINGLE/DOUBLE.

Writes PNGs under `--out-dir` (gitignored, same convention as every other
generated analysis product in this project) plus a small
`convolution_plots_manifest.json`.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (backend must be selected first)

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from analysis import omp_scaling_report
from analysis import wp10_summary
from harness import omp_scaling

_PRECISION_COLORS = {"SINGLE": "tab:orange", "DOUBLE": "tab:blue"}
_MODE_STYLES = {"pair": "-", "conv": "--"}


def _mode_of(aggregate_id):
    if "__pair__" in aggregate_id:
        return "pair"
    if "__conv__" in aggregate_id:
        return "conv"
    return None


def load_gpu_modes_by_cell(aggregates):
    """{(case_id, variant_id, size_id): {precision: {mode: aggregate}}} for
    every GPU steady_step_seconds aggregate that carries a pair/conv tag
    (aggregates without one -- e.g. an ordinary WP-10-style campaign with
    no convolution arm at all -- are simply absent here, not errored)."""
    out = {}
    for a in aggregates:
        if a["backend"] != "GPU" or a["metric"] != "steady_step_seconds":
            continue
        mode = _mode_of(a["aggregate_id"])
        if mode is None:
            continue
        key = (a["case_id"], a["variant_id"], a["size_id"])
        out.setdefault(key, {}).setdefault(a["requested_precision"], {})[mode] = a
    return out


def cpu_best_by_cell(aggregates):
    """{(case_id, variant_id, size_id): cpu_best_aggregate} via
    harness.omp_scaling.select_cpu_best over each cell's CPU thread sweep."""
    cells = {}
    for a in aggregates:
        if a["backend"] != "CPU" or a["metric"] != "steady_step_seconds":
            continue
        key = (a["case_id"], a["variant_id"], a["size_id"])
        cells.setdefault(key, []).append(a)
    return {
        key: omp_scaling.select_cpu_best(aggs)
        for key, aggs in cells.items()
        if aggs
    }


def plot_step_time_comparison(gpu_modes_by_cell, cpu_best, case_id, variant_id, size_ids, output_dir):
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    slug = f"{case_id}__{variant_id}"

    bar_specs = [("CPU-BEST", None, None)]
    for precision in ("SINGLE", "DOUBLE"):
        for mode in ("pair", "conv"):
            bar_specs.append((f"GPU {precision} {mode}", precision, mode))

    x = range(len(size_ids))
    width = 0.8 / len(bar_specs)
    figure, axis = plt.subplots(figsize=(7, 4.5))
    any_data = False
    for spec_index, (label, precision, mode) in enumerate(bar_specs):
        values = []
        for size_id in size_ids:
            key = (case_id, variant_id, size_id)
            if precision is None:
                agg = cpu_best.get(key)
            else:
                agg = gpu_modes_by_cell.get(key, {}).get(precision, {}).get(mode)
            values.append(agg["median"] if agg else None)
        if all(v is None for v in values):
            continue
        any_data = True
        color = "black" if precision is None else _PRECISION_COLORS[precision]
        hatch = None if mode != "conv" else "//"
        positions = [xi + spec_index * width for xi in x]
        axis.bar(
            positions, [v or 0 for v in values], width=width, label=label,
            color=color, alpha=0.5 if mode == "conv" else 1.0, hatch=hatch,
        )
    if not any_data:
        plt.close(figure)
        return None
    axis.set_yscale("log")
    axis.set_xticks([xi + 0.4 for xi in x])
    axis.set_xticklabels(size_ids)
    axis.set_ylabel("steady step time (s)")
    axis.set_title(f"{slug}: step time, CPU-BEST vs GPU pair/conv")
    axis.legend(fontsize="x-small", ncol=2)
    figure.tight_layout()
    path = output_dir / f"step_time_comparison__{slug}.png"
    figure.savefig(path)
    plt.close(figure)
    return path


def plot_convolution_speedup(gpu_modes_by_cell, cpu_best, case_id, variant_id, size_ids, output_dir):
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    slug = f"{case_id}__{variant_id}"

    rows = []
    for size_id in size_ids:
        key = (case_id, variant_id, size_id)
        cpu_agg = cpu_best.get(key)
        for precision in ("SINGLE", "DOUBLE"):
            modes = gpu_modes_by_cell.get(key, {}).get(precision, {})
            pair_agg, conv_agg = modes.get("pair"), modes.get("conv")
            if not pair_agg or not conv_agg:
                continue
            rows.append({
                "size_id": size_id, "precision": precision,
                "conv_vs_pair": pair_agg["median"] / conv_agg["median"],
                "conv_vs_cpu_best": (cpu_agg["median"] / conv_agg["median"]) if cpu_agg else None,
            })
    if not rows:
        return None

    figure, axis = plt.subplots(figsize=(7, 4.5))
    x = range(len(rows))
    axis.bar(x, [row["conv_vs_pair"] for row in rows], width=0.35, label="conv vs GPU pair", color="tab:green")
    if all(row["conv_vs_cpu_best"] is not None for row in rows):
        axis.bar(
            [xi + 0.35 for xi in x], [row["conv_vs_cpu_best"] for row in rows],
            width=0.35, label="conv vs CPU-BEST", color="tab:purple",
        )
    axis.axhline(1.0, linestyle=":", color="gray", linewidth=1)
    axis.set_xticks([xi + 0.175 for xi in x])
    axis.set_xticklabels([f"{row['size_id']}\n{row['precision']}" for row in rows], fontsize="small")
    axis.set_ylabel("speedup (x)")
    axis.set_title(f"{slug}: GPU convolution speedup")
    axis.legend(fontsize="small")
    figure.tight_layout()
    path = output_dir / f"convolution_speedup__{slug}.png"
    figure.savefig(path)
    plt.close(figure)
    return path


def plot_omp_scaling_with_conv_baselines(omp_summary, gpu_modes_by_cell, case_id, variant_id, output_dir):
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    written = []

    for size in omp_summary["sizes"]:
        threads = [row["omp_threads"] for row in size["rows"]]
        t_p = [row["t_p"] for row in size["rows"]]
        speedup = [row["speedup"] for row in size["rows"]]
        efficiency = [row["efficiency"] for row in size["rows"]]
        modes = gpu_modes_by_cell.get((case_id, variant_id, size["size_id"]), {})

        figure, axes = plt.subplots(1, 3, figsize=(15, 4))
        axes[0].plot(threads, t_p, marker="o", label="CPU", color="black")
        for precision in ("DOUBLE", "SINGLE"):
            for mode in ("pair", "conv"):
                agg = modes.get(precision, {}).get(mode)
                if agg is None:
                    continue
                axes[0].axhline(
                    agg["median"], linestyle=_MODE_STYLES[mode], color=_PRECISION_COLORS[precision],
                    label=f"GPU {precision} {mode}",
                )
        axes[0].set_yscale("log")
        axes[0].set_xlabel("OMP threads")
        axes[0].set_ylabel("t_p (s)")
        axes[0].set_title(f"{size['size_id']}: step time vs threads (+ GPU pair/conv baselines)")
        axes[0].legend(fontsize="x-small")

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
        path = output_dir / f"omp_scaling_with_conv__{size['size_id']}.png"
        figure.savefig(path)
        plt.close(figure)
        written.append(path)

    return written


def generate_all_plots(results_dir, out_dir):
    aggregates = wp10_summary.load_aggregates(results_dir)
    records = wp10_summary.load_raw_records(results_dir)
    wm_by_cell = wp10_summary.workload_metadata_by_cell(records)

    gpu_modes_by_cell = load_gpu_modes_by_cell(aggregates)
    cpu_best = cpu_best_by_cell(aggregates)
    out_dir = pathlib.Path(out_dir)
    manifest = {"comparison": {}, "speedup": {}, "omp_with_conv": {}, "skipped": []}

    case_variants = sorted({(a["case_id"], a["variant_id"]) for a in aggregates})
    for case_id, variant_id in case_variants:
        slug = f"{case_id}__{variant_id}"
        size_ids = sorted({
            a["size_id"] for a in aggregates
            if a["case_id"] == case_id and a["variant_id"] == variant_id
        })

        try:
            path = plot_step_time_comparison(gpu_modes_by_cell, cpu_best, case_id, variant_id, size_ids, out_dir)
            if path:
                manifest["comparison"][slug] = str(path)
        except Exception as exc:  # noqa: BLE001
            manifest["skipped"].append({"slug": slug, "stage": "comparison", "error": str(exc)})

        try:
            path = plot_convolution_speedup(gpu_modes_by_cell, cpu_best, case_id, variant_id, size_ids, out_dir)
            if path:
                manifest["speedup"][slug] = str(path)
        except Exception as exc:  # noqa: BLE001
            manifest["skipped"].append({"slug": slug, "stage": "speedup", "error": str(exc)})

        families = []
        for size_id in size_ids:
            cpu_aggs = [
                a for a in aggregates
                if a["case_id"] == case_id and a["variant_id"] == variant_id and a["size_id"] == size_id
                and a["backend"] == "CPU" and a["metric"] == "steady_step_seconds"
            ]
            wm = wm_by_cell.get((case_id, variant_id, size_id))
            if cpu_aggs and wm:
                families.append({"size_id": size_id, "natom": wm["natom"], "aggregates": cpu_aggs})
        if families:
            try:
                omp_summary = omp_scaling_report.build_scaling_summary(families)
                paths = plot_omp_scaling_with_conv_baselines(
                    omp_summary, gpu_modes_by_cell, case_id, variant_id, out_dir / "omp_scaling" / slug
                )
                manifest["omp_with_conv"][slug] = [str(p) for p in paths]
            except Exception as exc:  # noqa: BLE001
                manifest["skipped"].append({"slug": slug, "stage": "omp_with_conv", "error": str(exc)})

    manifest_path = out_dir / "convolution_plots_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    manifest = generate_all_plots(args.results_dir, args.out_dir)
    print(f"wrote {len(manifest['comparison'])} step-time-comparison PNGs")
    print(f"wrote {len(manifest['speedup'])} convolution-speedup PNGs")
    n_omp = sum(len(v) for v in manifest["omp_with_conv"].values())
    print(f"wrote {n_omp} OMP-scaling-with-conv-baseline PNGs")
    if manifest["skipped"]:
        print(f"skipped {len(manifest['skipped'])} plots:")
        for entry in manifest["skipped"]:
            print(" ", entry)
