"""WP-10 campaign summary: OMP scaling, GPU crossover, precision penalty,
workload-scaling comparison and setup-vs-steady-state economics, computed
from wp10_2026-08-28's aggregate records. Consumes records; produces one
big JSON plus a human-readable text summary. Writes nothing back into the
aggregate set itself.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from analysis import campaign_summary as campaign_summary_mod
from analysis import crossover as crossover_mod
from analysis import throughput as throughput_mod
from harness import omp_scaling
from harness import precision_metrics


def load_aggregates(results_dir):
    agg_dir = pathlib.Path(results_dir) / "aggregates"
    return [json.loads(p.read_text()) for p in agg_dir.glob("*.json")]


def load_raw_records(results_dir):
    records = []
    for p in pathlib.Path(results_dir).glob("*.json"):
        if p.name == "driver_log.jsonl":
            continue
        records.append(json.loads(p.read_text()))
    return records


def workload_metadata_by_cell(records):
    """{(case_id, variant_id, size_id): {natom, directed_interactions,
    fft_grid_points, workload_class}} from any COMPLETED raw record of that
    cell (these fields are identical across every sample/config of one
    cell -- resolved once per cell by the driver)."""
    out = {}
    for r in records:
        if r.get("run_status") != "COMPLETED":
            continue
        key = (r["case_id"], r["variant_id"], r["size_id"])
        if key not in out:
            out[key] = {
                "natom": r["natom"],
                "directed_interactions": r.get("directed_interactions"),
                "fft_grid_points": r.get("fft_grid_points"),
                "workload_class": r["workload_class"],
            }
    return out


def cells_from_aggregates(aggregates):
    cells = set()
    for a in aggregates:
        cells.add((a["case_id"], a["variant_id"], a["size_id"]))
    return sorted(cells)


def per_cell_omp_scaling(aggregates, case_id, variant_id, size_id):
    cpu_steady = [
        a for a in aggregates
        if a["case_id"] == case_id and a["variant_id"] == variant_id and a["size_id"] == size_id
        and a["backend"] == "CPU" and a["metric"] == "steady_step_seconds"
    ]
    if not cpu_steady:
        return None
    table = omp_scaling.compute_omp_speedup_table(cpu_steady)
    cpu_1t = omp_scaling.select_cpu_1t(cpu_steady)
    cpu_best = omp_scaling.select_cpu_best(cpu_steady)
    return {"table": table, "cpu_1t": cpu_1t, "cpu_best": cpu_best}


def per_cell_gpu(aggregates, case_id, variant_id, size_id):
    out = {}
    for precision in ("SINGLE", "DOUBLE"):
        matches = [
            a for a in aggregates
            if a["case_id"] == case_id and a["variant_id"] == variant_id and a["size_id"] == size_id
            and a["backend"] == "GPU" and a["metric"] == "steady_step_seconds"
            and a["requested_precision"] == precision
        ]
        if matches:
            out[precision] = matches[0]
    return out


def per_cell_setup(aggregates, case_id, variant_id, size_id, backend, **extra):
    matches = [
        a for a in aggregates
        if a["case_id"] == case_id and a["variant_id"] == variant_id and a["size_id"] == size_id
        and a["backend"] == backend and a["metric"] == "setup_seconds"
        and all(a.get(k) == v for k, v in extra.items())
    ]
    return matches[0] if matches else None


def build_full_summary(results_dir):
    aggregates = load_aggregates(results_dir)
    records = load_raw_records(results_dir)
    wm_by_cell = workload_metadata_by_cell(records)
    cells = cells_from_aggregates(aggregates)

    per_cell = {}
    for case_id, variant_id, size_id in cells:
        wm = wm_by_cell.get((case_id, variant_id, size_id), {})
        omp = per_cell_omp_scaling(aggregates, case_id, variant_id, size_id)
        gpu = per_cell_gpu(aggregates, case_id, variant_id, size_id)

        setup_cpu_best = None
        setup_gpu = {}
        if omp is not None:
            setup_cpu_best = per_cell_setup(
                aggregates, case_id, variant_id, size_id, "CPU",
                omp_threads=omp["cpu_best"]["omp_threads"],
            )
        for precision, gpu_agg in gpu.items():
            setup_gpu[precision] = per_cell_setup(
                aggregates, case_id, variant_id, size_id, "GPU", requested_precision=precision,
            )

        crossovers = {}
        for precision, gpu_agg in gpu.items():
            if omp is None:
                continue
            s_bestcpu = crossover_mod.compute_gpu_speedup(
                omp["cpu_best"], gpu_agg, label="S_GPU_BESTCPU"
            )
            s_1t = crossover_mod.compute_gpu_speedup(
                omp["cpu_1t"], gpu_agg, label="S_GPU_1T"
            )
            crossovers[precision] = {"vs_cpu_best": s_bestcpu, "vs_cpu_1t": s_1t}

        throughput = {}
        if wm.get("workload_class") and omp is not None:
            throughput["cpu_best"] = throughput_mod.compute_throughput(
                natom=wm["natom"], t_step=omp["cpu_best"]["median"], workload_class=wm["workload_class"],
                directed_interactions=wm.get("directed_interactions"), fft_grid_points=wm.get("fft_grid_points"),
            )
        for precision, gpu_agg in gpu.items():
            if wm.get("workload_class"):
                throughput[f"gpu_{precision}"] = throughput_mod.compute_throughput(
                    natom=wm["natom"], t_step=gpu_agg["median"], workload_class=wm["workload_class"],
                    directed_interactions=wm.get("directed_interactions"), fft_grid_points=wm.get("fft_grid_points"),
                )

        precision_penalty = None
        if "SINGLE" in gpu and "DOUBLE" in gpu:
            precision_penalty = precision_metrics.compute_r_gpu_32_64([gpu["SINGLE"], gpu["DOUBLE"]])

        per_cell[f"{case_id}::{variant_id}::{size_id}"] = {
            "case_id": case_id, "variant_id": variant_id, "size_id": size_id,
            "workload_metadata": wm,
            "omp_scaling": omp,
            "gpu": gpu,
            "setup_cpu_best": setup_cpu_best,
            "setup_gpu": setup_gpu,
            "crossovers": crossovers,
            "throughput": throughput,
            "precision_penalty": precision_penalty,
        }

    # Per (case, variant) crossover-threshold summaries across sizes (section G/F).
    case_variant_summaries = {}
    by_case_variant = {}
    for key, cell in per_cell.items():
        cv_key = (cell["case_id"], cell["variant_id"])
        by_case_variant.setdefault(cv_key, []).append(cell)

    for (case_id, variant_id), cells_list in by_case_variant.items():
        entries = []
        for cell in cells_list:
            if cell["omp_scaling"] is None:
                continue
            for precision, gpu_agg in cell["gpu"].items():
                entries.append({
                    "size_id": cell["size_id"],
                    "natom": cell["workload_metadata"]["natom"],
                    "directed_interactions": cell["workload_metadata"].get("directed_interactions"),
                    "fft_grid_points": cell["workload_metadata"].get("fft_grid_points"),
                    "workload_class": cell["workload_metadata"]["workload_class"],
                    "requested_precision": precision,
                    "gpu_aggregate": gpu_agg,
                    "cpu_best_aggregate": cell["omp_scaling"]["cpu_best"],
                    "cpu_1t_aggregate": cell["omp_scaling"]["cpu_1t"],
                })
        if entries:
            try:
                case_variant_summaries[f"{case_id}::{variant_id}"] = campaign_summary_mod.build_case_summary(entries)
            except campaign_summary_mod.CampaignSummaryError as exc:
                case_variant_summaries[f"{case_id}::{variant_id}"] = {"error": str(exc)}

    return {"per_cell": per_cell, "case_variant_summaries": case_variant_summaries}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    summary = build_full_summary(args.results_dir)
    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"wrote {out_path}, {len(summary['per_cell'])} cells, "
          f"{len(summary['case_variant_summaries'])} case/variant summaries")
