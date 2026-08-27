"""Machine-readable GPU-crossover campaign summary (WP-07 section G).

`build_case_summary` is the single entry point: given one case/variant's
measured cells across sizes and GPU precisions, it produces a plain dict
covering everything section G asks for --

* CPU-BEST (per size, from `harness.omp_scaling`-selected aggregates the
  caller passes in);
* crossover thresholds (`analysis.crossover`, 1.0x/1.25x/2.0x/5.0x, never
  extrapolated);
* throughput (`analysis.throughput`, per size and precision);
* precision penalty (`harness.precision_metrics.compute_r_gpu_32_64`, per size
  where both SINGLE and DOUBLE were measured);
* measured range (`natom` min/max actually tested, per precision);
* quality status (`authoritative` and the union of quality flags across every
  aggregate the summary draws on).

`write_summary_json`/`write_summary_csv` render it to disk. The CSV is the
flat per-(size, precision) table (one row per measured cell); crossover
thresholds and precision penalties, which are per-precision/per-size-pair
rather than per-row facts, get their own small CSVs
(`write_crossover_csv`/`write_precision_penalty_csv`) rather than being forced
into the same table.

Callers assemble one `entries` list per (campaign, case, variant, machine):
each entry describes one measured (size, requested_precision) GPU cell --
`size_id`, `natom`, `workload_class`, `directed_interactions`/
`fft_grid_points` where applicable, `requested_precision`, and the three
`steady_step_seconds` aggregates that cell needs (`gpu_aggregate`,
`cpu_best_aggregate`, `cpu_1t_aggregate` -- the same CPU-BEST/CPU-1T
aggregates `harness.omp_scaling.determine_p_best`/`select_cpu_1t` select for
that size). This module never selects CPU-BEST/CPU-1T itself and never runs a
measurement; it only combines records callers already have.
"""

from __future__ import annotations

import csv
import json
import pathlib

from analysis import crossover as crossover_mod
from analysis import throughput as throughput_mod
from harness import precision_metrics

_REQUIRED_ENTRY_FIELDS = (
    "size_id", "natom", "workload_class", "requested_precision",
    "gpu_aggregate", "cpu_best_aggregate", "cpu_1t_aggregate",
)

_CASE_IDENTITY_FIELDS = ("campaign_id", "case_id", "variant_id", "machine_id")


class CampaignSummaryError(ValueError):
    """A campaign summary could not be built from the given entries."""


def _require_entry_fields(entry):
    missing = [
        field for field in _REQUIRED_ENTRY_FIELDS
        if entry.get(field) is None and field != "natom"
    ]
    if "natom" not in entry:
        missing.append("natom")
    if missing:
        raise CampaignSummaryError(f"entry is missing required field(s): {missing}")


def _case_identity(entry):
    gpu = entry["gpu_aggregate"]
    return {field: gpu.get(field) for field in _CASE_IDENTITY_FIELDS}


def _require_one_case_identity(entries):
    reference = _case_identity(entries[0])
    for entry in entries[1:]:
        identity = _case_identity(entry)
        if identity != reference:
            raise CampaignSummaryError(
                f"entries do not share one (campaign, case, variant, machine); "
                f"{reference} vs {identity}"
            )
    return reference


def _quality_status(*aggregates):
    flags = set()
    for aggregate in aggregates:
        flags.update(aggregate.get("quality_flags_union", []) or [])
    authoritative = all(bool(aggregate.get("authoritative")) for aggregate in aggregates)
    return {
        "authoritative": authoritative,
        "quality_flags_union": sorted(flags),
    }


def _build_entry_row(entry):
    _require_entry_fields(entry)
    gpu = entry["gpu_aggregate"]
    cpu_best = entry["cpu_best_aggregate"]
    cpu_1t = entry["cpu_1t_aggregate"]

    s_bestcpu = crossover_mod.compute_gpu_speedup(cpu_best, gpu, label="S_GPU_BESTCPU")
    s_1t = crossover_mod.compute_gpu_speedup(cpu_1t, gpu, label="S_GPU_1T")

    tp = throughput_mod.compute_throughput(
        natom=entry["natom"], t_step=gpu["median"], workload_class=entry["workload_class"],
        directed_interactions=entry.get("directed_interactions"),
        fft_grid_points=entry.get("fft_grid_points"),
    )

    return {
        "size_id": entry["size_id"],
        "natom": entry["natom"],
        "directed_interactions": entry.get("directed_interactions"),
        "fft_grid_points": entry.get("fft_grid_points"),
        "requested_precision": entry["requested_precision"],
        "workload_class": entry["workload_class"],
        "s_gpu_bestcpu": s_bestcpu["speedup"],
        "s_gpu_1t": s_1t["speedup"],
        "t_gpu_seconds": gpu["median"],
        "t_gpu_mad": gpu.get("mad"),
        "t_cpu_best_seconds": cpu_best["median"],
        "t_cpu_best_mad": cpu_best.get("mad"),
        "t_cpu_1t_seconds": cpu_1t["median"],
        "t_cpu_1t_mad": cpu_1t.get("mad"),
        "cpu_best_omp_threads": cpu_best.get("omp_threads"),
        "throughput": tp,
        "gpu_aggregate_id": gpu["aggregate_id"],
        "cpu_best_aggregate_id": cpu_best["aggregate_id"],
        "cpu_1t_aggregate_id": cpu_1t["aggregate_id"],
        "quality_status": _quality_status(gpu, cpu_best, cpu_1t),
    }


def _measured_range(rows):
    natoms = [row["natom"] for row in rows]
    return {"min_natom": min(natoms), "max_natom": max(natoms)}


def _precision_penalties(entries_by_precision):
    """`R_GPU_32_64` at every size where both SINGLE and DOUBLE were measured."""
    single_by_size = {e["size_id"]: e for e in entries_by_precision.get("SINGLE", [])}
    double_by_size = {e["size_id"]: e for e in entries_by_precision.get("DOUBLE", [])}
    penalties = []
    for size_id in sorted(set(single_by_size) & set(double_by_size)):
        pair = [single_by_size[size_id]["gpu_aggregate"], double_by_size[size_id]["gpu_aggregate"]]
        result = precision_metrics.compute_r_gpu_32_64(pair)
        penalties.append({"size_id": size_id, **result})
    return penalties


def build_case_summary(entries, *, thresholds=crossover_mod.DEFAULT_CROSSOVER_THRESHOLDS):
    """Build the section-G machine-readable summary for one case/variant.

    Returns a dict: `identity` (campaign/case/variant/machine), and one entry
    under `by_precision[requested_precision]` per GPU precision present,
    holding that precision's `rows` (per-size table), `measured_range`,
    `crossover` (per threshold) and `quality_status`; plus a top-level
    `precision_penalty` list (`R_GPU_32_64` at every size measured in both
    SINGLE and DOUBLE).
    """
    entries = list(entries)
    if not entries:
        raise CampaignSummaryError("no entries given")
    identity = _require_one_case_identity(entries)

    entries_by_precision = {}
    for entry in entries:
        entries_by_precision.setdefault(entry["requested_precision"], []).append(entry)

    by_precision = {}
    for precision, precision_entries in entries_by_precision.items():
        rows = [_build_entry_row(entry) for entry in precision_entries]
        rows.sort(key=lambda row: row["natom"])

        curve = crossover_mod.build_speedup_curve([
            {"natom": row["natom"], "speedup": row["s_gpu_bestcpu"], "size_id": row["size_id"]}
            for row in rows
        ])
        crossovers = {
            threshold: crossover_mod.find_crossover(curve, threshold) for threshold in thresholds
        }

        by_precision[precision] = {
            "rows": rows,
            "measured_range": _measured_range(rows),
            "crossover": {str(t): result for t, result in crossovers.items()},
            "quality_status": _quality_status(
                *(e["gpu_aggregate"] for e in precision_entries),
                *(e["cpu_best_aggregate"] for e in precision_entries),
                *(e["cpu_1t_aggregate"] for e in precision_entries),
            ),
        }

    return {
        "identity": identity,
        "by_precision": by_precision,
        "precision_penalty": _precision_penalties(entries_by_precision),
    }


def write_summary_json(summary, path):
    """Write `build_case_summary`'s output as indented, sorted-key JSON."""
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return path


_ROW_CSV_FIELDS = (
    "requested_precision", "size_id", "natom", "workload_class",
    "s_gpu_bestcpu", "s_gpu_1t", "t_gpu_seconds", "t_cpu_best_seconds", "t_cpu_1t_seconds",
    "cpu_best_omp_threads", "spin_steps_per_second", "directed_interaction_visits_per_second",
    "fft_grid_points_per_second", "authoritative",
)


def write_summary_csv(summary, path):
    """Flat per-(size, precision) table: one row per measured cell.

    Crossover thresholds and precision penalties are not per-row facts (they
    are per-precision and per-size-pair respectively), so they are not forced
    into this table -- see `write_crossover_csv`/`write_precision_penalty_csv`.
    """
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_ROW_CSV_FIELDS)
        writer.writeheader()
        for block in summary["by_precision"].values():
            for row in block["rows"]:
                throughput = row["throughput"]
                writer.writerow({
                    "requested_precision": row["requested_precision"],
                    "size_id": row["size_id"],
                    "natom": row["natom"],
                    "workload_class": row["workload_class"],
                    "s_gpu_bestcpu": row["s_gpu_bestcpu"],
                    "s_gpu_1t": row["s_gpu_1t"],
                    "t_gpu_seconds": row["t_gpu_seconds"],
                    "t_cpu_best_seconds": row["t_cpu_best_seconds"],
                    "t_cpu_1t_seconds": row["t_cpu_1t_seconds"],
                    "cpu_best_omp_threads": row["cpu_best_omp_threads"],
                    "spin_steps_per_second": throughput.get("spin_steps_per_second"),
                    "directed_interaction_visits_per_second":
                        throughput.get("directed_interaction_visits_per_second"),
                    "fft_grid_points_per_second": throughput.get("fft_grid_points_per_second"),
                    "authoritative": row["quality_status"]["authoritative"],
                })
    return path


_CROSSOVER_CSV_FIELDS = (
    "requested_precision", "threshold", "status", "crossover_natom", "interpolated",
    "tested_range_min_natom", "tested_range_max_natom",
)


def write_crossover_csv(summary, path):
    """One row per (precision, threshold): the section-B crossover table."""
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_CROSSOVER_CSV_FIELDS)
        writer.writeheader()
        for precision, block in sorted(summary["by_precision"].items()):
            crossovers = sorted(block["crossover"].values(), key=lambda result: result["threshold"])
            for result in crossovers:
                writer.writerow({
                    "requested_precision": precision,
                    "threshold": result["threshold"],
                    "status": result["status"],
                    "crossover_natom": result["crossover_natom"],
                    "interpolated": result["interpolated"],
                    "tested_range_min_natom": result["tested_range"][0],
                    "tested_range_max_natom": result["tested_range"][1],
                })
    return path


_PENALTY_CSV_FIELDS = ("size_id", "r_gpu_32_64", "t_single", "t_double")


def write_precision_penalty_csv(summary, path):
    """One row per size where both SINGLE and DOUBLE were measured."""
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_PENALTY_CSV_FIELDS)
        writer.writeheader()
        for penalty in summary["precision_penalty"]:
            writer.writerow({
                "size_id": penalty["size_id"],
                "r_gpu_32_64": penalty["r_gpu_32_64"],
                "t_single": penalty["t_single"],
                "t_double": penalty["t_double"],
            })
    return path
