"""WP-10 campaign driver.

The one piece benchmarks/harness/README.md says is still missing: something
that actually hands a loaded campaign manifest's resolved cells to
`runner.run_sample`/`gpu_campaign` and writes aggregates. Everything else
(case/campaign resolution, OpenMP sweep resolution, provenance capture,
steady-state fitting, aggregation, GPU memory/precision handling) is WP-01
through WP-09 infrastructure this module only orchestrates -- it invents no
new measurement semantics.

Per (case, variant, size) cell:

1. Resolve workload metadata once (`runner.resolve_workload_metadata`),
   using the CPU binary regardless of which backends the cell measures
   (natom/directed_interactions/fft_grid_points are workload properties,
   independent of backend) -- then immediately delete that probe's work
   directory. B04_dhcpNd's long-range neighbour list makes the
   `do_prnstruct=1` probe's `struct.<simid>.out` genuinely large (observed
   >9 GiB and still growing at 131072 atoms on this host); every probe's
   temporary output is removed as soon as its metadata has been extracted,
   not just B04's.
2. Pick a fixed multi-nstep triple from the measured fixed-cost/natom
   relationship (see NSTEP_TIERS) -- not `steady_state.calibrate_step_span`,
   to avoid spending additional pilot runs recalibrating per cell; this is a
   disclosed methodology choice (docs/BENCHMARK_PRODUCTION_CHARACTERIZATION
   report, methodology section), not an oversight.
3. For CPU: every thread count in the campaign's resolved
   FULL_PHYSICAL_SWEEP, `sample_count` independent `fit_multi_nstep` fits
   each, aggregated via `aggregate.build_fit_aggregate`.
4. For GPU: every requested precision, a GPU-memory preflight
   (`gpu_memory.classify_gpu_memory_availability`; unavailable sizes get
   `gpu_campaign.build_unavailable_memory_record` instead of being
   attempted), then the same fit-aggregate pattern with GPU contamination
   snapshots bracketing every sample folded into the aggregate's
   `quality_flags_union` (provenance.py's documented "post-hoc quality
   pass" pattern -- raw sample records are never rewritten in place).

`environment_quality_mode: STRICT` is honored at the aggregate/report layer
(an aggregate's own `authoritative` field, already `numerical_valid and
environment_valid and not quality_flags_union`) rather than via
`runner.run_sample`'s raise-based per-sample gate: a single flagged sample
mid-campaign logs and continues rather than aborting an unattended
multi-hour run, while flagged data is still never silently promoted to
authoritative evidence -- the substantive STRICT guarantee is unchanged.
Every generated run work directory is deleted immediately after its sample
is recorded; the schema-valid JSON record in `results_dir` is the retained
evidence, not the working directory.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import shutil
import sys
import time
import traceback

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from harness import aggregate as aggregate_mod
from harness import campaigns as campaigns_mod
from harness import cases as cases_mod
from harness import cpu_provenance
from harness import gpu_campaign
from harness import gpu_memory
from harness import gpu_provenance
from harness import omp_sweep
from harness import omp_topology
from harness import provenance as provenance_mod
from harness import runner
from harness import steady_state

# Fixed multi-nstep triples, chosen from pilot calibration on this host
# (docs/BENCHMARK_PRODUCTION_CHARACTERIZATION_2026-08-28.md methodology
# section): fixed/setup cost was measured to scale ~linearly with
# max(natom, fft_grid_points) at ~0.2ms/atom for the symmetry-expanded
# medium-range case (B01_bccFe), and per-step marginal cost is small
# relative to it at every tested size, so a wide nstep spread stays cheap
# relative to the dominant fixed cost while keeping the regression well
# resolved above run-to-run jitter.
NSTEP_TIERS = (
    (20000, (100, 200, 400)),
    (100000, (40, 80, 160)),
    (float("inf"), (15, 30, 60)),
)

# Per-case GPU overrides beyond the universal gpu_mode=1 backend dispatch
# (found the hard way: the first campaign attempt lost 100% of B02/B03 GPU
# samples and 100% of B05 CPU+GPU samples to these two gaps).
#
# B02_skyrmion2D/B03_skyrmion3D templates default to `skyno T` (skyrmion
# number via triangulation); the GPU measurement path unconditionally
# throws on triangulation (source/gpu_files/measurement/gpuMeasurement.cpp,
# per B02_skyrmion2D/README.md "Backend dispatch") -- exactly the
# "nonzero exit code 1; ... 'Simulation finished' missing" failure observed
# on every B02/B03 GPU sample. `skyno` is in both cases' own
# allowed_input_overrides for exactly this reason; CPU keeps the template's
# default triangulation method (skyno stays unoverridden there).
#
# B05_dipoleFFT's `dipole_on` variant only sets `do_dip=1` (correct for the
# CPU arm); a GPU dipole sample needs all three of `gpu_mode=1`, `do_dip=0`
# and `gpu_dipole_mode=OPEN_FFT` set together (B05_dipoleFFT/README.md
# "GPU dispatch") -- otherwise the run either hard-errors ("OPEN_FFT
# requires do_dip=0") or, worse, silently reports zero dipole field.
GPU_EXTRA_OVERRIDES = {
    "B02_skyrmion2D": {"skyno": "Y"},
    "B03_skyrmion3D": {"skyno": "Y"},
    "B05_dipoleFFT": {"do_dip": 0, "gpu_dipole_mode": "OPEN_FFT"},
}

# dipole_backend (schema-required as a non-null string whenever
# workload_class is FFT_DIPOLE/NEIGHBOR_LIST_PLUS_FFT_DIPOLE): the first
# campaign attempt never set this, so schema_validate.validate_record
# rejected -- and runner.run_sample never persists -- every single
# B05_dipoleFFT raw sample, CPU and GPU alike. "BRUTE_FORCE" names the
# validated CPU do_dip=1 finite direct point-dipole sum
# (cases/B05_dipoleFFT/README.md: "the only WP11-validated CPU mode");
# "FFT_CUFFT" matches the vocabulary already used by
# tests/fixtures/valid/cuda_dipole_fft_plus_neighbours.json for the GPU
# OPEN_FFT/cuFFT path. Every other case has no dipole term, so its context
# keeps the schema's ordinary `null`.
DIPOLE_BACKEND = {
    "B05_dipoleFFT": {"CPU": "BRUTE_FORCE", "GPU": "FFT_CUFFT"},
}

RUN_TIMEOUT_SECONDS = 900

CELL_FIELDS = (
    "campaign_id", "case_id", "variant_id", "size_id", "machine_id",
    "build_id", "backend", "gpu_backend", "omp_threads", "requested_precision",
    "measurement_profile",
)


def pick_nstep(scale):
    for threshold, nsteps in NSTEP_TIERS:
        if scale <= threshold:
            return list(nsteps)
    return list(NSTEP_TIERS[-1][1])


class DriverLog:
    def __init__(self, path):
        self.path = pathlib.Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = open(self.path, "a")

    def event(self, kind, **fields):
        record = {"t": time.time(), "kind": kind, **fields}
        line = json.dumps(record, sort_keys=True, default=str)
        self._fh.write(line + "\n")
        self._fh.flush()
        print(line)


def cleanup_dir(path):
    shutil.rmtree(path, ignore_errors=True)


def preserve_failure_evidence(run_dir, results_dir, run_id):
    """Copy a failed run's stdout/stderr (small text files, not the full
    copied-template work directory) to results_dir/failures/<run_id>/ before
    that work directory is deleted -- the first campaign attempt deleted
    every failure's only diagnostic evidence unconditionally, which cost
    real time re-diagnosing two systematic bugs from case documentation
    alone rather than the actual failing stderr.
    """
    run_dir = pathlib.Path(run_dir)
    if not run_dir.exists():
        return
    dest = pathlib.Path(results_dir) / "failures" / run_id
    dest.mkdir(parents=True, exist_ok=True)
    for name in ("stdout.log", "stderr.log", "timing.json"):
        src = run_dir / name
        if src.is_file():
            shutil.copy2(src, dest / name)


def resolve_metadata(case, variant_id, size_id, work_root, probe_binary, log):
    probe_run_id = f"probe__{case.id}__{variant_id}__{size_id}"
    probe_dir = pathlib.Path(work_root) / probe_run_id
    if probe_dir.exists():
        cleanup_dir(probe_dir)
    wm = runner.resolve_workload_metadata(
        case, variant_id, size_id, work_root, probe_binary,
        timeout_seconds=RUN_TIMEOUT_SECONDS, run_id=probe_run_id,
    )
    cleanup_dir(probe_dir)
    log.event("workload_metadata", case_id=case.id, variant_id=variant_id, size_id=size_id,
              natom=wm["natom"], directed_interactions=wm.get("directed_interactions"),
              fft_grid_points=wm.get("fft_grid_points"))
    return wm


def run_fit_family(
    *, case, variant_id, size_id, nsteps, sample_count, campaign_id, work_root, results_dir,
    binary_path, measurement_profile, workload_metadata, base_context, env, extra_overrides,
    cell_tag, log, gpu_bracket=None,
):
    """Run `sample_count` independent multi-nstep fits for one (backend,
    thread-count-or-precision) configuration of one cell. Returns
    (fits, fit_run_ids, completed_records, contamination_flags) -- the last
    is the union of any GPU contamination flags observed across every
    sample's before/after bracket (empty for CPU).
    """
    fits = []
    fit_run_ids = []
    completed_records = []
    contamination_flags = set()

    for sample_index in range(sample_count):
        points = []
        sample_run_ids = []
        ok = True
        for nstep in nsteps:
            run_id = f"{campaign_id}__{cell_tag}__s{sample_index}__n{nstep}"
            before = gpu_bracket() if gpu_bracket else None
            try:
                result = runner.run_sample(
                    case, variant_id, size_id, nstep=nstep, run_id=run_id,
                    campaign_id=campaign_id, sample_index=sample_index,
                    work_root=work_root, results_dir=results_dir, binary_path=binary_path,
                    measurement_profile=measurement_profile, workload_metadata=workload_metadata,
                    context=base_context, env=env, timeout_seconds=RUN_TIMEOUT_SECONDS,
                    extra_overrides=extra_overrides, require_clean_environment=False,
                )
            except Exception as exc:  # noqa: BLE001 -- record and continue, never abort the campaign
                log.event("sample_error", run_id=run_id, error=str(exc), traceback=traceback.format_exc())
                preserve_failure_evidence(pathlib.Path(work_root) / run_id, results_dir, run_id)
                cleanup_dir(pathlib.Path(work_root) / run_id)
                ok = False
                break
            after = gpu_bracket() if gpu_bracket else None
            if before is not None and after is not None:
                contamination_flags |= provenance_mod.classify_gpu_environment(
                    "CUDA", before, after,
                )
            record = result.record
            if record["run_status"] != "COMPLETED":
                log.event("sample_not_completed", run_id=run_id, run_status=record["run_status"],
                          notes=record.get("notes"))
                preserve_failure_evidence(result.run_directory.path, results_dir, run_id)
                cleanup_dir(result.run_directory.path)
                ok = False
                break
            cleanup_dir(result.run_directory.path)
            points.append((nstep, record["process_wall_seconds"]))
            sample_run_ids.append(run_id)
            completed_records.append(record)
        if not ok:
            continue
        try:
            fit = steady_state.fit_multi_nstep(points)
        except steady_state.SteadyStateError as exc:
            log.event("fit_error", cell_tag=cell_tag, sample_index=sample_index, error=str(exc))
            continue
        fits.append(fit)
        fit_run_ids.append(sample_run_ids[0])

    return fits, fit_run_ids, completed_records, contamination_flags


def write_aggregate(fits, fit_run_ids, completed_records, contamination_flags, *, aggregate_id, results_dir, log):
    if not fits:
        log.event("no_fits", aggregate_id=aggregate_id)
        return None
    cell = {field: completed_records[-1][field] for field in CELL_FIELDS}
    quality_flags_union = set(contamination_flags)
    for record in completed_records:
        quality_flags_union.update(record.get("quality_flags", []))
    all_numerical_valid = all(r.get("numerical_valid") for r in completed_records)
    all_environment_valid = all(r.get("environment_valid") for r in completed_records) and not contamination_flags
    agg = aggregate_mod.build_fit_aggregate(
        fits, metric="steady_step_seconds", run_ids=fit_run_ids, cell=cell, aggregate_id=aggregate_id,
        all_numerical_valid=all_numerical_valid, all_environment_valid=all_environment_valid,
        quality_flags_union=quality_flags_union,
    )
    out_dir = pathlib.Path(results_dir) / "aggregates"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{aggregate_id}.json").write_text(json.dumps(agg, indent=2, sort_keys=True) + "\n")
    log.event("aggregate_written", aggregate_id=aggregate_id, median=agg["median"], mad=agg["mad"],
              sample_count=agg["sample_count"], authoritative=agg["authoritative"])
    return agg


def run_campaign(args):
    campaign = campaigns_mod.load_campaign_manifest(args.campaign)
    cases_by_id = cases_mod.discover_cases(args.cases_dir)
    log = DriverLog(pathlib.Path(args.results_dir) / "driver_log.jsonl")
    log.event("campaign_start", campaign_id=campaign.id, tier=campaign.tier)

    topology = omp_topology.capture_host_topology()
    cpu_sweep = campaigns_mod.resolve_cpu_sweep(campaign, topology) if campaign.has_cpu else None
    cpu_run_configs = omp_sweep.resolve_run_configurations(cpu_sweep, topology) if cpu_sweep else []

    cpu_static_context = None
    if campaign.has_cpu:
        cpu_static_context = provenance_mod.build_static_context(
            binary_path=args.cpu_binary, build_dir=args.cpu_build_dir, machine_id=args.machine_id,
        )
        log.event("cpu_build_context", binary_checksum=cpu_static_context["binary_checksum"],
                  git_commit=cpu_static_context["git_commit"], compiler=cpu_static_context["compiler"])

    gpu_binaries = {"SINGLE": args.gpu_fp32_binary, "DOUBLE": args.gpu_fp64_binary}
    gpu_build_dirs = {"SINGLE": args.gpu_fp32_build_dir, "DOUBLE": args.gpu_fp64_build_dir}
    gpu_static_contexts = {}
    if campaign.has_gpu:
        for precision in campaign.gpu_precisions:
            ctx = provenance_mod.build_static_context(
                binary_path=gpu_binaries[precision], build_dir=gpu_build_dirs[precision],
                machine_id=args.machine_id, gpu_index=campaign.gpu_index,
            )
            gpu_static_contexts[precision] = ctx
            log.event("gpu_build_context", precision=precision, binary_checksum=ctx["binary_checksum"],
                      gpu_model=ctx["gpu_model"])

    cells = campaigns_mod.resolve_cells(campaign, cases_by_id)
    unique_cells = []
    seen = set()
    for cell in cells:
        key = (cell["case_id"], cell["variant_id"], cell["size_id"])
        if key not in seen:
            seen.add(key)
            unique_cells.append(key)

    for case_id, variant_id, size_id in unique_cells:
        case = cases_by_id[case_id]
        log.event("cell_start", case_id=case_id, variant_id=variant_id, size_id=size_id)
        try:
            wm = resolve_metadata(case, variant_id, size_id, args.work_root, args.cpu_binary, log)
        except Exception as exc:  # noqa: BLE001
            log.event("cell_metadata_failed", case_id=case_id, variant_id=variant_id, size_id=size_id,
                      error=str(exc), traceback=traceback.format_exc())
            continue

        scale = max(wm["natom"], wm.get("fft_grid_points") or 0)
        nsteps = pick_nstep(scale)
        log.event("nstep_chosen", case_id=case_id, variant_id=variant_id, size_id=size_id,
                  scale=scale, nsteps=nsteps)

        if campaign.has_cpu:
            for run_config in cpu_run_configs:
                threads = run_config.threads
                context = dict(cpu_static_context)
                omp_env_fields = cpu_provenance.capture_omp_environment(run_config.env)
                context.update({
                    "omp_threads": omp_env_fields["omp_threads"],
                    "omp_places": omp_env_fields["omp_places"],
                    "omp_proc_bind": omp_env_fields["omp_proc_bind"],
                    "omp_dynamic": omp_env_fields["omp_dynamic"],
                    "dipole_backend": DIPOLE_BACKEND.get(case_id, {}).get("CPU"),
                })
                cell_tag = f"CPU__{case_id}__{variant_id}__{size_id}__t{threads}"
                fits, fit_run_ids, completed, _ = run_fit_family(
                    case=case, variant_id=variant_id, size_id=size_id, nsteps=nsteps,
                    sample_count=campaign.sample_count, campaign_id=campaign.id,
                    work_root=args.work_root, results_dir=args.results_dir, binary_path=args.cpu_binary,
                    measurement_profile=campaign.measurement_profile, workload_metadata=wm,
                    base_context=context, env=run_config.env, extra_overrides={}, cell_tag=cell_tag, log=log,
                )
                write_aggregate(
                    fits, fit_run_ids, completed, set(),
                    aggregate_id=f"{campaign.id}__{cell_tag}__steady_step_seconds",
                    results_dir=args.results_dir, log=log,
                )

        if campaign.has_gpu:
            devices = gpu_provenance.query_cuda_devices()
            device = next((d for d in devices if str(d.get("index")) == str(campaign.gpu_index)), None)
            for precision in campaign.gpu_precisions:
                context = dict(gpu_static_contexts[precision])
                context["dipole_backend"] = DIPOLE_BACKEND.get(case_id, {}).get("GPU")
                effective_gpu_precision = context["effective_gpu_precision"]
                if device is not None:
                    classification = gpu_memory.classify_gpu_memory_availability(
                        natom=wm["natom"], directed_interactions=wm.get("directed_interactions"),
                        effective_gpu_precision=effective_gpu_precision,
                        device_memory_total_mib=device["memory_total_mib"],
                        device_memory_used_mib=device["memory_used_mib"],
                    )
                    if classification["classification"] == gpu_memory.UNAVAILABLE_MEMORY:
                        run_id = f"{campaign.id}__GPU__{case_id}__{variant_id}__{size_id}__{precision}__unavailable"
                        record = gpu_campaign.build_unavailable_memory_record(
                            run_id=run_id, campaign_id=campaign.id, case=case, variant_id=variant_id,
                            size_id=size_id, sample_index=0, machine_id=args.machine_id, gpu_backend="CUDA",
                            requested_precision=precision, workload_metadata=wm, context=context,
                            memory_classification=classification, device=device,
                            temperature=0.0, timestep=0.0, nstep=nsteps[0],
                            measurement_profile=campaign.measurement_profile,
                        )
                        gpu_campaign.write_record(record, args.results_dir)
                        log.event("gpu_unavailable_memory", case_id=case_id, variant_id=variant_id,
                                  size_id=size_id, precision=precision, classification=classification)
                        continue

                env = gpu_campaign.resolve_gpu_host_env()
                gpu_index = campaign.gpu_index

                def bracket():
                    return provenance_mod.capture_gpu_snapshot("CUDA", gpu_index)

                gpu_extra_overrides = {"gpu_mode": 1, **GPU_EXTRA_OVERRIDES.get(case_id, {})}
                cell_tag = f"GPU__{case_id}__{variant_id}__{size_id}__{precision}"
                fits, fit_run_ids, completed, contamination = run_fit_family(
                    case=case, variant_id=variant_id, size_id=size_id, nsteps=nsteps,
                    sample_count=campaign.sample_count, campaign_id=campaign.id,
                    work_root=args.work_root, results_dir=args.results_dir,
                    binary_path=gpu_binaries[precision], measurement_profile=campaign.measurement_profile,
                    workload_metadata=wm, base_context=context, env=env,
                    extra_overrides=gpu_extra_overrides, cell_tag=cell_tag, log=log, gpu_bracket=bracket,
                )
                write_aggregate(
                    fits, fit_run_ids, completed, contamination,
                    aggregate_id=f"{campaign.id}__{cell_tag}__steady_step_seconds",
                    results_dir=args.results_dir, log=log,
                )

        log.event("cell_done", case_id=case_id, variant_id=variant_id, size_id=size_id)

    log.event("campaign_done")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", required=True)
    parser.add_argument("--cases-dir", default="cases")
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--work-root", required=True)
    parser.add_argument("--machine-id", required=True)
    parser.add_argument("--cpu-binary", required=True)
    parser.add_argument("--cpu-build-dir", required=True)
    parser.add_argument("--gpu-fp32-binary")
    parser.add_argument("--gpu-fp32-build-dir")
    parser.add_argument("--gpu-fp64-binary")
    parser.add_argument("--gpu-fp64-build-dir")
    return parser.parse_args(argv)


if __name__ == "__main__":
    run_campaign(parse_args())
