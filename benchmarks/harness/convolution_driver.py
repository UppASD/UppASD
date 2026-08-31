"""GPU convolution-vs-pair-interaction campaign driver.

Adds `do_gpu_convolution` as a third axis to any campaign manifest's GPU
cells, alongside the precision axis `wp10_driver.py` already drives: for
every GPU (case, variant, size, precision) cell this also attempts a
convolution ("conv") run, in addition to the existing sparse-neighbour-list
("pair") run, whenever the case's own template passes
`harness.gpu_convolution.check_convolution_gate`. A case that does not
pass (this project's own cases: everything except B04_dhcpNd -- B01/B03
lack `do_reduced Y`, B02/B05 are not fully periodic) is logged and skipped
for the conv arm only; its pair-interaction run still proceeds normally.

This is deliberately usable two ways:

1. Point it at `benchmarks/campaigns/dhcpnd_convolution_2026-08-31.yaml`
   for the specific ask (dhcp Nd at 40x40x40 and 100x100x3, CPU OMP
   scaling, GPU SINGLE/DOUBLE pair *and* convolution).
2. Point it at any other campaign manifest (e.g.
   `benchmarks/campaigns/wp10_2026-08-28.yaml`) to add the convolution
   dimension to an otherwise-ordinary campaign -- every GPU cell whose
   case passes the gate gets a conv arm alongside its existing pair arm,
   for direct comparison against already-collected pair-interaction data
   sharing the same `campaign_id`.

CPU cells are unaffected by any of this (`do_gpu_convolution` only means
something on the GPU path) and are driven exactly as `wp10_driver.py`
already does -- this module imports and reuses its CPU-sweep, aggregation
and per-case GPU-override machinery directly rather than duplicating it,
so a fix or convention change there does not need to be kept in sync by
hand. The one genuinely new execution path is the conv arm's own sample
loop, `run_fit_family_with_verification`, which the plain
`wp10_driver.run_fit_family` cannot be reused for: `do_gpu_convolution Y`
falls back to the sparse kernel *silently* on a case/build that does not
satisfy its runtime gate (a real exit code, a real "Simulation finished",
no error -- `gpuHamiltonianCalculations.cpp:530-542`), so every conv
sample's own stdout must be inspected for the activation marker before its
timing is trusted as convolution evidence, not just accepted on ordinary
run-status validity.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import traceback

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from harness import campaigns as campaigns_mod
from harness import cases as cases_mod
from harness import cpu_provenance
from harness import gpu_campaign
from harness import gpu_convolution
from harness import gpu_memory
from harness import gpu_provenance
from harness import omp_sweep
from harness import omp_topology
from harness import provenance as provenance_mod
from harness import runner
from harness import steady_state
from harness import wp10_driver


def run_fit_family_with_verification(
    *, case, variant_id, size_id, nsteps, sample_count, campaign_id, work_root, results_dir,
    binary_path, measurement_profile, workload_metadata, base_context, env, extra_overrides,
    cell_tag, log, gpu_bracket=None,
):
    """`wp10_driver.run_fit_family`, plus one extra step per sample: read
    that sample's own stdout back before cleanup and require
    `gpu_convolution.verify_convolution_active` before trusting its timing.
    A sample that ran to a normal `COMPLETED` status but silently fell back
    to the sparse kernel is logged (`convolution_fallback_detected`,
    preserving stdout/stderr for inspection the same way an ordinary
    failure does) and excluded from the fit exactly like a failed sample --
    it measured the wrong thing, cleanly.

    Also forces `do_prnstruct=0` on every sample -- see
    `wp10_driver.run_fit_family`'s docstring for why: B04_dhcpNd's own
    template defaults to `do_prnstruct 2` (cheap, natom-proportional), but
    this function must never depend on a case's particular default, since
    nothing here is the one-time metadata probe that legitimately needs
    diagnostic output.
    """
    fits = []
    fit_run_ids = []
    completed_records = []
    contamination_flags = set()
    extra_overrides = {"do_prnstruct": 0, **(extra_overrides or {})}

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
                    context=base_context, env=env, timeout_seconds=wp10_driver.RUN_TIMEOUT_SECONDS,
                    extra_overrides=extra_overrides, require_clean_environment=False,
                )
            except Exception as exc:  # noqa: BLE001 -- record and continue, never abort the campaign
                log.event("sample_error", run_id=run_id, error=str(exc), traceback=traceback.format_exc())
                wp10_driver.preserve_failure_evidence(
                    pathlib.Path(work_root) / run_id, results_dir, run_id
                )
                wp10_driver.cleanup_dir(pathlib.Path(work_root) / run_id)
                ok = False
                break
            after = gpu_bracket() if gpu_bracket else None
            if before is not None and after is not None:
                contamination_flags |= provenance_mod.classify_gpu_environment("CUDA", before, after)

            stdout_text = (result.run_directory.path / "stdout.log").read_text()
            record = result.record
            if record["run_status"] == "COMPLETED" and not gpu_convolution.verify_convolution_active(stdout_text):
                reason = gpu_convolution.fallback_reason(stdout_text)
                log.event(
                    "convolution_fallback_detected", run_id=run_id,
                    reason=reason or "no activation marker found in stdout",
                )
                wp10_driver.preserve_failure_evidence(result.run_directory.path, results_dir, run_id)
                wp10_driver.cleanup_dir(result.run_directory.path)
                ok = False
                break

            if record["run_status"] != "COMPLETED":
                log.event("sample_not_completed", run_id=run_id, run_status=record["run_status"],
                          notes=record.get("notes"))
                wp10_driver.preserve_failure_evidence(result.run_directory.path, results_dir, run_id)
                wp10_driver.cleanup_dir(result.run_directory.path)
                ok = False
                break

            wp10_driver.cleanup_dir(result.run_directory.path)
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


def run_campaign(args):
    campaign = campaigns_mod.load_campaign_manifest(args.campaign)
    cases_by_id = cases_mod.discover_cases(args.cases_dir)
    log = wp10_driver.DriverLog(pathlib.Path(args.results_dir) / "driver_log.jsonl")
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
            wm = wp10_driver.resolve_metadata(case, variant_id, size_id, args.work_root, args.cpu_binary, log)
        except Exception as exc:  # noqa: BLE001
            log.event("cell_metadata_failed", case_id=case_id, variant_id=variant_id, size_id=size_id,
                      error=str(exc), traceback=traceback.format_exc())
            continue

        scale = max(wm["natom"], wm.get("fft_grid_points") or 0)
        nsteps = wp10_driver.pick_nstep(scale)
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
                    "dipole_backend": wp10_driver.DIPOLE_BACKEND.get(case_id, {}).get("CPU"),
                })
                cell_tag = f"CPU__{case_id}__{variant_id}__{size_id}__t{threads}"
                fits, fit_run_ids, completed, _ = wp10_driver.run_fit_family(
                    case=case, variant_id=variant_id, size_id=size_id, nsteps=nsteps,
                    sample_count=campaign.sample_count, campaign_id=campaign.id,
                    work_root=args.work_root, results_dir=args.results_dir, binary_path=args.cpu_binary,
                    measurement_profile=campaign.measurement_profile, workload_metadata=wm,
                    base_context=context, env=run_config.env, extra_overrides={}, cell_tag=cell_tag, log=log,
                )
                wp10_driver.write_aggregate(
                    fits, fit_run_ids, completed, set(),
                    aggregate_id=f"{campaign.id}__{cell_tag}__steady_step_seconds",
                    results_dir=args.results_dir, log=log,
                )

        if campaign.has_gpu:
            conv_eligible, conv_reason = gpu_convolution.check_convolution_gate(case)
            log.event("convolution_gate", case_id=case_id, eligible=conv_eligible, reason=conv_reason)

            devices = gpu_provenance.query_cuda_devices()
            device = next((d for d in devices if str(d.get("index")) == str(campaign.gpu_index)), None)
            for precision in campaign.gpu_precisions:
                context = dict(gpu_static_contexts[precision])
                context["dipole_backend"] = wp10_driver.DIPOLE_BACKEND.get(case_id, {}).get("GPU")
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

                base_gpu_overrides = {"gpu_mode": 1, **wp10_driver.GPU_EXTRA_OVERRIDES.get(case_id, {})}

                # Pair (sparse neighbour-list) baseline -- identical to
                # wp10_driver.py's own GPU arm.
                pair_tag = f"GPU__{case_id}__{variant_id}__{size_id}__{precision}__pair"
                fits, fit_run_ids, completed, contamination = wp10_driver.run_fit_family(
                    case=case, variant_id=variant_id, size_id=size_id, nsteps=nsteps,
                    sample_count=campaign.sample_count, campaign_id=campaign.id,
                    work_root=args.work_root, results_dir=args.results_dir,
                    binary_path=gpu_binaries[precision], measurement_profile=campaign.measurement_profile,
                    workload_metadata=wm, base_context=context, env=env,
                    extra_overrides=base_gpu_overrides, cell_tag=pair_tag, log=log, gpu_bracket=bracket,
                )
                wp10_driver.write_aggregate(
                    fits, fit_run_ids, completed, contamination,
                    aggregate_id=f"{campaign.id}__{pair_tag}__steady_step_seconds",
                    results_dir=args.results_dir, log=log,
                )

                # Convolution arm -- only where the case's own template
                # passes the gate; otherwise explicitly skipped, never
                # silently attempted.
                if not conv_eligible:
                    log.event("convolution_skipped", case_id=case_id, variant_id=variant_id,
                              size_id=size_id, precision=precision, reason=conv_reason)
                    continue

                conv_overrides = {**base_gpu_overrides, "do_gpu_convolution": "Y"}
                conv_tag = f"GPU__{case_id}__{variant_id}__{size_id}__{precision}__conv"
                fits, fit_run_ids, completed, contamination = run_fit_family_with_verification(
                    case=case, variant_id=variant_id, size_id=size_id, nsteps=nsteps,
                    sample_count=campaign.sample_count, campaign_id=campaign.id,
                    work_root=args.work_root, results_dir=args.results_dir,
                    binary_path=gpu_binaries[precision], measurement_profile=campaign.measurement_profile,
                    workload_metadata=wm, base_context=context, env=env,
                    extra_overrides=conv_overrides, cell_tag=conv_tag, log=log, gpu_bracket=bracket,
                )
                wp10_driver.write_aggregate(
                    fits, fit_run_ids, completed, contamination,
                    aggregate_id=f"{campaign.id}__{conv_tag}__steady_step_seconds",
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
