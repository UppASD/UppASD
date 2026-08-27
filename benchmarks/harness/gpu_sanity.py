"""T=0 CPU/GPU physical-consistency sanity check (WP-06 section G).

`omp_sanity.compare_restart_moments` already compares the real
`restart.<simid>.out` final moment state of two T=0 runs for a *gross*
physics failure (NaN/Inf, a diverged magnitude, a flipped moment direction)
without caring which run-configuration dimension differs between them --
`sanity_check_thread_counts` is only its OpenMP-specific guard. This module
supplies the CPU/GPU-appropriate guard instead (same cell, opposite
backend), plus a looser default tolerance for a GPU SINGLE-precision
comparison, since fp32 accumulation over many timesteps legitimately drifts
further from an fp64 CPU reference than two fp64 runs differing only in
reduction order.

This is not a substitute for UppASD's scientific validation suite (section
G: "This is not intended to replace the scientific validation suite") --
it exists to catch a benchmark build/run that is silently producing garbage
before its timing is trusted.
"""

from __future__ import annotations

from harness import omp_sanity

DEFAULT_MAGNITUDE_RTOL_DOUBLE = omp_sanity.DEFAULT_MAGNITUDE_RTOL
# An order of magnitude looser than the DOUBLE-vs-DOUBLE bound: documented,
# not derived -- fp32 GPU state accumulates real rounding error a same
# precision comparison would not have, and this only needs to catch a gross
# divergence, not bound it tightly.
DEFAULT_MAGNITUDE_RTOL_SINGLE = 1e-2
DEFAULT_MIN_COSINE_SIMILARITY = omp_sanity.DEFAULT_MIN_COSINE_SIMILARITY


class GpuSanityError(omp_sanity.OmpSanityError):
    """A CPU/GPU sanity comparison's preconditions were not met."""


def default_magnitude_rtol(effective_gpu_precision):
    """The magnitude tolerance appropriate to the GPU sample's own audited
    precision (`precision_audit.effective_gpu_precision`)."""
    if effective_gpu_precision == "SINGLE":
        return DEFAULT_MAGNITUDE_RTOL_SINGLE
    return DEFAULT_MAGNITUDE_RTOL_DOUBLE


def sanity_check_cpu_vs_gpu(cpu_record, cpu_restart_path, gpu_record, gpu_restart_path, **tolerance_kwargs):
    """`compare_restart_moments` guarded by the T=0/same-cell/CPU-vs-GPU
    preconditions. Both samples must be `COMPLETED` and report
    `temperature == 0`, and must agree on `case_id`/`variant_id`/`size_id` --
    this check is only meaningful for the same deterministic cell run on two
    different backends.
    """
    for record in (cpu_record, gpu_record):
        if record.get("run_status") != "COMPLETED":
            raise GpuSanityError(
                f"both samples must be COMPLETED to compare their final state; "
                f"got run_status={record.get('run_status')!r} for run_id={record.get('run_id')!r}"
            )
        if record.get("temperature") != 0:
            raise GpuSanityError(
                f"CPU/GPU sanity check requires a deterministic T=0 variant; "
                f"run_id={record.get('run_id')!r} has temperature={record.get('temperature')!r}"
            )
    if cpu_record.get("backend") != "CPU":
        raise GpuSanityError(f"cpu_record must be backend=CPU, got {cpu_record.get('backend')!r}")
    if gpu_record.get("backend") != "GPU":
        raise GpuSanityError(f"gpu_record must be backend=GPU, got {gpu_record.get('backend')!r}")

    identity_fields = ("case_id", "variant_id", "size_id")
    mismatched = {field for field in identity_fields if cpu_record.get(field) != gpu_record.get(field)}
    if mismatched:
        raise GpuSanityError(f"samples are not the same cell; differing fields: {sorted(mismatched)}")

    tolerance_kwargs.setdefault(
        "magnitude_rtol", default_magnitude_rtol(gpu_record.get("effective_gpu_precision"))
    )
    tolerance_kwargs.setdefault("min_cosine_similarity", DEFAULT_MIN_COSINE_SIMILARITY)
    return omp_sanity.compare_restart_moments(cpu_restart_path, gpu_restart_path, **tolerance_kwargs)
