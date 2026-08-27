"""GPU memory-limit classification (WP-06 section E).

"If a case does not fit GPU memory: classify = unavailable_memory. Record
requested system size and available GPU memory. Do not silently reduce the
problem size."

`estimate_gpu_memory_required_mib` is a deliberately conservative *order of
magnitude* estimate, not a precise device memory model -- getting an exact
number would require either instrumenting the real allocator or a dry-run
probe, neither of which this work package builds. It is grounded in the
component audit (`precision_audit.py`): the GPU device state
(`gpuStructures.hpp`) holds several per-atom 3-vectors at the effective GPU
`real` width (moments, fields, predictor/corrector staging buffers for
Depondt) plus a per-directed-interaction neighbour table (index + coupling
constant). `PER_ATOM_REAL_COMPONENTS`/`BYTES_PER_INTERACTION_INDEX` are
named, overridable constants precisely so a caller with better evidence (a
real allocation trace) can supply a tighter estimate rather than trusting
this one blindly.

`classify_gpu_memory_availability` applies a `safety_margin` on top of the
estimate before comparing against real queried device memory
(`gpu_provenance.query_cuda_devices`' `memory_total_mib`/`memory_used_mib`)
-- the margin exists so a case is never attempted right at the edge of a
crash, not to make the estimate itself more accurate.
"""

from __future__ import annotations

UNAVAILABLE_MEMORY = "unavailable_memory"

# Conservative per-atom state: moments (emom, emomM: 2 x 3 components),
# effective/thermal/torque fields and Depondt predictor/corrector staging
# buffers -- documented as ~8 three-component vectors per atom, per
# `gpuStructures.hpp`'s deviceLattice layout audited in precision_audit.py.
PER_ATOM_REAL_COMPONENTS = 24
# Directed-interaction neighbour table: one interaction index (always int32
# on the GPU regardless of build precision) plus one coupling constant at
# the effective GPU real width.
BYTES_PER_INTERACTION_INDEX = 4

DEFAULT_SAFETY_MARGIN = 1.5


class GpuMemoryError(ValueError):
    """A memory-availability classification could not be computed from the
    given inputs."""


def _real_bytes(effective_gpu_precision):
    if effective_gpu_precision == "SINGLE":
        return 4
    if effective_gpu_precision == "DOUBLE":
        return 8
    raise GpuMemoryError(
        f"cannot estimate GPU memory without an audited effective_gpu_precision "
        f"(got {effective_gpu_precision!r}); see precision_audit.effective_gpu_precision"
    )


def estimate_gpu_memory_required_mib(natom, directed_interactions, effective_gpu_precision):
    """Order-of-magnitude estimate of device memory one GPU sample needs.

    `directed_interactions` may be `None` (pure FFT/dipole workloads carry no
    neighbour-list metadata) and is then simply not counted -- an FFT grid's
    memory footprint is a different, unaudited quantity this function does
    not attempt to estimate.
    """
    if natom <= 0:
        raise GpuMemoryError(f"natom must be positive, got {natom}")
    real_bytes = _real_bytes(effective_gpu_precision)
    state_bytes = natom * PER_ATOM_REAL_COMPONENTS * real_bytes
    interaction_bytes = 0
    if directed_interactions:
        interaction_bytes = directed_interactions * (BYTES_PER_INTERACTION_INDEX + real_bytes)
    return (state_bytes + interaction_bytes) / (1024 ** 2)


def classify_gpu_memory_availability(
    *, natom, directed_interactions, effective_gpu_precision,
    device_memory_total_mib, device_memory_used_mib=0, safety_margin=DEFAULT_SAFETY_MARGIN,
):
    """Decide whether one (size, precision) cell is expected to fit the GPU.

    Returns a dict with `available` (bool), `classification`
    (`"available"` or `UNAVAILABLE_MEMORY`), `estimated_required_mib` and
    `available_mib` (real device headroom: total minus already-used) -- the
    exact numbers a caller records per section E ("record requested system
    size and available GPU memory").
    """
    if device_memory_total_mib is None:
        raise GpuMemoryError("device_memory_total_mib is required (device could not be queried)")
    estimated_required_mib = estimate_gpu_memory_required_mib(
        natom, directed_interactions, effective_gpu_precision
    )
    available_mib = device_memory_total_mib - (device_memory_used_mib or 0)
    available = (estimated_required_mib * safety_margin) <= available_mib
    return {
        "available": available,
        "classification": "available" if available else UNAVAILABLE_MEMORY,
        "estimated_required_mib": estimated_required_mib,
        "available_mib": available_mib,
        "safety_margin": safety_margin,
    }


# CUDA out-of-memory error signatures actually emitted by the runtime/driver
# or by cuFFT, as real evidence a FAILED sample was a memory failure rather
# than a generic bug -- a secondary, post-hoc confirmation signal, never a
# substitute for the preflight classify_gpu_memory_availability check above
# (which exists specifically so a doomed run is never attempted at all).
_OOM_SIGNATURES = (
    "out of memory",
    "cudaerrormemoryallocation",
    "cuda_error_out_of_memory",
    "cufft_alloc_failed",
    "cannot allocate memory",
)


def detect_gpu_oom(stdout, stderr):
    """Whether a completed process's captured output shows a real CUDA
    out-of-memory signature. Case-insensitive substring match against known
    runtime/driver/cuFFT phrasing; never a guess from a bare nonzero exit
    code alone.
    """
    combined = f"{stdout or ''}\n{stderr or ''}".lower()
    return any(signature in combined for signature in _OOM_SIGNATURES)
