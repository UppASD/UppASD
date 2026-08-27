"""Host (CPU) memory-limit classification -- the B04_dhcpNd counterpart of
`gpu_memory.py`'s section E (WP-06).

`B04_dhcpNd` (WP-08D) is the first admitted case whose interaction list is
large enough that the *host*, not only the GPU, can genuinely run out of
memory: mean directed-interactions/atom for this case is ~1338 (vs 6-96 for
B01-B03), so its interaction-list footprint dominates total process memory
the same way it does on the GPU. Blueprint B04 section D requires the same
"classify = unavailable_memory, record requested size and memory context,
never a code failure" contract `gpu_memory.py` already implements for the
GPU -- this module gives the host side of it.

Unlike `gpu_memory.estimate_gpu_memory_required_mib` (an *architectural*
order-of-magnitude estimate, because no CUDA build prints its own real
allocation total), `BYTES_PER_DIRECTED_INTERACTION` here is an empirical
*measurement*: UppASD's own production memory accountant
(`source/Tools/profiling.f90::memocc`, enabled by `do_meminfo 1`, already
the B04_dhcpNd template's own baseline) reports a real peak allocation byte
count at every CPU run. Fitting a line through two real `B04_dhcpNd`
`do_meminfo` runs (see `benchmarks/cases/B04_dhcpNd/README.md` section D):

    16x16x16  natom=16,384  directed_interactions=21,921,792  peak=265,807,837 B
    20x20x20  natom=32,000  directed_interactions=42,816,000  peak=517,975,005 B

gives ~12.07 bytes/directed-interaction plus a ~1.2 MB fixed overhead,
negligible at this case's scale. No `effective_precision` parameter is
needed (unlike the GPU estimator): CPU is DOUBLE in every build regardless
of `UPPASD_PRECISION` (`dblprec` is unconditional -- see
`harness/precision_audit.py` and `benchmarks/cases/README.md`'s
BENCH-WP06 finding), so there is only one CPU memory law to fit.
"""

from __future__ import annotations

import pathlib

UNAVAILABLE_MEMORY = "unavailable_memory"

# Empirically measured, not architecturally estimated -- see module
# docstring. Grounded in two real `do_meminfo 1` peak-allocation
# measurements on B04_dhcpNd, the case whose long-range interaction cloud
# motivated adding this module in the first place.
BYTES_PER_DIRECTED_INTERACTION = 12.07
FIXED_OVERHEAD_BYTES = 1_243_000

DEFAULT_SAFETY_MARGIN = 1.5

_MEMINFO_PATH = pathlib.Path("/proc/meminfo")


class HostMemoryError(ValueError):
    """A host memory-availability classification could not be computed."""


def estimate_host_memory_required_mib(natom, directed_interactions):
    """Predicted CPU peak allocation (MiB) for one run at this workload size.

    Linear in `directed_interactions` (see module docstring); `natom` is
    accepted for interface symmetry with `gpu_memory.estimate_gpu_memory_required_mib`
    and to reject an obviously-invalid call, but does not enter the formula
    directly -- for every case admitted so far, the per-atom state (moments,
    fields) is a small fraction of an interaction-list-dominated workload's
    total footprint (confirmed directly for B04_dhcpNd in its README).
    """
    if natom <= 0:
        raise HostMemoryError(f"natom must be positive, got {natom}")
    if directed_interactions is None or directed_interactions < 0:
        raise HostMemoryError(
            f"directed_interactions must be a non-negative count, got {directed_interactions!r}"
        )
    required_bytes = directed_interactions * BYTES_PER_DIRECTED_INTERACTION + FIXED_OVERHEAD_BYTES
    return required_bytes / (1024 ** 2)


def query_host_memory_mib(meminfo_path=None):
    """Real host memory from `/proc/meminfo`: `{"total_mib", "available_mib"}`.

    `MemAvailable` (not `MemFree`) is used for headroom: it already accounts
    for reclaimable page/buffer cache, the same quantity the kernel itself
    uses to decide whether an allocation can succeed -- matching this
    project's shared-host-contention practice of never trusting a naive
    "free" number (see the CGP-05 host-quota finding).
    """
    path = pathlib.Path(meminfo_path) if meminfo_path is not None else _MEMINFO_PATH
    values = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] in ("MemTotal:", "MemAvailable:"):
            values[parts[0].rstrip(":")] = int(parts[1])  # kB
    if "MemTotal" not in values or "MemAvailable" not in values:
        raise HostMemoryError(f"{path}: missing MemTotal/MemAvailable")
    return {
        "total_mib": values["MemTotal"] / 1024,
        "available_mib": values["MemAvailable"] / 1024,
    }


def classify_host_memory_availability(
    *, natom, directed_interactions, host_memory_available_mib, safety_margin=DEFAULT_SAFETY_MARGIN,
):
    """Decide whether one (case, size) cell is expected to fit host memory.

    Same shape as `gpu_memory.classify_gpu_memory_availability`: returns
    `available` (bool), `classification` (`"available"` or
    `UNAVAILABLE_MEMORY`), `estimated_required_mib` and
    `available_mib` -- the numbers blueprint B04 section D requires
    recording ("requested size and memory context") rather than silently
    reducing `size_id`.
    """
    if host_memory_available_mib is None:
        raise HostMemoryError("host_memory_available_mib is required (host could not be queried)")
    estimated_required_mib = estimate_host_memory_required_mib(natom, directed_interactions)
    available = (estimated_required_mib * safety_margin) <= host_memory_available_mib
    return {
        "available": available,
        "classification": "available" if available else UNAVAILABLE_MEMORY,
        "estimated_required_mib": estimated_required_mib,
        "available_mib": host_memory_available_mib,
        "safety_margin": safety_margin,
    }


# Real signatures a host-memory failure can actually produce: UppASD's own
# allocation accountant (`source/Tools/profiling.f90::memocc`, the
# `stat=i_stat` check every `allocate(...)` in this codebase reports through)
# prints "problem of allocation of array" and calls Fortran `stop` on a
# failed allocate; gfortran's own runtime prints a distinct message for an
# allocate that fails via a different path; the Linux OOM killer terminates
# the process with SIGKILL (visible to a `subprocess`-based caller as
# returncode -9, or 137 under a shell wrapper). Never a guess from a bare
# nonzero exit code alone -- the same evidence-based discipline
# `gpu_memory.detect_gpu_oom` uses.
_HOST_OOM_SIGNATURES = (
    "problem of allocation of array",
    "insufficient virtual memory",
    "out of memory",
    "cannot allocate memory",
    "std::bad_alloc",
    "killed",
)


def detect_host_oom(returncode, stdout, stderr):
    """Whether a completed CPU process's exit status/output shows a real
    host out-of-memory signature (SIGKILL exit code, or a known allocation-
    failure message from UppASD's own accountant, gfortran, or the OS)."""
    if returncode in (-9, 137):
        return True
    combined = f"{stdout or ''}\n{stderr or ''}".lower()
    return any(signature in combined for signature in _HOST_OOM_SIGNATURES)
