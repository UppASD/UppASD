"""Effective numerical precision audit (WP-06 section A/B).

The blueprint's warning is not hypothetical here: `UPPASD_PRECISION=SINGLE`
does **not** make UppASD's CPU arithmetic single precision. Source
inspection establishes a clean split:

* The Fortran (CPU) tree defines its working kind unconditionally --
  ``dblprec = selected_real_kind(15, 307)`` in
  ``source/Parameters/parameters.f90:8`` -- with no preprocessor guard. No
  ``.f90``/``.F90`` file anywhere in ``source/`` tests the ``SINGLE_PREC``
  macro. The CPU moment state (``source/System/momentdata.f90:10-16``), the
  Hamiltonian (``source/Hamiltonian/*.f90``, all ``real(dblprec)``) and the
  Depondt integrator (``source/Evolution/depondt.f90:62-85``) are therefore
  **double precision in every build**, including a ``UPPASD_PRECISION=SINGLE``
  configuration.
* ``UPPASD_PRECISION=SINGLE`` instead defines the ``SINGLE_PREC`` compiler
  macro on the shared ``asdlib`` target (``CMakeLists.txt:286-288``), which
  only the GPU (CUDA/HIP) C++ translation units under ``source/gpu_files/``
  test. The GPU spin state, Hamiltonian and Depondt kernels all key their
  ``real`` typedef off it (``gpu_files/real_type.h:50-58``,
  ``gpuStructures.hpp:213-223``, ``gpuHamiltonianCalculations.hpp:26-65``,
  ``gpuDepondtIntegrator.cpp:21-61``), and the thermal-field RNG switches
  between ``curandGenerateNormalDouble``/``curandGenerateNormal`` on the same
  macro (``gpu_wrappers.h:124-131``).
* The Fortran<->C++ boundary is deliberately fixed at double regardless of
  build precision: ``fortranData.hpp:8-10`` -- "The Fortran ABI always
  supplies real(dblprec) storage. Keep this type independent from the
  precision selected for device-side real." GPU energy/measurement
  reductions accumulate in the device's ``real`` (float in a SINGLE build)
  and are explicitly widened to ``double`` before crossing back into Fortran
  (``measurement/gpuMeasurement.cpp:835,972-973``).

The practical consequence: there is no such thing as an "effective CPU
SINGLE" mode in this codebase. `effective_cpu_precision` is DOUBLE for every
supported build; only `effective_gpu_precision` responds to
`requested_precision`. See ``benchmarks/PRECISION_AUDIT.md`` for the full
component table and citations this module encodes.
"""

from __future__ import annotations

# One row per WP-06 section A component. `double_build`/`single_build` are
# the effective precision under UPPASD_PRECISION=DOUBLE/SINGLE respectively;
# `citation` is the file:line evidence backing the row. This table is
# documentation encoded as data -- nothing in the harness parses or iterates
# it, it is read by humans (and `test_precision_audit.py`, which checks it
# stays internally consistent with the functions below).
COMPONENT_AUDIT = (
    {
        "component": "CPU magnetic moments / spin state",
        "double_build": "DOUBLE",
        "single_build": "DOUBLE",
        "citation": "source/Parameters/parameters.f90:8; source/System/momentdata.f90:10-16",
    },
    {
        "component": "CPU Hamiltonian evaluation",
        "double_build": "DOUBLE",
        "single_build": "DOUBLE",
        "citation": "source/Hamiltonian/*.f90 (real(dblprec) throughout)",
    },
    {
        "component": "CPU Depondt integrator",
        "double_build": "DOUBLE",
        "single_build": "DOUBLE",
        "citation": "source/Evolution/depondt.f90:62-85",
    },
    {
        "component": "GPU spin state (device emom/emomM/mmom)",
        "double_build": "DOUBLE",
        "single_build": "SINGLE",
        "citation": "gpu_files/real_type.h:50-58; gpuStructures.hpp:213-223",
    },
    {
        "component": "GPU Hamiltonian kernels",
        "double_build": "DOUBLE",
        "single_build": "SINGLE",
        "citation": "gpuHamiltonianCalculations.hpp:26-65",
    },
    {
        "component": "GPU Depondt kernels",
        "double_build": "DOUBLE",
        "single_build": "SINGLE",
        "citation": "gpuDepondtIntegrator.cpp:21-61",
    },
    {
        "component": "Thermal field (curand)",
        "double_build": "DOUBLE",
        "single_build": "SINGLE",
        "citation": "gpu_wrappers.h:124-131 (curandGenerateNormalDouble vs curandGenerateNormal); "
                     "CPU-side RNG stays real(dblprec) unconditionally",
    },
    {
        "component": "Energy / measurement reductions",
        "double_build": "DOUBLE",
        "single_build": "SINGLE partial sums, widened to DOUBLE before returning to Fortran",
        "citation": "measurement/kernels.hpp:36,49,145 (real sum); "
                     "measurement/gpuMeasurement.cpp:835,972-973 (static_cast<double>(...)); "
                     "CPU-side energy stays real(dblprec): source/Hamiltonian/energy.f90:49,108",
    },
)

# Components 4-6 (GPU spin state, GPU Hamiltonian, GPU Depondt) all key off
# the same SINGLE_PREC-gated `real` typedef, so effective_gpu_precision is a
# single value per (requested_precision, gpu_backend) pair, not per
# component. Energy reductions (component 8) are the one exception -- their
# *device-side working precision* still follows SINGLE_PREC, it is only the
# final host-bound value that is widened -- so they are not modelled
# separately here either.

CPU_EFFECTIVE_PRECISION = "DOUBLE"
"""The CPU (Fortran) numerical path's precision in every supported build.
Grounded in COMPONENT_AUDIT rows 1-3: `dblprec` is defined unconditionally
and no `.f90`/`.F90` file tests `SINGLE_PREC`."""

# Only CUDA has been audited against real source (the GPU real-typedef and
# curand-selection sites cited above). HIP shares the same gpu_files/ headers
# textually, but per the blueprint's "no HIP performance claim without real
# hardware" rule this module does not claim to have *verified* HIP's
# effective precision -- it stays UNKNOWN, same as gpu_provenance.py's
# unaudited HIP device metadata.
_AUDITED_GPU_BACKENDS = ("CUDA",)


def resolve_precision_support_state(requested_precision, backend, gpu_backend):
    """Whether this build's precision behaviour is audited well enough to
    classify a comparison (WP-06 section B), from real backend identity
    rather than a blanket default.

    `source_provenance.build_source_context` always reports `"unaudited"`
    for a real (non-MIXED) build -- correctly, since that module (section A:
    git/build identity) has no precision-audit knowledge of its own. This is
    the missing link `provenance.build_static_context` (which *does* import
    this module) is expected to apply on top: CPU is `"supported"`
    unconditionally (`CPU_EFFECTIVE_PRECISION` is unconditional); a GPU
    backend is `"supported"` iff it is in `_AUDITED_GPU_BACKENDS` (CUDA
    today), else `"unaudited"` (HIP -- no HIP performance claim without real
    hardware, per the blueprint). Without this, `classify_comparison_
    precision_class` -- which trusts whatever `precision_support_state` it is
    given -- would report every real CUDA/CPU comparison as `UNAUDITED`
    despite this project's own component-level audit (`PRECISION_AUDIT.md`)
    establishing otherwise; found exactly that way running WP-10's first full
    campaign, since no earlier work package exercised the whole
    provenance -> classification pipeline end to end.
    """
    if requested_precision == "MIXED":
        return "unsupported"
    if backend == "CPU":
        return "supported"
    if backend == "GPU":
        return "supported" if gpu_backend in _AUDITED_GPU_BACKENDS else "unaudited"
    raise ValueError(f"unrecognized backend {backend!r}")


def effective_cpu_precision(requested_precision):
    """CPU-path effective precision for a build requesting `requested_precision`.

    Always `DOUBLE` for a real (buildable) configuration: `UPPASD_PRECISION`
    never changes what `dblprec` resolves to. `MIXED` cannot be built at all
    (CMakeLists.txt fatal-errors at configure time), so there is nothing to
    audit for it; `UNKNOWN` is returned rather than a fabricated `DOUBLE`.
    """
    if requested_precision == "MIXED":
        return "UNKNOWN"
    if requested_precision not in ("SINGLE", "DOUBLE"):
        raise ValueError(f"unrecognized requested_precision {requested_precision!r}")
    return CPU_EFFECTIVE_PRECISION


def effective_gpu_precision(requested_precision, gpu_backend):
    """GPU-path effective precision for `requested_precision` on `gpu_backend`.

    `None` for a CPU-only run (no GPU path exists to have a precision).
    Mirrors `requested_precision` exactly for an audited backend (CUDA):
    `SINGLE_PREC` gates the device `real` typedef straight off
    `UPPASD_PRECISION`, per COMPONENT_AUDIT rows 4-6. `UNKNOWN` for an
    unaudited backend (HIP) or for MIXED (unbuildable).
    """
    if gpu_backend is None:
        return None
    if gpu_backend not in ("CUDA", "HIP"):
        raise ValueError(f"unrecognized gpu_backend {gpu_backend!r}")
    if requested_precision == "MIXED":
        return "UNKNOWN"
    if requested_precision not in ("SINGLE", "DOUBLE"):
        raise ValueError(f"unrecognized requested_precision {requested_precision!r}")
    if gpu_backend not in _AUDITED_GPU_BACKENDS:
        return "UNKNOWN"
    return requested_precision


def classify_comparison_precision_class(
    *, backend, gpu_backend, precision_support_state, effective_cpu_precision, effective_gpu_precision, run_status,
):
    """WP-06 section B: PRECISION_MATCHED vs PRODUCTION_CONFIGURATION vs
    UNAUDITED vs null (VALIDITY.md's four comparison classes), decided from
    real audited precision rather than asserted.

    Ordering matters:

    1. `unsupported` (MIXED) is always null -- the schema requires this and
       an unbuildable configuration can never enter a comparison.
    2. A sample that did not `COMPLETED` carries no execution evidence to
       back any comparison claim, regardless of what the build supports.
    3. `unaudited` precision support (e.g. HIP, or a non-authoritative
       developer context) reports UNAUDITED rather than fabricating a class.
    4. For a `supported`, `COMPLETED` build: a CPU record is
       PRODUCTION_CONFIGURATION -- it is real production evidence, but
       "matched" is a claim about a *pairing*, not a standalone CPU record.
       A GPU record is PRECISION_MATCHED only when its own
       effective_cpu_precision (the host side) and effective_gpu_precision
       genuinely agree -- which, since effective_cpu_precision is always
       DOUBLE (see `effective_cpu_precision`), only ever happens for a GPU
       DOUBLE build. A GPU SINGLE build is therefore always
       PRODUCTION_CONFIGURATION against the always-DOUBLE CPU host, never
       PRECISION_MATCHED.
    """
    if precision_support_state == "unsupported":
        return None
    if run_status != "COMPLETED":
        return None
    if precision_support_state == "unaudited":
        return "UNAUDITED"
    if precision_support_state != "supported":
        raise ValueError(f"unrecognized precision_support_state {precision_support_state!r}")

    if backend == "CPU":
        return "PRODUCTION_CONFIGURATION"
    if backend != "GPU":
        raise ValueError(f"unrecognized backend {backend!r}")

    if effective_gpu_precision in ("SINGLE", "DOUBLE") and effective_gpu_precision == effective_cpu_precision:
        return "PRECISION_MATCHED"
    if effective_gpu_precision in ("SINGLE", "DOUBLE"):
        return "PRODUCTION_CONFIGURATION"
    return "UNAUDITED"
