"""GPU lattice-convolution Hamiltonian evaluator support.

`do_gpu_convolution Y` switches the GPU exchange/DM/tensor Hamiltonian
evaluator from the sparse neighbour-list kernel (the template default for
every case in this project) to an FFT convolution over the device lattice
(`source/gpu_files/gpuLatticeConvolutionHamiltonian.{hpp,cpp}`). It is
entirely separate from the dipole FFT path (`gpu_dipole_mode`/`do_dip`,
`harness/gpu_campaign.py`) -- this one applies to the ordinary exchange
Hamiltonian any `NEIGHBOR_LIST`-class case already builds.

Two things this module exists for, both because the underlying feature has
a *silent* failure mode (`source/gpu_files/gpuHamiltonianCalculations.cpp:
530-542`): a case/build that does not satisfy the gate falls back to the
sparse kernel with a log line, not an error, so a successful exit code is
not evidence the convolution path actually ran.

1. `check_convolution_gate` -- decide, from the case's own *template*
   (never from the case manifest's declared metadata, which does not carry
   `BC`/`do_reduced`), whether this case can even ask for
   `do_gpu_convolution Y` and expect it to engage. Static, offline, no
   process spawned -- meant to be checked before ever running anything.
2. `verify_convolution_active` -- decide, from one sample's actual captured
   stdout, whether the request was honoured. A campaign driver that
   requests convolution must call this before trusting that sample's
   timing as convolution evidence; if the marker is absent, that sample
   measured the sparse kernel mislabelled as convolution and must not be
   folded into a convolution aggregate.

Precision: no separate gating exists (`source/gpu_files/gpuFftWrapper.hpp`'s
`GpuFftComplex`/`GPUFFT_R2C` aliases switch on the same `SINGLE_PREC` every
other GPU code path already uses) -- SINGLE and DOUBLE are both fully
supported, unlike the dipole `OPEN_FFT` path's HIP-only fp32 restriction.
"""

from __future__ import annotations

import pathlib

CONVOLUTION_ACTIVE_MARKER = "device lattice convolution active"
CONVOLUTION_DISABLED_MARKER = "device lattice convolution requested but disabled"

_INPUT_FILENAME = "inpsd.dat"


class GpuConvolutionError(ValueError):
    """A convolution-gate or activation check could not be completed."""


def _read_keyword_all_tokens(template_dir, keyword):
    """Every token following the first line whose token (case-insensitive)
    equals `keyword`, or `None` if the keyword is absent. Mirrors
    `cases.read_keyword`'s case-insensitive matching, but reads a template
    directory directly rather than a generated run directory (this check
    must run before any run directory exists)."""
    inpsd_path = pathlib.Path(template_dir) / _INPUT_FILENAME
    for line in inpsd_path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if parts and parts[0].lower() == keyword.lower():
            return parts[1:]
    return None


def check_convolution_gate(case):
    """Whether `case`'s own template satisfies
    `gpuHamiltonianCalculations.cpp:530-542`'s hard gate for
    `do_gpu_convolution Y` to actually engage (not fall back silently).

    Checks only what is verifiable statically from the template file, not
    from a generated run's actual `Natom`/`NH` (the `Natom == N1*N2*N3*NA`
    and `NH == NA` per-run checks; the latter is "the common case" to get
    wrong per this project's own GPU implementation checklist -- forgetting
    `do_reduced Y`): all-periodic boundary conditions (`BC P P P`) and
    `do_reduced Y`. Returns `(eligible, reason)`; `reason` explains a
    negative result, or is `None` for a positive one.
    """
    bc_tokens = _read_keyword_all_tokens(case.template_dir, "BC")
    if bc_tokens is None:
        return False, "template declares no BC keyword"
    if [token.upper() for token in bc_tokens[:3]] != ["P", "P", "P"]:
        return False, f"BC is {' '.join(bc_tokens[:3])!r}, not fully periodic (P P P required)"

    do_reduced_tokens = _read_keyword_all_tokens(case.template_dir, "do_reduced")
    if do_reduced_tokens is None or do_reduced_tokens[0].upper() != "Y":
        return False, (
            f"do_reduced is {do_reduced_tokens[0] if do_reduced_tokens else 'unset'!r}, "
            "not 'Y' (required; a case with do_reduced N always has NH != NA)"
        )

    return True, None


def verify_convolution_active(stdout_text):
    """Whether a sample that requested `do_gpu_convolution Y` actually ran
    the convolution kernel, from its own captured stdout.

    Returns `True` only when the activation marker is present. `False`
    (never raises) when the disabled-fallback marker is present, or when
    neither marker appears (e.g. a run that failed before reaching
    Hamiltonian setup) -- a caller must treat either as "this sample is not
    convolution evidence", not attempt to distinguish the two further.
    """
    if CONVOLUTION_ACTIVE_MARKER in stdout_text:
        return True
    return False


def fallback_reason(stdout_text):
    """The disabled-fallback banner's own text, if present, else `None`.
    Purely diagnostic -- for logging why `verify_convolution_active` was
    `False`, never for deciding validity itself."""
    for line in stdout_text.splitlines():
        if CONVOLUTION_DISABLED_MARKER in line:
            return line.strip()
    return None
