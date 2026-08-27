"""T=0 OpenMP thread-count correctness sanity check (WP-05 section G).

Compares the real production `restart.<simid>.out` output (the file
`runner.assess_validity` already requires to exist for a `COMPLETED` sample)
of two otherwise-identical T=0 runs that differ only in OpenMP thread count.
The purpose is catching a *gross* physics failure a thread-count change
introduced -- a diverged/garbage final state, a NaN/Inf moment, a flipped
moment direction -- not proving bitwise reproducibility, which OpenMP
reduction-ordering makes unreasonable to demand (blueprint section 9: "Do not
demand bitwise equality where OpenMP reduction ordering makes that
unreasonable").

Tolerance, documented: two independent, deliberately loose checks on the
final-iteration moment state, chosen to comfortably absorb legitimate
reduction-order floating-point differences in a damped, non-chaotic T=0
relaxation while still catching a real divergence:

* per-atom moment magnitude relative difference <= 1e-3 (0.1%);
* per-atom moment direction cosine similarity >= 0.999.

Never used to compare two different `variant_id`s or a finite-temperature run
(where run-to-run RNG differences already make even same-thread-count states
diverge -- see WP-04's finite-T GPU nondeterminism note): callers must supply
two T=0 samples of the same case/variant/size.
"""

from __future__ import annotations

import math
import pathlib
import re

_HEADER_NATOM_RE = re.compile(r"^#\s*Number of atoms:\s*(\d+)", re.MULTILINE)
# Numeric field value, or the NaN/Infinity spellings a Fortran formatted
# write can legitimately produce for a genuinely broken floating-point
# result -- this check needs to recognize those as data, not silently skip
# the row (which would otherwise hide exactly the failure it exists to catch).
_VALUE_TOKEN = r"[+-]?(?:\d[0-9.dDeE+-]*|NaN|Inf(?:inity)?)"
_DATA_ROW_RE = re.compile(
    rf"^\s*(\d+)\s+(\d+)\s+(\d+)\s+({_VALUE_TOKEN})\s+({_VALUE_TOKEN})\s+({_VALUE_TOKEN})\s+({_VALUE_TOKEN})\s*$",
    re.IGNORECASE,
)

DEFAULT_MAGNITUDE_RTOL = 1e-3
DEFAULT_MIN_COSINE_SIMILARITY = 0.999


class OmpSanityError(ValueError):
    """A restart file could not be parsed, or the comparison inputs are
    not a legal T=0 thread-count pair."""


class MomentRecord:
    __slots__ = ("iteration", "ensemble", "iatom", "mmom", "mx", "my", "mz")

    def __init__(self, iteration, ensemble, iatom, mmom, mx, my, mz):
        self.iteration = iteration
        self.ensemble = ensemble
        self.iatom = iatom
        self.mmom = mmom
        self.mx = mx
        self.my = my
        self.mz = mz


class SanityResult:
    """Outcome of comparing two T=0 thread-count runs' final moment state."""

    def __init__(self, *, consistent, max_magnitude_rel_diff, min_cosine_similarity, reasons):
        self.consistent = consistent
        self.max_magnitude_rel_diff = max_magnitude_rel_diff
        self.min_cosine_similarity = min_cosine_similarity
        self.reasons = list(reasons)

    def __repr__(self):
        return (
            f"SanityResult(consistent={self.consistent!r}, "
            f"max_magnitude_rel_diff={self.max_magnitude_rel_diff!r}, "
            f"min_cosine_similarity={self.min_cosine_similarity!r})"
        )


def _parse_fortran_float(token):
    return float(token.replace("D", "E").replace("d", "e"))


def parse_restart_moments(path):
    """Parse a real `restart.<simid>.out` file's final-iteration moment state.

    Returns `(records, declared_natom)` where `records` maps
    `(ensemble, iatom) -> MomentRecord` for the highest `iter` value present
    in the file. Raises :class:`OmpSanityError` if the file does not look
    like production restart output -- this never fabricates a state to
    compare.
    """
    path = pathlib.Path(path)
    text = path.read_text()

    header_match = _HEADER_NATOM_RE.search(text)
    if header_match is None:
        raise OmpSanityError(f"{path}: could not find a 'Number of atoms' header")
    declared_natom = int(header_match.group(1))

    rows = []
    for line in text.splitlines():
        match = _DATA_ROW_RE.match(line)
        if match is None:
            continue
        iteration, ensemble, iatom, mmom, mx, my, mz = match.groups()
        rows.append(MomentRecord(
            int(iteration), int(ensemble), int(iatom),
            _parse_fortran_float(mmom), _parse_fortran_float(mx),
            _parse_fortran_float(my), _parse_fortran_float(mz),
        ))
    if not rows:
        raise OmpSanityError(f"{path}: no moment data rows found")

    last_iteration = max(row.iteration for row in rows)
    final = {
        (row.ensemble, row.iatom): row for row in rows if row.iteration == last_iteration
    }
    return final, declared_natom


def compare_restart_moments(
    path_a, path_b, *,
    magnitude_rtol=DEFAULT_MAGNITUDE_RTOL,
    min_cosine_similarity=DEFAULT_MIN_COSINE_SIMILARITY,
):
    """Compare the final moment state of two restart files for a gross
    physics failure. See module docstring for the tolerance rationale.
    """
    final_a, natom_a = parse_restart_moments(path_a)
    final_b, natom_b = parse_restart_moments(path_b)

    reasons = []
    if natom_a != natom_b:
        reasons.append(f"declared atom count differs: {natom_a} (a) vs {natom_b} (b)")
    if set(final_a) != set(final_b):
        reasons.append("final moment records do not cover the same (ensemble, atom) keys")

    shared_keys = sorted(set(final_a) & set(final_b))
    if not shared_keys:
        reasons.append("no shared (ensemble, atom) records to compare")

    max_rel_diff = 0.0
    min_cosine = 1.0
    for key in shared_keys:
        record_a, record_b = final_a[key], final_b[key]
        components = (record_a.mmom, record_a.mx, record_a.my, record_a.mz,
                      record_b.mmom, record_b.mx, record_b.my, record_b.mz)
        if any(math.isnan(value) or math.isinf(value) for value in components):
            reasons.append(f"NaN/Inf moment component at (ensemble, atom)={key}")
            continue

        denom = max(abs(record_a.mmom), abs(record_b.mmom), 1e-12)
        rel_diff = abs(record_a.mmom - record_b.mmom) / denom
        max_rel_diff = max(max_rel_diff, rel_diff)

        norm_a = math.sqrt(record_a.mx ** 2 + record_a.my ** 2 + record_a.mz ** 2) or 1e-12
        norm_b = math.sqrt(record_b.mx ** 2 + record_b.my ** 2 + record_b.mz ** 2) or 1e-12
        cosine = (
            record_a.mx * record_b.mx + record_a.my * record_b.my + record_a.mz * record_b.mz
        ) / (norm_a * norm_b)
        min_cosine = min(min_cosine, cosine)

    if max_rel_diff > magnitude_rtol:
        reasons.append(
            f"moment magnitude relative difference {max_rel_diff:.3g} exceeds tolerance {magnitude_rtol:.3g}"
        )
    if min_cosine < min_cosine_similarity:
        reasons.append(
            f"moment direction cosine similarity {min_cosine:.6f} below tolerance {min_cosine_similarity:.6f}"
        )

    return SanityResult(
        consistent=not reasons,
        max_magnitude_rel_diff=max_rel_diff,
        min_cosine_similarity=min_cosine,
        reasons=reasons,
    )


def sanity_check_thread_counts(record_a, restart_path_a, record_b, restart_path_b, **tolerance_kwargs):
    """`compare_restart_moments` guarded by the T=0/same-cell preconditions.

    `record_a`/`record_b` are raw sample records (`runner.run_sample`'s
    output): both must be `COMPLETED`, both must report `temperature == 0`,
    and they must agree on `case_id`/`variant_id`/`size_id` -- this check is
    only meaningful for two thread counts of the *same* deterministic cell.
    """
    for record in (record_a, record_b):
        if record.get("run_status") != "COMPLETED":
            raise OmpSanityError(
                f"both samples must be COMPLETED to compare their final state; "
                f"got run_status={record.get('run_status')!r} for run_id={record.get('run_id')!r}"
            )
        if record.get("temperature") != 0:
            raise OmpSanityError(
                f"OpenMP thread-count sanity check requires a deterministic T=0 variant; "
                f"run_id={record.get('run_id')!r} has temperature={record.get('temperature')!r}"
            )
    identity_fields = ("case_id", "variant_id", "size_id")
    mismatched = {field for field in identity_fields if record_a.get(field) != record_b.get(field)}
    if mismatched:
        raise OmpSanityError(
            f"samples are not the same cell; differing fields: {sorted(mismatched)}"
        )
    return compare_restart_moments(restart_path_a, restart_path_b, **tolerance_kwargs)
