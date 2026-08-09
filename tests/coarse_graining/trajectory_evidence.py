"""RCG-04C trajectory parsing, validation, and comparison infrastructure.

This module implements the reusable "canonical moving-evidence record"
parsing/comparison layer specified by ``docs/RCG-04_MOVING_E2E_EVIDENCE.md``
section 5 (RCG-04A) and required by the RCG-04C prompt
(``docs/RCG-04_MOVING_E2E_PROMPT_PACK.md`` section 7). It does not run the
UppASD executable, does not construct fixtures, does not implement a
physical oracle, and does not select or accept any final numerical
tolerance -- comparison routines here report error statistics only. Later
slices (RCG-04D onward) own tolerance selection and moving-dynamics claims.

Production-output audit (before writing this module; see
``docs/RCG-04_MOVING_E2E_EVIDENCE.md`` section 2 for the fuller version and
this docstring for what it means for the design below):

- ``moment.<simid>.out`` (``do_tottraj='Y'``, written by
  ``source/Measurement/prn_trajectories.f90:prn_tottraj`` ->
  ``source/System/restart.f90:prn_mag_conf_iter`` with ``type='M'``) and
  ``restart.<simid>.out`` (``type='R'``, written once at the end of a run)
  share **exactly the same per-record text format**: a 7-line header
  (2 rule lines, file type, simulation type, atom count, ensemble count,
  column header) followed by ``mstep ens iatom |Mom| Mx My Mz`` data rows
  (Fortran format ``i8,i8,i8,2x,es16.8,es16.8,es16.8,es16.8``). The restart
  file is therefore parseable as a trajectory of length exactly one by the
  same parser, with no special-casing -- this module treats it that way.
  ``Mx/My/Mz`` are UppASD's unit moment-direction vector (``emom``, per its
  own doc comment); ``|Mom|`` (``mmom``) is the separately stored magnitude.
- ``moment.<simid>.out`` is written identically regardless of CPU or GPU
  backend: the GPU measurement path calls back into the same Fortran
  ``print_trajectories``/``prn_mag_conf_iter`` routine with synced-back
  ``emom``/``mmom`` (``source/chelper.f90:409-434``). This is genuinely
  backend-neutral production output, reused here rather than adding a new
  diagnostic, per the RCG-04C prompt's "prefer existing ordinary production
  output" instruction.
- The per-event adaptive transition log (``AdaptiveCG: transition
  step=... block=... old_state=... new_state=... accepted=... reason=...
  outcome=... energies_j=<before> <after> <jump>
  polarization_ratio=...``, ``source/CoarseGraining/adaptivecgproduction.f90:
  print_transition_events``, line format at 1450-1461) and the per-step
  resolution-state history (``AdaptiveCG: resolution_state label=step
  step=<n> values=<csv>``, emitted at ``cg_diagnostics>=2`` for every
  synchronization step, not only ``label=initial``/``label=final``) are
  **CPU-only**. They exist and are already exercised by at least one
  tracked fixture (``e2e/adaptive_mixed`` sets ``cg_diagnostics 2``), but
  ``tests/coarse_graining/run_production_e2e.py`` today only regexes the
  aggregate ``accepted_transitions=``/``rejected_transitions=`` counts and
  the final ``resolution_state``/``final_state`` line -- it never parses an
  individual transition event or an intermediate resolution sample. This
  module adds that parser (``parse_transition_events``,
  ``parse_resolution_state_history``) without touching production code.
- **Backend gap found and documented, not fixed, per governing rule** ("if a
  backend genuinely cannot emit an equivalent state record, document the
  gap before changing production code"): the GPU diagnostic path
  (``source/gpu_files/gpuSimulation.cpp:807-980``) prints only
  ``initial``/``final`` snapshots. It has no per-step resolution-state or
  per-event transition log at all, and its one aggregate summary line
  hardcodes ``rejected_transitions=0`` (line ~949) rather than tracking
  rejections. ``parse_transition_events``/``parse_resolution_state_history``
  below will therefore return an empty/single-sample result on GPU stdout
  today; this is a real capability gap for RCG-04G/RCG-04I to consider, not
  a parser defect.
- Named per-step energies/fields are **not** currently available on either
  backend: ``print_adaptive_cg_summary`` (CPU) and its GPU mirror are each
  called exactly once, after the run completes (confirmed by grepping for
  their sole call sites: ``source/uppasd.f90:523`` and
  ``gpuSimulation.cpp``'s ``release()``). ``parse_energy_field_series``
  below is written to accept and merge an arbitrary number of samples (in
  case a later slice adds per-step emission) but will only ever produce one
  sample against today's production output -- this is documented here, not
  silently hidden, matching the RCG-04A audit's identical finding.

None of the three findings above required or received a production-code
change in this slice: every parser and comparison routine below operates
purely on already-emitted stdout text and already-written trajectory files.
"""

from __future__ import annotations

import cmath
import math
import re
from dataclasses import dataclass
from pathlib import Path

from moving_state_generator import Geometry, manifest_json  # noqa: F401  (re-exported for callers)

MODULE_VERSION = "rcg04c-trajectory-evidence-v1"


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------


class TrajectoryError(ValueError):
    """Base class for every error this module raises."""


class TrajectoryParseError(TrajectoryError):
    """Base class for strict spin-trajectory parser failures."""


class TruncatedTrajectoryFileError(TrajectoryParseError):
    pass


class MalformedHeaderError(TrajectoryParseError):
    pass


class NonMonotonicStepError(TrajectoryParseError):
    """A later record's step is not strictly greater than an earlier one."""


class MissingStepError(TrajectoryParseError):
    """A step consistent with the inferred sampling cadence never appears."""


class MissingStepDataError(TrajectoryParseError):
    """A step's records do not cover every declared (ensemble, atom) pair."""


class DuplicateRecordError(TrajectoryParseError):
    """The same (step, ensemble, atom) key appears more than once."""


class InconsistentAtomCountError(TrajectoryParseError):
    """A record references an atom/ensemble index outside the declared range."""


class NonFiniteValueError(TrajectoryParseError):
    pass


class NonUnitDirectionError(TrajectoryParseError):
    """A direction vector's norm is outside the accepted normalization budget."""


class AmbiguousSimulationIdentifierError(TrajectoryError):
    """More than one candidate output file matches an unqualified pattern."""


class MissingOutputFileError(TrajectoryError):
    pass


class TrajectoryKeyMismatchError(TrajectoryError):
    """Two trajectories being compared do not cover the same (step, ens, atom) keys."""


class TrajectoryStageError(TrajectoryError):
    """A derived-metric routine was given input outside its stated contract."""


# ---------------------------------------------------------------------------
# Ambiguity-safe output-file lookup
# ---------------------------------------------------------------------------


def find_unique_output(directory: Path, pattern: str) -> Path:
    """Return the single file in ``directory`` matching ``pattern``.

    Replaces the ambiguity-prone ``next(case.glob(pattern))`` idiom (used by
    ``run_production_e2e.py``'s pre-RCG-04C ``restart_state`` helper, which
    silently picks the first glob match): a directory with zero matches or
    more than one match is a data-provenance defect, not a "pick one and
    continue" situation, per the RCG-04C prompt's explicit requirement to
    detect and reject "ambiguous simulation identifiers".
    """
    matches = sorted(directory.glob(pattern))
    if not matches:
        raise MissingOutputFileError(
            f"no file matching {pattern!r} in {directory}"
        )
    if len(matches) > 1:
        raise AmbiguousSimulationIdentifierError(
            f"{len(matches)} files match {pattern!r} in {directory}, expected exactly "
            f"one (ambiguous simulation identifier): {[m.name for m in matches]}"
        )
    return matches[0]


# ---------------------------------------------------------------------------
# Strict spin-trajectory parser
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SpinRecord:
    moment: float
    direction: tuple[float, float, float]


@dataclass(frozen=True)
class TrajectoryStep:
    step: int
    # keyed by (ensemble, atom), both 1-based, matching the file's own columns
    records: dict[tuple[int, int], SpinRecord]


@dataclass(frozen=True)
class Trajectory:
    """A complete, validated spin trajectory (length 1 for a restart file)."""

    natom: int
    mensemble: int
    file_type: str  # 'M' (moment/trajectory) or 'R' (restart)
    simulation_mode: str
    steps: tuple[TrajectoryStep, ...]
    source: str

    def step_numbers(self) -> tuple[int, ...]:
        return tuple(s.step for s in self.steps)

    def step(self, step_number: int) -> TrajectoryStep:
        for s in self.steps:
            if s.step == step_number:
                return s
        raise KeyError(f"step {step_number} not present in trajectory from {self.source}")

    def direction(self, step_number: int, ensemble: int, atom: int) -> tuple[float, float, float]:
        return self.step(step_number).records[(ensemble, atom)].direction


_HEADER_RULE = re.compile(r"^#{10,}$")
_HEADER_FIELD = re.compile(r"^# (File type|Simulation type|Number of atoms|Number of ensembles):\s*(.+)$")
# The numeric alternative also accepts NaN/Inf (case-insensitive): ordinary
# production ``es16.8`` output never emits these, but the strict parser must
# still recognize and then reject them via NonFiniteValueError (checked by
# math.isfinite below), not bounce them off the syntax check as a truncated
# row -- a defect that silently produced NaN must fail loudly with the
# specific "non-finite value" diagnostic, not a generic parse failure.
_FORTRAN_NUMBER = r"[+-]?(?:\d+\.?\d*(?:[EeDd][+-]?\d+)?|\.\d+(?:[EeDd][+-]?\d+)?|nan|inf(?:inity)?)"
_DATA_LINE = re.compile(
    r"^\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+"
    rf"({_FORTRAN_NUMBER})\s+"
    rf"({_FORTRAN_NUMBER})\s+"
    rf"({_FORTRAN_NUMBER})\s+"
    rf"({_FORTRAN_NUMBER})\s*$",
    re.IGNORECASE,
)


def _parse_fortran_float(text: str) -> float:
    # Fortran ``D`` exponent marker (double precision) is legal in this
    # format's write statements (``es16.8``); Python's float() only accepts
    # 'e'/'E'.
    return float(text.replace("D", "E").replace("d", "e"))


def parse_mag_conf_text(
    text: str,
    *,
    source: str = "<text>",
    normalization_tolerance: float = 1.0e-6,
    require_monotonic_cadence: bool = True,
) -> Trajectory:
    """Strictly parse a ``moment.<simid>.out``/``restart.<simid>.out`` body.

    Rejects (raising a specific subclass of :class:`TrajectoryParseError`):
    a truncated/malformed header or data row, non-finite values, a
    direction vector outside the normalization budget, duplicate
    ``(step, ensemble, atom)`` records, a step whose records do not cover
    every declared ``(ensemble, atom)`` pair, an atom/ensemble index outside
    the declared range, non-monotonically increasing step numbers (records
    "reordered" across step blocks), and a step gap inconsistent with the
    cadence inferred from the first two observed steps.

    Record order *within* a step is not significant: every record is
    indexed by its own explicit ``(step, ensemble, atom)`` columns, never by
    position, so a step's rows may appear in any (ensemble, atom) order.
    """
    lines = text.splitlines()
    if len(lines) < 7:
        raise TruncatedTrajectoryFileError(
            f"{source}: expected a 7-line header, found only {len(lines)} line(s)"
        )
    if not (_HEADER_RULE.match(lines[0]) and _HEADER_RULE.match(lines[5])):
        raise MalformedHeaderError(f"{source}: header rule lines (1 and 6) are missing or malformed")

    header_fields: dict[str, str] = {}
    for line in lines[1:5]:
        match = _HEADER_FIELD.match(line)
        if not match:
            raise MalformedHeaderError(f"{source}: malformed header field line: {line!r}")
        header_fields[match.group(1)] = match.group(2).strip()

    missing_fields = {"File type", "Simulation type", "Number of atoms", "Number of ensembles"} - header_fields.keys()
    if missing_fields:
        raise MalformedHeaderError(f"{source}: header is missing field(s): {sorted(missing_fields)}")

    file_type = header_fields["File type"]
    if file_type not in ("M", "R"):
        raise MalformedHeaderError(f"{source}: unrecognized File type {file_type!r} (expected 'M' or 'R')")
    simulation_mode = header_fields["Simulation type"]
    try:
        natom = int(header_fields["Number of atoms"])
        mensemble = int(header_fields["Number of ensembles"])
    except ValueError as exc:
        raise MalformedHeaderError(f"{source}: non-integer atom/ensemble count in header") from exc
    if natom <= 0 or mensemble <= 0:
        raise MalformedHeaderError(
            f"{source}: header declares non-positive natom={natom} or mensemble={mensemble}"
        )

    if not lines[6].lstrip().startswith("#"):
        raise MalformedHeaderError(
            f"{source}: line 7 (column-header line, e.g. '#iter...') does not start with '#': "
            f"{lines[6]!r}"
        )

    data_lines = lines[7:]
    if not data_lines or not any(line.strip() for line in data_lines):
        raise TruncatedTrajectoryFileError(f"{source}: header present but no data rows follow it")

    expected_keys = {(ens, atom) for ens in range(1, mensemble + 1) for atom in range(1, natom + 1)}

    step_order: list[int] = []
    step_records: dict[int, dict[tuple[int, int], SpinRecord]] = {}
    seen_keys: dict[int, set[tuple[int, int]]] = {}

    for line_number, raw_line in enumerate(data_lines, start=8):
        if not raw_line.strip():
            continue
        match = _DATA_LINE.match(raw_line)
        if not match:
            raise TruncatedTrajectoryFileError(
                f"{source}:{line_number}: data row does not match the expected "
                f"'mstep ens iatom |Mom| Mx My Mz' format: {raw_line!r}"
            )
        mstep = int(match.group(1))
        ens = int(match.group(2))
        atom = int(match.group(3))
        try:
            moment = _parse_fortran_float(match.group(4))
            direction = (
                _parse_fortran_float(match.group(5)),
                _parse_fortran_float(match.group(6)),
                _parse_fortran_float(match.group(7)),
            )
        except ValueError as exc:
            raise TruncatedTrajectoryFileError(
                f"{source}:{line_number}: could not parse numeric fields: {raw_line!r}"
            ) from exc

        for value in (moment, *direction):
            if not math.isfinite(value):
                raise NonFiniteValueError(
                    f"{source}:{line_number}: non-finite value in record "
                    f"step={mstep} ens={ens} atom={atom}: {raw_line!r}"
                )

        if not (1 <= ens <= mensemble and 1 <= atom <= natom):
            raise InconsistentAtomCountError(
                f"{source}:{line_number}: record references ens={ens} atom={atom}, "
                f"outside the declared range ens in [1,{mensemble}], atom in [1,{natom}]"
            )

        norm = math.sqrt(sum(c * c for c in direction))
        if abs(norm - 1.0) > normalization_tolerance:
            raise NonUnitDirectionError(
                f"{source}:{line_number}: direction vector norm {norm} deviates from 1.0 "
                f"by more than {normalization_tolerance} at step={mstep} ens={ens} atom={atom}"
            )

        if mstep not in step_records:
            if step_order and mstep <= step_order[-1]:
                raise NonMonotonicStepError(
                    f"{source}:{line_number}: step {mstep} is not strictly greater than "
                    f"the previous step {step_order[-1]} -- records for a step must not "
                    "be interleaved or reordered across step blocks"
                )
            step_order.append(mstep)
            step_records[mstep] = {}
            seen_keys[mstep] = set()

        key = (ens, atom)
        if key in seen_keys[mstep]:
            raise DuplicateRecordError(
                f"{source}:{line_number}: duplicate record for step={mstep} ens={ens} atom={atom}"
            )
        seen_keys[mstep].add(key)
        step_records[mstep][key] = SpinRecord(moment=moment, direction=direction)

    for mstep in step_order:
        missing = expected_keys - seen_keys[mstep]
        if missing:
            raise MissingStepDataError(
                f"{source}: step {mstep} is missing {len(missing)} of "
                f"{len(expected_keys)} declared (ensemble, atom) record(s): "
                f"{sorted(missing)[:5]}{'...' if len(missing) > 5 else ''}"
            )

    if require_monotonic_cadence and len(step_order) >= 3:
        cadence = step_order[1] - step_order[0]
        if cadence <= 0:
            raise NonMonotonicStepError(f"{source}: inferred step cadence {cadence} is not positive")
        for previous, current in zip(step_order, step_order[1:]):
            if current - previous != cadence:
                raise MissingStepError(
                    f"{source}: step {previous} -> {current} gap ({current - previous}) does not "
                    f"match the cadence ({cadence}) inferred from the first two steps "
                    f"({step_order[0]} -> {step_order[1]}); a sampled step appears to be missing"
                )

    steps = tuple(
        TrajectoryStep(step=mstep, records=step_records[mstep]) for mstep in step_order
    )
    return Trajectory(
        natom=natom, mensemble=mensemble, file_type=file_type,
        simulation_mode=simulation_mode, steps=steps, source=source,
    )


def parse_mag_conf_file(path: Path, **kwargs) -> Trajectory:
    return parse_mag_conf_text(path.read_text(), source=str(path), **kwargs)


def load_restart_state(case_dir: Path, **kwargs) -> Trajectory:
    """Load the unique ``restart.*.out`` in ``case_dir`` as a length-1 trajectory."""
    path = find_unique_output(case_dir, "restart.*.out")
    trajectory = parse_mag_conf_file(path, **kwargs)
    if trajectory.file_type != "R":
        raise MalformedHeaderError(f"{path}: File type {trajectory.file_type!r}, expected 'R'")
    if len(trajectory.steps) != 1:
        raise TrajectoryStageError(
            f"{path}: a restart file must contain exactly one state snapshot, found "
            f"{len(trajectory.steps)}"
        )
    return trajectory


def load_moment_trajectory(case_dir: Path, *, simid: str | None = None, **kwargs) -> Trajectory:
    """Load the unique ``moment.<simid>.out`` in ``case_dir`` as a full trajectory."""
    pattern = f"moment.{simid}.out" if simid else "moment.*.out"
    path = find_unique_output(case_dir, pattern)
    trajectory = parse_mag_conf_file(path, **kwargs)
    if trajectory.file_type != "M":
        raise MalformedHeaderError(f"{path}: File type {trajectory.file_type!r}, expected 'M'")
    return trajectory


# ---------------------------------------------------------------------------
# Comparison and derived metrics (parsing only; no acceptance policy)
# ---------------------------------------------------------------------------


def _matching_step_keys(a: Trajectory, b: Trajectory) -> list[tuple[int, int, int]]:
    """(step, ens, atom) triples present in both trajectories; requires an exact key-set match."""
    a_keys = {(s.step, ens, atom) for s in a.steps for (ens, atom) in s.records}
    b_keys = {(s.step, ens, atom) for s in b.steps for (ens, atom) in s.records}
    if a_keys != b_keys:
        only_a = sorted(a_keys - b_keys)[:5]
        only_b = sorted(b_keys - a_keys)[:5]
        raise TrajectoryKeyMismatchError(
            f"trajectories {a.source!r} and {b.source!r} do not cover the same "
            f"(step, ensemble, atom) keys: {len(a_keys - b_keys)} only in the first "
            f"(e.g. {only_a}), {len(b_keys - a_keys)} only in the second (e.g. {only_b})"
        )
    return sorted(a_keys)


def _record(trajectory: Trajectory, step: int, ens: int, atom: int) -> SpinRecord:
    return trajectory.step(step).records[(ens, atom)]


@dataclass(frozen=True)
class ComponentErrorSummary:
    max_abs_error: tuple[float, float, float]
    rms_error: tuple[float, float, float]
    max_abs_error_overall: float
    rms_error_overall: float
    sample_count: int

    def as_dict(self) -> dict:
        return {
            "max_abs_error_xyz": list(self.max_abs_error),
            "rms_error_xyz": list(self.rms_error),
            "max_abs_error_overall": self.max_abs_error_overall,
            "rms_error_overall": self.rms_error_overall,
            "sample_count": self.sample_count,
        }


def component_trajectory_error(reference: Trajectory, candidate: Trajectory) -> ComponentErrorSummary:
    """Per-component maximum and RMS spin-direction error over every sampled step."""
    keys = _matching_step_keys(reference, candidate)
    if not keys:
        raise TrajectoryStageError("cannot compute a trajectory error over zero matched records")
    max_abs = [0.0, 0.0, 0.0]
    sum_sq = [0.0, 0.0, 0.0]
    max_overall = 0.0
    sum_sq_overall = 0.0
    for step, ens, atom in keys:
        ref = _record(reference, step, ens, atom).direction
        cand = _record(candidate, step, ens, atom).direction
        for c in range(3):
            err = abs(ref[c] - cand[c])
            max_abs[c] = max(max_abs[c], err)
            sum_sq[c] += err * err
            max_overall = max(max_overall, err)
            sum_sq_overall += err * err
    n = len(keys)
    return ComponentErrorSummary(
        max_abs_error=tuple(max_abs),
        rms_error=tuple(math.sqrt(s / n) for s in sum_sq),
        max_abs_error_overall=max_overall,
        rms_error_overall=math.sqrt(sum_sq_overall / (3 * n)),
        sample_count=n,
    )


def _stable_angle_between(u: tuple[float, float, float], v: tuple[float, float, float]) -> float:
    """Angle between two (near-)unit vectors, stable near parallel/antiparallel.

    Uses ``angle = 2*atan2(|u-v|, |u+v|)`` rather than
    ``acos(u.v)``: ``acos`` has vanishing derivative (catastrophic
    cancellation in finite precision) exactly where trajectory comparisons
    most need precision -- near-parallel (small-angle dynamics) and
    near-antiparallel (a spin that has flipped) configurations.
    """
    diff = math.sqrt(sum((u[c] - v[c]) ** 2 for c in range(3)))
    summ = math.sqrt(sum((u[c] + v[c]) ** 2 for c in range(3)))
    return 2.0 * math.atan2(diff, summ)


@dataclass(frozen=True)
class AngularErrorSummary:
    max_radians: float
    rms_radians: float
    sample_count: int

    def as_dict(self) -> dict:
        return {
            "max_radians": self.max_radians,
            "max_degrees": math.degrees(self.max_radians),
            "rms_radians": self.rms_radians,
            "rms_degrees": math.degrees(self.rms_radians),
            "sample_count": self.sample_count,
        }


def angular_trajectory_error(reference: Trajectory, candidate: Trajectory) -> AngularErrorSummary:
    keys = _matching_step_keys(reference, candidate)
    if not keys:
        raise TrajectoryStageError("cannot compute an angular error over zero matched records")
    max_angle = 0.0
    sum_sq = 0.0
    for step, ens, atom in keys:
        ref = _record(reference, step, ens, atom).direction
        cand = _record(candidate, step, ens, atom).direction
        angle = _stable_angle_between(ref, cand)
        max_angle = max(max_angle, angle)
        sum_sq += angle * angle
    n = len(keys)
    return AngularErrorSummary(max_radians=max_angle, rms_radians=math.sqrt(sum_sq / n), sample_count=n)


@dataclass(frozen=True)
class DisplacementSummary:
    """Per-(ensemble, atom) displacement of a single trajectory from its own initial state."""

    final_displacement: dict[tuple[int, int], float]
    max_over_time_displacement: dict[tuple[int, int], float]
    max_final_displacement: float
    max_over_time_displacement_overall: float

    def as_dict(self) -> dict:
        return {
            "max_final_displacement": self.max_final_displacement,
            "max_over_time_displacement_overall": self.max_over_time_displacement_overall,
            "per_atom_final_displacement": {
                f"{ens}:{atom}": value for (ens, atom), value in self.final_displacement.items()
            },
            "per_atom_max_over_time_displacement": {
                f"{ens}:{atom}": value for (ens, atom), value in self.max_over_time_displacement.items()
            },
        }


def spin_displacement(trajectory: Trajectory) -> DisplacementSummary:
    """Nonzero-evolution evidence: initial-to-final and max-over-time per-spin displacement.

    Uses the stable angular distance (radians) as the displacement metric,
    so a spin that reverses through 180 degrees is not undercounted the way
    a naive chord-length or component-difference metric would understate it
    at large angles.
    """
    if len(trajectory.steps) < 2:
        raise TrajectoryStageError(
            f"{trajectory.source}: displacement requires at least 2 steps, found "
            f"{len(trajectory.steps)}"
        )
    initial = trajectory.steps[0]
    final_displacement: dict[tuple[int, int], float] = {}
    max_over_time: dict[tuple[int, int], float] = {}
    for key, initial_record in initial.records.items():
        max_over_time[key] = 0.0
        for step in trajectory.steps:
            record = step.records.get(key)
            if record is None:
                raise TrajectoryKeyMismatchError(
                    f"{trajectory.source}: step {step.step} is missing key {key} present at "
                    f"the initial step {initial.step}"
                )
            distance = _stable_angle_between(initial_record.direction, record.direction)
            max_over_time[key] = max(max_over_time[key], distance)
        final_record = trajectory.steps[-1].records[key]
        final_displacement[key] = _stable_angle_between(initial_record.direction, final_record.direction)
    return DisplacementSummary(
        final_displacement=final_displacement,
        max_over_time_displacement=max_over_time,
        max_final_displacement=max(final_displacement.values()),
        max_over_time_displacement_overall=max(max_over_time.values()),
    )


# ---------------------------------------------------------------------------
# Named energy and field-checksum series (retaining term identity)
# ---------------------------------------------------------------------------

_ENERGY_TERM_KEYS = (
    "atomistic_bilinear", "atomistic_onsite", "coarse_exchange", "coarse_spiralization",
    "coarse_anisotropy", "coarse_external", "coarse_dipole",
)
_FIELD_CHECKSUM_KEYS = ("atom_sum", "atom_norm2", "coarse_sum", "coarse_norm2")

_FLOAT_RE = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
_ENERGY_LINE_RE = re.compile(r"last_energy_j |last_total_energy_j=")
_FIELD_LINE_RE = re.compile(r"last_field_checksums_t ")


@dataclass(frozen=True)
class EnergyFieldSample:
    """One production diagnostic emission: named energies plus field checksums.

    ``step`` is ``None`` when the emitting diagnostic does not itself carry
    a step number (today's ``print_adaptive_cg_summary``, called once at the
    end of the run -- see the module docstring's production-output audit).
    """

    sample_index: int
    step: int | None
    energies_j: dict[str, float]
    field_checksums_t: dict[str, float]

    def as_dict(self) -> dict:
        return {
            "sample_index": self.sample_index,
            "step": self.step,
            "energies_j": dict(self.energies_j),
            "field_checksums_t": dict(self.field_checksums_t),
        }


def parse_energy_field_series(stdout: str) -> list[EnergyFieldSample]:
    """Parse every named-energy/field-checksum diagnostic emission, in order.

    Handles both the CPU ``AdaptiveCG:`` and GPU ``Gpu: AdaptiveCG`` line
    prefixes. Each production emission spans several adjacent printed
    lines (``last_energy_j`` term line(s), then ``last_total_energy_j=``,
    then ``last_field_checksums_t``); this parser groups consecutive lines
    belonging to one emission by that fixed adjacency, starting a new
    sample whenever a ``last_energy_j`` line is seen. It does not assume
    there is only one emission: if a later slice extends the diagnostic to
    print per-step, this parser returns one sample per step without
    modification (exercised by a synthetic multi-sample test).
    """
    samples: list[EnergyFieldSample] = []
    current_energies: dict[str, float] | None = None
    current_fields: dict[str, float] | None = None

    def flush() -> None:
        nonlocal current_energies, current_fields
        if current_energies is not None:
            samples.append(
                EnergyFieldSample(
                    sample_index=len(samples), step=None,
                    energies_j=current_energies, field_checksums_t=current_fields or {},
                )
            )
        current_energies, current_fields = None, None

    for line in stdout.splitlines():
        if _ENERGY_LINE_RE.search(line):
            if "last_energy_j" in line and current_energies is not None and current_fields is not None:
                # A new emission started before the previous one saw its field line.
                flush()
            if current_energies is None:
                current_energies = {}
            for key in _ENERGY_TERM_KEYS + ("last_total_energy_j",):
                # es24.16 is a fixed-width field: '=' may be followed by
                # padding spaces before the (possibly signed) number, e.g.
                # 'atomistic_bilinear= -1.33...' or 'coarse_dipole=  0.0...'.
                found = re.search(rf"\b{re.escape(key)}=\s*({_FLOAT_RE})", line)
                if found:
                    out_key = "total" if key == "last_total_energy_j" else key
                    current_energies[out_key] = float(found.group(1))
        elif _FIELD_LINE_RE.search(line):
            if current_fields is None:
                current_fields = {}
            for key in _FIELD_CHECKSUM_KEYS:
                found = re.search(rf"\b{re.escape(key)}=\s*({_FLOAT_RE})", line)
                if found:
                    current_fields[key] = float(found.group(1))
            flush()

    flush()
    return samples


def compare_energy_field_series(
    reference: list[EnergyFieldSample], candidate: list[EnergyFieldSample]
) -> dict:
    """Term-by-term energy/field-checksum comparison, retaining term identity.

    Compares samples pairwise by index (both series must be produced with
    the same sampling policy by the caller; this routine does not invent a
    step-alignment heuristic). Reports absolute and relative differences
    per term; it does not decide pass/fail.
    """
    if len(reference) != len(candidate):
        raise TrajectoryStageError(
            f"reference has {len(reference)} energy/field sample(s), candidate has "
            f"{len(candidate)}; cannot compare series of different length"
        )
    per_sample = []
    for ref_sample, cand_sample in zip(reference, candidate):
        term_errors = {}
        for key in set(ref_sample.energies_j) | set(cand_sample.energies_j):
            if key not in ref_sample.energies_j or key not in cand_sample.energies_j:
                raise TrajectoryKeyMismatchError(f"energy term {key!r} present in only one series")
            ref_val, cand_val = ref_sample.energies_j[key], cand_sample.energies_j[key]
            term_errors[key] = {
                "reference": ref_val, "candidate": cand_val,
                "abs_error": abs(ref_val - cand_val),
                "rel_error": abs(ref_val - cand_val) / abs(ref_val) if ref_val != 0.0 else math.inf,
            }
        field_errors = {}
        for key in set(ref_sample.field_checksums_t) | set(cand_sample.field_checksums_t):
            if key not in ref_sample.field_checksums_t or key not in cand_sample.field_checksums_t:
                raise TrajectoryKeyMismatchError(f"field checksum {key!r} present in only one series")
            ref_val, cand_val = ref_sample.field_checksums_t[key], cand_sample.field_checksums_t[key]
            field_errors[key] = {
                "reference": ref_val, "candidate": cand_val,
                "abs_error": abs(ref_val - cand_val),
                "rel_error": abs(ref_val - cand_val) / abs(ref_val) if ref_val != 0.0 else math.inf,
            }
        per_sample.append({"energies_j": term_errors, "field_checksums_t": field_errors})
    return {"sample_count": len(per_sample), "samples": per_sample}


# ---------------------------------------------------------------------------
# Adaptive transition history and resolution-state history (RCG-03 diagnostics)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TransitionEvent:
    step: int
    block: int
    old_state: int
    new_state: int
    accepted: bool
    reason: str
    outcome: str
    energy_before_j: float
    energy_after_j: float
    energy_jump_j: float
    polarization_ratio: float  # sentinel -1.0 when not measured, per production convention

    def as_dict(self) -> dict:
        return {
            "step": self.step, "block": self.block, "old_state": self.old_state,
            "new_state": self.new_state, "accepted": self.accepted, "reason": self.reason,
            "outcome": self.outcome, "energy_before_j": self.energy_before_j,
            "energy_after_j": self.energy_after_j, "energy_jump_j": self.energy_jump_j,
            "polarization_ratio": self.polarization_ratio,
        }


_TRANSITION_RE = re.compile(
    r"AdaptiveCG: transition step=(-?\d+) block=(-?\d+) old_state=(-?\d+) new_state=(-?\d+) "
    # es24.16 is a fixed-width field: values shorter than 24 characters are
    # left-padded with spaces, so an arbitrary amount of whitespace can
    # precede each of the three energies and the polarization ratio.
    rf"accepted=([TF]) reason=(\S+) outcome=(\S+) energies_j=\s*({_FLOAT_RE})\s+({_FLOAT_RE})\s+({_FLOAT_RE})\s+"
    rf"polarization_ratio=\s*({_FLOAT_RE})"
)


def parse_transition_events(stdout: str) -> list[TransitionEvent]:
    """Parse every ``AdaptiveCG: transition ...`` line, in emission order.

    Matches the exact field set ``print_transition_events`` writes
    (``source/CoarseGraining/adaptivecgproduction.f90``, format string at
    line ~1450): step, block, old/new state, accepted flag, reason,
    outcome, before/after/jump energy, and polarization ratio. Returns an
    empty list on stdout that contains none (including all GPU stdout
    today -- see the module docstring's documented backend gap).
    """
    events = []
    for match in _TRANSITION_RE.finditer(stdout):
        events.append(
            TransitionEvent(
                step=int(match.group(1)), block=int(match.group(2)),
                old_state=int(match.group(3)), new_state=int(match.group(4)),
                accepted=(match.group(5) == "T"), reason=match.group(6), outcome=match.group(7),
                energy_before_j=float(match.group(8)), energy_after_j=float(match.group(9)),
                energy_jump_j=float(match.group(10)), polarization_ratio=float(match.group(11)),
            )
        )
    return events


@dataclass(frozen=True)
class ResolutionSample:
    label: str
    step: int
    values: tuple[int, ...]

    def as_dict(self) -> dict:
        return {"label": self.label, "step": self.step, "values": list(self.values)}


_RESOLUTION_RE = re.compile(
    r"AdaptiveCG: resolution_state label=(\S+) step=(-?\d+) values=([0-9,\-]*)"
)


def parse_resolution_state_history(stdout: str) -> list[ResolutionSample]:
    """Parse every ``AdaptiveCG: resolution_state label=... step=... values=...`` line.

    Unlike ``run_production_e2e.py``'s existing ``final_state()`` helper
    (which regexes only the *last* occurrence, i.e. ``label=final``), this
    returns every sample in emission order -- including ``label=initial``
    and, at ``cg_diagnostics>=2``, every intermediate ``label=step`` sample
    -- giving a genuine per-step resolution-state time series wherever
    production already emits one.
    """
    samples = []
    for match in _RESOLUTION_RE.finditer(stdout):
        values_text = match.group(3)
        values = tuple(int(v) for v in values_text.split(",")) if values_text else ()
        samples.append(ResolutionSample(label=match.group(1), step=int(match.group(2)), values=values))
    return samples


# ---------------------------------------------------------------------------
# Conical-mode amplitude, phase, and fitted frequency
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ConicalModeSample:
    step: int
    amplitude: float
    phase_radians: float  # wrapped, (-pi, pi]


@dataclass(frozen=True)
class FrequencyFit:
    angular_frequency: float  # radians per unit of x (step index, or time if provided)
    phase_intercept: float
    residual_rms: float
    sample_count: int

    def as_dict(self) -> dict:
        return {
            "angular_frequency": self.angular_frequency,
            "phase_intercept": self.phase_intercept,
            "residual_rms": self.residual_rms,
            "sample_count": self.sample_count,
        }


def conical_mode_series(
    trajectory: Trajectory, *, axis: tuple[float, float, float],
    in_plane_reference: tuple[float, float, float], phase_by_atom: dict[int, float],
    ensemble: int = 1,
) -> list[ConicalModeSample]:
    """Per-step complex order parameter of a conical-spiral mode.

    For each atom ``a`` with known modulation phase ``phase_by_atom[a]``
    (``q.x_a``, e.g. as produced by
    ``moving_state_generator.conical_spiral_state``'s own per-atom phase
    convention -- callers own that mapping, decoupling this metric from any
    one generator), projects the spin onto the orthonormal in-plane frame
    ``(e1, e2) = (in_plane_reference_perp, axis x e1)`` and forms
    ``z_a = (S_a.e1 + i S_a.e2) * exp(-i * phase_by_atom[a])``. The
    per-step order parameter is the mean of ``z_a`` over every atom with a
    supplied phase; its magnitude is the mode amplitude and its argument is
    the mode phase.
    """
    axis_unit = _normalize(axis)
    e1 = _orthogonalize(in_plane_reference, axis_unit)
    e2 = _cross(axis_unit, e1)

    samples = []
    for step in trajectory.steps:
        total = complex(0.0, 0.0)
        count = 0
        for atom, atom_phase in phase_by_atom.items():
            record = step.records.get((ensemble, atom))
            if record is None:
                continue
            d = record.direction
            u = d[0] * e1[0] + d[1] * e1[1] + d[2] * e1[2]
            v = d[0] * e2[0] + d[1] * e2[1] + d[2] * e2[2]
            total += complex(u, v) * cmath.exp(-1j * atom_phase)
            count += 1
        if count == 0:
            raise TrajectoryStageError(f"no atoms from phase_by_atom found at step {step.step}")
        mean = total / count
        samples.append(ConicalModeSample(step=step.step, amplitude=abs(mean), phase_radians=cmath.phase(mean)))
    return samples


def _unwrap_phases(phases: list[float]) -> list[float]:
    if not phases:
        return []
    unwrapped = [phases[0]]
    for phase in phases[1:]:
        delta = phase - phases[len(unwrapped) - 1]
        delta -= 2.0 * math.pi * round(delta / (2.0 * math.pi))
        unwrapped.append(unwrapped[-1] + delta)
    return unwrapped


def fit_conical_mode_frequency(
    samples: list[ConicalModeSample], *, x_by_step: dict[int, float] | None = None
) -> FrequencyFit:
    """Least-squares angular frequency from a conical-mode phase series.

    ``x_by_step`` maps step number to a physical x-coordinate (typically
    simulated time); if omitted, the step index itself is used, so the
    returned frequency is in radians per step.
    """
    if len(samples) < 2:
        raise TrajectoryStageError("frequency fitting requires at least 2 samples")
    xs = [x_by_step[s.step] if x_by_step is not None else float(s.step) for s in samples]
    unwrapped = _unwrap_phases([s.phase_radians for s in samples])
    n = len(xs)
    mean_x = sum(xs) / n
    mean_y = sum(unwrapped) / n
    cov = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, unwrapped))
    var = sum((x - mean_x) ** 2 for x in xs)
    if var == 0.0:
        raise TrajectoryStageError("all sample x-coordinates are identical; cannot fit a slope")
    slope = cov / var
    intercept = mean_y - slope * mean_x
    residual_sq = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(xs, unwrapped))
    return FrequencyFit(
        angular_frequency=slope, phase_intercept=intercept,
        residual_rms=math.sqrt(residual_sq / n), sample_count=n,
    )


# ---------------------------------------------------------------------------
# Signed chirality
# ---------------------------------------------------------------------------


def signed_chirality(
    step: TrajectoryStep, *, directed_bonds: list[tuple[int, int]],
    axis: tuple[float, float, float], ensemble: int = 1,
) -> float:
    """Mean signed chirality over a list of directed bonds.

    ``chi = mean_over_bonds[ axis . (S_i x S_j) ]`` for each directed bond
    ``(i, j)`` in ``directed_bonds``. This uses the *same* triple-product
    orientation as the accepted RCG-02 DMI convention
    (``docs/RCG-02_DMI_HANDEDNESS_EVIDENCE.md``: positive
    ``D_12 . (M_1 x M_2)`` is called right-handed, and positive `D_zx`
    raises that right-handed dimer), so a caller comparing this value's
    sign against a DMI-sign expectation is using one consistent handedness
    definition throughout, not inventing a second one. The bond direction
    (which atom is ``i`` and which is ``j``) is the caller's documented
    convention -- typically consecutive atoms along the propagation axis,
    from the geometry the trajectory came from.
    """
    if not directed_bonds:
        raise TrajectoryStageError("signed_chirality requires at least one directed bond")
    axis_unit = _normalize(axis)
    total = 0.0
    for i, j in directed_bonds:
        si = step.records[(ensemble, i)].direction
        sj = step.records[(ensemble, j)].direction
        cross = _cross(si, sj)
        total += axis_unit[0] * cross[0] + axis_unit[1] * cross[1] + axis_unit[2] * cross[2]
    return total / len(directed_bonds)


def axis_chain_bonds(geometry: Geometry, *, axis_cell_index: int, i0: int = 1,
                      other_cell_indices: tuple[int, int] = (0, 0)) -> list[tuple[int, int]]:
    """Directed nearest-neighbour bonds along one periodic axis, increasing index.

    A small geometry-walking helper shared by the chirality and domain-wall
    metrics above/below: both need an ordered chain of atoms along one
    ``ncell`` axis with the other two cell indices and the basis site held
    fixed, matching ``moving_state_generator.Geometry``'s atom-index
    convention exactly (no separate ordering is invented here).
    """
    if axis_cell_index not in (0, 1, 2):
        raise TrajectoryStageError(f"axis_cell_index must be 0, 1, or 2, got {axis_cell_index}")
    axis_length = (geometry.n1, geometry.n2, geometry.n3)[axis_cell_index]
    other = list(other_cell_indices)
    atoms = []
    for n in range(axis_length):
        cell = [0, 0, 0]
        cell[axis_cell_index] = n
        remaining = [c for c in range(3) if c != axis_cell_index]
        for slot, c in zip(other, remaining):
            cell[c] = slot
        atoms.append(geometry.atom_index(i0, cell[0], cell[1], cell[2]))
    return [(atoms[n], atoms[(n + 1) % axis_length]) for n in range(axis_length)]


# ---------------------------------------------------------------------------
# Domain-wall centre(s), displacement, and periodic crossing events
# ---------------------------------------------------------------------------


def _normalize(vector: tuple[float, float, float]) -> tuple[float, float, float]:
    norm = math.sqrt(sum(c * c for c in vector))
    if norm <= 0.0:
        raise TrajectoryStageError(f"cannot normalize a zero-norm vector {vector!r}")
    return (vector[0] / norm, vector[1] / norm, vector[2] / norm)


def _cross(u: tuple[float, float, float], v: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def _orthogonalize(vector: tuple[float, float, float], axis_unit: tuple[float, float, float]
                    ) -> tuple[float, float, float]:
    parallel = sum(vector[c] * axis_unit[c] for c in range(3))
    perp = tuple(vector[c] - parallel * axis_unit[c] for c in range(3))
    norm = math.sqrt(sum(c * c for c in perp))
    if norm < 1.0e-8:
        raise TrajectoryStageError(f"in_plane_reference {vector!r} is parallel to axis {axis_unit!r}")
    return tuple(c / norm for c in perp)


def domain_wall_centers(
    step: TrajectoryStep, geometry: Geometry, *, axis_cell_index: int, easy_axis: tuple[float, float, float],
    i0: int = 1, other_cell_indices: tuple[int, int] = (0, 0), ensemble: int = 1,
) -> list[float]:
    """Interpolated zero-crossing positions (in cell units) of the easy-axis component.

    Walks the periodic chain from :func:`axis_chain_bonds`, evaluates each
    atom's projection onto ``easy_axis``, and linearly interpolates the
    fractional cell position of every sign change (including the periodic
    wrap between the last and first atom), matching
    ``moving_state_generator.domain_wall_pair_state``'s two-kink
    up-down-up construction without depending on that generator's
    parameters.
    """
    axis_length = (geometry.n1, geometry.n2, geometry.n3)[axis_cell_index]
    bonds = axis_chain_bonds(
        geometry, axis_cell_index=axis_cell_index, i0=i0, other_cell_indices=other_cell_indices
    )
    axis_unit = _normalize(easy_axis)
    atoms = [i for i, _ in bonds]  # ordered, one full periodic loop
    projections = []
    for n, atom in enumerate(atoms):
        d = step.records[(ensemble, atom)].direction
        projections.append(axis_unit[0] * d[0] + axis_unit[1] * d[1] + axis_unit[2] * d[2])

    centers = []
    for n in range(axis_length):
        p0, p1 = projections[n], projections[(n + 1) % axis_length]
        if p0 == 0.0:
            centers.append(float(n))
        elif (p0 < 0.0) != (p1 < 0.0):
            # Linear interpolation of the zero crossing between cell n and n+1.
            fraction = p0 / (p0 - p1)
            centers.append(n + fraction)
    return centers


@dataclass(frozen=True)
class WallCrossingSummary:
    unwrapped_centers_by_step: dict[int, list[float]]
    displacement: float | None  # None when the tracked wall count changes across the trajectory
    crossing_events: list[tuple[int, int, int]]  # (step, from_block, to_block)

    def as_dict(self) -> dict:
        return {
            "unwrapped_centers_by_step": {
                str(step): centers for step, centers in self.unwrapped_centers_by_step.items()
            },
            "displacement": self.displacement,
            "crossing_events": [list(event) for event in self.crossing_events],
        }


def track_wall_crossings(
    center_series: list[tuple[int, list[float]]], *, axis_length_cells: float, block_length_cells: float,
) -> WallCrossingSummary:
    """Unwrap a per-step wall-centre series and detect physical block-boundary crossings.

    ``center_series`` is ``[(step, [center_0, center_1, ...]), ...]`` as
    produced by repeated calls to :func:`domain_wall_centers`, one call per
    step, in increasing step order (the caller's responsibility -- this
    routine does not re-sort, since silently reordering a caller's series
    would hide the same "reordered without explicit index handling" hazard
    the trajectory parser above rejects). Centres are unwrapped per-track
    using nearest-instant matching (each step's centre list must have the
    same length): a jump greater than half the periodic length is assumed
    to be wraparound, not a discontinuous physical jump, and is unwrapped
    accordingly -- the same convention :func:`_unwrap_phases` uses for
    phase.
    """
    if not center_series:
        raise TrajectoryStageError("track_wall_crossings requires at least one sample")
    steps = [s for s, _ in center_series]
    if steps != sorted(steps):
        raise TrajectoryStageError("center_series must be provided in increasing step order")
    wall_count = len(center_series[0][1])
    if wall_count == 0:
        raise TrajectoryStageError("first sample has zero tracked wall centres")

    unwrapped_by_step: dict[int, list[float]] = {}
    tracks: list[list[float]] = [[] for _ in range(wall_count)]
    for step, centers in center_series:
        if len(centers) != wall_count:
            unwrapped_by_step[step] = list(centers)
            for track in tracks:
                track.append(math.nan)
            continue
        for w, center in enumerate(centers):
            if not tracks[w]:
                tracks[w].append(center)
            else:
                previous = tracks[w][-1]
                delta = center - previous
                delta -= axis_length_cells * round(delta / axis_length_cells)
                tracks[w].append(previous + delta)
        unwrapped_by_step[step] = [tracks[w][-1] for w in range(wall_count)]

    consistent = all(len(centers) == wall_count for _, centers in center_series)
    displacement = None
    if consistent:
        displacement = max(
            abs(unwrapped_by_step[steps[-1]][w] - unwrapped_by_step[steps[0]][w])
            for w in range(wall_count)
        )

    crossing_events: list[tuple[int, int, int]] = []
    if consistent:
        for w in range(wall_count):
            previous_block = math.floor((tracks[w][0] % axis_length_cells) / block_length_cells)
            for step, position in zip(steps[1:], tracks[w][1:]):
                current_block = math.floor((position % axis_length_cells) / block_length_cells)
                if current_block != previous_block:
                    crossing_events.append((step, int(previous_block), int(current_block)))
                    previous_block = current_block

    return WallCrossingSummary(
        unwrapped_centers_by_step=unwrapped_by_step, displacement=displacement,
        crossing_events=crossing_events,
    )


# ---------------------------------------------------------------------------
# Machine-readable summaries
# ---------------------------------------------------------------------------


def trajectory_summary(trajectory: Trajectory, *, provenance: dict | None = None) -> dict:
    """Machine-readable record matching ``docs/RCG-04_MOVING_E2E_EVIDENCE.md`` section 5.

    Includes the per-step-per-atom trajectory itself (schema item 2) and a
    ``provenance``/``module_version`` slot (item 1) for the caller to fill
    with commit/build/command evidence -- this module has no access to that
    context, so it is accepted as an opaque pass-through dict rather than
    guessed.
    """
    return {
        "module_version": MODULE_VERSION,
        "provenance": provenance or {},
        "source": trajectory.source,
        "file_type": trajectory.file_type,
        "simulation_mode": trajectory.simulation_mode,
        "natom": trajectory.natom,
        "mensemble": trajectory.mensemble,
        "step_numbers": list(trajectory.step_numbers()),
        "steps": [
            {
                "step": s.step,
                "records": {
                    f"{ens}:{atom}": {"moment": rec.moment, "direction": list(rec.direction)}
                    for (ens, atom), rec in s.records.items()
                },
            }
            for s in trajectory.steps
        ],
    }
