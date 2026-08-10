#!/usr/bin/env python3
"""RCG-04G: E2E-MOVING-ADAPTIVE.

Demonstrates accepted adaptive-resolution transitions during genuine
domain-wall motion, including resolution-boundary crossing and the RCG-03
polarization-safety contract.

Fixtures: ``e2e/moving_wall_feature_off`` (plain, AdaptiveCG-disabled
physics reference) and ``e2e/moving_wall_adaptive`` (``cg_mask_mode
ADAPTIVE``), both consuming a byte-identical ``restart_seed.out`` generated
by ``moving_state_generator.domain_wall_pair_state`` (RCG-04B), a
deterministic periodic two-kink Neel wall pair, plus a new
``kfile_cg_wall`` uniaxial anisotropy (K1=-0.5 per basis site, easy axis
matching the wall's), needed so the exchange/anisotropy competition gives
the wall pair a genuine, finite width rather than collapsing/diffusing
under exchange alone. See each fixture's ``README.md`` for full provenance.

No external field is used (or usable: ``adaptivecgproduction.f90``
explicitly rejects nonzero ``hfield``, see the RCG-04G evidence section).
The drive is the pair's own intrinsic exchange/anisotropy relaxation
dynamics under damped LLG, the same category of mechanism RCG-04D/E/F use
(a deliberately nonstationary initial condition relaxing under the real
Hamiltonian/integrator, not an externally applied torque) -- see the
evidence document for the observed direction (the confined "down" domain
between the two walls expands, i.e. the walls move apart, consistent with
two same-chirality Neel kinks repelling) and its physical justification.

RCG-04G found and fixed a real, previously-undiscovered production defect
while building this fixture: ``BlockTopology``'s ``atom_to_block`` map used
a block-major traversal counter instead of the canonical global atom index,
silently mislabelling which physical atom belongs to which spatial block
for any block spanning more than one cell along any axis (i.e. essentially
every fixture in this suite). See the evidence document section 15 for the
full writeup; the fix is ``source/CoarseGraining/blocktopology.f90``, not
part of this file.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import math
import subprocess
from pathlib import Path

from fixture_dependencies import MOVING_WALL_ADAPTIVE_CASE, MOVING_WALL_FEATURE_OFF_CASE
from moving_state_generator import Geometry, domain_wall_pair_state
import torque_oracle as orc
import trajectory_evidence as te

ROOT = Path(__file__).with_name("e2e")

GEOMETRY = Geometry(na=2, n1=24, n2=2, n3=2, basis=((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)))
GENERATOR_PARAMS = dict(
    axis_cell_index=0, wall_centers_cells=(9.0, 15.0), width_cells=1.0,
    easy_axis=(0.0, 0.0, 1.0), wall_type="NEEL", chirality=1, cant_deg=2.0,
    moment_magnitude=2.23, separation_margin_widths=4.0,
)
ANISOTROPY_AXIS = (0.0, 0.0, 1.0)
ANISOTROPY_K1 = {1: -0.5, 2: -0.5}
AXIS_LENGTH_CELLS = 24.0
BLOCK_LENGTH_CELLS = 1.0
OFF_SIMID = "cg105gof"
ADAPTIVE_SIMID = "cg105gad"

# First synchronization step; any accepted/rejected transition at step=1 is
# the initial ownership setup (misalignment/polarization evaluated on the
# as-generated state before any dynamics has run), not motion-driven -- see
# the evidence document's classification. A transition strictly after this
# step, and after WALL_MOTION_ONSET_STEP below, is the motion-driven claim.
INITIAL_SETUP_STEP = 1

# Nontriviality floors (independent oracle, exchange+anisotropy combined;
# observed on this fixture: max_torque~0.178, rms_torque~0.073,
# max_field_misalignment_deg~0.80 -- see torque_oracle.combined_torque_report).
MIN_MAX_TORQUE = 0.02
MIN_RMS_TORQUE = 0.02
MIN_FIELD_MISALIGNMENT_DEG = 0.1
MIN_DISPLACEMENT_RAD = 0.02

# Provisional fp64 tolerances (Human review pending; final cross-precision
# budgets are RCG-04I's responsibility). Observed on this fixture: angular
# max~0.0766 rad, rms~0.0159 rad; restart max_abs comparable; independent
# exchange-energy max abs error~0.408 J (reduced units) against an energy
# scale of ~2405 (relative ~1.7e-4). Budgets below give roughly 2x headroom,
# following RCG-04E/F's precedent of absorbing a real, documented physical
# approximation (the coarse operator's own precession-rate mismatch,
# RCG-04E's open item) rather than fitting it away.
MAX_ANGULAR_ERROR_RAD = 0.15
MAX_COMPONENT_ERROR = 0.15
MAX_RESTART_ERROR = 0.15
MAX_ENERGY_ERROR_J_REDUCED = 1.0

# Wall must displace by at least this many cells (net, initial-to-final;
# the pair's damped-oscillatory dynamics -- see the evidence document --
# means net displacement is smaller than the peak mid-run excursion, so
# this floor is deliberately looser than BLOCK_LENGTH_CELLS) and cross at
# least one physical block boundary during the run (observed net
# displacement on this fixture: ~0.176-0.174 cells).
MIN_WALL_DISPLACEMENT_CELLS = 0.1


class NegativeControlDidNotFailError(AssertionError):
    """Raised when a deliberately broken trajectory/claim passed the check anyway."""


def run_case(binary: Path, case: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(binary)], cwd=case, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120, check=False,
    )


def verify_byte_identical_initial_state(off_dir: Path, adaptive_dir: Path) -> None:
    off_bytes = (off_dir / "restart_seed.out").read_bytes()
    adaptive_bytes = (adaptive_dir / "restart_seed.out").read_bytes()
    assert off_bytes == adaptive_bytes, (
        "feature-off and adaptive restart_seed.out are not byte-identical"
    )
    assert hashlib.sha256(off_bytes).hexdigest() == hashlib.sha256(adaptive_bytes).hexdigest()
    off_kfile = (off_dir / "kfile_cg_wall").read_bytes()
    adaptive_kfile = (adaptive_dir / "kfile_cg_wall").read_bytes()
    assert off_kfile == adaptive_kfile, "feature-off and adaptive kfile_cg_wall are not byte-identical"


def verify_normalized_inputs_equal(off_dir: Path, adaptive_dir: Path) -> None:
    ignored_keys = {"simid", "do_adaptive_cg", "block_size_x", "block_size_y",
                     "block_size_z", "cg_operator", "cg_mask_mode", "cg_selector",
                     "cg_refine_threshold", "cg_coarsen_threshold", "cg_polarization_threshold",
                     "cg_update_interval", "cg_minimum_dwell_updates", "cg_buffer_blocks",
                     "cg_reconstruction", "cg_energy_jump_limit", "cg_diagnostics"}

    def normalized(path: Path) -> dict[str, str]:
        records: dict[str, str] = {}
        for line in path.read_text().splitlines():
            fields = line.split()
            if not fields or fields[0] in ("cell",) or fields[0].replace(".", "").replace(
                    "-", "").isdigit():
                continue
            key = fields[0]
            if key not in ignored_keys:
                records[key] = " ".join(fields[1:])
        return records

    off_records = normalized(off_dir / "inpsd.dat")
    adaptive_records = normalized(adaptive_dir / "inpsd.dat")
    mismatched = {
        key: (off_records.get(key), adaptive_records.get(key))
        for key in set(off_records) | set(adaptive_records)
        if off_records.get(key) != adaptive_records.get(key)
    }
    assert not mismatched, (
        "feature-off/adaptive inpsd.dat differ in a physically/numerically relevant "
        f"key beyond the intended AdaptiveCG block: {mismatched}"
    )


def compute_independent_nontriviality() -> orc.TorqueReport:
    state = domain_wall_pair_state(GEOMETRY, **GENERATOR_PARAMS)
    shells = orc.parse_jfile_shells((ROOT / "jfile").read_text())
    return orc.combined_torque_report(
        GEOMETRY, shells, state.direction_by_atom,
        anisotropy_axis=ANISOTROPY_AXIS, anisotropy_k1=ANISOTROPY_K1,
    )


def independent_exchange_energy_series(trajectory: te.Trajectory,
                                        shells: tuple[orc.ExchangeShell, ...],
                                        ) -> dict[int, float]:
    """Independent per-step exchange energy (not either production diagnostic).

    Exchange only, same convention/limitation as RCG-04D's
    ``independent_energy_series``: the anisotropy contribution is not
    included here (this oracle's calibrated convention, see
    torque_oracle.py, was derived for the exchange term only). Production's
    own ``last_energy_j`` (which includes ``atomistic_onsite``, the
    anisotropy term) is checked separately, as supplementary evidence, per
    RCG-04C's documented "final step only" limitation.
    """
    energies: dict[int, float] = {}
    for step in trajectory.step_numbers():
        directions = {
            atom: record.direction
            for (ensemble, atom), record in trajectory.step(step).records.items()
        }
        per_atom = orc.exchange_energy_per_atom(GEOMETRY, shells, directions)
        energies[step] = sum(per_atom.values())
    return energies


def wall_center_series(trajectory: te.Trajectory) -> list[tuple[int, list[float]]]:
    series = []
    for step in trajectory.step_numbers():
        centers = te.domain_wall_centers(
            trajectory.step(step), GEOMETRY, axis_cell_index=0, easy_axis=ANISOTROPY_AXIS,
        )
        series.append((step, centers))
    return series


def _mutate_freeze(trajectory: te.Trajectory) -> te.Trajectory:
    initial = trajectory.steps[0]
    frozen_steps = tuple(
        dataclasses.replace(step, records=initial.records) for step in trajectory.steps
    )
    return dataclasses.replace(trajectory, steps=frozen_steps)


def _mutate_drop_step(trajectory: te.Trajectory) -> te.Trajectory:
    assert len(trajectory.steps) > 2
    kept = trajectory.steps[:1] + trajectory.steps[2:]
    return dataclasses.replace(trajectory, steps=kept)


def _mutate_perturb_one_component(trajectory: te.Trajectory) -> te.Trajectory:
    steps = list(trajectory.steps)
    final = steps[-1]
    key = next(iter(final.records))
    record = final.records[key]
    corrupted_direction = (record.direction[0] + 0.5, record.direction[1], record.direction[2])
    new_records = dict(final.records)
    new_records[key] = dataclasses.replace(record, direction=corrupted_direction)
    steps[-1] = dataclasses.replace(final, records=new_records)
    return dataclasses.replace(trajectory, steps=tuple(steps))


def run_trajectory_negative_controls(off_traj: te.Trajectory, adaptive_traj: te.Trajectory) -> None:
    """Prove the trajectory-parity assertions above are defect-sensitive, not vacuous.

    Every mutation here is applied to an in-memory, already-parsed
    ``Trajectory`` object; no tracked file is ever written or needs
    restoring. The selector-disabled and stationary-wall controls (RCG-04G's
    own required negative controls, distinct from this generic trajectory
    defect-sensitivity check) are documented separately in the evidence
    document as disposable production runs, following the RCG-04F precedent
    for a source/input-level (not in-memory) mutation.
    """
    frozen = _mutate_freeze(adaptive_traj)
    frozen_error = te.angular_trajectory_error(off_traj, frozen)
    if frozen_error.max_radians <= MAX_ANGULAR_ERROR_RAD:
        raise NegativeControlDidNotFailError(
            f"frozen-evolution negative control did not fail: max_radians="
            f"{frozen_error.max_radians} <= budget {MAX_ANGULAR_ERROR_RAD}"
        )
    print(
        f"negative control [freeze]: max_radians={frozen_error.max_radians:.6g} "
        f"> budget {MAX_ANGULAR_ERROR_RAD:.3g} -- failed as expected"
    )

    dropped = _mutate_drop_step(adaptive_traj)
    try:
        te.angular_trajectory_error(off_traj, dropped)
    except te.TrajectoryKeyMismatchError:
        print("negative control [drop-step]: TrajectoryKeyMismatchError raised as expected")
    else:
        raise NegativeControlDidNotFailError(
            "drop-step negative control did not fail: mismatched keys silently accepted"
        )

    perturbed = _mutate_perturb_one_component(adaptive_traj)
    perturbed_error = te.angular_trajectory_error(off_traj, perturbed)
    if perturbed_error.max_radians <= MAX_ANGULAR_ERROR_RAD:
        raise NegativeControlDidNotFailError(
            f"single-component-perturbation negative control did not fail: max_radians="
            f"{perturbed_error.max_radians} <= budget {MAX_ANGULAR_ERROR_RAD}"
        )
    print(
        f"negative control [perturb-one-component]: max_radians={perturbed_error.max_radians:.6g} "
        f"> budget {MAX_ANGULAR_ERROR_RAD:.3g} -- failed as expected"
    )

    # Boundary-crossing claim must itself be defect-sensitive: a trajectory
    # with the wall held exactly at its initial position (no crossing) must
    # not satisfy track_wall_crossings' displacement/crossing evidence.
    stationary_series = [(s, wall_center_series(adaptive_traj)[0][1]) for s, _ in
                          wall_center_series(adaptive_traj)]
    stationary_summary = te.track_wall_crossings(
        stationary_series, axis_length_cells=AXIS_LENGTH_CELLS, block_length_cells=BLOCK_LENGTH_CELLS,
    )
    if stationary_summary.crossing_events:
        raise NegativeControlDidNotFailError(
            "stationary-wall-centers negative control did not fail: a trajectory "
            "held at its initial centers reported a crossing event"
        )
    print("negative control [stationary-centers]: zero crossing_events reported, as expected")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    args = parser.parse_args()
    binary = args.binary.resolve()

    off_dir = ROOT / MOVING_WALL_FEATURE_OFF_CASE
    adaptive_dir = ROOT / MOVING_WALL_ADAPTIVE_CASE

    verify_byte_identical_initial_state(off_dir, adaptive_dir)
    verify_normalized_inputs_equal(off_dir, adaptive_dir)
    print("input equivalence: restart_seed.out/kfile_cg_wall byte-identical; inpsd.dat records "
          "equal apart from the intended AdaptiveCG block")

    report = compute_independent_nontriviality()
    assert report.max_torque > MIN_MAX_TORQUE, (
        f"oracle max_torque={report.max_torque} does not exceed nontriviality floor {MIN_MAX_TORQUE}"
    )
    assert report.rms_torque > MIN_RMS_TORQUE
    assert report.max_field_misalignment_deg > MIN_FIELD_MISALIGNMENT_DEG
    print(
        f"independent oracle (exchange+anisotropy): max_torque={report.max_torque:.6g} "
        f"rms_torque={report.rms_torque:.6g} "
        f"max_field_misalignment_deg={report.max_field_misalignment_deg:.6g} "
        "-- state is demonstrably not a stationary special case"
    )

    off_result = run_case(binary, off_dir)
    assert off_result.returncode == 0, off_result.stdout
    assert "AdaptiveCG:" not in off_result.stdout, off_result.stdout

    adaptive_result = run_case(binary, adaptive_dir)
    assert adaptive_result.returncode == 0, adaptive_result.stdout
    assert "AdaptiveCG: capability accepted" in adaptive_result.stdout, adaptive_result.stdout

    off_traj = te.load_moment_trajectory(off_dir, simid=OFF_SIMID)
    adaptive_traj = te.load_moment_trajectory(adaptive_dir, simid=ADAPTIVE_SIMID)
    assert off_traj.step_numbers() == adaptive_traj.step_numbers()

    off_displacement = te.spin_displacement(off_traj)
    adaptive_displacement = te.spin_displacement(adaptive_traj)
    assert off_displacement.max_final_displacement > MIN_DISPLACEMENT_RAD
    assert adaptive_displacement.max_final_displacement > MIN_DISPLACEMENT_RAD
    print(
        f"spin displacement: feature-off={off_displacement.max_final_displacement:.6g} rad, "
        f"adaptive={adaptive_displacement.max_final_displacement:.6g} rad"
    )

    component = te.component_trajectory_error(off_traj, adaptive_traj)
    angular = te.angular_trajectory_error(off_traj, adaptive_traj)
    assert component.max_abs_error_overall <= MAX_COMPONENT_ERROR, component.as_dict()
    assert angular.max_radians <= MAX_ANGULAR_ERROR_RAD, angular.as_dict()
    print(
        f"trajectory parity vs all-atomistic reference: component max_abs="
        f"{component.max_abs_error_overall:.6g}; angular max={angular.max_radians:.6g} rad "
        f"({math.degrees(angular.max_radians):.6g} deg) rms={angular.rms_radians:.6g} rad"
    )

    off_restart = te.load_restart_state(off_dir)
    adaptive_restart = te.load_restart_state(adaptive_dir)
    restart_component = te.component_trajectory_error(off_restart, adaptive_restart)
    assert restart_component.max_abs_error_overall <= MAX_RESTART_ERROR, restart_component.as_dict()
    print(f"restart-state parity: max_abs={restart_component.max_abs_error_overall:.6g}")

    shells = orc.parse_jfile_shells((ROOT / "jfile").read_text())
    off_energy = independent_exchange_energy_series(off_traj, shells)
    adaptive_energy = independent_exchange_energy_series(adaptive_traj, shells)
    assert set(off_energy) == set(adaptive_energy)
    max_energy_error = max(abs(off_energy[s] - adaptive_energy[s]) for s in off_energy)
    assert max_energy_error <= MAX_ENERGY_ERROR_J_REDUCED, (
        f"independent-oracle exchange energy max abs error {max_energy_error} exceeds "
        f"budget {MAX_ENERGY_ERROR_J_REDUCED}"
    )
    print(f"named energy (independent oracle, exchange only): max abs error={max_energy_error:.6g} "
          f"(energy scale ~{abs(off_energy[0]):.4g})")

    adaptive_series = te.parse_energy_field_series(adaptive_result.stdout)
    assert len(adaptive_series) == 1
    assert adaptive_series[0].energies_j["atomistic_onsite"] != 0.0, (
        "expected a nonzero anisotropy (atomistic_onsite) energy term"
    )
    print(f"production last_energy_j (final step only, RCG-04C limitation): "
          f"atomistic_onsite={adaptive_series[0].energies_j['atomistic_onsite']:.6g} J")

    # --- Wall-centre tracking with periodic unwrapping (independent of the
    # production resolution-state diagnostic; parsed purely from the
    # per-step, per-atom trajectory files) ---
    off_series = wall_center_series(off_traj)
    adaptive_series_centers = wall_center_series(adaptive_traj)
    off_crossings = te.track_wall_crossings(
        off_series, axis_length_cells=AXIS_LENGTH_CELLS, block_length_cells=BLOCK_LENGTH_CELLS,
    )
    adaptive_crossings = te.track_wall_crossings(
        adaptive_series_centers, axis_length_cells=AXIS_LENGTH_CELLS, block_length_cells=BLOCK_LENGTH_CELLS,
    )
    assert off_crossings.displacement is not None and off_crossings.displacement > MIN_WALL_DISPLACEMENT_CELLS, (
        f"feature-off wall displacement {off_crossings.displacement} does not exceed floor "
        f"{MIN_WALL_DISPLACEMENT_CELLS} cells"
    )
    assert adaptive_crossings.displacement is not None and (
        adaptive_crossings.displacement > MIN_WALL_DISPLACEMENT_CELLS
    ), (
        f"adaptive wall displacement {adaptive_crossings.displacement} does not exceed floor "
        f"{MIN_WALL_DISPLACEMENT_CELLS} cells"
    )
    assert off_crossings.crossing_events, "feature-off reference shows no block-boundary crossing"
    assert adaptive_crossings.crossing_events, "adaptive run shows no block-boundary crossing"
    first_crossing_step = min(step for step, _, _ in adaptive_crossings.crossing_events)
    assert first_crossing_step > INITIAL_SETUP_STEP, (
        f"first crossing at step={first_crossing_step} is not after the initial setup step "
        f"{INITIAL_SETUP_STEP} -- would be a startup artifact, not motion"
    )
    print(
        f"wall tracking: feature-off displacement={off_crossings.displacement:.4g} cells, "
        f"crossings={off_crossings.crossing_events}; adaptive displacement="
        f"{adaptive_crossings.displacement:.4g} cells, crossings={adaptive_crossings.crossing_events} "
        f"(first at step={first_crossing_step}, after the step={INITIAL_SETUP_STEP} setup evaluation)"
    )

    # --- Transition history: full per-event record, classified ---
    events = te.parse_transition_events(adaptive_result.stdout)
    assert events, "no transition events recorded"
    setup_events = [e for e in events if e.step == INITIAL_SETUP_STEP]
    motion_events = [e for e in events if e.step > INITIAL_SETUP_STEP]
    assert setup_events, "expected initial-ownership coarsening at the first synchronization step"
    assert all(e.accepted and e.reason == "coarsen-request" for e in setup_events), (
        "initial-setup transitions were expected to be accepted far-block coarsening only"
    )
    assert motion_events, (
        "expected at least one accepted/rejected transition strictly after the initial setup "
        "step and after real wall motion has begun -- none were recorded"
    )
    accepted_motion_events = [e for e in motion_events if e.accepted]
    assert accepted_motion_events, (
        "expected at least one ACCEPTED transition strictly after the initial setup step, "
        "spatially tied to the moving wall; only rejections (if any) were recorded"
    )
    for event in accepted_motion_events:
        assert event.step > first_crossing_step, (
            f"accepted motion-driven transition at step={event.step} occurs before the first "
            f"wall crossing at step={first_crossing_step} -- not demonstrably motion-driven"
        )
    print(
        f"transition history: {len(setup_events)} initial-setup coarsen events (step="
        f"{INITIAL_SETUP_STEP}); {len(motion_events)} motion-window events, "
        f"{len(accepted_motion_events)} accepted, all strictly after the first wall crossing "
        f"(step={first_crossing_step}): "
        + "; ".join(
            f"step={e.step} block={e.block} {e.old_state}->{e.new_state} reason={e.reason}"
            for e in accepted_motion_events
        )
    )
    if not any(e.reason == "refine-request" and e.accepted for e in events):
        print(
            "REVIEWED LIMITATION (not a defect): no accepted refine-request (coarse-to-atomistic) "
            "transition occurs in this fixture -- the wall's total excursion (crossing exactly "
            "one block boundary each) never reaches back into a block that coarsened during "
            "initial setup. See the evidence document for the parameter search that did produce "
            "genuine refine-REQUESTS (rejected by ALIGNED reconstruction, a distinct, itself "
            "informative outcome) at a tighter geometry, not adopted here to keep this fixture's "
            "accepted-case evidence clean and unambiguous."
        )

    # --- RCG-03 polarization safety: every transition's ratio, and every
    # accepted coarsen must be safely polarized ---
    coarsen_events = [e for e in events if e.reason == "coarsen-request" and e.accepted]
    assert all(e.polarization_ratio > 0.9 for e in coarsen_events), (
        "an accepted coarsen-request had a polarization_ratio at/under the "
        "cg_polarization_threshold=0.9 gate -- would itself be a safety-contract violation"
    )
    print(
        f"RCG-03 polarization safety: all {len(coarsen_events)} accepted coarsen-request "
        f"polarization ratios exceed the 0.9 gate (min observed "
        f"{min(e.polarization_ratio for e in coarsen_events):.6g})"
    )

    # --- Defect sensitivity ---
    run_trajectory_negative_controls(off_traj, adaptive_traj)

    print("RCG-04G E2E-MOVING-ADAPTIVE passed")


if __name__ == "__main__":
    main()
