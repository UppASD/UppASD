#!/usr/bin/env python3
"""RCG-04C tests for trajectory_evidence.py: parser, metrics, and negative controls."""

from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path

from moving_state_generator import Geometry
import trajectory_evidence as te


RULE = "#" * 80


def header(file_type: str, natom: int, mensemble: int, mode: str = "S") -> str:
    return "\n".join([
        RULE,
        f"# File type: {file_type}",
        f"# Simulation type: {mode}",
        f"# Number of atoms: {natom}",
        f"# Number of ensembles: {mensemble}",
        RULE,
        "#iter     ens   iatom            |Mom|             M_x             M_y             M_z",
    ])


def row(step: int, ens: int, atom: int, moment: float, direction: tuple[float, float, float]) -> str:
    mx, my, mz = direction
    return f"{step:8d}{ens:8d}{atom:8d}  {moment:16.8e}{mx:16.8e}{my:16.8e}{mz:16.8e}"


def uniform_step_rows(step: int, natom: int, mensemble: int,
                       direction: tuple[float, float, float] = (1.0, 0.0, 0.0),
                       moment: float = 2.23, order: str = "forward") -> list[str]:
    keys = [(ens, atom) for ens in range(1, mensemble + 1) for atom in range(1, natom + 1)]
    if order == "reverse":
        keys = list(reversed(keys))
    return [row(step, ens, atom, moment, direction) for ens, atom in keys]


def make_text(file_type: str, natom: int, mensemble: int, step_rows: list[list[str]],
              mode: str = "S") -> str:
    lines = [header(file_type, natom, mensemble, mode)]
    for rows in step_rows:
        lines.extend(rows)
    return "\n".join(lines) + "\n"


class BasicParsingTests(unittest.TestCase):
    def test_restart_style_single_step(self):
        text = make_text("R", 2, 1, [uniform_step_rows(0, 2, 1)])
        trajectory = te.parse_mag_conf_text(text, source="synthetic")
        self.assertEqual(trajectory.file_type, "R")
        self.assertEqual(trajectory.natom, 2)
        self.assertEqual(trajectory.mensemble, 1)
        self.assertEqual(trajectory.step_numbers(), (0,))
        self.assertAlmostEqual(trajectory.direction(0, 1, 1)[0], 1.0)

    def test_moment_style_multi_step(self):
        text = make_text("M", 3, 1, [
            uniform_step_rows(0, 3, 1),
            uniform_step_rows(10, 3, 1),
            uniform_step_rows(20, 3, 1),
        ])
        trajectory = te.parse_mag_conf_text(text, source="synthetic")
        self.assertEqual(trajectory.step_numbers(), (0, 10, 20))
        self.assertEqual(len(trajectory.steps[0].records), 3)

    def test_multi_ensemble(self):
        text = make_text("M", 2, 2, [uniform_step_rows(0, 2, 2)])
        trajectory = te.parse_mag_conf_text(text, source="synthetic")
        self.assertEqual(len(trajectory.steps[0].records), 4)
        self.assertIn((2, 2), trajectory.steps[0].records)

    def test_moment_magnitude_preserved(self):
        text = make_text("M", 1, 1, [uniform_step_rows(0, 1, 1, moment=3.14159)])
        trajectory = te.parse_mag_conf_text(text, source="synthetic")
        self.assertAlmostEqual(trajectory.steps[0].records[(1, 1)].moment, 3.14159, places=6)

    def test_fortran_d_exponent_accepted(self):
        text = header("M", 1, 1) + "\n" + "       0       1       1  2.23000000D+00  1.00000000D+00  0.00000000D+00  0.00000000D+00\n"
        trajectory = te.parse_mag_conf_text(text, source="synthetic")
        self.assertAlmostEqual(trajectory.steps[0].records[(1, 1)].moment, 2.23)


class ReorderingIsHandledByExplicitIndexTests(unittest.TestCase):
    def test_within_step_reordering_gives_identical_trajectory(self):
        forward = make_text("M", 4, 1, [uniform_step_rows(0, 4, 1, direction=(0.0, 1.0, 0.0), order="forward")])
        reverse = make_text("M", 4, 1, [uniform_step_rows(0, 4, 1, direction=(0.0, 1.0, 0.0), order="reverse")])
        self.assertNotEqual(forward, reverse)
        t_forward = te.parse_mag_conf_text(forward, source="forward")
        t_reverse = te.parse_mag_conf_text(reverse, source="reverse")
        self.assertEqual(t_forward.steps[0].records, t_reverse.steps[0].records)


class CorruptionNegativeControlTests(unittest.TestCase):
    """One test per corruption category the RCG-04C prompt requires rejecting."""

    def test_truncated_short_header(self):
        with self.assertRaises(te.TruncatedTrajectoryFileError):
            te.parse_mag_conf_text("\n".join([RULE, "# File type: M"]), source="x")

    def test_truncated_no_data_rows(self):
        text = header("M", 1, 1) + "\n"
        with self.assertRaises(te.TruncatedTrajectoryFileError):
            te.parse_mag_conf_text(text, source="x")

    def test_truncated_partial_data_row(self):
        text = header("M", 1, 1) + "\n" + "       0       1       1  2.23000000E+00  1.00000000E+00\n"
        with self.assertRaises(te.TruncatedTrajectoryFileError):
            te.parse_mag_conf_text(text, source="x")

    def test_duplicate_atom_record(self):
        rows = uniform_step_rows(0, 2, 1) + [row(0, 1, 1, 2.23, (0.0, 1.0, 0.0))]
        text = make_text("M", 2, 1, [rows])
        with self.assertRaises(te.DuplicateRecordError):
            te.parse_mag_conf_text(text, source="x")

    def test_missing_step_data(self):
        rows = uniform_step_rows(0, 4, 1)[:-1]  # drop the last atom's record
        text = make_text("M", 4, 1, [rows])
        with self.assertRaises(te.MissingStepDataError):
            te.parse_mag_conf_text(text, source="x")

    def test_inconsistent_atom_count_out_of_range(self):
        rows = uniform_step_rows(0, 2, 1) + [row(0, 1, 5, 2.23, (1.0, 0.0, 0.0))]
        text = make_text("M", 2, 1, [rows])
        with self.assertRaises(te.InconsistentAtomCountError):
            te.parse_mag_conf_text(text, source="x")

    def test_non_finite_value_rejected(self):
        text = header("M", 1, 1) + "\n" + \
            f"{0:8d}{1:8d}{1:8d}  {'NaN':>16}{1.0:16.8e}{0.0:16.8e}{0.0:16.8e}\n"
        with self.assertRaises(te.NonFiniteValueError):
            te.parse_mag_conf_text(text, source="x")

    def test_non_unit_direction_rejected(self):
        text = make_text("M", 1, 1, [[row(0, 1, 1, 2.23, (2.0, 0.0, 0.0))]])
        with self.assertRaises(te.NonUnitDirectionError):
            te.parse_mag_conf_text(text, source="x")

    def test_non_monotonic_steps_interleaved(self):
        # Step 10's block appears, then a step-0 record is appended afterward:
        # positional/append-only parsing would silently treat this as valid.
        rows_10 = uniform_step_rows(10, 2, 1)
        rows_0 = uniform_step_rows(0, 2, 1)
        text = make_text("M", 2, 1, [rows_10, rows_0])
        with self.assertRaises(te.NonMonotonicStepError):
            te.parse_mag_conf_text(text, source="x")

    def test_missing_step_gap_inconsistent_with_cadence(self):
        # Cadence 10 inferred from steps 0->10; step 30 (not 20) is next: a
        # sampled step was silently dropped.
        text = make_text("M", 2, 1, [
            uniform_step_rows(0, 2, 1), uniform_step_rows(10, 2, 1), uniform_step_rows(30, 2, 1),
        ])
        with self.assertRaises(te.MissingStepError):
            te.parse_mag_conf_text(text, source="x")

    def test_ambiguous_simulation_identifier(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "restart.caseA.out").write_text("x")
            (directory / "restart.caseB.out").write_text("x")
            with self.assertRaises(te.AmbiguousSimulationIdentifierError):
                te.find_unique_output(directory, "restart.*.out")

    def test_missing_output_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(te.MissingOutputFileError):
                te.find_unique_output(Path(tmp), "restart.*.out")

    def test_unique_output_file_is_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "restart.only.out").write_text("x")
            found = te.find_unique_output(directory, "restart.*.out")
            self.assertEqual(found.name, "restart.only.out")


class LoadHelperTests(unittest.TestCase):
    def test_load_restart_state_requires_type_r(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            text = make_text("M", 1, 1, [uniform_step_rows(0, 1, 1)])
            (directory / "restart.x.out").write_text(text)
            with self.assertRaises(te.MalformedHeaderError):
                te.load_restart_state(directory)

    def test_load_restart_state_rejects_multi_step(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            text = make_text("R", 1, 1, [uniform_step_rows(0, 1, 1), uniform_step_rows(10, 1, 1)])
            (directory / "restart.x.out").write_text(text)
            with self.assertRaises(te.TrajectoryStageError):
                te.load_restart_state(directory)

    def test_load_moment_trajectory_requires_type_m(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            text = make_text("R", 1, 1, [uniform_step_rows(0, 1, 1)])
            (directory / "moment.x.out").write_text(text)
            with self.assertRaises(te.MalformedHeaderError):
                te.load_moment_trajectory(directory)

    def test_load_moment_trajectory_with_explicit_simid(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            text = make_text("M", 1, 1, [uniform_step_rows(0, 1, 1)])
            (directory / "moment.abc.out").write_text(text)
            (directory / "moment.def.out").write_text(text)
            trajectory = te.load_moment_trajectory(directory, simid="abc")
            self.assertEqual(trajectory.step_numbers(), (0,))
            with self.assertRaises(te.AmbiguousSimulationIdentifierError):
                te.load_moment_trajectory(directory)


class ComponentAndAngularErrorTests(unittest.TestCase):
    def _trajectory(self, direction):
        text = make_text("M", 1, 1, [uniform_step_rows(0, 1, 1, direction=direction)])
        return te.parse_mag_conf_text(text, source="t")

    def test_identical_trajectories_have_zero_error(self):
        a = self._trajectory((1.0, 0.0, 0.0))
        b = self._trajectory((1.0, 0.0, 0.0))
        summary = te.component_trajectory_error(a, b)
        self.assertEqual(summary.max_abs_error_overall, 0.0)
        angular = te.angular_trajectory_error(a, b)
        self.assertEqual(angular.max_radians, 0.0)

    def test_perpendicular_vectors_have_known_error(self):
        a = self._trajectory((1.0, 0.0, 0.0))
        b = self._trajectory((0.0, 1.0, 0.0))
        summary = te.component_trajectory_error(a, b)
        self.assertAlmostEqual(summary.max_abs_error_overall, 1.0)
        angular = te.angular_trajectory_error(a, b)
        self.assertAlmostEqual(angular.max_radians, math.pi / 2.0, places=10)

    def test_antiparallel_vectors_stable_at_pi(self):
        # acos(u.v) for u.v very close to -1 loses precision; the stable
        # atan2 formula must still return pi to high precision.
        a = self._trajectory((1.0, 0.0, 0.0))
        b = self._trajectory((-1.0, 0.0, 0.0))
        angular = te.angular_trajectory_error(a, b)
        self.assertAlmostEqual(angular.max_radians, math.pi, places=12)

    def test_near_parallel_vectors_stable_near_zero(self):
        epsilon = 1.0e-9
        u = (math.sqrt(1.0 - epsilon**2), epsilon, 0.0)
        a = self._trajectory((1.0, 0.0, 0.0))
        b = self._trajectory(u)
        angular = te.angular_trajectory_error(a, b)
        self.assertGreater(angular.max_radians, 0.0)
        self.assertAlmostEqual(angular.max_radians, epsilon, delta=1.0e-6)

    def test_key_mismatch_is_rejected(self):
        a_text = make_text("M", 2, 1, [uniform_step_rows(0, 2, 1)])
        b_text = make_text("M", 1, 1, [uniform_step_rows(0, 1, 1)])
        a = te.parse_mag_conf_text(a_text, source="a")
        b = te.parse_mag_conf_text(b_text, source="b")
        with self.assertRaises(te.TrajectoryKeyMismatchError):
            te.component_trajectory_error(a, b)


class DisplacementTests(unittest.TestCase):
    def test_rotating_spin_has_expected_final_and_max_displacement(self):
        # A single spin rotating in the xy-plane by 30 degrees per sampled
        # step, for 4 steps: 0, 30, 60, 90 degrees from the initial state.
        step_rows = []
        for n, angle_deg in enumerate((0, 30, 60, 90)):
            angle = math.radians(angle_deg)
            step_rows.append(uniform_step_rows(n * 10, 1, 1, direction=(math.cos(angle), math.sin(angle), 0.0)))
        text = make_text("M", 1, 1, step_rows)
        trajectory = te.parse_mag_conf_text(text, source="t")
        summary = te.spin_displacement(trajectory)
        self.assertAlmostEqual(summary.max_final_displacement, math.radians(90), places=6)
        self.assertAlmostEqual(summary.max_over_time_displacement_overall, math.radians(90), places=6)

    def test_stationary_trajectory_has_zero_displacement(self):
        text = make_text("M", 1, 1, [uniform_step_rows(0, 1, 1), uniform_step_rows(10, 1, 1)])
        trajectory = te.parse_mag_conf_text(text, source="t")
        summary = te.spin_displacement(trajectory)
        self.assertEqual(summary.max_final_displacement, 0.0)

    def test_displacement_requires_multiple_steps(self):
        text = make_text("M", 1, 1, [uniform_step_rows(0, 1, 1)])
        trajectory = te.parse_mag_conf_text(text, source="t")
        with self.assertRaises(te.TrajectoryStageError):
            te.spin_displacement(trajectory)


class EnergyFieldSeriesTests(unittest.TestCase):
    # Fixed-width es24.16 padding, exactly as real production stdout emits
    # it (a negative value gets one leading space before its '-', a
    # positive value gets two): this is the format
    # ``print_adaptive_cg_summary`` actually writes (verified against a
    # real production run in this slice's evidence), not a simplified
    # no-space synthetic shorthand -- a prior version of this test used
    # a no-space variant and passed against a parser bug that failed on
    # real production output.
    SAMPLE_STDOUT = (
        "AdaptiveCG: last_energy_j atomistic_bilinear= -1.2340000000000000E-20"
        " atomistic_onsite=  2.0000000000000000E-21 coarse_exchange=  0.0000000000000000E+00"
        " coarse_spiralization=  0.0000000000000000E+00\n"
        "AdaptiveCG: last_energy_j coarse_anisotropy=  0.0000000000000000E+00"
        " coarse_external=  0.0000000000000000E+00 coarse_dipole=  0.0000000000000000E+00\n"
        "AdaptiveCG: last_total_energy_j= -1.1340000000000000E-20\n"
        "AdaptiveCG: last_field_checksums_t atom_sum=  1.0000000000000000E+00"
        " atom_norm2=  2.0000000000000000E+00 coarse_sum=  0.0000000000000000E+00"
        " coarse_norm2=  0.0000000000000000E+00\n"
    )

    def test_single_emission_parsed(self):
        samples = te.parse_energy_field_series(self.SAMPLE_STDOUT)
        self.assertEqual(len(samples), 1)
        sample = samples[0]
        self.assertAlmostEqual(sample.energies_j["atomistic_bilinear"], -1.234e-20)
        self.assertAlmostEqual(sample.energies_j["coarse_dipole"], 0.0)
        self.assertAlmostEqual(sample.energies_j["total"], -1.134e-20)
        self.assertAlmostEqual(sample.field_checksums_t["atom_sum"], 1.0)

    def test_gpu_prefix_parsed_identically(self):
        gpu_stdout = self.SAMPLE_STDOUT.replace("AdaptiveCG:", "Gpu: AdaptiveCG")
        samples = te.parse_energy_field_series(gpu_stdout)
        self.assertEqual(len(samples), 1)

    def test_multi_emission_forward_compatible(self):
        doubled = self.SAMPLE_STDOUT + self.SAMPLE_STDOUT
        samples = te.parse_energy_field_series(doubled)
        self.assertEqual(len(samples), 2)

    def test_compare_identical_series_zero_error(self):
        samples = te.parse_energy_field_series(self.SAMPLE_STDOUT)
        comparison = te.compare_energy_field_series(samples, samples)
        for term in comparison["samples"][0]["energies_j"].values():
            self.assertEqual(term["abs_error"], 0.0)

    def test_compare_mismatched_length_rejected(self):
        one = te.parse_energy_field_series(self.SAMPLE_STDOUT)
        two = te.parse_energy_field_series(self.SAMPLE_STDOUT + self.SAMPLE_STDOUT)
        with self.assertRaises(te.TrajectoryStageError):
            te.compare_energy_field_series(one, two)


class TransitionEventParsingTests(unittest.TestCase):
    # Two trailing spaces before 'polarization_ratio=', matching real
    # production output (verified in this slice's evidence): a prior
    # version of this fixture used a single space and passed against a
    # parser bug that failed to match real production stdout.
    LINE = (
        "AdaptiveCG: transition step=12 block=3 old_state=0 new_state=2 accepted=T "
        "reason=misalignment outcome=applied "
        "energies_j=  1.2340000000000000E-21  1.2500000000000000E-21  1.6000000000000000E-23  "
        "polarization_ratio=  5.0000000000000000E-01"
    )
    REJECTED_LINE = (
        "AdaptiveCG: transition step=13 block=4 old_state=2 new_state=0 accepted=F "
        "reason=polarization-unsafe outcome=rejected "
        "energies_j=  1.0000000000000000E-21  1.0000000000000000E-21  0.0000000000000000E+00  "
        "polarization_ratio=  4.0000000000000000E-01"
    )

    def test_single_event_fields(self):
        events = te.parse_transition_events(self.LINE)
        self.assertEqual(len(events), 1)
        event = events[0]
        self.assertEqual(event.step, 12)
        self.assertEqual(event.block, 3)
        self.assertEqual(event.old_state, 0)
        self.assertEqual(event.new_state, 2)
        self.assertTrue(event.accepted)
        self.assertEqual(event.reason, "misalignment")
        self.assertEqual(event.outcome, "applied")
        self.assertAlmostEqual(event.energy_jump_j, 1.6e-23)
        self.assertAlmostEqual(event.polarization_ratio, 0.5)

    def test_rejected_event_accepted_false(self):
        events = te.parse_transition_events(self.REJECTED_LINE)
        self.assertFalse(events[0].accepted)
        self.assertEqual(events[0].reason, "polarization-unsafe")

    def test_multiple_events_in_order(self):
        stdout = self.LINE + "\n" + self.REJECTED_LINE + "\n"
        events = te.parse_transition_events(stdout)
        self.assertEqual([e.step for e in events], [12, 13])

    def test_no_events_returns_empty_list(self):
        self.assertEqual(te.parse_transition_events("nothing relevant here"), [])


class ResolutionStateHistoryTests(unittest.TestCase):
    def test_multiple_labels_and_steps_in_order(self):
        stdout = "\n".join([
            "AdaptiveCG: resolution_state label=initial step=0 values=2,2,0,0",
            "AdaptiveCG: resolution_state label=step step=5 values=2,1,0,0",
            "AdaptiveCG: resolution_state label=step step=10 values=2,2,0,0",
            "AdaptiveCG: resolution_state label=final step=10 values=2,2,0,0",
        ])
        samples = te.parse_resolution_state_history(stdout)
        self.assertEqual(len(samples), 4)
        self.assertEqual([s.label for s in samples], ["initial", "step", "step", "final"])
        self.assertEqual(samples[1].values, (2, 1, 0, 0))

    def test_only_final_is_a_subset_of_full_history(self):
        stdout = "AdaptiveCG: resolution_state label=final step=2 values=0,0\n"
        samples = te.parse_resolution_state_history(stdout)
        self.assertEqual(len(samples), 1)
        self.assertEqual(samples[0].label, "final")


class ConicalModeTests(unittest.TestCase):
    def _geometry(self, n1: int) -> Geometry:
        return Geometry(na=1, n1=n1, n2=1, n3=1, basis=((0.0, 0.0, 0.0),))

    def test_amplitude_and_frequency_recovered_from_synthetic_rotation(self):
        n1 = 12
        geometry = self._geometry(n1)
        cone_theta = math.radians(45.0)
        turns = 1
        omega_true = 0.1  # radians per step
        phase_by_atom = {}
        for atom, i0, i1, i2, i3 in geometry.iter_atoms():
            phase_by_atom[atom] = 2.0 * math.pi * turns * i1 / n1

        step_rows = []
        steps = list(range(0, 50, 5))
        for step in steps:
            rows = []
            for atom, i0, i1, i2, i3 in geometry.iter_atoms():
                local_phase = phase_by_atom[atom] + omega_true * step
                direction = (
                    math.sin(cone_theta) * math.cos(local_phase),
                    math.sin(cone_theta) * math.sin(local_phase),
                    math.cos(cone_theta),
                )
                rows.append(row(step, 1, atom, 2.23, direction))
            step_rows.append(rows)
        text = make_text("M", geometry.natom, 1, step_rows)
        trajectory = te.parse_mag_conf_text(text, source="t")

        samples = te.conical_mode_series(
            trajectory, axis=(0.0, 0.0, 1.0), in_plane_reference=(1.0, 0.0, 0.0),
            phase_by_atom=phase_by_atom,
        )
        self.assertEqual(len(samples), len(steps))
        for sample in samples:
            self.assertAlmostEqual(sample.amplitude, math.sin(cone_theta), places=6)

        fit = te.fit_conical_mode_frequency(samples)
        self.assertAlmostEqual(fit.angular_frequency, omega_true, places=6)
        self.assertLess(fit.residual_rms, 1.0e-8)

    def test_stationary_mode_has_zero_frequency(self):
        geometry = self._geometry(4)
        phase_by_atom = {atom: 0.0 for atom, *_ in geometry.iter_atoms()}
        step_rows = [uniform_step_rows(n * 10, geometry.natom, 1, direction=(1.0, 0.0, 0.0)) for n in range(3)]
        text = make_text("M", geometry.natom, 1, step_rows)
        trajectory = te.parse_mag_conf_text(text, source="t")
        samples = te.conical_mode_series(
            trajectory, axis=(0.0, 0.0, 1.0), in_plane_reference=(0.0, 1.0, 0.0),
            phase_by_atom=phase_by_atom,
        )
        fit = te.fit_conical_mode_frequency(samples)
        self.assertAlmostEqual(fit.angular_frequency, 0.0, places=10)


class ChiralityTests(unittest.TestCase):
    def _geometry(self, n1: int) -> Geometry:
        return Geometry(na=1, n1=n1, n2=1, n3=1, basis=((0.0, 0.0, 0.0),))

    def _chain_step(self, geometry, turns, cone_theta_deg=90.0):
        theta = math.radians(cone_theta_deg)
        rows = []
        for atom, i0, i1, i2, i3 in geometry.iter_atoms():
            phase = 2.0 * math.pi * turns * i1 / geometry.n1
            direction = (math.sin(theta) * math.cos(phase), math.sin(theta) * math.sin(phase), math.cos(theta))
            rows.append(row(0, 1, atom, 2.23, direction))
        text = make_text("M", geometry.natom, 1, [rows])
        return te.parse_mag_conf_text(text, source="t").steps[0]

    def test_right_and_left_handed_chains_have_opposite_sign(self):
        geometry = self._geometry(8)
        bonds = te.axis_chain_bonds(geometry, axis_cell_index=0)
        right = self._chain_step(geometry, turns=1)
        left = self._chain_step(geometry, turns=-1)
        chi_right = te.signed_chirality(right, directed_bonds=bonds, axis=(0.0, 0.0, 1.0))
        chi_left = te.signed_chirality(left, directed_bonds=bonds, axis=(0.0, 0.0, 1.0))
        self.assertGreater(chi_right, 0.0)
        self.assertLess(chi_left, 0.0)
        self.assertAlmostEqual(chi_right, -chi_left, places=10)

    def test_reversed_axis_flips_sign(self):
        geometry = self._geometry(8)
        bonds = te.axis_chain_bonds(geometry, axis_cell_index=0)
        state = self._chain_step(geometry, turns=1)
        chi_plus = te.signed_chirality(state, directed_bonds=bonds, axis=(0.0, 0.0, 1.0))
        chi_minus = te.signed_chirality(state, directed_bonds=bonds, axis=(0.0, 0.0, -1.0))
        self.assertAlmostEqual(chi_plus, -chi_minus, places=10)

    def test_requires_at_least_one_bond(self):
        geometry = self._geometry(4)
        state = self._chain_step(geometry, turns=1)
        with self.assertRaises(te.TrajectoryStageError):
            te.signed_chirality(state, directed_bonds=[], axis=(0.0, 0.0, 1.0))


class DomainWallTests(unittest.TestCase):
    def _geometry(self, n1: int) -> Geometry:
        return Geometry(na=1, n1=n1, n2=1, n3=1, basis=((0.0, 0.0, 0.0),))

    def _wall_pair_step(self, geometry, centers, width, step=0):
        rows = []
        for atom, i0, i1, i2, i3 in geometry.iter_atoms():
            u = float(i1)
            c1, c2 = centers
            theta = 2.0 * math.atan(math.exp((u - c1) / width)) - 2.0 * math.atan(math.exp((u - c2) / width))
            direction = (math.sin(theta), 0.0, math.cos(theta))
            rows.append(row(step, 1, atom, 2.23, direction))
        text = make_text("M", geometry.natom, 1, [rows])
        return te.parse_mag_conf_text(text, source="t").steps[0]

    def test_two_wall_centers_detected_near_expected_positions(self):
        geometry = self._geometry(40)
        step = self._wall_pair_step(geometry, centers=(10.0, 30.0), width=1.5)
        centers = te.domain_wall_centers(step, geometry, axis_cell_index=0, easy_axis=(0.0, 0.0, 1.0))
        self.assertEqual(len(centers), 2)
        self.assertAlmostEqual(min(centers), 10.0, delta=0.2)
        self.assertAlmostEqual(max(centers), 30.0, delta=0.2)

    def test_wall_crossing_detected_when_wall_moves_across_block_boundary(self):
        series = [
            (0, [9.6, 30.0]),
            (5, [10.4, 30.0]),
            (10, [11.6, 30.0]),
        ]
        summary = te.track_wall_crossings(series, axis_length_cells=40.0, block_length_cells=1.0)
        self.assertEqual(summary.displacement, 2.0)
        self.assertTrue(any(step == 10 for step, _, _ in summary.crossing_events))

    def test_periodic_unwrap_does_not_report_spurious_large_jump(self):
        # A wall drifting past the periodic boundary (39.5 -> 0.5) should be
        # unwrapped to a small displacement, not a ~39-cell jump.
        series = [(0, [39.5]), (1, [0.5])]
        summary = te.track_wall_crossings(series, axis_length_cells=40.0, block_length_cells=1.0)
        self.assertAlmostEqual(summary.displacement, 1.0, places=6)

    def test_center_series_must_be_increasing_step_order(self):
        series = [(10, [1.0]), (0, [1.0])]
        with self.assertRaises(te.TrajectoryStageError):
            te.track_wall_crossings(series, axis_length_cells=40.0, block_length_cells=1.0)


class SummaryTests(unittest.TestCase):
    def test_trajectory_summary_is_json_serializable(self):
        text = make_text("M", 2, 1, [uniform_step_rows(0, 2, 1), uniform_step_rows(10, 2, 1)])
        trajectory = te.parse_mag_conf_text(text, source="t")
        summary = te.trajectory_summary(trajectory, provenance={"commit": "deadbeef"})
        serialized = json.dumps(summary)
        reloaded = json.loads(serialized)
        self.assertEqual(reloaded["natom"], 2)
        self.assertEqual(reloaded["provenance"]["commit"], "deadbeef")
        self.assertEqual(reloaded["step_numbers"], [0, 10])


if __name__ == "__main__":
    unittest.main(verbosity=2)
