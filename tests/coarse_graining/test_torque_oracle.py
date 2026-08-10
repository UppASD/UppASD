#!/usr/bin/env python3
"""Unit tests for the RCG-04D independent torque oracle (torque_oracle.py).

Run directly: ``python3 test_torque_oracle.py -v``. Also registered as the
CTest ``coarse-graining-torque-oracle`` (host-only Python, no UppASD/CUDA/HIP
dependency). The real-production calibration this module's docstring
describes is *not* re-run here (that requires the compiled binary and is a
one-time, disposable calibration step, not a repeatable unit test); this
file instead pins the exact calibrated numbers as regression constants, and
separately proves the oracle's analytic building blocks (sign convention,
degenerate/nondegenerate limits, minimum-image safety) from first principles.
"""

from __future__ import annotations

import pathlib
import unittest

import torque_oracle as orc
from moving_state_generator import Geometry, chiral_partner_pair, conical_spiral_state

E2E_JFILE = pathlib.Path(__file__).with_name("e2e") / "jfile"


def _host_geometry() -> Geometry:
    return Geometry(na=2, n1=6, n2=2, n3=2, basis=((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)))


class ParseJfileTests(unittest.TestCase):
    def test_parses_real_e2e_jfile(self) -> None:
        shells = orc.parse_jfile_shells(E2E_JFILE.read_text())
        self.assertEqual(len(shells), 4)
        distances = sorted(round(shell.distance, 6) for shell in shells)
        self.assertEqual(distances, [0.866025, 1.0, 1.414214, 1.658312])

    def test_rejects_wrong_field_count(self) -> None:
        with self.assertRaises(orc.MalformedJfileError):
            orc.parse_jfile_shells("1 1 0.5 0.5 0.5\n")

    def test_rejects_non_numeric_field(self) -> None:
        with self.assertRaises(orc.MalformedJfileError):
            orc.parse_jfile_shells("1 1 x 0.5 0.5 1.0\n")

    def test_rejects_empty_text(self) -> None:
        with self.assertRaises(orc.MalformedJfileError):
            orc.parse_jfile_shells("\n\n")

    def test_ignores_blank_lines(self) -> None:
        shells = orc.parse_jfile_shells("\n1 1 1.0 0.0 0.0 0.5\n\n")
        self.assertEqual(len(shells), 1)


class MinimumImageSafetyTests(unittest.TestCase):
    def test_real_e2e_geometry_is_safe(self) -> None:
        shells = orc.parse_jfile_shells(E2E_JFILE.read_text())
        orc.validate_minimum_image_safe(_host_geometry(), shells)  # no raise

    def test_oversized_shell_component_raises(self) -> None:
        geometry = _host_geometry()  # box (6, 2, 2): half-lengths (3, 1, 1)
        unsafe = (orc.ExchangeShell(dx=0.0, dy=1.5, dz=0.0, jij=1.0),)
        with self.assertRaises(orc.UnsafeMinimumImageError):
            orc.validate_minimum_image_safe(geometry, unsafe)

    def test_exactly_half_box_is_safe(self) -> None:
        geometry = _host_geometry()
        boundary = (orc.ExchangeShell(dx=0.0, dy=1.0, dz=0.0, jij=1.0),)
        orc.validate_minimum_image_safe(geometry, boundary)  # no raise


class GeometricBondsCalibrationTests(unittest.TestCase):
    """Pins the neighbour-list structure validated against real production output.

    See torque_oracle.py's module docstring for the two real ``sd.f95`` runs
    (uniform-aligned and sublattice-flip momfiles) whose printed
    ``localenergy.<simid>.out`` exchange energy these constants reproduce
    exactly (9 significant figures); this test only checks that this
    module's own geometric construction still gives that same topology, not
    that production still agrees with it (which was checked once,
    out-of-band, when the oracle was built).
    """

    def setUp(self) -> None:
        self.geometry = _host_geometry()
        self.shells = orc.parse_jfile_shells(E2E_JFILE.read_text())

    def test_every_atom_has_25_neighbours(self) -> None:
        coordination = orc.coordination_numbers(self.geometry, self.shells)
        self.assertEqual(len(coordination), 48)
        self.assertEqual(set(coordination.values()), {25})

    def test_atom_one_shell_breakdown(self) -> None:
        bonds = orc.build_geometric_bonds(self.geometry, self.shells)
        by_distance: dict[float, int] = {}
        distance_by_j = {round(shell.jij, 6): round(shell.distance, 6) for shell in self.shells}
        for _neighbour, jij in bonds[1]:
            key = distance_by_j[round(jij, 6)]
            by_distance[key] = by_distance.get(key, 0) + 1
        self.assertEqual(by_distance, {0.866025: 8, 1.0: 4, 1.414214: 5, 1.658312: 8})

    def test_uniform_aligned_energy_matches_real_production_run(self) -> None:
        uniform = {atom: (1.0, 0.0, 0.0) for atom, *_ in self.geometry.iter_atoms()}
        energies = orc.exchange_energy_per_atom(self.geometry, self.shells, uniform)
        self.assertEqual(len(energies), 48)
        for atom_energy in energies.values():
            self.assertAlmostEqual(atom_energy, -12.72518322837978, places=9)

    def test_sublattice_flip_energy_matches_real_production_run(self) -> None:
        flipped = {
            atom: ((0.0, 0.0, 1.0) if i0 == 1 else (1.0, 0.0, 0.0))
            for atom, i0, i1, i2, i3 in self.geometry.iter_atoms()
        }
        energies = orc.exchange_energy_per_atom(self.geometry, self.shells, flipped)
        for atom_energy in energies.values():
            self.assertAlmostEqual(atom_energy, -2.729371179633701, places=9)


class AnalyticLimitTests(unittest.TestCase):
    """Hand-derived two-atom toy models, independent of any jfile/production run."""

    def _toy_geometry(self) -> Geometry:
        return Geometry(na=1, n1=2, n2=1, n3=1, basis=((0.0, 0.0, 0.0),))

    def test_two_atom_perpendicular_spins_known_cross_product(self) -> None:
        geometry = self._toy_geometry()
        shells = (orc.ExchangeShell(dx=1.0, dy=0.0, dz=0.0, jij=1.0),)
        directions = {1: (0.0, 0.0, 1.0), 2: (1.0, 0.0, 0.0)}
        report = orc.initial_torque_report(geometry, shells, directions)
        # B_1 = J * S_2 = (1,0,0); torque_1 = S_1 x B_1 = (0,0,1)x(1,0,0) = (0,1,0)
        self.assertAlmostEqual(report.field_by_atom[1][0], 1.0, places=12)
        self.assertAlmostEqual(report.torque_by_atom[1][1], 1.0, places=12)
        self.assertAlmostEqual(orc._norm(report.torque_by_atom[1]), 1.0, places=12)
        # B_2 = J * S_1 = (0,0,1); torque_2 = S_2 x B_2 = (1,0,0)x(0,0,1) = (0,-1,0)
        self.assertAlmostEqual(report.torque_by_atom[2][1], -1.0, places=12)

    def test_two_atom_aligned_spins_zero_torque(self) -> None:
        geometry = self._toy_geometry()
        shells = (orc.ExchangeShell(dx=1.0, dy=0.0, dz=0.0, jij=1.0),)
        directions = {1: (0.0, 0.0, 1.0), 2: (0.0, 0.0, 1.0)}
        report = orc.initial_torque_report(geometry, shells, directions)
        self.assertAlmostEqual(report.max_torque, 0.0, places=12)
        self.assertAlmostEqual(report.rms_torque, 0.0, places=12)

    def test_two_atom_antiparallel_spins_zero_torque(self) -> None:
        # S x (-S) = 0 identically: antiparallel is also a stationary special case.
        geometry = self._toy_geometry()
        shells = (orc.ExchangeShell(dx=1.0, dy=0.0, dz=0.0, jij=1.0),)
        directions = {1: (0.0, 0.0, 1.0), 2: (0.0, 0.0, -1.0)}
        report = orc.initial_torque_report(geometry, shells, directions)
        self.assertAlmostEqual(report.max_torque, 0.0, places=12)

    def test_no_matching_shell_gives_zero_field_and_torque(self) -> None:
        geometry = self._toy_geometry()
        shells = (orc.ExchangeShell(dx=0.25, dy=0.0, dz=0.0, jij=5.0),)  # doesn't match distance 1
        directions = {1: (0.0, 0.0, 1.0), 2: (1.0, 0.0, 0.0)}
        report = orc.initial_torque_report(geometry, shells, directions)
        self.assertEqual(report.field_by_atom[1], (0.0, 0.0, 0.0))
        self.assertEqual(report.max_torque, 0.0)


class ConicalSpiralDegeneracyTests(unittest.TestCase):
    """Closes the RCG-04A/C open item on ``initmag_spin_spiral``'s zero torque.

    Uses the *same* generator, geometry, and jfile RCG-04D's own moving
    fixture uses, calling ``conical_spiral_state`` with
    ``degeneracy_tolerance_deg=0.0`` to let an exactly-degenerate cone angle
    through (0.0 < 0.0 is false, so ``_require_nondegenerate_cone`` does not
    raise) while still using the generator's own rotation math, rather than
    hand-building a direction map.
    """

    def setUp(self) -> None:
        self.geometry = _host_geometry()
        self.shells = orc.parse_jfile_shells(E2E_JFILE.read_text())

    def test_uniform_state_zero_torque(self) -> None:
        state = conical_spiral_state(
            self.geometry, cone_angle_deg=0.0, turns=1, degeneracy_tolerance_deg=0.0,
        )
        report = orc.initial_torque_report(self.geometry, self.shells, state.direction_by_atom)
        self.assertAlmostEqual(report.max_torque, 0.0, places=9)
        self.assertAlmostEqual(report.rms_torque, 0.0, places=9)

    def test_planar_spiral_zero_torque(self) -> None:
        """theta=90 degrees: the RCG-04A-flagged ``initmag_spin_spiral`` construction."""
        state = conical_spiral_state(
            self.geometry, cone_angle_deg=90.0, turns=1, degeneracy_tolerance_deg=0.0,
        )
        report = orc.initial_torque_report(self.geometry, self.shells, state.direction_by_atom)
        self.assertAlmostEqual(report.max_torque, 0.0, places=9)
        self.assertAlmostEqual(report.rms_torque, 0.0, places=9)

    def test_generic_cone_angle_nonzero_torque(self) -> None:
        """The RCG-04D fixture's actual parameters: genuinely nonstationary."""
        state = conical_spiral_state(self.geometry, cone_angle_deg=40.0, turns=1)
        report = orc.initial_torque_report(self.geometry, self.shells, state.direction_by_atom)
        self.assertGreater(report.max_torque, 1.0e-3)
        self.assertGreater(report.rms_torque, 1.0e-3)

    def test_zero_turns_is_rejected_by_generator_before_reaching_oracle(self) -> None:
        with self.assertRaises(Exception):
            conical_spiral_state(
                self.geometry, cone_angle_deg=40.0, turns=0, degeneracy_tolerance_deg=0.0,
            )


class TorqueReportSerializationTests(unittest.TestCase):
    def test_as_dict_is_json_serializable(self) -> None:
        import json

        geometry = Geometry(na=1, n1=2, n2=1, n3=1, basis=((0.0, 0.0, 0.0),))
        shells = (orc.ExchangeShell(dx=1.0, dy=0.0, dz=0.0, jij=1.0),)
        directions = {1: (0.0, 0.0, 1.0), 2: (1.0, 0.0, 0.0)}
        report = orc.initial_torque_report(geometry, shells, directions)
        json.dumps(report.as_dict())  # must not raise


class DmiOracleTests(unittest.TestCase):
    """RCG-04H: independent DMI field/energy oracle, cross-checked against
    RCG-02's own accepted worked dimer example
    (docs/RCG-02_DMI_HANDEDNESS_EVIDENCE.md), independent of the harness
    this oracle will gate."""

    def test_rcg02_dimer_cross_check(self) -> None:
        # RCG-02: M_1=+x, M_2=+y, D_12=+D*z, D_21=-D*z -> E_reduced=D,
        # B_1=-D*x, B_2=-D*y.
        geometry = Geometry(na=1, n1=2, n2=1, n3=1, basis=((0.0, 0.0, 0.0),))
        d = 2.5
        dmi_bonds = (
            orc.DmiBond(dx=1.0, dy=0.0, dz=0.0, d=(0.0, 0.0, d)),
            orc.DmiBond(dx=-1.0, dy=0.0, dz=0.0, d=(0.0, 0.0, -d)),
        )
        directions = {1: (1.0, 0.0, 0.0), 2: (0.0, 1.0, 0.0)}
        field = orc.dmi_field(geometry, dmi_bonds, directions)
        self.assertAlmostEqual(field[1][0], -d, places=12)
        self.assertAlmostEqual(field[1][1], 0.0, places=12)
        self.assertAlmostEqual(field[1][2], 0.0, places=12)
        self.assertAlmostEqual(field[2][1], -d, places=12)
        energy = orc.dmi_energy_reduced(geometry, dmi_bonds, directions)
        self.assertAlmostEqual(energy, d, places=12)

    def test_negate_dmi_bonds_reverses_field_and_energy(self) -> None:
        geometry = Geometry(na=1, n1=2, n2=1, n3=1, basis=((0.0, 0.0, 0.0),))
        dmi_bonds = (
            orc.DmiBond(dx=1.0, dy=0.0, dz=0.0, d=(0.0, 0.0, 1.0)),
            orc.DmiBond(dx=-1.0, dy=0.0, dz=0.0, d=(0.0, 0.0, -1.0)),
        )
        directions = {1: (1.0, 0.0, 0.0), 2: (0.0, 1.0, 0.0)}
        energy = orc.dmi_energy_reduced(geometry, dmi_bonds, directions)
        reversed_bonds = orc.negate_dmi_bonds(dmi_bonds)
        reversed_energy = orc.dmi_energy_reduced(geometry, reversed_bonds, directions)
        self.assertAlmostEqual(reversed_energy, -energy, places=12)

    def test_planar_chiral_partner_energy_ordering_matches_rcg02(self) -> None:
        """Positive D_zx must favour -q over +q, matching RCG-02's own stated
        result and the accepted DMI-HYBRID-CROSSING operator fixture --
        independent of the run_moving_dmi_chiral.py harness this gates."""
        geometry = _host_geometry()
        d = 0.02
        dmi_bonds = (
            orc.DmiBond(dx=1.0, dy=0.0, dz=0.0, d=(0.0, 0.0, d)),
            orc.DmiBond(dx=-1.0, dy=0.0, dz=0.0, d=(0.0, 0.0, -d)),
        )
        plus, minus = chiral_partner_pair(
            geometry, cone_angle_deg=40.0, turns=1, axis=(0.0, 0.0, 1.0),
        )
        energy_plus = orc.dmi_energy_reduced(geometry, dmi_bonds, plus.direction_by_atom)
        energy_minus = orc.dmi_energy_reduced(geometry, dmi_bonds, minus.direction_by_atom)
        self.assertLess(energy_minus, energy_plus)
        # Reversed operator must flip the ordering, not merely change magnitude.
        reversed_bonds = orc.negate_dmi_bonds(dmi_bonds)
        reversed_energy_plus = orc.dmi_energy_reduced(geometry, reversed_bonds, plus.direction_by_atom)
        reversed_energy_minus = orc.dmi_energy_reduced(geometry, reversed_bonds, minus.direction_by_atom)
        self.assertGreater(reversed_energy_minus, reversed_energy_plus)

    def test_parse_dmfile_bonds_roundtrip(self) -> None:
        text = "1 1 1.0 0.0 0.0 0.0 0.0 0.02\n1 1 -1.0 0.0 0.0 0.0 0.0 -0.02\n"
        bonds = orc.parse_dmfile_bonds(text)
        self.assertEqual(len(bonds), 2)
        self.assertEqual(bonds[0].d, (0.0, 0.0, 0.02))
        self.assertEqual(bonds[1].d, (0.0, 0.0, -0.02))

    def test_parse_dmfile_bonds_rejects_wrong_field_count(self) -> None:
        with self.assertRaises(orc.MalformedJfileError):
            orc.parse_dmfile_bonds("1 1 1.0 0.0 0.0 0.02\n")

    def test_parse_dmfile_bonds_rejects_empty_text(self) -> None:
        with self.assertRaises(orc.MalformedJfileError):
            orc.parse_dmfile_bonds("\n\n")

    def test_anisotropy_energy_reduced_uniaxial(self) -> None:
        geometry = Geometry(na=1, n1=1, n2=1, n3=1, basis=((0.0, 0.0, 0.0),))
        directions = {1: (1.0, 0.0, 0.0)}
        energy = orc.anisotropy_energy_reduced(
            geometry, directions, axis=(1.0, 0.0, 0.0), k1_by_site={1: -0.002},
        )
        self.assertAlmostEqual(energy, -0.002, places=12)

    def test_dmi_anisotropy_torque_report_nonzero_for_displaced_state(self) -> None:
        geometry = _host_geometry()
        shells = orc.parse_jfile_shells(E2E_JFILE.read_text())
        d = 0.02
        dmi_bonds = (
            orc.DmiBond(dx=1.0, dy=0.0, dz=0.0, d=(0.0, 0.0, d)),
            orc.DmiBond(dx=-1.0, dy=0.0, dz=0.0, d=(0.0, 0.0, -d)),
        )
        state = conical_spiral_state(geometry, cone_angle_deg=40.0, turns=1)
        report = orc.dmi_anisotropy_torque_report(
            geometry, shells, dmi_bonds, state.direction_by_atom,
            anisotropy_axis=(1.0, 0.0, 0.0), anisotropy_k1={1: -0.002, 2: -0.003},
        )
        self.assertGreater(report.max_torque, 1.0e-3)
        self.assertGreater(report.rms_torque, 1.0e-3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
