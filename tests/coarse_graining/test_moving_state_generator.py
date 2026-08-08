#!/usr/bin/env python3
"""RCG-04B focused tests for the deterministic moving-state generators.

These are generator-level unit tests only: they check properties of the
generated momfile/restart text and per-atom direction arrays directly, in
Python, independent of the UppASD executable. No production dynamics is run
or claimed here (see ``docs/RCG-04_MOVING_E2E_EVIDENCE.md``, RCG-04B).
"""
from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import moving_state_generator as gen  # noqa: E402


# Matches tests/coarse_graining/e2e/posfile: two basis atoms at (0,0,0) and
# (0.5,0.5,0.5), the geometry shared by every existing e2e fixture.
STANDARD_BASIS = ((0.0, 0.0, 0.0), (0.5, 0.5, 0.5))


def standard_geometry(n1: int = 6, n2: int = 2, n3: int = 2) -> gen.Geometry:
    return gen.Geometry(na=2, n1=n1, n2=n2, n3=n3, basis=STANDARD_BASIS)


class GeometryAtomOrderingTests(unittest.TestCase):
    def test_atom_indices_are_a_permutation_of_1_to_natom(self) -> None:
        geometry = standard_geometry()
        indices = [atom for atom, *_ in geometry.iter_atoms()]
        self.assertEqual(sorted(indices), list(range(1, geometry.natom + 1)))

    def test_iteration_order_matches_basis_fastest_then_i1_i2_i3(self) -> None:
        geometry = gen.Geometry(na=2, n1=3, n2=1, n3=1, basis=STANDARD_BASIS)
        produced = [(atom, i0, i1, i2, i3) for atom, i0, i1, i2, i3 in geometry.iter_atoms()]
        # i0=1,i1=0 -> 1 ; i0=2,i1=0 -> 2 ; i0=1,i1=1 -> 3 ; i0=2,i1=1 -> 4 ; ...
        expected = [
            (1, 1, 0, 0, 0), (2, 2, 0, 0, 0),
            (3, 1, 1, 0, 0), (4, 2, 1, 0, 0),
            (5, 1, 2, 0, 0), (6, 2, 2, 0, 0),
        ]
        self.assertEqual(produced, expected)

    def test_atom_index_formula_matches_magnetizationinit(self) -> None:
        geometry = standard_geometry()
        # i = I0 + I1*NA + I2*N1*NA + I3*N2*N1*NA (I0 1-based, I1/I2/I3 0-based)
        i0, i1, i2, i3 = 2, 3, 1, 1
        expected = i0 + i1 * geometry.na + i2 * geometry.n1 * geometry.na + \
            i3 * geometry.n2 * geometry.n1 * geometry.na
        self.assertEqual(geometry.atom_index(i0, i1, i2, i3), expected)

    def test_geometry_rejects_mismatched_basis_length(self) -> None:
        with self.assertRaises(gen.MalformedGeneratorInputError):
            gen.Geometry(na=2, n1=6, n2=2, n3=2, basis=STANDARD_BASIS[:1])

    def test_geometry_rejects_nonpositive_extent(self) -> None:
        with self.assertRaises(gen.MalformedGeneratorInputError):
            gen.Geometry(na=2, n1=0, n2=2, n3=2, basis=STANDARD_BASIS)


class ConicalSpiralTests(unittest.TestCase):
    def test_every_direction_is_normalized(self) -> None:
        geometry = standard_geometry()
        state = gen.conical_spiral_state(geometry, cone_angle_deg=45.0, turns=1.0)
        for atom, vector in state.direction_by_atom.items():
            norm = math.sqrt(sum(c * c for c in vector))
            self.assertAlmostEqual(norm, 1.0, places=12, msg=f"atom {atom} not normalized")

    def test_momfile_covers_every_basis_site_once(self) -> None:
        geometry = standard_geometry()
        state = gen.conical_spiral_state(geometry, cone_angle_deg=45.0, turns=1.0)
        lines = state.momfile_text.strip("\n").split("\n")
        self.assertEqual(len(lines), geometry.na)
        sites = [int(line.split()[0]) for line in lines]
        self.assertEqual(sites, list(range(1, geometry.na + 1)))

    def test_axis_component_equals_cos_theta_everywhere(self) -> None:
        geometry = standard_geometry()
        state = gen.conical_spiral_state(geometry, cone_angle_deg=30.0, turns=2.0,
                                          axis=(0.0, 0.0, 1.0))
        expected = math.cos(math.radians(30.0))
        for vector in state.direction_by_atom.values():
            self.assertAlmostEqual(vector[2], expected, places=10)

    def test_periodic_closure_for_integer_turns(self) -> None:
        # An integer number of turns across the full supercell along the
        # modulation axis must bring the phase back to its start-of-cell
        # value at the next periodic image; check this indirectly by
        # confirming the in-plane azimuth advances by exactly
        # turns*2*pi/n1 per unit cell step (a linear, not wrapped-random,
        # progression), so n1 consecutive steps close on 2*pi*turns.
        geometry = standard_geometry(n1=6)
        state = gen.conical_spiral_state(geometry, cone_angle_deg=60.0, turns=1.0,
                                          modulation_cell_axis=0)
        azimuths = []
        for i1 in range(geometry.n1):
            atom = geometry.atom_index(1, i1, 0, 0)
            x, y, _ = state.direction_by_atom[atom]
            azimuths.append(math.atan2(y, x))
        # Unwrap and check linear advance of turns*2*pi/n1 per step.
        unwrapped = [azimuths[0]]
        for value in azimuths[1:]:
            previous = unwrapped[-1]
            delta = value - (previous % (2 * math.pi))
            while delta > math.pi:
                delta -= 2 * math.pi
            while delta < -math.pi:
                delta += 2 * math.pi
            unwrapped.append(previous + delta)
        step = 2 * math.pi * 1.0 / geometry.n1
        for index in range(1, len(unwrapped)):
            self.assertAlmostEqual(unwrapped[index] - unwrapped[index - 1], step, places=8)
        total_advance = unwrapped[-1] - unwrapped[0] + step  # + one more step wraps to image 0
        self.assertAlmostEqual(total_advance, 2 * math.pi * 1.0, places=8)

    def test_rejects_planar_cone_angle(self) -> None:
        geometry = standard_geometry()
        with self.assertRaises(gen.DegenerateStateError):
            gen.conical_spiral_state(geometry, cone_angle_deg=90.0, turns=1.0)

    def test_rejects_near_planar_cone_angle_within_tolerance(self) -> None:
        geometry = standard_geometry()
        with self.assertRaises(gen.DegenerateStateError):
            gen.conical_spiral_state(geometry, cone_angle_deg=87.0, turns=1.0,
                                      degeneracy_tolerance_deg=5.0)

    def test_rejects_collinear_cone_angle(self) -> None:
        geometry = standard_geometry()
        with self.assertRaises(gen.DegenerateStateError):
            gen.conical_spiral_state(geometry, cone_angle_deg=0.0, turns=1.0)
        with self.assertRaises(gen.DegenerateStateError):
            gen.conical_spiral_state(geometry, cone_angle_deg=180.0, turns=1.0)

    def test_rejects_zero_turns(self) -> None:
        geometry = standard_geometry()
        with self.assertRaises(gen.DegenerateStateError):
            gen.conical_spiral_state(geometry, cone_angle_deg=45.0, turns=0.0)

    def test_accepts_angle_just_outside_tolerance(self) -> None:
        geometry = standard_geometry()
        # Should not raise: 80 degrees is 10 degrees from the 90-degree
        # degenerate case, outside the default 5-degree tolerance.
        gen.conical_spiral_state(geometry, cone_angle_deg=80.0, turns=1.0)

    def test_torque_scale_hint_vanishes_only_at_rejected_angles(self) -> None:
        geometry = standard_geometry()
        state = gen.conical_spiral_state(geometry, cone_angle_deg=45.0, turns=1.0)
        self.assertGreater(state.initial_torque_scale_hint, 0.0)

    def test_rejects_malformed_modulation_axis(self) -> None:
        geometry = standard_geometry()
        with self.assertRaises(gen.MalformedGeneratorInputError):
            gen.conical_spiral_state(geometry, cone_angle_deg=45.0, turns=1.0,
                                      modulation_cell_axis=3)

    def test_rejects_nonpositive_moment_magnitude(self) -> None:
        geometry = standard_geometry()
        with self.assertRaises(gen.MalformedGeneratorInputError):
            gen.conical_spiral_state(geometry, cone_angle_deg=45.0, turns=1.0,
                                      moment_magnitude=0.0)

    def test_arbitrary_axis_matches_z_axis_case_by_symmetry(self) -> None:
        # Rotating the whole construction (axis and modulation direction)
        # by a fixed 3D rotation should leave the cone angle (axis-aligned
        # component) invariant; spot-check with axis=(0,0,1) vs a permuted
        # axis=(1,0,0) using modulation along the matching cell direction.
        geometry = standard_geometry()
        state_z = gen.conical_spiral_state(geometry, cone_angle_deg=40.0, turns=1.0,
                                            axis=(0.0, 0.0, 1.0), modulation_cell_axis=0)
        state_x = gen.conical_spiral_state(geometry, cone_angle_deg=40.0, turns=1.0,
                                            axis=(1.0, 0.0, 0.0), modulation_cell_axis=0)
        for atom in state_z.direction_by_atom:
            axis_component_z = state_z.direction_by_atom[atom][2]
            axis_component_x = state_x.direction_by_atom[atom][0]
            self.assertAlmostEqual(axis_component_z, axis_component_x, places=10)


class ChiralPartnerPairTests(unittest.TestCase):
    def test_partners_share_identical_cone_profile(self) -> None:
        geometry = standard_geometry()
        plus, minus = gen.chiral_partner_pair(geometry, cone_angle_deg=35.0, turns=1.0,
                                               axis=(0.0, 0.0, 1.0))
        for atom in plus.direction_by_atom:
            self.assertAlmostEqual(
                plus.direction_by_atom[atom][2], minus.direction_by_atom[atom][2], places=10
            )

    def test_partners_wind_in_opposite_directions(self) -> None:
        geometry = standard_geometry()
        plus, minus = gen.chiral_partner_pair(geometry, cone_angle_deg=35.0, turns=1.0,
                                               axis=(0.0, 0.0, 1.0), modulation_cell_axis=0)
        atom_a = geometry.atom_index(1, 0, 0, 0)
        atom_b = geometry.atom_index(1, 1, 0, 0)

        def azimuth_delta(state: gen.MomfileSpiralState) -> float:
            xa, ya, _ = state.direction_by_atom[atom_a]
            xb, yb, _ = state.direction_by_atom[atom_b]
            delta = math.atan2(yb, xb) - math.atan2(ya, xa)
            while delta > math.pi:
                delta -= 2 * math.pi
            while delta < -math.pi:
                delta += 2 * math.pi
            return delta

        self.assertGreater(azimuth_delta(plus), 0.0)
        self.assertLess(azimuth_delta(minus), 0.0)

    def test_partners_are_tagged_in_manifest(self) -> None:
        geometry = standard_geometry()
        plus, minus = gen.chiral_partner_pair(geometry, cone_angle_deg=35.0, turns=1.0)
        self.assertEqual(plus.manifest["extra"]["chirality"], "+q")
        self.assertEqual(minus.manifest["extra"]["chirality"], "-q")

    def test_does_not_depend_on_any_dmi_file(self) -> None:
        # Structural guard: neither generator function accepts a DMI-file
        # argument at all, so the +-q choice cannot be derived from one.
        import inspect
        for name in ("conical_spiral_state", "chiral_partner_pair"):
            parameters = inspect.signature(getattr(gen, name)).parameters
            self.assertNotIn("dmfile", parameters)
            self.assertNotIn("dmi_file", parameters)


class DomainWallPairTests(unittest.TestCase):
    def _geometry(self) -> gen.Geometry:
        return standard_geometry(n1=24, n2=2, n3=2)

    def test_every_direction_is_normalized(self) -> None:
        geometry = self._geometry()
        state = gen.domain_wall_pair_state(
            geometry, axis_cell_index=0, wall_centers_cells=(6.0, 18.0), width_cells=1.0,
        )
        for atom, vector in state.direction_by_atom.items():
            norm = math.sqrt(sum(c * c for c in vector))
            self.assertAlmostEqual(norm, 1.0, places=10, msg=f"atom {atom} not normalized")

    def test_restart_text_has_exactly_natom_rows_after_header(self) -> None:
        geometry = self._geometry()
        state = gen.domain_wall_pair_state(
            geometry, axis_cell_index=0, wall_centers_cells=(6.0, 18.0), width_cells=1.0,
        )
        lines = state.restart_text.strip("\n").split("\n")
        header, rows = lines[:7], lines[7:]
        self.assertTrue(header[0].startswith("#" * 10))
        self.assertEqual(len(rows), geometry.natom)

    def test_exactly_two_wall_crossings_along_periodic_axis(self) -> None:
        geometry = self._geometry()
        state = gen.domain_wall_pair_state(
            geometry, axis_cell_index=0, wall_centers_cells=(6.0, 18.0), width_cells=1.0,
        )
        signs = []
        for i1 in range(geometry.n1):
            atom = geometry.atom_index(1, i1, 0, 0)
            signs.append(1 if state.direction_by_atom[atom][2] >= 0.0 else -1)
        crossings = sum(1 for i in range(len(signs)) if signs[i] != signs[i - 1])
        self.assertEqual(crossings, 2)

    def test_wall_crossings_are_near_requested_centres(self) -> None:
        geometry = self._geometry()
        centres = (6.0, 18.0)
        state = gen.domain_wall_pair_state(
            geometry, axis_cell_index=0, wall_centers_cells=centres, width_cells=1.0,
        )
        crossing_positions = []
        previous_sign = None
        for i1 in range(geometry.n1 + 1):
            atom = geometry.atom_index(1, i1 % geometry.n1, 0, 0)
            sign = 1 if state.direction_by_atom[atom][2] >= 0.0 else -1
            if previous_sign is not None and sign != previous_sign:
                crossing_positions.append(i1)
            previous_sign = sign
        self.assertEqual(len(crossing_positions), 2)
        for observed, expected in zip(sorted(crossing_positions), sorted(centres)):
            self.assertLess(abs(observed - expected), 1.5)

    def test_periodic_closure_up_state_matches_across_boundary(self) -> None:
        # The two-kink profile closes only up to an exponentially small
        # residual in exp(-separation/width) (see _wall_polar_angle's
        # docstring); with centres (6, 18) on a 24-cell axis and width 1,
        # the atom nearest each boundary is ~5-6 widths from its nearest
        # wall, so a residual of order exp(-5) ~= 0.007 (in radians of
        # polar angle, i.e. a comparable deviation in the z-component near
        # +1) is expected and is not evidence of a construction error.
        geometry = self._geometry()
        state = gen.domain_wall_pair_state(
            geometry, axis_cell_index=0, wall_centers_cells=(6.0, 18.0), width_cells=1.0,
        )
        near_start = state.direction_by_atom[geometry.atom_index(1, 0, 0, 0)]
        near_end = state.direction_by_atom[geometry.atom_index(1, geometry.n1 - 1, 0, 0)]
        for component in range(3):
            self.assertLess(abs(near_start[component] - near_end[component]), 0.05)

    def test_neel_and_bloch_walls_rotate_in_different_planes(self) -> None:
        geometry = self._geometry()
        neel = gen.domain_wall_pair_state(
            geometry, axis_cell_index=0, wall_centers_cells=(6.0, 18.0), width_cells=1.0,
            wall_type="NEEL", cant_deg=0.0,
        )
        bloch = gen.domain_wall_pair_state(
            geometry, axis_cell_index=0, wall_centers_cells=(6.0, 18.0), width_cells=1.0,
            wall_type="BLOCH", cant_deg=0.0,
        )
        # At the first wall centre, a Neel wall's core moment should have a
        # nonzero x-component (rotation in the x-z plane) and negligible
        # y-component; Bloch is the opposite.
        core_atom = geometry.atom_index(1, 6, 0, 0)
        nx, ny, _ = neel.direction_by_atom[core_atom]
        bx, by, _ = bloch.direction_by_atom[core_atom]
        self.assertGreater(abs(nx), abs(ny))
        self.assertGreater(abs(by), abs(bx))

    def test_chirality_flips_winding_sense(self) -> None:
        geometry = self._geometry()
        plus = gen.domain_wall_pair_state(
            geometry, axis_cell_index=0, wall_centers_cells=(6.0, 18.0), width_cells=1.0,
            chirality=1, cant_deg=0.0,
        )
        minus = gen.domain_wall_pair_state(
            geometry, axis_cell_index=0, wall_centers_cells=(6.0, 18.0), width_cells=1.0,
            chirality=-1, cant_deg=0.0,
        )
        core_atom = geometry.atom_index(1, 6, 0, 0)
        vp = plus.direction_by_atom[core_atom]
        vm = minus.direction_by_atom[core_atom]
        self.assertAlmostEqual(vp[0], -vm[0], places=8)
        self.assertAlmostEqual(vp[1], -vm[1], places=8)
        self.assertAlmostEqual(vp[2], vm[2], places=8)

    def test_cant_breaks_exact_symmetry_at_wall_core(self) -> None:
        geometry = self._geometry()
        uncanted = gen.domain_wall_pair_state(
            geometry, axis_cell_index=0, wall_centers_cells=(6.0, 18.0), width_cells=1.0,
            wall_type="NEEL", cant_deg=0.0,
        )
        canted = gen.domain_wall_pair_state(
            geometry, axis_cell_index=0, wall_centers_cells=(6.0, 18.0), width_cells=1.0,
            wall_type="NEEL", cant_deg=2.0,
        )
        core_atom = geometry.atom_index(1, 6, 0, 0)
        uncanted_y = uncanted.direction_by_atom[core_atom][1]
        canted_y = canted.direction_by_atom[core_atom][1]
        self.assertAlmostEqual(uncanted_y, 0.0, places=8)
        self.assertNotAlmostEqual(canted_y, 0.0, places=4)

    def test_rejects_walls_too_close_together(self) -> None:
        geometry = self._geometry()
        with self.assertRaises(gen.DegenerateStateError):
            gen.domain_wall_pair_state(
                geometry, axis_cell_index=0, wall_centers_cells=(11.5, 12.5), width_cells=1.0,
            )

    def test_rejects_walls_too_close_to_periodic_boundary(self) -> None:
        geometry = self._geometry()
        with self.assertRaises(gen.DegenerateStateError):
            gen.domain_wall_pair_state(
                geometry, axis_cell_index=0, wall_centers_cells=(0.5, 12.0), width_cells=1.0,
            )

    def test_rejects_unordered_or_out_of_range_centres(self) -> None:
        geometry = self._geometry()
        with self.assertRaises(gen.MalformedGeneratorInputError):
            gen.domain_wall_pair_state(
                geometry, axis_cell_index=0, wall_centers_cells=(18.0, 6.0), width_cells=1.0,
            )
        with self.assertRaises(gen.MalformedGeneratorInputError):
            gen.domain_wall_pair_state(
                geometry, axis_cell_index=0, wall_centers_cells=(6.0, 999.0), width_cells=1.0,
            )

    def test_rejects_nonpositive_width(self) -> None:
        geometry = self._geometry()
        with self.assertRaises(gen.MalformedGeneratorInputError):
            gen.domain_wall_pair_state(
                geometry, axis_cell_index=0, wall_centers_cells=(6.0, 18.0), width_cells=0.0,
            )

    def test_rejects_invalid_wall_type_and_chirality(self) -> None:
        geometry = self._geometry()
        with self.assertRaises(gen.MalformedGeneratorInputError):
            gen.domain_wall_pair_state(
                geometry, axis_cell_index=0, wall_centers_cells=(6.0, 18.0), width_cells=1.0,
                wall_type="ICE",
            )
        with self.assertRaises(gen.MalformedGeneratorInputError):
            gen.domain_wall_pair_state(
                geometry, axis_cell_index=0, wall_centers_cells=(6.0, 18.0), width_cells=1.0,
                chirality=0,
            )


class DeterministicRegenerationTests(unittest.TestCase):
    def test_conical_spiral_is_byte_stable_across_calls(self) -> None:
        geometry = standard_geometry()
        first = gen.conical_spiral_state(geometry, cone_angle_deg=45.0, turns=1.0)
        second = gen.conical_spiral_state(geometry, cone_angle_deg=45.0, turns=1.0)
        self.assertEqual(first.momfile_text, second.momfile_text)
        self.assertEqual(first.manifest["content_sha256"], second.manifest["content_sha256"])
        self.assertEqual(first.direction_by_atom, second.direction_by_atom)

    def test_domain_wall_pair_is_byte_stable_across_calls(self) -> None:
        geometry = standard_geometry(n1=24)
        first = gen.domain_wall_pair_state(
            geometry, axis_cell_index=0, wall_centers_cells=(6.0, 18.0), width_cells=1.0,
        )
        second = gen.domain_wall_pair_state(
            geometry, axis_cell_index=0, wall_centers_cells=(6.0, 18.0), width_cells=1.0,
        )
        self.assertEqual(first.restart_text, second.restart_text)
        self.assertEqual(first.manifest["content_sha256"], second.manifest["content_sha256"])

    def test_changing_a_parameter_changes_the_hash(self) -> None:
        geometry = standard_geometry()
        base = gen.conical_spiral_state(geometry, cone_angle_deg=45.0, turns=1.0)
        changed = gen.conical_spiral_state(geometry, cone_angle_deg=46.0, turns=1.0)
        self.assertNotEqual(base.manifest["content_sha256"], changed.manifest["content_sha256"])

    def test_manifest_json_is_deterministic_and_parses(self) -> None:
        import json
        geometry = standard_geometry()
        state = gen.conical_spiral_state(geometry, cone_angle_deg=45.0, turns=1.0)
        first = gen.manifest_json(state.manifest)
        second = gen.manifest_json(state.manifest)
        self.assertEqual(first, second)
        parsed = json.loads(first)
        self.assertEqual(parsed["generator_name"], "conical_spiral_state")
        self.assertEqual(parsed["generator_version"], gen.GENERATOR_VERSION)
        self.assertIn("content_sha256", parsed)

    def test_content_hash_is_a_pure_function_of_text(self) -> None:
        self.assertEqual(gen.content_hash("abc"), gen.content_hash("abc"))
        self.assertNotEqual(gen.content_hash("abc"), gen.content_hash("abd"))


if __name__ == "__main__":
    unittest.main()
