#!/usr/bin/env python3
"""Unit tests for static_topology_oracle.py (RCG-04F)."""

from __future__ import annotations

import math
import unittest
from pathlib import Path

import static_topology_oracle as sto
import torque_oracle as orc
from moving_state_generator import Geometry, manifest_json

E2E_ROOT = Path(__file__).with_name("e2e")

GEOMETRY = Geometry(na=2, n1=24, n2=2, n3=2, basis=((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)))
JFILE_TEXT = (
    "1 1 0.5 0.5 0.5 1.33767484769984\n"
    "1 1 1.0 0.0 0.0 0.75703576545650\n"
    "1 1 1.0 1.0 0.0 -0.05975437643846\n"
    "1 1 1.5 0.5 0.5 -0.08819834160658\n"
)
SHELLS = orc.parse_jfile_shells(JFILE_TEXT)


class BufferWidthTests(unittest.TestCase):
    def test_max_radius_is_largest_shell(self):
        self.assertAlmostEqual(sto.max_interaction_radius_m(SHELLS), 2.75 ** 0.5, places=12)

    def test_buffer_width_block_size_1(self):
        self.assertEqual(sto.buffer_width_blocks_x(SHELLS, 1), 2)

    def test_buffer_width_block_size_2(self):
        self.assertEqual(sto.buffer_width_blocks_x(SHELLS, 2), 1)

    def test_buffer_width_adds_cg_buffer_blocks(self):
        self.assertEqual(sto.buffer_width_blocks_x(SHELLS, 1, cg_buffer_blocks=3), 5)

    def test_rejects_nonpositive_block_size(self):
        with self.assertRaises(ValueError):
            sto.buffer_width_blocks_x(SHELLS, 0)

    def test_rejects_negative_buffer_blocks(self):
        with self.assertRaises(ValueError):
            sto.buffer_width_blocks_x(SHELLS, 1, cg_buffer_blocks=-1)


class ExpectedTopologyTests(unittest.TestCase):
    def test_bs1_partition_matches_hand_derivation(self):
        topology = sto.compute_expected_topology(
            GEOMETRY, SHELLS, block_size_x=1, block_size_y=2, block_size_z=2,
            fine_block_ids=set(range(1, 7)),
        )
        self.assertEqual(topology.block_counts(), {"fine": 6, "interface": 4, "coarse": 14})
        self.assertEqual(topology.atom_counts(), {"fine": 48, "interface": 32, "coarse": 112})
        self.assertEqual(topology.block_ids_in_state(sto.FINE), [1, 2, 3, 4, 5, 6])
        self.assertEqual(topology.block_ids_in_state(sto.BUFFER), [7, 8, 23, 24])

    def test_bs2_partition_matches_bs1_atom_counts(self):
        """bs2 is constructed to cover the exact same physical atom split as bs1
        (RCG-04F's spatial-refinement pair) -- checked, not merely asserted by
        construction, so a future edit to either fixture's FINE ids would be
        caught here before it silently broke the refinement comparison."""
        topology = sto.compute_expected_topology(
            GEOMETRY, SHELLS, block_size_x=2, block_size_y=2, block_size_z=2,
            fine_block_ids=set(range(1, 4)),
        )
        self.assertEqual(topology.block_counts(), {"fine": 3, "interface": 2, "coarse": 7})
        self.assertEqual(topology.atom_counts(), {"fine": 48, "interface": 32, "coarse": 112})

    def test_shifted_partition_preserves_counts_different_ids(self):
        topology = sto.compute_expected_topology(
            GEOMETRY, SHELLS, block_size_x=1, block_size_y=2, block_size_z=2,
            fine_block_ids=set(range(13, 19)),
        )
        self.assertEqual(topology.block_counts(), {"fine": 6, "interface": 4, "coarse": 14})
        self.assertEqual(topology.block_ids_in_state(sto.FINE), [13, 14, 15, 16, 17, 18])
        self.assertEqual(topology.block_ids_in_state(sto.BUFFER), [11, 12, 19, 20])

    def test_rejects_nondivisible_block_size(self):
        with self.assertRaises(ValueError):
            sto.compute_expected_topology(
                GEOMETRY, SHELLS, block_size_x=5, block_size_y=2, block_size_z=2,
                fine_block_ids={1},
            )

    def test_genuine_3d_block_grid_is_now_supported(self):
        """RCG-05B: a block grid that does not fully span y/z (``block_size_y=1``
        on this GEOMETRY's ``n2=2``, giving a genuine ``grid_y=2``) used to raise
        ``UnsupportedTopologyError``; it now computes a real 2D partition.
        Every number below is independently hand-derived (not merely re-printed
        from the module): ``block_vectors=diag(1,1,2)``, so
        ``inverse=diag(1,1,0.5)``, giving ``width=(ceil(1.6583),ceil(1.6583),
        ceil(0.8292))=(2,2,1)``; block 2 (grid coord (1,0,0)) has periodic delta
        (1,0,0) <= (2,2,1) -> BUFFER; block 3 (coord (2,0,0)) has periodic delta
        (2,0,0) <= (2,2,1) -> still BUFFER (width_x=2, not 1); block 25 (coord
        (0,1,0), the seed's periodic y-neighbour) has periodic delta (0,1,0) <=
        (2,2,1) -> BUFFER."""
        topology = sto.compute_expected_topology(
            GEOMETRY, SHELLS, block_size_x=1, block_size_y=1, block_size_z=2,
            fine_block_ids={1},
        )
        self.assertEqual(topology.block_grid, (24, 2, 1))
        self.assertEqual(topology.buffer_width, (2, 2, 1))
        self.assertEqual(topology.block_counts(), {"fine": 1, "interface": 9, "coarse": 38})
        self.assertEqual(topology.block_state_by_id[1], sto.FINE)
        self.assertEqual(topology.block_state_by_id[2], sto.BUFFER)
        self.assertEqual(topology.block_state_by_id[3], sto.BUFFER)
        self.assertEqual(topology.block_state_by_id[25], sto.BUFFER)

    def test_rejects_nondivisible_block_size_y(self):
        """The divisibility check now applies on every axis, not only x."""
        with self.assertRaises(ValueError):
            sto.compute_expected_topology(
                GEOMETRY, SHELLS, block_size_x=1, block_size_y=3, block_size_z=2,
                fine_block_ids={1},
            )

    def test_rejects_empty_fine_ids(self):
        with self.assertRaises(ValueError):
            sto.compute_expected_topology(
                GEOMETRY, SHELLS, block_size_x=1, block_size_y=2, block_size_z=2,
                fine_block_ids=set(),
            )

    def test_rejects_out_of_range_fine_id(self):
        with self.assertRaises(ValueError):
            sto.compute_expected_topology(
                GEOMETRY, SHELLS, block_size_x=1, block_size_y=2, block_size_z=2,
                fine_block_ids={25},
            )

    def test_atom_block_mapping_matches_floor_division(self):
        topology = sto.compute_expected_topology(
            GEOMETRY, SHELLS, block_size_x=2, block_size_y=2, block_size_z=2,
            fine_block_ids={1},
        )
        # atom 1 (basis 1, cell (0,0,0)) -> block 1; the next basis-1 atom
        # along x (cell (2,0,0), the first cell of block 2) -> block 2.
        atom_cell0 = GEOMETRY.atom_index(1, 0, 0, 0)
        atom_cell2 = GEOMETRY.atom_index(1, 2, 0, 0)
        self.assertEqual(topology.atom_block_by_atom[atom_cell0], 1)
        self.assertEqual(topology.atom_block_by_atom[atom_cell2], 2)

    def test_distance_from_boundary_zero_on_boundary_blocks(self):
        topology = sto.compute_expected_topology(
            GEOMETRY, SHELLS, block_size_x=1, block_size_y=2, block_size_z=2,
            fine_block_ids=set(range(1, 7)),
        )
        # Block 6 is FINE and adjacent to block 7 (BUFFER): distance 1, not 0.
        # Block 3/4 are interior FINE blocks, farthest from any BUFFER/COARSE block.
        self.assertEqual(topology.distance_from_boundary_blocks[6], 1)
        self.assertGreaterEqual(topology.distance_from_boundary_blocks[3], 2)
        # Block 9 is the first COARSE block after the BUFFER window {7,8}: distance 1.
        self.assertEqual(topology.distance_from_boundary_blocks[9], 1)

    def test_resolution_state_values_ordered_by_block_id(self):
        topology = sto.compute_expected_topology(
            GEOMETRY, SHELLS, block_size_x=2, block_size_y=2, block_size_z=2,
            fine_block_ids={1},
        )
        values = topology.resolution_state_values()
        self.assertEqual(len(values), topology.nblocks_x)
        self.assertEqual(values[0], sto.FINE)


class MaskFileTests(unittest.TestCase):
    def test_write_mask_file_round_trips_fine_ids(self):
        topology = sto.compute_expected_topology(
            GEOMETRY, SHELLS, block_size_x=1, block_size_y=2, block_size_z=2,
            fine_block_ids={3, 1, 2},
        )
        text = sto.write_mask_file(topology)
        lines = [line for line in text.splitlines() if not line.startswith("#")]
        self.assertEqual(sorted(lines), ["1 FINE", "2 FINE", "3 FINE"])


class InterfaceBondCountTests(unittest.TestCase):
    def test_bs1_and_bs2_report_the_same_physical_bond_count(self):
        """Same physical fine/interface/coarse atom partition -> same boundary
        bonds, independent of how finely that partition is discretised into
        blocks (a regression check on the topology construction itself)."""
        t1 = sto.compute_expected_topology(
            GEOMETRY, SHELLS, block_size_x=1, block_size_y=2, block_size_z=2,
            fine_block_ids=set(range(1, 7)),
        )
        t2 = sto.compute_expected_topology(
            GEOMETRY, SHELLS, block_size_x=2, block_size_y=2, block_size_z=2,
            fine_block_ids=set(range(1, 4)),
        )
        count1 = sto.interface_bond_count(GEOMETRY, SHELLS, t1)
        count2 = sto.interface_bond_count(GEOMETRY, SHELLS, t2)
        self.assertGreater(count1, 0)
        self.assertEqual(count1, count2)

    def test_all_fine_topology_has_no_interface_bonds(self):
        topology = sto.compute_expected_topology(
            GEOMETRY, SHELLS, block_size_x=1, block_size_y=2, block_size_z=2,
            fine_block_ids=set(range(1, 25)),
        )
        self.assertEqual(sto.interface_bond_count(GEOMETRY, SHELLS, topology), 0)


# RCG-05B: general (anisotropic/skew) buffer-width metric, genuine 3D block
# grid, skew periodic-image handling, and the new block-geometry generators.

_RCG05A_SHELLS = (orc.ExchangeShell(dx=4.0, dy=0.0, dz=0.0, jij=0.01),)
_SKEW_CELL = ((1.0, 0.0, 0.0), (0.5, 1.0, 0.0), (0.0, 0.0, 1.0))
_SKEW_SHELLS = (orc.ExchangeShell(dx=0.5, dy=1.0, dz=0.0, jij=0.01),)  # exactly C2


class GeneralBufferWidthTests(unittest.TestCase):
    def test_orthogonal_anisotropic_block_shape_matches_rcg05a_hardware_run(self):
        """Same parameters RCG-05A reproduced through the real CPU and CUDA
        executables (``docs/RCG-05_GEOMETRY_OWNERSHIP_EVIDENCE.md`` SS2.3):
        ``block_size=(1,2,4)``, radius 4.0 -> CPU printed
        ``cpu_buffer_width_blocks=4 2 1``."""
        self.assertEqual(
            sto.buffer_width_blocks(_RCG05A_SHELLS, sto._IDENTITY_CELL, (1, 2, 4)), (4, 2, 1),
        )

    def test_skew_cell_widths_hand_derived(self):
        """cell=((1,0,0),(0.5,1,0),(0,0,1)), radius=|C2|=sqrt(1.25)=1.118034.
        block_size=(1,2,2): block_vectors columns (1,0,0)/(1,2,0)/(0,0,2),
        determinant 4; inverse rows (1,0,0)/(-0.5,0.5,0)/(0,0,0.5); row norms
        1, 0.7071, 0.5 -> widths ceil(1.118),ceil(0.7906),ceil(0.559) = 2,1,1.
        block_size=(2,2,2): scales every row norm by 1/2 -> ceil(0.559),
        ceil(0.3953),ceil(0.2795) = 1,1,1. Independently computed by a
        standalone script before this module's generalization was written,
        not re-derived from the module under test."""
        self.assertEqual(
            sto.buffer_width_blocks(_SKEW_SHELLS, _SKEW_CELL, (1, 2, 2)), (2, 1, 1),
        )
        self.assertEqual(
            sto.buffer_width_blocks(_SKEW_SHELLS, _SKEW_CELL, (2, 2, 2)), (1, 1, 1),
        )

    def test_backward_compatible_with_buffer_width_blocks_x(self):
        """``buffer_width_blocks_x`` must equal axis-0 of the general formula on
        an identity cell with ``block_shape=(bs,1,1)``, for every shells/block
        size this file's original ``BufferWidthTests`` already covers."""
        for block_size_x, cg in ((1, 0), (2, 0), (1, 3)):
            with self.subTest(block_size_x=block_size_x, cg=cg):
                self.assertEqual(
                    sto.buffer_width_blocks_x(SHELLS, block_size_x, cg_buffer_blocks=cg),
                    sto.buffer_width_blocks(SHELLS, sto._IDENTITY_CELL, (block_size_x, 1, 1), cg)[0],
                )

    def test_rejects_nonpositive_block_shape(self):
        with self.assertRaises(ValueError):
            sto.buffer_width_blocks(SHELLS, sto._IDENTITY_CELL, (1, 0, 1))

    def test_rejects_singular_cell(self):
        singular = ((1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 0.0, 1.0))  # C1, C2 collinear
        with self.assertRaises(sto.DegenerateGeometryError):
            sto.buffer_width_blocks(SHELLS, singular, (1, 1, 1))


class Full3DOwnershipTests(unittest.TestCase):
    """Genuine 3D grid (grid_x, grid_y, grid_z all > 1), hand-derived, not merely
    printed from the module under test: na=1, n1=n2=n3=5, identity cell, a
    single distance-1.0 shell, block_size=(1,1,1) -> width=(1,1,1). A block is
    atomistic iff every axis's periodic delta from the seed is <= 1, i.e. its
    grid coordinate lies in ``{0,1,4}`` on every axis: a 3x3x3=27-block box
    (1 FINE + 26 BUFFER) out of 5^3=125, leaving 98 COARSE."""

    @classmethod
    def setUpClass(cls):
        cls.geometry = Geometry(na=1, n1=5, n2=5, n3=5, basis=((0.0, 0.0, 0.0),))
        cls.shells = (orc.ExchangeShell(dx=1.0, dy=0.0, dz=0.0, jij=1.0),)
        cls.topology = sto.compute_expected_topology(
            cls.geometry, cls.shells, block_size_x=1, block_size_y=1, block_size_z=1,
            fine_block_ids={1},
        )

    def test_width_and_counts(self):
        self.assertEqual(self.topology.buffer_width, (1, 1, 1))
        self.assertEqual(self.topology.block_grid, (5, 5, 5))
        self.assertEqual(self.topology.block_counts(), {"fine": 1, "interface": 26, "coarse": 98})

    def test_specific_block_states_by_hand_derived_coordinate(self):
        cases = {
            (0, 0, 0): sto.FINE,      # the seed itself
            (1, 0, 0): sto.BUFFER,    # periodic delta (1,0,0), within the box
            (4, 0, 0): sto.BUFFER,    # periodic wrap: delta 4 -> min(4,1)=1
            (2, 0, 0): sto.COARSE,    # periodic delta 2 > width 1 on axis 0
            (1, 1, 1): sto.BUFFER,    # within the box on every axis
            (2, 2, 2): sto.COARSE,    # outside the box on every axis
        }
        for coord, expected_state in cases.items():
            block_id = sto._block_id(coord, self.topology.block_grid)
            with self.subTest(coord=coord):
                self.assertEqual(self.topology.block_state_by_id[block_id], expected_state)

    def test_atom_block_mapping_uses_all_three_axes(self):
        """With block_size=(1,1,1) every cell is its own block, so atom-to-block
        must equal the coordinate-to-block-id formula directly on every axis,
        not only axis 0 (the previous 1D-only formula's coverage)."""
        for i1, i2, i3 in ((0, 0, 0), (1, 2, 3), (4, 4, 4)):
            atom = self.geometry.atom_index(1, i1, i2, i3)
            expected_block = sto._block_id((i1, i2, i3), self.topology.block_grid)
            self.assertEqual(self.topology.atom_block_by_atom[atom], expected_block)

    def test_distance_from_boundary_is_chebyshev_over_all_three_axes(self):
        """Block (1,1,1) is BUFFER, immediately adjacent (Chebyshev distance 1,
        in every one of its 26 neighbours) to at least one COARSE block, e.g.
        (2,1,1) (periodic delta (1,0,0) from (1,1,1) -> coarse since (2,1,1)'s
        own periodic delta from the seed is (2,1,1), failing axis 0)."""
        block_111 = sto._block_id((1, 1, 1), self.topology.block_grid)
        self.assertEqual(self.topology.distance_from_boundary_blocks[block_111], 1)


class SkewCellOwnershipTests(unittest.TestCase):
    """A genuinely non-orthogonal cell, at the same scale as RCG-05B's own
    ``skew_cell_fixture`` default, with counts independently computed by a
    standalone script before this test was written (not merely re-derived
    from the module under test -- see the RCG-05B evidence document)."""

    def test_skew_partition_counts_match_independent_computation(self):
        geometry = Geometry(na=2, n1=8, n2=8, n3=8,
                             basis=((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)), cell_vectors=_SKEW_CELL)
        topology = sto.compute_expected_topology(
            geometry, _SKEW_SHELLS, block_size_x=1, block_size_y=2, block_size_z=2,
            fine_block_ids={1}, cell_vectors=_SKEW_CELL,
        )
        self.assertEqual(topology.block_grid, (8, 4, 4))
        self.assertEqual(topology.buffer_width, (2, 1, 1))
        self.assertEqual(topology.block_counts(), {"fine": 1, "interface": 44, "coarse": 83})


class LatticeMinimumImageTests(unittest.TestCase):
    """Hand-derived skew periodic-image reduction: displacement = (integer
    lattice translation) + (small residual) must reduce to exactly the
    residual, for a non-orthogonal cell where the old per-axis orthogonal
    wrap (``torque_oracle._minimum_image``) would give the wrong answer."""

    def test_reduces_exactly_c2_plus_residual(self):
        displacement = (0.6, 1.05, 0.0)  # = 1*C2 + (0.1, 0.05, 0.0)
        reduced = sto._lattice_minimum_image(displacement, _SKEW_CELL)
        for got, want in zip(reduced, (0.1, 0.05, 0.0)):
            self.assertAlmostEqual(got, want, places=10)

    def test_reduces_exactly_2c1_minus_c3_plus_residual(self):
        displacement = (2.2, -0.1, -0.95)  # = 2*C1 - 1*C3 + (0.2, -0.1, 0.05)
        reduced = sto._lattice_minimum_image(displacement, _SKEW_CELL)
        for got, want in zip(reduced, (0.2, -0.1, 0.05)):
            self.assertAlmostEqual(got, want, places=10)

    def test_identity_cell_matches_plain_per_axis_wrap_for_small_displacements(self):
        displacement = (0.3, -0.2, 0.1)
        reduced = sto._lattice_minimum_image(displacement, sto._IDENTITY_CELL)
        for got, want in zip(reduced, displacement):
            self.assertAlmostEqual(got, want, places=10)

    def test_rejects_singular_cell(self):
        singular = ((1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 0.0, 1.0))
        with self.assertRaises(sto.DegenerateGeometryError):
            sto._lattice_minimum_image((0.1, 0.1, 0.1), singular)


class GeneralBondBuilderRegressionTests(unittest.TestCase):
    """``build_geometric_bonds_general``, on an identity cell, must reproduce
    ``torque_oracle.build_geometric_bonds`` exactly on the same
    production-calibrated GEOMETRY/SHELLS this file already uses (see
    ``torque_oracle.py``'s module docstring for that calibration) -- a
    stronger regression than a small synthetic case, since GEOMETRY/SHELLS
    here were validated against real production exchange-energy output."""

    def test_matches_orthogonal_build_geometric_bonds_exactly(self):
        general = sto.build_geometric_bonds_general(GEOMETRY, SHELLS, sto._IDENTITY_CELL)
        orthogonal = orc.build_geometric_bonds(GEOMETRY, SHELLS)
        self.assertEqual(general, orthogonal)

    def test_interface_bond_count_general_matches_interface_bond_count(self):
        topology = sto.compute_expected_topology(
            GEOMETRY, SHELLS, block_size_x=1, block_size_y=2, block_size_z=2,
            fine_block_ids=set(range(1, 7)),
        )
        self.assertEqual(
            sto.interface_bond_count_general(GEOMETRY, SHELLS, topology, sto._IDENTITY_CELL),
            sto.interface_bond_count(GEOMETRY, SHELLS, topology),
        )


def _reference_1d_topology(n1: int, block_size_x: int, shells, fine_block_ids: set[int]):
    """Independent re-implementation of RCG-04F's *original* (pre-RCG-05B)
    1D-only algorithm, transcribed directly from the removed code (not
    calling ``static_topology_oracle`` at all), used only as a regression
    oracle-of-the-oracle against the new general 3D implementation."""
    radius = max(shell.distance for shell in shells)
    fractional = radius / block_size_x
    width = math.ceil(max(0.0, fractional - sto._EPS_GUARD))
    nblocks_x = n1 // block_size_x
    fine_mask = [False] * nblocks_x
    for bid in fine_block_ids:
        fine_mask[bid - 1] = True
    atomistic_block = list(fine_mask)
    for seed_bx in range(nblocks_x):
        if not fine_mask[seed_bx]:
            continue
        for bx in range(nblocks_x):
            delta = abs(bx - seed_bx)
            periodic_delta = min(delta, nblocks_x - delta)
            if periodic_delta <= width:
                atomistic_block[bx] = True
    state_by_id = {}
    for bx in range(nblocks_x):
        bid = bx + 1
        if fine_mask[bx]:
            state_by_id[bid] = sto.FINE
        elif atomistic_block[bx]:
            state_by_id[bid] = sto.BUFFER
        else:
            state_by_id[bid] = sto.COARSE
    return nblocks_x, width, state_by_id


class RealFixtureGeometryRegressionTests(unittest.TestCase):
    """Regression: the generalized 3D ``compute_expected_topology``, called
    with every ``(ncell_x, block_size_x)`` combination an actual tracked
    RCG-04 ``tests/coarse_graining/e2e/*/inpsd.dat`` fixture uses (every one
    of which has ``ncell _ 2 2``/``block_size_y=block_size_z=2``, confirmed
    by grep across the tracked tree before writing this test), must agree
    exactly with ``_reference_1d_topology``'s independent transcription of
    the original algorithm this slice replaced. Uses ``SHELLS``, parsed from
    this file's ``JFILE_TEXT``, which is confirmed byte-identical to the
    tracked ``e2e/jfile`` every one of these fixtures actually reads."""

    # (ncell_x, block_size_x), grep'd from every e2e/*/inpsd.dat with a
    # block_size_x record and a divisible ncell_x (excluding
    # invalid_partial_block's deliberately-nondivisible ncell 5/block_size 2
    # negative control, which is covered by test_rejects_nondivisible_block_size).
    REAL_COMBOS = ((10, 1), (24, 1), (24, 2), (24, 4), (24, 8), (6, 1), (6, 2))

    def test_matches_reference_1d_algorithm_on_every_real_fixture_geometry(self):
        with open(E2E_ROOT / "jfile", encoding="utf-8") as handle:
            tracked_jfile_text = handle.read()
        self.assertEqual(tracked_jfile_text, JFILE_TEXT)

        for ncell_x, block_size_x in self.REAL_COMBOS:
            geometry = Geometry(na=2, n1=ncell_x, n2=2, n3=2,
                                 basis=((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)))
            fine_block_ids = {1}
            expected_nblocks, expected_width, expected_state = _reference_1d_topology(
                ncell_x, block_size_x, SHELLS, fine_block_ids,
            )
            topology = sto.compute_expected_topology(
                geometry, SHELLS, block_size_x=block_size_x, block_size_y=2, block_size_z=2,
                fine_block_ids=fine_block_ids,
            )
            with self.subTest(ncell_x=ncell_x, block_size_x=block_size_x):
                self.assertEqual(topology.block_grid, (expected_nblocks, 1, 1))
                self.assertEqual(topology.nblocks_x, expected_nblocks)
                self.assertEqual(topology.buffer_width_x, expected_width)
                self.assertEqual(topology.block_state_by_id, expected_state)


class BlockGeometryGeneratorTests(unittest.TestCase):
    def test_unequal_width_orthogonal_fixture_all_states_nonempty(self):
        fixture = sto.unequal_width_orthogonal_fixture(
            block_size=(1, 2, 4), ncell=(12, 16, 20), bond=(4.0, 0.0, 0.0, 0.01),
        )
        self.assertEqual(fixture.expected_topology.buffer_width, (4, 2, 1))
        counts = fixture.expected_topology.block_counts()
        self.assertTrue(all(count > 0 for count in counts.values()), counts)

    def test_skew_cell_fixture_all_states_nonempty(self):
        fixture = sto.skew_cell_fixture(
            cell_vectors=_SKEW_CELL, block_size=(1, 2, 2), ncell=(8, 8, 8),
            bond=(0.5, 1.0, 0.0, 0.01),
        )
        self.assertEqual(fixture.expected_topology.buffer_width, (2, 1, 1))
        counts = fixture.expected_topology.block_counts()
        self.assertTrue(all(count > 0 for count in counts.values()), counts)

    def test_deterministic_regeneration(self):
        first = sto.unequal_width_orthogonal_fixture(
            block_size=(1, 2, 4), ncell=(12, 16, 20), bond=(4.0, 0.0, 0.0, 0.01),
        )
        second = sto.unequal_width_orthogonal_fixture(
            block_size=(1, 2, 4), ncell=(12, 16, 20), bond=(4.0, 0.0, 0.0, 0.01),
        )
        self.assertEqual(first.jfile_text, second.jfile_text)
        self.assertEqual(first.mask_text, second.mask_text)
        self.assertEqual(first.manifest["content_sha256"], second.manifest["content_sha256"])
        self.assertEqual(manifest_json(first.manifest), manifest_json(second.manifest))

    def test_changing_a_parameter_changes_the_hash(self):
        baseline = sto.unequal_width_orthogonal_fixture(
            block_size=(1, 2, 4), ncell=(12, 16, 20), bond=(4.0, 0.0, 0.0, 0.01),
        )
        changed = sto.unequal_width_orthogonal_fixture(
            block_size=(1, 2, 4), ncell=(12, 16, 20), bond=(4.0, 0.0, 0.0, 0.02),
        )
        self.assertNotEqual(baseline.manifest["content_sha256"], changed.manifest["content_sha256"])

    def test_rejects_isotropic_block_size(self):
        with self.assertRaises(ValueError):
            sto.unequal_width_orthogonal_fixture(
                block_size=(2, 2, 2), ncell=(12, 16, 20), bond=(4.0, 0.0, 0.0, 0.01),
            )

    def test_rejects_nondivisible_ncell(self):
        with self.assertRaises(ValueError):
            sto.unequal_width_orthogonal_fixture(
                block_size=(1, 2, 4), ncell=(12, 15, 20), bond=(4.0, 0.0, 0.0, 0.01),
            )

    def test_rejects_orthogonal_cell_for_skew_fixture(self):
        with self.assertRaises(ValueError):
            sto.skew_cell_fixture(
                cell_vectors=sto._IDENTITY_CELL, block_size=(1, 2, 2), ncell=(8, 8, 8),
                bond=(1.0, 0.0, 0.0, 0.01),
            )

    def test_rejects_host_too_small_for_nonempty_coarse(self):
        """A host too small to leave every state nonempty must fail clearly,
        not silently return a degenerate (e.g. all-buffer) fixture."""
        with self.assertRaises(ValueError):
            sto.unequal_width_orthogonal_fixture(
                block_size=(1, 2, 4), ncell=(4, 4, 4), bond=(4.0, 0.0, 0.0, 0.01),
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
