#!/usr/bin/env python3
"""Unit tests for static_topology_oracle.py (RCG-04F)."""

from __future__ import annotations

import unittest

import static_topology_oracle as sto
import torque_oracle as orc
from moving_state_generator import Geometry

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

    def test_rejects_block_grid_not_spanning_yz(self):
        with self.assertRaises(sto.UnsupportedTopologyError):
            sto.compute_expected_topology(
                GEOMETRY, SHELLS, block_size_x=1, block_size_y=1, block_size_z=2,
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


if __name__ == "__main__":
    unittest.main(verbosity=2)
