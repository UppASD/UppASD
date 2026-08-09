#!/usr/bin/env python3
"""Focused tests for coarse_torque_oracle.py (RCG-04E)."""

from __future__ import annotations

import math
import unittest

import coarse_torque_oracle as cto
from moving_state_generator import Geometry, conical_spiral_state
from torque_oracle import ExchangeShell


class BlockAverageDirectionsTests(unittest.TestCase):
    def test_uniform_state_has_no_neighbor_misalignment(self) -> None:
        geometry = Geometry(na=1, n1=4, n2=1, n3=1, basis=((0.0, 0.0, 0.0),))
        directions = {atom: (1.0, 0.0, 0.0) for atom, *_ in geometry.iter_atoms()}
        _raw, normalized, raw_norm = cto.block_average_directions(geometry, 2, directions)
        self.assertEqual(len(normalized), 2)
        for norm in raw_norm.values():
            self.assertAlmostEqual(norm, 1.0, places=12)
        # every block points the same way -> zero neighbor angle, zero laplacian, zero torque
        report = cto.coarse_chain_nontriviality(
            normalized, raw_norm, block_length_x=2.0, stiffness_xx=1.0,
        )
        self.assertAlmostEqual(report.max_neighbor_angle_deg, 0.0, places=10)
        self.assertAlmostEqual(report.max_torque_proxy, 0.0, places=10)

    def test_conical_spiral_gives_nonuniform_block_average(self) -> None:
        geometry = Geometry(na=2, n1=24, n2=2, n3=2, basis=((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)))
        state = conical_spiral_state(
            geometry, cone_angle_deg=40.0, turns=1, axis=(0.0, 0.0, 1.0),
            phase_deg=0.0, modulation_cell_axis=0, moment_magnitude=2.23, landeg=2.0,
        )
        _raw, normalized, raw_norm = cto.block_average_directions(geometry, 4, state.direction_by_atom)
        self.assertEqual(len(normalized), 6)
        # neighbouring blocks are rotated relative to each other, not collinear
        angle = math.degrees(math.acos(
            max(-1.0, min(1.0, sum(a * b for a, b in zip(normalized[0], normalized[1]))))
        ))
        self.assertGreater(angle, 1.0)
        for norm in raw_norm.values():
            self.assertGreater(norm, 0.9)  # little intra-block averaging loss at this block size

    def test_rejects_block_size_not_dividing_n1(self) -> None:
        geometry = Geometry(na=1, n1=6, n2=1, n3=1, basis=((0.0, 0.0, 0.0),))
        directions = {atom: (1.0, 0.0, 0.0) for atom, *_ in geometry.iter_atoms()}
        with self.assertRaises(ValueError):
            cto.block_average_directions(geometry, 4, directions)


class StiffnessEstimateTests(unittest.TestCase):
    def test_single_shell_matches_hand_derivation(self) -> None:
        # A toy 1D-like chain: one shell at dx=1 (dy=dz=0), J=2.0. Two atoms
        # per cell so the geometry stays non-degenerate; only the x-shell
        # contributes to D_xx.
        geometry = Geometry(na=1, n1=8, n2=1, n3=1, basis=((0.0, 0.0, 0.0),))
        shells = (
            ExchangeShell(dx=1.0, dy=0.0, dz=0.0, jij=2.0),
        )
        dxx = cto.estimate_unregularized_stiffness_xx(geometry, shells, cell_volume=1.0)
        # Two bonds (+x, -x) each contribute J*dx**2 = 2.0*1.0 = 2.0 -> sum=4.0, /(2*V)=2.0
        self.assertAlmostEqual(dxx, 2.0, places=10)


class NonTrivialityFloorTests(unittest.TestCase):
    def test_wide_fixture_block_sizes_are_all_nontrivial(self) -> None:
        geometry = Geometry(na=2, n1=24, n2=2, n3=2, basis=((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)))
        state = conical_spiral_state(
            geometry, cone_angle_deg=40.0, turns=1, axis=(0.0, 0.0, 1.0),
            phase_deg=0.0, modulation_cell_axis=0, moment_magnitude=2.23, landeg=2.0,
        )
        shells = (
            ExchangeShell(dx=0.5, dy=0.5, dz=0.5, jij=1.33767484769984),
            ExchangeShell(dx=1.0, dy=0.0, dz=0.0, jij=0.75703576545650),
            ExchangeShell(dx=1.0, dy=1.0, dz=0.0, jij=-0.05975437643846),
            ExchangeShell(dx=1.5, dy=0.5, dz=0.5, jij=-0.08819834160658),
        )
        dxx = cto.estimate_unregularized_stiffness_xx(geometry, shells)
        self.assertGreater(dxx, 0.0)
        for block_size_x in (1, 2, 4, 8):
            _raw, normalized, raw_norm = cto.block_average_directions(
                geometry, block_size_x, state.direction_by_atom,
            )
            report = cto.coarse_chain_nontriviality(
                normalized, raw_norm, block_length_x=float(block_size_x), stiffness_xx=dxx,
            )
            self.assertGreater(report.max_torque_proxy, 1.0e-3, msg=f"block_size_x={block_size_x}")
            self.assertGreater(report.max_neighbor_angle_deg, 1.0, msg=f"block_size_x={block_size_x}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
