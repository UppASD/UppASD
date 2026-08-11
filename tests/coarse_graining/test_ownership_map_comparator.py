#!/usr/bin/env python3
"""Unit tests for ownership_map_comparator.py (RCG-05C).

Fast, host-only, no executable required: the CPU/CUDA stdout strings below
are verbatim excerpts (only the ``AdaptiveCG:``/``Gpu: AdaptiveCG:`` lines
the parsers actually match) captured from a real run of this repository's
tracked ``e2e/ownership_aniso_buffer`` fixture -- CPU via
``build_rcg05b_cpu/bin/sd.f95``, CUDA via ``build_rcg05a_cuda/bin/sd.f95.cuda``
-- not synthesized by hand, so the parsing/comparison logic is exercised
against production's actual output format. The real, fresh, full build/test
evidence lives in ``docs/RCG-05_GEOMETRY_OWNERSHIP_EVIDENCE.md`` (RCG-05C
section); this file only proves the comparator's own logic is correct.
"""

from __future__ import annotations

import unittest

import ownership_map_comparator as omc
import static_topology_oracle as sto
from moving_state_generator import Geometry
from torque_oracle import parse_jfile_shells

# Real captured stdout excerpts, ``e2e/ownership_aniso_buffer`` (block_size
# 1/2/3, ncell 6 10 9, single FINE seed at block 1, single dx=2.0 shell).
CPU_STDOUT = (
    "AdaptiveCG: resolution_state label=initial step=0 values="
    "2,1,1,0,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,1,1,0,1,1,"
    "1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,0,1,1,"
    "0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1\n"
    "AdaptiveCG: resolution_state label=final step=1 values="
    "2,1,1,0,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,2,1,1,0,1,1,"
    "1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,2,1,1,0,1,1,1,1,1,0,1,1,"
    "0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1\n"
)
GPU_STDOUT = (
    "Gpu: AdaptiveCG final_state values="
    "2,1,1,1,1,1,1,1,1,1,1,1,2,2,1,1,2,2,2,2,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1,1,"
    "1,1,1,1,1,1,2,2,1,1,2,2,2,2,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,"
    "2,2,1,1,2,2,2,2,1,1,2,2,1,1,1,1,1,1\n"
)

GEOMETRY = Geometry(na=2, n1=6, n2=10, n3=9, basis=((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)))
JFILE_TEXT = "1 1 2.0 0.0 0.0 1.0\n"
SHELLS = parse_jfile_shells(JFILE_TEXT)
BLOCK_SIZE = dict(block_size_x=1, block_size_y=2, block_size_z=3)


class ParseFinalOwnershipMapTests(unittest.TestCase):
    def test_parses_cpu_final_line(self):
        ownership = omc.parse_final_ownership_map(CPU_STDOUT, "CPU")
        self.assertEqual(ownership.nblocks, 90)
        self.assertEqual(ownership.block_counts(), {"fine": 3, "interface": 42, "coarse": 45})
        self.assertEqual(ownership.fine_block_ids(), {1, 31, 61})

    def test_parses_gpu_final_state_line(self):
        ownership = omc.parse_final_ownership_map(GPU_STDOUT, "CUDA")
        self.assertEqual(ownership.nblocks, 90)
        self.assertEqual(ownership.block_counts(), {"fine": 25, "interface": 65, "coarse": 0})

    def test_uses_the_last_occurrence_not_the_first(self):
        # CPU_STDOUT contains an ``initial`` sample before ``final``; the
        # parser must return the *final* one (run_production_e2e.final_state
        # already guarantees this -- regression-checked here too).
        ownership = omc.parse_final_ownership_map(CPU_STDOUT, "CPU")
        self.assertEqual(ownership.block_state_by_id[31], sto.FINE)

    def test_missing_final_state_raises(self):
        with self.assertRaises(AssertionError):
            omc.parse_final_ownership_map("no adaptive diagnostics here", "CPU")


class CompareOwnershipMapsTests(unittest.TestCase):
    def test_identical_maps_match(self):
        a = omc.OwnershipMap(source="CPU", block_state_by_id={1: 2, 2: 0})
        b = omc.OwnershipMap(source="CUDA", block_state_by_id={1: 2, 2: 0})
        comparison = omc.compare_ownership_maps(a, b, label="identity")
        self.assertTrue(comparison.matches)
        self.assertEqual(comparison.mismatched_block_ids, ())

    def test_mismatch_is_reported_by_block_id_and_state_never_only_a_count(self):
        a = omc.OwnershipMap(source="CPU", block_state_by_id={1: 2, 2: 0, 3: 0})
        b = omc.OwnershipMap(source="CUDA", block_state_by_id={1: 2, 2: 1, 3: 0})
        comparison = omc.compare_ownership_maps(a, b, label="one-block-diff")
        self.assertFalse(comparison.matches)
        self.assertEqual(comparison.mismatched_block_ids, (2,))
        self.assertEqual(comparison.mismatch_detail, {2: (0, 1)})

    def test_rejects_maps_with_different_block_id_sets(self):
        a = omc.OwnershipMap(source="CPU", block_state_by_id={1: 2})
        b = omc.OwnershipMap(source="CUDA", block_state_by_id={1: 2, 2: 0})
        with self.assertRaises(omc.MapShapeMismatchError):
            omc.compare_ownership_maps(a, b, label="shape-mismatch")

    def test_real_cpu_vs_gpu_maps_disagree_on_47_of_90_blocks(self):
        cpu = omc.parse_final_ownership_map(CPU_STDOUT, "CPU")
        gpu = omc.parse_final_ownership_map(GPU_STDOUT, "CUDA")
        comparison = omc.compare_ownership_maps(cpu, gpu, label="cpu-vs-cuda")
        self.assertFalse(comparison.matches)
        self.assertEqual(len(comparison.mismatched_block_ids), 47)


class SelfConsistencyCheckTests(unittest.TestCase):
    """The RCG-05C mechanism that isolates the buffer-width defect: does a
    backend's own reported map match the correct (per-axis) or the
    isotropic (scalar-collapsed) re-dilation of its own FINE seed set."""

    def test_cpu_matches_the_correct_oracle_not_the_isotropic_one(self):
        cpu = omc.parse_final_ownership_map(CPU_STDOUT, "CPU")
        result = omc.self_consistency_check(cpu, GEOMETRY, SHELLS, **BLOCK_SIZE)
        self.assertEqual(result.fine_seed_ids, frozenset({1, 31, 61}))
        self.assertEqual(result.correct_buffer_width, (2, 1, 1))
        self.assertTrue(result.matches_correct_oracle)
        self.assertFalse(result.matches_isotropic_oracle)

    def test_gpu_matches_the_isotropic_oracle_not_the_correct_one(self):
        gpu = omc.parse_final_ownership_map(GPU_STDOUT, "CUDA")
        result = omc.self_consistency_check(gpu, GEOMETRY, SHELLS, **BLOCK_SIZE)
        self.assertEqual(result.isotropic_buffer_width, (2, 2, 2))
        self.assertFalse(result.matches_correct_oracle)
        self.assertTrue(result.matches_isotropic_oracle)

    def test_a_fabricated_correctly_dilated_map_matches_only_the_correct_oracle(self):
        # Negative control: a map that *is* the correct dilation of a
        # different seed set must not spuriously match the isotropic oracle
        # too (i.e. the two oracle checks are not vacuously both true).
        topology = sto.compute_expected_topology(
            GEOMETRY, SHELLS, fine_block_ids={5}, **BLOCK_SIZE,
        )
        fabricated = omc.OwnershipMap(source="synthetic", block_state_by_id=dict(topology.block_state_by_id))
        result = omc.self_consistency_check(fabricated, GEOMETRY, SHELLS, **BLOCK_SIZE)
        self.assertTrue(result.matches_correct_oracle)
        self.assertFalse(result.matches_isotropic_oracle)


class PeriodicWrapAxesTests(unittest.TestCase):
    def test_wrap_is_exercised_on_every_axis_of_the_real_fixture(self):
        cpu = omc.parse_final_ownership_map(CPU_STDOUT, "CPU")
        topology = sto.compute_expected_topology(
            GEOMETRY, SHELLS, fine_block_ids=cpu.fine_block_ids(), **BLOCK_SIZE,
        )
        self.assertEqual(omc.periodic_wrap_axes(topology), (True, True, True))

    def test_a_single_axis_with_no_wrap_headroom_is_not_reported_as_exercised(self):
        # A block_size_z equal to n3 collapses the z grid to one block, so
        # the z axis can never require a wrapped (as opposed to raw, both
        # identically zero) reduction -- must be reported False, not True.
        geometry = Geometry(na=2, n1=6, n2=10, n3=3, basis=((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)))
        topology = sto.compute_expected_topology(
            geometry, SHELLS, block_size_x=1, block_size_y=2, block_size_z=3, fine_block_ids={1},
        )
        wrap = omc.periodic_wrap_axes(topology)
        self.assertFalse(wrap[2])


class BondCoverageTests(unittest.TestCase):
    def test_matching_maps_report_zero_disagreeing_endpoints(self):
        cpu = omc.parse_final_ownership_map(CPU_STDOUT, "CPU")
        topology = sto.compute_expected_topology(
            GEOMETRY, SHELLS, fine_block_ids=cpu.fine_block_ids(), **BLOCK_SIZE,
        )
        result = omc.bond_coverage(GEOMETRY, SHELLS, topology.atom_block_by_atom, cpu, cpu)
        self.assertGreater(result.total_bonds, 0)
        self.assertEqual(result.disagreeing_bond_endpoints, ())

    def test_real_cpu_vs_gpu_maps_disagree_on_576_bond_endpoints(self):
        cpu = omc.parse_final_ownership_map(CPU_STDOUT, "CPU")
        gpu = omc.parse_final_ownership_map(GPU_STDOUT, "CUDA")
        topology = sto.compute_expected_topology(
            GEOMETRY, SHELLS, fine_block_ids=cpu.fine_block_ids(), **BLOCK_SIZE,
        )
        result = omc.bond_coverage(GEOMETRY, SHELLS, topology.atom_block_by_atom, cpu, gpu)
        self.assertEqual(result.interface_bonds_a, 576)
        self.assertEqual(result.interface_bonds_b, 0)
        self.assertEqual(len(result.disagreeing_bond_endpoints), 576)


class EvidenceRecordTests(unittest.TestCase):
    def test_evidence_record_is_json_serializable_and_complete(self):
        import json

        cpu = omc.parse_final_ownership_map(CPU_STDOUT, "CPU")
        gpu = omc.parse_final_ownership_map(GPU_STDOUT, "CUDA")
        comparison = omc.compare_ownership_maps(cpu, gpu, label="cpu-vs-cuda")
        cpu_sc = omc.self_consistency_check(cpu, GEOMETRY, SHELLS, **BLOCK_SIZE)
        gpu_sc = omc.self_consistency_check(gpu, GEOMETRY, SHELLS, **BLOCK_SIZE)
        topology = sto.compute_expected_topology(
            GEOMETRY, SHELLS, fine_block_ids=cpu.fine_block_ids(), **BLOCK_SIZE,
        )
        wrap = omc.periodic_wrap_axes(topology)
        coverage = omc.bond_coverage(GEOMETRY, SHELLS, topology.atom_block_by_atom, cpu, gpu)
        record = omc.evidence_record(
            fixture_name="ownership_aniso_buffer", geometry=GEOMETRY,
            block_size=(1, 2, 3), cell_vectors=sto._IDENTITY_CELL, cg_buffer_blocks=0,
            maps={"CPU": cpu, "CUDA": gpu}, comparisons=[comparison],
            self_consistency=[cpu_sc, gpu_sc], wrap_axes=wrap, coverage=coverage,
        )
        serialized = json.dumps(record)  # must not raise
        self.assertIn("ownership_aniso_buffer", serialized)
        self.assertFalse(record["comparisons"][0]["matches"])
        self.assertTrue(record["self_consistency"][0]["matches_correct_oracle"])
        self.assertTrue(record["self_consistency"][1]["matches_isotropic_oracle"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
