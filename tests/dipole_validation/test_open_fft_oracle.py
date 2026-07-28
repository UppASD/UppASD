#!/usr/bin/env python3
"""Host-only acceptance tests for Luna's independent OPEN_FFT oracle."""
from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import open_fft_oracle as oracle  # noqa: E402


GOLDENS = json.loads((HERE / "open_fft_goldens.json").read_text())


def case_from_record(record: dict) -> oracle.FiniteCase:
    moments = tuple(
        tuple(tuple(tuple(vector) for vector in cell) for cell in ensemble)
        for ensemble in record["moments_mu_b"]
    )
    return oracle.FiniteCase(
        record["name"], tuple(record["grid"]), tuple(record["primitive_vectors"]),
        tuple(tuple(offset) for offset in record["basis_offsets"]), moments, record["alat_m"],
    )


def assert_nested_close(test: unittest.TestCase, actual, expected, places: int = 13) -> None:
    if isinstance(expected, list):
        test.assertEqual(len(actual), len(expected))
        for got, want in zip(actual, expected):
            assert_nested_close(test, got, want, places)
    else:
        test.assertAlmostEqual(actual, expected, places=places)


class OpenFftOracleTests(unittest.TestCase):
    def test_frozen_goldens_match_the_independent_evaluator(self) -> None:
        self.assertEqual(GOLDENS["oracle"], "independent finite open-boundary point-dipole sum")
        self.assertEqual(len(GOLDENS["cases"]), 10)
        for record in GOLDENS["cases"]:
            result = oracle.evaluate(case_from_record(record))
            assert_nested_close(self, result["fields_dimensionless"], record["fields_dimensionless"])
            assert_nested_close(self, result["fields_t"], record["fields_t"])
            assert_nested_close(self, result["total_energy_j"], record["total_energy_j"])
            assert_nested_close(self, result["energy_per_atom_mry"], record["energy_per_atom_mry"])

    def test_pair_sign_and_exact_self(self) -> None:
        cases = {case.name: case for case in oracle.deterministic_cases()}
        axial = oracle.evaluate(cases["pair_2x1x1_axial"])["fields_dimensionless"][0]
        transverse = oracle.evaluate(cases["pair_2x1x1_transverse"])["fields_dimensionless"][0]
        self.assertEqual(axial, (((2.0, 0.0, 0.0),), ((2.0, 0.0, 0.0),)))
        self.assertEqual(transverse, (((0.0, 0.0, -1.0),), ((0.0, 0.0, -1.0),)))
        single = oracle.evaluate(cases["single_1x1x1"])
        self.assertEqual(single["fields_dimensionless"], ((((0.0, 0.0, 0.0),),),))
        self.assertEqual(single["total_energy_j"], (-0.0,))

    def test_units_are_applied_after_dimensionless_tensor(self) -> None:
        case = next(case for case in oracle.deterministic_cases() if case.name == "skew_2x2x1")
        result = oracle.evaluate(case)
        scale = oracle.physical_prefactor(case.alat_m)
        for ensemble in range(len(case.moments_mu_b)):
            for cell in range(case.active_cells):
                for basis in range(case.basis):
                    for component in range(3):
                        self.assertEqual(
                            result["fields_t"][ensemble][cell][basis][component],
                            result["fields_dimensionless"][ensemble][cell][basis][component] * scale,
                        )

    def test_padded_storage_has_zero_gap_and_active_source_extent(self) -> None:
        case = next(case for case in oracle.deterministic_cases() if case.name == "nonuniform_2x3x1")
        padded = (4, 6, 1)
        kernel = oracle.embed_kernel(case, padded)
        fft_cells = padded[0] * padded[1] * padded[2]
        used_q = {(dx % padded[0], dy % padded[1], dz % padded[2])
                  for dx in range(-1, 2) for dy in range(-2, 3) for dz in range(0, 1)}
        for qz in range(padded[2]):
            for qy in range(padded[1]):
                for qx in range(padded[0]):
                    if (qx, qy, qz) not in used_q:
                        qcell = oracle.fft_cell(padded, (qx, qy, qz))
                        self.assertTrue(all(value == 0.0 for value in kernel[qcell::fft_cells]))
        packed = oracle.embed_active_moments(case, 0, padded)
        self.assertEqual(len(packed), fft_cells * 3 * case.basis)
        active_q = {oracle.fft_cell(padded, cell) for cell in oracle.cell_coordinates(case.grid)}
        for qcell in range(fft_cells):
            if qcell not in active_q:
                self.assertTrue(all(value == 0.0 for value in packed[qcell::fft_cells]))

    def test_same_embedded_kernel_matches_finite_sum_without_prefactor(self) -> None:
        for case in oracle.deterministic_cases():
            expected = oracle.dimensionless_fields(case)
            for ensemble in range(len(case.moments_mu_b)):
                actual = oracle.padded_direct_convolution(case, ensemble)
                for cell in range(case.active_cells):
                    for basis in range(case.basis):
                        for component in range(3):
                            self.assertAlmostEqual(
                                actual[cell][basis][component],
                                expected[ensemble][cell][basis][component],
                                places=12,
                                msg=f"{case.name} ensemble={ensemble} cell={cell} basis={basis}",
                            )

    def test_delta_sources_at_every_active_corner_do_not_wrap(self) -> None:
        case = next(case for case in oracle.deterministic_cases() if case.name == "nonuniform_2x3x1")
        corners = [(x, y, 0) for x in (0, case.grid[0] - 1) for y in (0, case.grid[1] - 1)]
        for source_xyz in corners:
            source_cell = oracle.cell_index(case.grid, *source_xyz)
            moments = tuple(
                ((1.0, 0.0, 0.0),) if cell == source_cell else ((0.0, 0.0, 0.0),)
                for cell in range(case.active_cells)
            )
            impulse = replace(case, moments_mu_b=(moments,))
            expected = oracle.dimensionless_fields(impulse)[0]
            actual = oracle.padded_direct_convolution(impulse, 0)
            for cell in range(case.active_cells):
                self.assertEqual(actual[cell][0], expected[cell][0], msg=f"corner={source_xyz}, target={cell}")

    def test_no_opposite_face_wrap_and_g3_fastest_ordering(self) -> None:
        case = next(case for case in oracle.deterministic_cases() if case.name == "g3_gt_g2_2x1x3")
        self.assertGreater(case.grid[2], case.grid[1])
        # The source is at (0,0,0); target (1,0,2) must use d=(+1,0,+2),
        # not a periodic shortcut such as (-1,0,-1).
        source_moments = tuple(
            tuple(((1.0, 0.0, 0.0),) if cell == 0 else ((0.0, 0.0, 0.0),))
            for cell in range(case.active_cells)
        )
        impulse = replace(case, moments_mu_b=(source_moments,))
        got = oracle.dimensionless_fields(impulse)[0][5][0]
        displacement = (1.0, 0.0, 2.0)
        want = oracle.tensor_vector(oracle.point_tensor(displacement), (1.0, 0.0, 0.0))
        self.assertEqual(got, want)

    def test_energy_finite_difference_derivative(self) -> None:
        case = next(case for case in oracle.deterministic_cases() if case.name == "basis_2x1x1_nonuniform")
        epsilon = 1.0e-7
        def shifted(delta: float) -> oracle.FiniteCase:
            moments = [[list(vector) for vector in cell] for cell in case.moments_mu_b[0]]
            moments[0][1][0] += delta
            return replace(case, moments_mu_b=(tuple(tuple(tuple(vector) for vector in cell) for cell in moments),))
        plus = oracle.evaluate(shifted(epsilon))["total_energy_j"][0]
        minus = oracle.evaluate(shifted(-epsilon))["total_energy_j"][0]
        numerical = (plus - minus) / (2.0 * epsilon)
        base = oracle.evaluate(case)
        analytic = -oracle.MU_B_J_PER_T * oracle.physical_prefactor(case.alat_m) * \
            base["fields_dimensionless"][0][0][1][0]
        self.assertAlmostEqual(numerical, analytic, delta=max(abs(analytic) * 2.0e-10, 1.0e-40))

    def test_ensemble_channels_are_distinct(self) -> None:
        case = next(case for case in oracle.deterministic_cases() if case.name == "four_distinct_ensembles")
        result = oracle.evaluate(case)
        self.assertEqual(len(result["fields_t"]), 4)
        self.assertEqual(len(set(result["total_energy_j"])), 4)


if __name__ == "__main__":
    unittest.main()
