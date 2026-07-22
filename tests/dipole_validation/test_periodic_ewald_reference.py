#!/usr/bin/env python3
"""Self-checks for the independent periodic dipole Ewald evaluator."""
from __future__ import annotations

from periodic_ewald_reference import add, evaluate


def close(left: float, right: float, tolerance: float = 2e-10) -> None:
    if abs(left - right) > tolerance * max(1.0, abs(left), abs(right)):
        raise AssertionError(f"{left:.17e} != {right:.17e}")


def compare(first, second) -> None:
    close(first.energy, second.energy)
    for lhs, rhs in zip(first.fields, second.fields):
        for left, right in zip(lhs, rhs):
            close(left, right)


def compare_golden(result, energy: float, fields) -> None:
    close(result.energy, energy, tolerance=5e-13)
    for actual, expected in zip(result.fields, fields):
        for got, want in zip(actual, expected):
            close(got, want, tolerance=5e-13)


def main() -> int:
    cubic = ((3.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 3.0))
    positions = ((0.0, 0.0, 0.0), (0.7, 0.4, 1.1))
    moments = ((1.0, -0.2, 0.4), (-0.3, 0.8, 0.1))
    # The Ewald split parameter must not change a sufficiently converged sum.
    low_alpha = evaluate(positions, moments, cubic, alpha=0.7, real_images=(5, 5, 5), reciprocal_images=(5, 5, 5))
    high_alpha = evaluate(positions, moments, cubic, alpha=1.1, real_images=(5, 5, 5), reciprocal_images=(5, 5, 5))
    compare(low_alpha, high_alpha)

    # Translating every macrocell by the same lattice vector changes no field.
    shifted = tuple(add(position, cubic[0]) for position in positions)
    compare(high_alpha, evaluate(shifted, moments, cubic, alpha=1.1,
                                 real_images=(5, 5, 5), reciprocal_images=(5, 5, 5)))

    # Fixed independently reviewed values protect the oracle against an
    # accidental sign, self-term or reciprocal-tensor edit.  The first case is
    # precisely the block-one macrocell limit of one source in a cubic cell.
    cubic_golden = evaluate(((0.0, 0.0, 0.0),), ((1.0, -0.2, 0.4),), cubic,
                            alpha=0.9, real_images=(7, 7, 7), reciprocal_images=(7, 7, 7))
    compare_golden(cubic_golden, -0.09308422677303108,
                   ((0.1551403779550518, -0.031028075591010243, 0.06205615118202068),))

    skew = ((3.1, 0.0, 0.0), (0.4, 2.8, 0.0), (0.2, 0.3, 3.3))
    skew_positions = ((0.0, 0.0, 0.0), (0.9, 0.7, 1.2))
    skew_result = evaluate(skew_positions, moments, skew, alpha=0.9,
                           real_images=(7, 7, 7), reciprocal_images=(7, 7, 7))
    compare_golden(skew_result, -0.264979280195023,
                   ((0.22112113378570364, -0.0332487057527215, 0.12398892258605723),
                    (0.08604483604361035, 0.311214561767032, 0.29433917818832656)))

    # Ensembles are independent: a production batched PME path must agree with
    # a separate reference evaluation of every ensemble, including M=4.
    ensemble_moments = (
        moments,
        tuple((-x, -y, -z) for x, y, z in moments),
        ((0.4, 0.0, 1.0), (0.7, -0.5, 0.2)),
        ((-0.8, 0.3, 0.2), (0.1, 0.9, -0.4)),
    )
    ensemble_results = [evaluate(skew_positions, item, skew, alpha=0.9,
                                 real_images=(7, 7, 7), reciprocal_images=(7, 7, 7))
                        for item in ensemble_moments]
    compare(ensemble_results[0], skew_result)
    close(ensemble_results[0].energy, ensemble_results[1].energy)
    for original, inverted in zip(ensemble_results[0].fields, ensemble_results[1].fields):
        for left, right in zip(original, inverted):
            close(left, -right)  # global M -> -M makes B odd while E is even
    print("PASS periodic dipole Ewald reference self-checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
