#!/usr/bin/env python3
"""Self-checks for the independent periodic dipole Ewald evaluator."""
from __future__ import annotations

from periodic_ewald_reference import (
    add,
    apply_displacement_kernel,
    assert_reciprocal_identity,
    default_alpha,
    evaluate,
    evaluate_converged,
    finite_difference_field,
    generate_displacement_kernel,
    grid_positions,
    kernel_energy,
    reciprocal_alias_tensor,
    single_representative_reciprocal_tensor,
)


def close(left: float, right: float, tolerance: float = 2e-10) -> None:
    if abs(left - right) > tolerance * max(1.0, abs(left), abs(right)):
        raise AssertionError(f"{left:.17e} != {right:.17e}")


def compare(first, second, tolerance: float = 2e-10) -> None:
    close(first.energy, second.energy, tolerance)
    for lhs, rhs in zip(first.fields, second.fields):
        for left, right in zip(lhs, rhs):
            close(left, right, tolerance)


def compare_golden(result, energy: float, fields) -> None:
    close(result.energy, energy, tolerance=5e-13)
    for actual, expected in zip(result.fields, fields):
        for got, want in zip(actual, expected):
            close(got, want, tolerance=5e-13)


def compare_kernel_to_direct(grid, cell, basis, moments, *, alpha: float,
                             real_images=(7, 7, 7), reciprocal_images=(7, 7, 7)) -> None:
    positions = grid_positions(grid, cell, basis)
    direct = evaluate(positions, moments, cell, alpha=alpha,
                      real_images=real_images, reciprocal_images=reciprocal_images)
    kernel = generate_displacement_kernel(
        grid, cell, basis, alpha=alpha, real_images=real_images,
        reciprocal_images=reciprocal_images)
    fields = apply_displacement_kernel(kernel, grid, len(basis), moments)
    for actual, expected in zip(fields, direct.fields):
        for got, want in zip(actual, expected):
            close(got, want, tolerance=4e-11)
    close(kernel_energy(fields, moments), direct.energy, tolerance=4e-11)

    # K(d,a,b) = K(-d,b,a)^T.  The packed tensors are symmetric, so the
    # transpose is represented by the same six entries here.
    for (displacement, target_basis, source_basis), tensor in kernel.items():
        reverse = tuple((-displacement[axis]) % grid[axis] for axis in range(3))
        reciprocal = kernel[(reverse, source_basis, target_basis)]
        for left, right in zip(tensor, reciprocal):
            close(left, right, tolerance=4e-11)


def main() -> int:
    cubic = ((3.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 3.0))
    positions = ((0.0, 0.0, 0.0), (0.7, 0.4, 1.1))
    moments = ((1.0, -0.2, 0.4), (-0.3, 0.8, 0.1))

    # Geometry is checked independently for every cell used below.  In
    # particular, this catches accidentally using primitive vectors where the
    # full periodic supercell H is required.
    assert_reciprocal_identity(cubic)
    skew = ((3.1, 0.0, 0.0), (0.4, 2.8, 0.0), (0.2, 0.3, 3.3))
    assert_reciprocal_identity(skew)

    # The Ewald split parameter must not change a sufficiently converged sum.
    low_alpha = evaluate(positions, moments, cubic, alpha=0.7,
                         real_images=(5, 5, 5), reciprocal_images=(5, 5, 5))
    high_alpha = evaluate(positions, moments, cubic, alpha=1.1,
                          real_images=(5, 5, 5), reciprocal_images=(5, 5, 5))
    compare(low_alpha, high_alpha)

    # Automatic convergence checks expand real and reciprocal shells
    # independently and require two consecutive passing changes per side.
    converged = evaluate_converged(positions, moments, cubic, alpha=0.9,
                                   tolerance=1e-10, max_shell=10)
    assert converged.convergence is not None
    assert converged.convergence.real_images[0] >= 2
    assert converged.convergence.reciprocal_images[0] >= 2
    assert converged.convergence.real_shell_changes[-1] <= 1e-10
    assert converged.convergence.reciprocal_shell_changes[-1] <= 1e-10

    # Translating every macrocell by the same lattice vector changes no field.
    shifted = tuple(add(position, cubic[0]) for position in positions)
    compare(high_alpha, evaluate(shifted, moments, cubic, alpha=1.1,
                                real_images=(5, 5, 5), reciprocal_images=(5, 5, 5)))

    # Fixed independently reviewed values protect against an accidental sign,
    # self-term or reciprocal-tensor edit.  The first case is precisely the
    # block-one macrocell limit of one source in a cubic cell.
    cubic_golden = evaluate(((0.0, 0.0, 0.0),), ((1.0, -0.2, 0.4),), cubic,
                            alpha=0.9, real_images=(7, 7, 7), reciprocal_images=(7, 7, 7))
    compare_golden(cubic_golden, -0.09308422677303108,
                   ((0.1551403779550518, -0.031028075591010243, 0.06205615118202068),))

    skew_positions = ((0.0, 0.0, 0.0), (0.9, 0.7, 1.2))
    skew_result = evaluate(skew_positions, moments, skew, alpha=0.9,
                           real_images=(7, 7, 7), reciprocal_images=(7, 7, 7))
    compare_golden(skew_result, -0.264979280195023,
                   ((0.22112113378570364, -0.0332487057527215, 0.12398892258605723),
                    (0.08604483604361035, 0.311214561767032, 0.29433917818832656)))

    # Energy and field are one Hamiltonian: -dE/dM equals B.  Check every
    # component in a skew two-particle cell, not just a convenient orientation.
    for atom in range(len(moments)):
        for component in range(3):
            derivative = finite_difference_field(
                skew_positions, moments, skew, atom=atom, component=component,
                step=1e-5, alpha=0.9, real_images=(7, 7, 7), reciprocal_images=(7, 7, 7))
            close(derivative, skew_result.fields[atom][component], tolerance=2e-8)

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

    # Explicit small-grid displacement kernels cover the required indexing and
    # geometry matrix.  The direct Ewald evaluator is the cross-formulation
    # reference; no UppASD output is used.
    compare_kernel_to_direct(
        (1, 1, 1), cubic, ((0.0, 0.0, 0.0),),
        ((1.0, -0.2, 0.4),), alpha=0.9)
    auto_kernel = generate_displacement_kernel(
        (1, 1, 1), cubic, ((0.0, 0.0, 0.0),), tolerance=1e-10, max_shell=10)
    auto_fields = apply_displacement_kernel(auto_kernel, (1, 1, 1), 1,
                                             ((1.0, -0.2, 0.4),))
    auto_direct = evaluate(((0.0, 0.0, 0.0),), ((1.0, -0.2, 0.4),), cubic,
                           alpha=default_alpha(cubic), real_images=(7, 7, 7),
                           reciprocal_images=(7, 7, 7))
    for actual, expected in zip(auto_fields, auto_direct.fields):
        for got, want in zip(actual, expected):
            close(got, want, tolerance=1e-9)
    close(kernel_energy(auto_fields, ((1.0, -0.2, 0.4),)), auto_direct.energy,
          tolerance=1e-9)
    two_cells = ((6.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 3.0))
    compare_kernel_to_direct(
        (2, 1, 1), two_cells, ((0.0, 0.0, 0.0),),
        ((1.0, -0.2, 0.4), (-0.3, 0.8, 0.1)), alpha=0.8,
        real_images=(6, 6, 6), reciprocal_images=(6, 6, 6))
    non_cubic = ((6.0, 0.0, 0.0), (0.0, 9.0, 0.0), (0.0, 0.0, 3.0))
    non_cubic_moments = tuple((0.2 + 0.03 * index, -0.4 + 0.02 * index,
                               0.1 - 0.01 * index) for index in range(6))
    compare_kernel_to_direct(
        (2, 3, 1), non_cubic, ((0.0, 0.0, 0.0),), non_cubic_moments,
        alpha=0.8, real_images=(5, 5, 5), reciprocal_images=(5, 5, 5))
    skew_grid = ((6.0, 0.0, 0.0), (0.8, 3.0, 0.0), (0.3, 0.2, 3.4))
    compare_kernel_to_direct(
        (2, 1, 1), skew_grid, ((0.0, 0.0, 0.0),),
        ((1.0, -0.2, 0.4), (-0.3, 0.8, 0.1)), alpha=0.8,
        real_images=(6, 6, 6), reciprocal_images=(7, 7, 7))
    two_basis_cell = ((4.5, 0.0, 0.0), (0.4, 4.0, 0.0), (0.2, 0.3, 4.4))
    two_basis = ((0.0, 0.0, 0.0), (1.1, 0.7, 0.9))
    two_basis_moments = (
        (1.0, -0.2, 0.4), (-0.3, 0.8, 0.1),
    )
    compare_kernel_to_direct(
        (1, 1, 1), two_basis_cell, two_basis, two_basis_moments,
        alpha=0.8, real_images=(6, 6, 6), reciprocal_images=(7, 7, 7))

    # WP6 fixture: an L1_0-style two-basis lattice with deliberately skewed
    # tetragonal vectors.  This is independent of the C++/GPU fixture and
    # checks alpha, energy derivative, and basis-label invariance directly.
    l10_cell = ((7.6, 0.45, 0.2), (0.7, 3.7, 0.25), (0.3, 0.2, 5.3))
    l10_basis = ((0.0, 0.0, 0.0), (1.95, 1.65, 2.55))
    l10_positions = grid_positions((2, 1, 1), l10_cell, l10_basis)
    l10_moments = (
        (0.20, -0.40, 0.10), (-0.30, 0.80, 0.10),
        (0.26, -0.36, 0.08), (0.19, -0.08, -0.23),
    )
    l10_low = evaluate_converged(l10_positions, l10_moments, l10_cell, alpha=0.70,
                                 tolerance=1e-10, max_shell=16)
    l10_high = evaluate_converged(l10_positions, l10_moments, l10_cell, alpha=1.05,
                                  tolerance=1e-10, max_shell=16)
    compare(l10_low, l10_high, tolerance=3e-10)
    for atom in range(len(l10_moments)):
        for component in range(3):
            derivative = finite_difference_field(
                l10_positions, l10_moments, l10_cell, atom=atom, component=component,
                step=1e-5, alpha=1.05, real_images=l10_high.convergence.real_images,
                reciprocal_images=l10_high.convergence.reciprocal_images)
            close(derivative, l10_high.fields[atom][component], tolerance=3e-8)
    swapped_positions = (l10_positions[1], l10_positions[0], l10_positions[3], l10_positions[2])
    swapped_moments = (l10_moments[1], l10_moments[0], l10_moments[3], l10_moments[2])
    swapped_l10 = evaluate(swapped_positions, swapped_moments, l10_cell, alpha=1.05,
                           real_images=l10_high.convergence.real_images,
                           reciprocal_images=l10_high.convergence.reciprocal_images)
    close(swapped_l10.energy, l10_high.energy, tolerance=3e-10)
    for original, swapped in zip(l10_high.fields, (swapped_l10.fields[1], swapped_l10.fields[0],
                                                     swapped_l10.fields[3], swapped_l10.fields[2])):
        for left, right in zip(original, swapped):
            close(left, right, tolerance=3e-10)

    # Red regression for the current single-representative reciprocal builder:
    # q=0 on a 1x1x1 mesh is not the physical k=0 term.  Every nonzero
    # reciprocal vector aliases into that bin, so the full bin is nonzero.
    aliases = reciprocal_alias_tensor(
        (1, 1, 1), cubic, (0, 0, 0), alpha=0.9,
        reciprocal_images=(7, 7, 7))
    representative = single_representative_reciprocal_tensor(
        (1, 1, 1), cubic, (0, 0, 0), alpha=0.9)
    assert max(abs(value) for value in aliases) > 1e-3
    assert representative == (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # A non-cubic image extent is also a valid explicit convergence study.
    rectangular = evaluate(positions, moments, cubic, alpha=0.9,
                           real_images=(6, 5, 4), reciprocal_images=(6, 5, 4))
    compare(rectangular, evaluate(positions, moments, cubic, alpha=0.9,
                                  real_images=(7, 7, 7), reciprocal_images=(7, 7, 7)),
            tolerance=3e-9)
    print("PASS periodic dipole Ewald reference self-checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
