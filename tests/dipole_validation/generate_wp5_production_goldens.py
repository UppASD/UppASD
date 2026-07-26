#!/usr/bin/env python3
"""Regenerate the versioned WP5 production oracle fixture.

This imports only the independent periodic Ewald reference.  It never calls
UppASD, Builder A, directConvolution, or a GPU.  Review the emitted JSON
before replacing wp5_production_goldens_v1.json.
"""
from __future__ import annotations

import json

from periodic_ewald_reference import cell_vector, evaluate_converged


CASES = {
    "one_cell_uniform": ((1, 1, 1), ((3.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 3.0))),
    "three_cell_uniform": ((3, 1, 1), ((9.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 3.0))),
}
MOMENT = (1.0, -0.2, 0.4)


def positions(grid, cell):
    return tuple(cell_vector(cell, (x / grid[0], y / grid[1], z / grid[2]))
                 for z in range(grid[2]) for y in range(grid[1]) for x in range(grid[0]))


def main() -> None:
    primitive_grid, primitive_cell = CASES["one_cell_uniform"]
    primitive = evaluate_converged(positions(primitive_grid, primitive_cell), (MOMENT,), primitive_cell,
                                   tolerance=1e-13, max_shell=16)
    cases = {
        "one_cell_uniform": {
            "grid": primitive_grid,
            "cell": primitive_cell,
            "moments": (MOMENT,),
            "fields": primitive.fields,
            "energy_total": primitive.energy,
        }
    }
    # A uniform repetition of a primitive P/P/P cell represents the identical
    # infinite crystal.  Derive this N=3 oracle from the independently
    # converged primitive reference, rather than re-summing a finite Ewald
    # shell in a larger, symmetry-equivalent supercell.
    grid, cell = CASES["three_cell_uniform"]
    replicas = grid[0] * grid[1] * grid[2]
    cases["three_cell_uniform"] = {
        "grid": grid,
        "cell": cell,
        "moments": (MOMENT,) * replicas,
        "fields": primitive.fields * replicas,
        "energy_total": primitive.energy * replicas,
    }
    print(json.dumps({
        "schema": "uppasd-wp5-production-oracle-v1",
        "provenance": {
            "generator": "generate_wp5_production_goldens.py",
            "oracle": "periodic_ewald_reference.evaluate_converged",
            "tolerance": 1e-13,
            "max_shell": 16,
            "replicated_uniform_supercells": "derived from the primitive independent oracle by periodic replication",
            "convention": "dimensionless B = Hessian(1/r) m; E = -0.5 sum m dot B; tin-foil 3D periodic",
        },
        "cases": cases,
    }, indent=2))


if __name__ == "__main__":
    main()
