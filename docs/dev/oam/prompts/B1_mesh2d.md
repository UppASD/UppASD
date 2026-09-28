# B1 — `mesh2d` module  ·  Sonnet  ·  commit `[OAM-B1]`  ·  gate G1 follows

Implement Blueprint B section **B1**. Paste the shared preamble first. You may run this in parallel with Blueprint A; the files don't overlap.

## Key facts

- `delaunay_tri_tri` (`topology.f90` 329–386) wraps indices but not coordinates. On 41×41 the total mesh area is 6400 instead of 1600.
- The fix is minimum-image *local* vertex positions per triangle. **Never write to `coord`.**
- Orientation must be counter-clockwise by signed area. The current `abs()` is harmless for areas but flips the sign of gradients.
- Print exactly the C7 diagnostic line. The harness parses it.

## Acceptance

- `run_traj_checks.py` check C0 (four mesh cases) GREEN. The OAM checks are expected to FAIL at this point; say so.
- The new FEM-exactness test passes: linear ψ gives its exact gradient to 1e-12.
- Every existing consumer (chirality, Pontryagin/skyrmion number, mesh print) uses `mesh2d`. There is one triangulation build per run.
- The `tests/` regression suite passes, or the differences are listed for B2.

**Commit `[OAM-B1]`, report, stop.**
