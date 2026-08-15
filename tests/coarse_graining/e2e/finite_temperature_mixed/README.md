# RCG-10-F1 E2E-FINITE-T-MIXED: finite temperature with coarse blocks present

**Purpose:** hold open the accepted temperature boundary — *finite temperature
at the atomistic level, coarse blocks at T=0* — in the one configuration that
actually exercises both halves at once. Human capability decision, Anders
Bergman, 2026-08-15; recorded in `docs/CG-13_RELEASE_VALIDATION.md` under
"Temperature" and analysed in `docs/RCG-10_RELEASE_RECONCILIATION.md` §9.1.

**Why the existing `../finite_temperature` case is not enough:** it is all-fine
by construction (`initial_coarse=0`, `resolution_counts fine=6 interface=0
coarse=0`), so it exercises the thermal path but never the boundary. Before
this fixture existed, nothing in the suite ran `Temp > 0` with a coarse block
present, and the combination was reachable without any test observing it.

**Construction:** byte-identical to `../static_mixed` except for `simid`,
`temp 0.0` -> `temp 10.0`, and an explicit `tseed 4711`. Keeping everything
else fixed is the point: `../static_mixed` is the T=0 reference this case is
compared against, so the only difference between the two trajectories is the
thermal field. The static mask (`mask.dat`, one `FINE` seed) with
`cg_buffer_blocks 0` yields `fine=1 interface=4 coarse=1` on the 6-block grid.

**What the harness asserts** (`tests/coarse_graining/run_production_e2e.py`):

1. the run is accepted (`AdaptiveCG: capability accepted`) and exits 0;
2. `initial_coarse >= 1` — the case really is mixed, not silently all-fine;
3. `active_atoms > 0` — fine atoms really are being integrated;
4. its `trajectory_checksums direction_sum` differs from `../static_mixed`'s by
   more than `1e-6` (observed `48.051222` against `48.000000`).

**Negative control, executed 2026-08-15:** neutralizing the thermal field
collapses the warm trajectory onto the cold one and assertion 4 fails with
`mixed-resolution thermal field did not move the fine atoms: warm=48.0
cold=48.0`. That is the silent failure this fixture exists to catch — a thermal
field that stops reaching fine atoms once any block coarsens would otherwise
leave every other test in the suite green.

**`tseed` is explicit** so the trajectory is reproducible run to run
(bit-identical across repeated executions) rather than relying on a default
seed that could change.

**Scope:** this is a boundary and regression fixture, not a finite-temperature
physics oracle. The coarse blocks carry no thermal fluctuations, so no
finite-temperature *coarse* claim is made or tested here; see CG-13's
"Temperature" error-budget paragraph.
