# RCG-04G: E2E-MOVING-ADAPTIVE, ADAPTIVE mask mode

**Purpose:** accepted adaptive-resolution transitions during genuine
domain-wall motion, including physical block-boundary crossing and the
RCG-03 polarization-safety contract. See
`docs/RCG-04_MOVING_E2E_EVIDENCE.md` (RCG-04G section) for the full writeup,
including the atom-to-block indexing defect this fixture's construction
found and fixed in `source/CoarseGraining/blocktopology.f90`.

**Construction:** deterministic periodic two-kink Neel domain-wall pair
(`moving_state_generator.domain_wall_pair_state`, RCG-04B), restart-format
(`Initmag 4`), on the standard `ncell 24 2 2` host geometry
(`../posfile`/`../jfile`), plus a new `kfile_cg_wall` uniaxial anisotropy
(K1=-0.5 per basis site, easy axis (0,0,1), matching the wall's own easy
axis) -- much stronger than the shared `e2e/kfile_cg`'s K1=-0.002/-0.003,
chosen so the exchange/anisotropy competition gives the wall a genuine,
finite width comparable to the imposed `width_cells=1.0` rather than the
wall immediately blowing up toward a much wider natural equilibrium profile
(see the evidence document for the parameter search).

```python
import moving_state_generator as gen
geometry = gen.Geometry(na=2, n1=24, n2=2, n3=2,
                         basis=((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)))
state = gen.domain_wall_pair_state(
    geometry, axis_cell_index=0, wall_centers_cells=(9.0, 15.0),
    width_cells=1.0, easy_axis=(0.0, 0.0, 1.0), wall_type="NEEL",
    chirality=1, cant_deg=2.0, moment_magnitude=2.23, simid="cg105gw",
    separation_margin_widths=4.0,
)
```

`block_size_x 1` (block_size_y/z 2, matching n2/n3) gives 24 physical
blocks, one cell wide each -- fine enough that the wall's modest excursion
still crosses a real block boundary. `cg_operator PROJECTED`, `cg_mask_mode
ADAPTIVE`, `cg_selector MAX_ANGLE`, `cg_refine_threshold 0.20`,
`cg_coarsen_threshold 0.05` (non-saturated, unlike the pre-existing
`adaptive_mixed`/`parity_adaptive_*` fixtures -- see RCG-04A section 3.2),
`cg_polarization_threshold 0.9` (RCG-03, unweakened),
`cg_energy_jump_limit 1.0e-15`, `damping 0.05`, `timestep 1.0e-16`, `Nstep
900`, `do_tottraj Y`/`tottraj_step 25` (37 sampled steps).

**Drive mechanism (no external field):** `adaptivecgproduction.f90`'s
`validate_configuration` explicitly rejects any nonzero `hfield` ("external
or time-dependent fields are not supported by the first production CG
path") -- confirmed by reading the source before design, not assumed. The
drive is therefore the wall pair's own intrinsic exchange/anisotropy
relaxation under damped LLG, the same category of mechanism RCG-04D/E/F use
(a deliberately nonstationary initial condition relaxing under the real
Hamiltonian/integrator). Both walls in this construction share one global
`chirality` parameter (the generator does not support independently-signed
kinks), so the pair behaves as two same-handedness Neel solitons; the
observed dynamics (see below) is consistent with same-chirality kink-kink
**repulsion**, not attraction -- the confined "down" domain between the two
walls expands and both walls move outward.

**Observed dynamics** (identical in the `moving_wall_feature_off` physics
reference and this ADAPTIVE run, cross-validating the physical claim
independently of AdaptiveCG): both wall centres are essentially static for
~350 steps (initial local relaxation, not yet translational), then move
outward, crossing their nearest physical block boundary at step 400 (wall
at cell~9 crosses into block 8; wall at cell~15 crosses into block 15), and
continue a damped, oscillatory excursion (net final displacement ~0.17-0.18
cells; full parameter search in the evidence document shows this settles
into a longer-period oscillation, and a closer-spaced variant shows the
pair fully annihilating after ~1200 steps -- not used here, kept out of
scope to preserve a clean, non-ambiguous accepted case).

**Accepted transitions:** 11 blocks (1-6, 20-24, all comfortably away from
either wall) coarsen at the first synchronization step (`step=1`) --
initial-ownership setup, not a motion claim. One further block (13, deep in
the expanding confined domain's interior) coarsens at `step=788`, strictly
after the first wall crossing (`step=400`) -- the motion-driven accepted
transition this fixture's core claim rests on. All accepted coarsen events'
`polarization_ratio` exceed 0.9998, comfortably above the
`cg_polarization_threshold=0.9` gate.

**Reviewed limitation (not fabricated, not a defect):** no accepted
refine-request (coarse-to-atomistic) transition occurs in this specific
fixture -- the wall's excursion (crossing exactly one block boundary each)
never travels far enough back to threaten a block that coarsened during
initial setup. During parameter search (see the evidence document), a
tighter geometry (`ncell 16 2 2`, walls 6 cells apart) *did* produce genuine
refine-**requests** as the wall's larger excursion approached
already-coarsened blocks, but every one was rejected by `ALIGNED`
reconstruction (`outcome=reconstruction-rejected`, zero energy jump --
distinct from an energy-jump rejection), a separate, itself informative
finding not adopted for this fixture's accepted case (kept here as a
documented, reviewed limitation per the RCG-04G prompt's explicit
allowance, rather than re-tuned until an accepted refine appeared).

**Negative controls** (see `run_moving_adaptive_wall.py` for the always-run
in-memory trajectory mutations, and the evidence document for the two
disposable production-run controls this fixture's own claims specifically
require): (1) `cg_update_interval` raised above `Nstep` in a disposable
copy of this input reproduces zero transition events at all (confirming
the selector's periodic synchronization, not something else, produces the
observed transitions); (2) the identical construction at a much wider wall
separation (12 cells instead of 6) shows a net wall displacement of
~1.2e-4 cells over 3000 steps -- below this fixture's own
`MIN_WALL_DISPLACEMENT_CELLS` floor -- confirming the boundary-crossing
claim is not satisfied by a stationary/negligible-drive construction.
