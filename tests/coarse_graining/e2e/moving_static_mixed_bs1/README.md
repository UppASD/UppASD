# RCG-04F: E2E-MOVING-STATIC, block_size_x=1 (finest)

**Purpose:** genuine static fine/buffer(interface)/coarse decomposition of
the RCG-04E `../moving_all_fine_wide` conical-spiral moving state
(byte-identical `momfile`, verified by `run_moving_static_mixed.py`), used
to validate `StaticHybridOperator`'s ownership dispatch and interface
coupling against the accepted all-fine atomistic reference.

**Construction:** `cg_operator PROJECTED` (matching `../static_mixed`'s
convention for a genuine partial-fine mask, as opposed to the
`TENSOR`-labelled all-fine/all-coarse fixtures where the operator choice is
physically irrelevant), `cg_mask_mode STATIC`, `mask.dat` marks one-based
blocks 1-6 `FINE` (48 atoms) out of 24 total `block_size_x=1` blocks; every
other physical/integration parameter (`ncell 24 2 2`, `damping 0.05`,
`timestep 1.0e-16`, `Nstep 50`, `do_tottraj Y`/`tottraj_step 5`) is identical
to `../moving_all_fine_wide`.

**Expected topology (independently derived by
`static_topology_oracle.compute_expected_topology`, not read back from
production):** with this shared `jfile`'s maximum exchange-shell radius
(`sqrt(2.75)~1.6583`, the `(1.5,0.5,0.5)` shell) and `block_size_x=1`,
`buffer_width_blocks_x = ceil(1.6583/1 - 64*eps) + cg_buffer_blocks(=0) = 2`
blocks. Dilating the six-block FINE seed (ids 1-6) by 2 blocks each way
(periodic, `nblocks_x=24`) gives interface/buffer blocks {23,24,7,8} (4
blocks, 32 atoms) and leaves blocks {9..22} (14 blocks, 112 atoms) COARSE.
`run_moving_static_mixed.py` asserts the runtime `AdaptiveCG:
resolution_state`/`interface_atoms=` diagnostics agree with this expectation
before trusting any comparison against the all-fine reference.

**Refinement pair:** `../moving_static_mixed_bs2` uses `block_size_x=2`
with FINE seed blocks 1-3, chosen so both fixtures cover the *exact same*
physical atom partition (48 fine / 32 interface / 112 coarse atoms) at two
different block granularities -- a controlled spatial-refinement comparison
at fixed physical evolution, not a different experiment.

**Shifted-interface pair:** `../moving_static_mixed_bs1_shifted` uses the
same `block_size_x=1` but FINE seed blocks 13-18 (half the 24-cell period
away). See that fixture's `README.md` for why the uniform-pitch conical
spiral's translational symmetry predicts this should *not* meaningfully
change any reported error.
