# B2 mesh regression log

Baseline: `b60b762`, built in `/private/tmp/uppasd-base-build`.

After: current `[OAM-B1]` worktree, built in `/private/tmp/uppasd-b1-build`.

All runs used `OMP_NUM_THREADS=1`, `mseed 1`, `tseed 1`, and `Nstep 2`. Skyrmion and chirality were run separately because enabling both legacy measurement paths in one run triggers an existing measurement-plumbing bus error in the baseline executable. The skyrmion example was `examples/SpecialFeatures/SkyrmionLattice`.

The `*.inpsd.append` files are applied to the corresponding case input. The triang2D overlays additionally select a short spin-only run: they disable the SLD/lattice initial phase and unrelated spectra/correlation measurements. This is required because the unmodified `tests/triang2D` SLD initial phase crashes on the pinned baseline before measurements begin.

Commands, from a copied case directory containing the original case data:

```text
OMP_NUM_THREADS=1 /private/tmp/uppasd-base-build/bin/uppasd
OMP_NUM_THREADS=1 /private/tmp/uppasd-b1-build/bin/uppasd
```

The full result table and the gate-G1 finding are recorded in `../OAM_QUESTIONS.md`. No simulation outputs or binary artifacts are part of this commit.
