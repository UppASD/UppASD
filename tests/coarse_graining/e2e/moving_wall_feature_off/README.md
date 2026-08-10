# RCG-04G: E2E-MOVING-ADAPTIVE, physics reference (AdaptiveCG disabled)

**Purpose:** plain (`do_adaptive_cg` unset) production run of the identical
domain-wall-pair initial condition and Hamiltonian used by
`../moving_wall_adaptive`, providing the all-atomistic reference trajectory
`run_moving_adaptive_wall.py` compares the ADAPTIVE run against. See that
fixture's `README.md` and `docs/RCG-04_MOVING_E2E_EVIDENCE.md` (RCG-04G
section) for full construction details, drive-mechanism justification, and
observed dynamics -- not repeated here.

`restart_seed.out` and `kfile_cg_wall` are byte-identical to
`../moving_wall_adaptive`'s own copies (verified by content hash in
`run_moving_adaptive_wall.py`, not merely by convention). The only intended
`inpsd.dat` difference from `../moving_wall_adaptive` is the absence of
every `do_adaptive_cg`/`block_size_*`/`cg_*` key and `simid` (also verified
programmatically, via a normalized key-by-key comparison, not a visual
diff).
