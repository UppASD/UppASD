# RCG-04E: E2E-MOVING-ALL-COARSE, block_size_x=1 (finest sweep resolution)

**Purpose:** finest-resolution point in the RCG-04E all-coarse spatial
convergence sweep. Compared against `../moving_all_fine_wide` (the accepted
atomistic oracle) by `run_moving_all_coarse.py`. Full analysis:
`docs/RCG-04_MOVING_E2E_EVIDENCE.md` (RCG-04E section).

**Construction:** same `ncell 24 2 2` host geometry, same byte-identical
conical-spiral `momfile` (see `../moving_all_fine_wide/README.md`), same
Hamiltonian/integration parameters (exchange-only, `damping 0.05`,
`timestep 1.0e-16`, `Nstep 50`, `do_tottraj Y`/`tottraj_step 5`) as every
other fixture in this sweep. The only intentional difference across the
`moving_all_coarse_bs{1,2,4,8}` family is `block_size_x` and
`cg_static_mask_file mask.dat` (a comment-only file: every one-based block
id is omitted, so every block defaults to `COARSE` -- same convention as
`../static_all_coarse/mask.dat`).

**Block topology:** `block_size_x/y/z = 1 2 2` on `ncell 24 2 2` gives
`blocks_x = n1/block_size_x = 24`, `blocks_y = blocks_z = 1` (the state is
uniform in `y`/`z`, so fully coarsening those directions loses no
information -- only `block_size_x`, aligned with the modulation axis,
matters for this mode). Physical block length along `x`:
`block_length_x = block_size_x * alat = 1.0`. Modulation wavevector
`q = 2*pi*turns/n1 = 2*pi/24 = 0.261799` rad per unit length (`alat=1`).
`q * block_length_x = 0.261799` rad (`15.0` degrees) -- 24 blocks per full
spiral wavelength, comfortably inside the long-wave continuum-coarsening
regime.

**Reconstructed trajectory:** `reconstruct_coarse_atoms`
(`source/CoarseGraining/adaptivecgproduction.f90:1152`, called
unconditionally every step from `adaptive_cg_cpu_step`, before
`do_tottraj`'s per-step write) broadcasts each block's macrospin direction
to every atom it owns using the `ALIGNED` reconstruction scheme (production
default; no `cg_reconstruction` override here). This means the ordinary
`moment.<simid>.out` trajectory this fixture emits (`do_tottraj Y`) already
*is* the reconstructed per-atom trajectory the RCG-04E prompt requires --
no additional postprocessing step reconstructs it separately.
