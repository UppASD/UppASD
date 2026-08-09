# RCG-04E: E2E-MOVING-ALL-COARSE primary atomistic oracle

**Purpose:** the accepted atomistic reference trajectory that every
`../moving_all_coarse_bs*` fixture is compared against. Per the RCG-04E
governing rule ("use the accepted atomistic trajectory as the primary
oracle rather than presenting a misleading formula"), this fixture -- not a
from-scratch analytic coarse-dynamics formula -- carries the RCG-04E
accuracy claim. See `../moving_feature_off_wide/README.md` for why a fresh
parity re-check against feature-off is run at this geometry before this
fixture is trusted as that oracle.

**Construction:** AdaptiveCG-enabled (`cg_operator TENSOR`, `cg_mask_mode
STATIC`, no mask file -> every block fine), same `ncell 24 2 2`
(192 atoms), same conical-spiral `momfile`/`Initmag=8` construction as
`../moving_feature_off_wide` (see that fixture's `README.md` for the shared
generator call and byte-identity claim). `block_size_x/y/z 1 2 2` here is
irrelevant to the physics (every block is forced fine by `cg_mask_mode
STATIC` with no mask file, regardless of block shape) but is kept
numerically identical to `../moving_all_coarse_bs1`'s block shape for a
clean side-by-side reading.

**Why a new fixture instead of reusing `../moving_all_fine`:** RCG-04E's
governing rule requires the physical sample, mode wavevector, timestep,
total simulated time, and observable definitions to stay fixed while only
spatial resolution (block size) changes across the sweep. `../moving_all_fine`
uses `ncell 6 2 2` (48 atoms), which only admits block sizes up to `n1=6`
along the modulation axis -- too few divisors, and too short a wavelength in
units of the periodic box, to give three-plus block sizes that stay
comfortably inside the long-wave regime (see the RCG-04E evidence section's
`q * block_length` discussion). `ncell 24 2 2` keeps the same one-full-turn
conical spiral (`turns=1`) but spread over 24 cells instead of 6, giving
`initpropvec_x = 1/24` instead of `1/6` and admitting block sizes 1, 2, 4
(comfortably long-wave: 24, 12, 6 blocks per wavelength) plus 8 (marginal,
3 blocks per wavelength, included only as an explicitly out-of-regime data
point -- see the evidence document). This is the "carefully justified
long-wave variant" of the accepted RCG-04D conical mode the RCG-04E prompt
explicitly permits, not an unrelated new construction: same cone angle, same
axis, same `turns`, same Hamiltonian, same integration parameters: only the
periodic box length along the modulation axis changed, and it changed for a
documented spatial-resolution reason.
