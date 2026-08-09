# RCG-04E: E2E-MOVING-ALL-COARSE, block_size_x=4

Third resolution point in the RCG-04E all-coarse spatial convergence sweep;
see `../moving_all_coarse_bs1/README.md` for the full shared construction,
byte-identity, and reconstruction-mechanism notes (identical here except
`block_size_x`).

**Block topology:** `block_size_x/y/z = 4 2 2` on `ncell 24 2 2` gives
`blocks_x = 6`, `blocks_y = blocks_z = 1`. Physical block length along `x`:
`4.0`. `q * block_length_x = 1.047198` rad (`60.0` degrees) -- 6 blocks per
full spiral wavelength: the coarsest point still treated as inside the
long-wave regime in this sweep's interpretation (a diminishing but still
meaningful number of samples per wavelength).
