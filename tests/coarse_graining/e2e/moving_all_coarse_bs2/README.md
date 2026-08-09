# RCG-04E: E2E-MOVING-ALL-COARSE, block_size_x=2

Second resolution point in the RCG-04E all-coarse spatial convergence
sweep; see `../moving_all_coarse_bs1/README.md` for the full shared
construction, byte-identity, and reconstruction-mechanism notes (identical
here except `block_size_x`).

**Block topology:** `block_size_x/y/z = 2 2 2` on `ncell 24 2 2` gives
`blocks_x = 12`, `blocks_y = blocks_z = 1`. Physical block length along
`x`: `2.0`. `q * block_length_x = 0.523599` rad (`30.0` degrees) -- 12
blocks per full spiral wavelength, still comfortably inside the long-wave
regime.
