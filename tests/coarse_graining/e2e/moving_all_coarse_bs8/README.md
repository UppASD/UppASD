# RCG-04E: E2E-MOVING-ALL-COARSE, block_size_x=8 (deliberately out-of-regime)

Fourth point in the RCG-04E all-coarse sweep; see
`../moving_all_coarse_bs1/README.md` for the full shared construction,
byte-identity, and reconstruction-mechanism notes (identical here except
`block_size_x`).

**Block topology:** `block_size_x/y/z = 8 2 2` on `ncell 24 2 2` gives
`blocks_x = 3`, `blocks_y = blocks_z = 1`. Physical block length along `x`:
`8.0`. `q * block_length_x = 2.094395` rad (`120.0` degrees) -- only 3
blocks per full spiral wavelength.

**This resolution is deliberately excluded from any convergence-order
claim** (governing rule: "do not claim an asymptotic order... from a regime
outside long-wave validity"). Three samples per wavelength is close enough
to the discrete Nyquist limit (2 samples per wavelength) that aliasing/
discretization error, not the leading-order long-wave truncation error the
other three resolutions probe, is expected to dominate. It is retained as a
single additional data point to show *where* the coarse operator's
long-wave approximation breaks down, not to extend the fitted convergence
trend from `bs1`/`bs2`/`bs4`.
