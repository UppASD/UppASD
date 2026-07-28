# WP10.1 Luna acceptance package

This package owns the independent finite/open-boundary reference and red tests
for `OPEN_FFT`. It does not construct production GPU kernels or dispatch a GPU
Hamiltonian.

## Test-only Terra seam

`tests/dipole_validation/open_fft_test_seam.hpp` declares the exact observable
needed by `test_open_fft_layout_red.cpp`. Terra must provide its implementation
as a test-only adapter around the open runtime:

- `active_cells = G1*G2*G3`, `active_macros = active_cells*NA`, and
  `fft_cells = P1*P2*P3`;
- `packed_moments` and `padded_fields` in
  `fft_cell + fft_cells*field_batch(alpha,a,e)` order;
- `active_fields` in
  `alpha + 3*(macro + active_macros*ensemble)` order;
- a host copy of the explicitly zeroed padded moment buffer;
- dimensionless active fields and one dimensionless energy per ensemble;
- no Tesla prefactor, Ewald term, periodic builder, or production dispatch.

The adapter must accept the `Input` geometry and active moment vector without
requiring an `Nfft`-sized active source array. The red executable is deliberately
not in the current CMake targets because there is no implementation of this
test seam yet. Once Terra supplies it, compile and run:

```text
c++ -std=c++17 -Itests/dipole_validation \
  tests/dipole_validation/test_open_fft_layout_red.cpp \
  <Terra's test-seam implementation and GPU test libraries>
```

The expected first failure signatures against the current periodic-only class
are a missing `run` implementation, an active/padded count mismatch, nonzero
padding, a field at the opposite face from an active-corner delta, or a
dimensionless axial result other than `K_xx=2`.

## Independent derivation

For target `(g,a)` and source `(h,b)`, positions are constructed directly as

```text
x(g,a) = g1*C1 + g2*C2 + g3*C3 + Bas(a)
r      = x(g,a) - x(h,b)
K      = 3*r*r^T/|r|^5 - I/|r|^3
```

The evaluator loops over every ordered finite target/source pair exactly once,
skips only the same cell and basis, and rejects any other zero displacement.
The tensor is dimensionless. For moments in `mu_B`, the independent unit layer
uses

```text
field_prefactor = 1e-7 * 9.2740100783e-24 / alat_m^3
B_T             = field_prefactor * B0
E_J             = -1/2 * mu_B * field_prefactor * sum(m dot B0)
```

`MRY_J=2.179872325e-21` is used only for the reported energy-per-atom value.
The sign agrees with the legacy tensor algebra, but the legacy routine is not
used to calculate any expected value.

## Golden cases

`tests/dipole_validation/open_fft_goldens.json` contains frozen results for:

- `single_1x1x1`: exact self exclusion, zero field and energy;
- axial and transverse `2x1x1`: longitudinal factor/sign and transverse sign;
- nonuniform `2x3x1`: non-cubic indexing and permutation sensitivity;
- `2x1x3`: `G3>G2`, exposing dimension-order mistakes;
- a skew `2x2x1` cell: full Cartesian primitive-vector metric;
- `NA=2` with unequal offsets and nonuniform basis moments: phase/channel order;
- `M=4`: ensemble batch independence;
- finite `4x3x1` films magnetized in-plane and normal: shape response.

The cases intentionally overlap in no essential observable: each adds a
distinct geometry, channel, batch, or magnetostatic-shape constraint.

## Checks run

```text
python3 -m py_compile tests/dipole_validation/open_fft_oracle.py
python3 tests/dipole_validation/open_fft_oracle.py \
  --write-goldens tests/dipole_validation/open_fft_goldens.json
python3 tests/dipole_validation/test_open_fft_oracle.py -v
```

All host tests pass. The GPU red executable is not yet compiled because the
test seam is not present. No expected value in the oracle or goldens came from
CPU output, GPU output, `UppASD`, periodic Ewald code, `do_dip=3`, or the future
open builder.
