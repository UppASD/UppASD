## Verdict

WP4b acceptance: **PASS for the accepted CUDA fp64 `NA=1` scope.**

Production approval: **FAIL — WP5 production wiring is not complete.**

Terra reran the CUDA suite and compute-sanitizer on commit
`ab4c5518a7b22078f69900ba028073e280710468`; Luna reran the independent Python
oracle self-check on the same commit.

| Gate | Result | Evidence |
|---|---|---|
| Python independent oracle | PASS | `test_periodic_ewald_reference.py` passes |
| Host Builder A vs oracle | PASS | CTest passes; cubic, non-cubic, skew golden fields pass |
| Single `1/Ngrid` normalization | PASS by code audit | Sole scaling kernel at [gpuDipoleConvolution.cpp:529](</home/andersb/SD/UppASD_gpu_hip_cu/source/gpu_files/gpuDipoleConvolution.cpp:529>); raw C2R at line 505 has no later scale |
| Displacement sign / n1-fast indexing | PASS | Delta maximum field error `2.78e-17` |
| Basis/source ordering | PASS for WP4a plumbing | Basis-2 impulse matrix passes with zero error |
| CUDA delta suite | PASS | All required grid/basis/M=1/M=4 cases pass |
| CUDA periodic convolution | PASS | Maximum field/energy errors `1.78e-14` / `1.14e-13` |
| GPU vs independent oracle | PASS for NA=1 matrix | Maximum field/energy errors `3.89e-16` / `1.81e-16` |
| Alpha invariance | PASS | Field/energy errors `1.08e-15` / `6.11e-16` |
| Translation/sign-flip/zero/changed-moment/M=4 | PASS | All reported errors are `0` |
| Reciprocity/Hermitian residual | PASS | Maximum `4.16e-17` / `2.40e-17` |
| `1x1x1`, q=0 nonzero/correct | PASS host/oracle | Tensor diagonal ≈ `0.155140377955052`; representative-only reciprocal term is zero while alias sum is nonzero |
| Compute-sanitizer | PASS | `ERROR SUMMARY: 0 errors` |
| Production field/energy wiring | **NOT APPROVED / WP5** | Intentional guard at [gpuHamiltonianCalculations.cpp:578](</home/andersb/SD/UppASD_gpu_hip_cu/source/gpu_files/gpuHamiltonianCalculations.cpp:578>) |

CUDA standalone numerical output was:

- `1x1x1`: field/energy error `0 / 0`
- `2x1x1`: `6.94e-18 / 1.73e-18`
- `2x3x1`: `2.78e-17 / 5.55e-17`
- skew: `3.47e-18 / 0`
- two-basis, M=4: `1.78e-14 / 1.14e-13`

The independent-oracle matrix additionally reported maximum field/energy errors
of `3.8857805861880479e-16` / `1.8041124150158794e-16`.

## Prioritized defects

1. **WP5 — Production backend remains intentionally disabled.**  
   `GpuHamiltonianCalculations::initiate()` throws:

   > “GPU dipole mode was requested, but its field/energy operator is not yet available”

   The operator is never invoked from `heisge()`.

2. **P0 — No production field or energy integration.**  
   The Tesla prefactor is only staged at [gpuDipoleConvolution.cpp:324](</home/andersb/SD/UppASD_gpu_hip_cu/source/gpu_files/gpuDipoleConvolution.cpp:324>). It is not applied to fields, and no dipole contribution is added to `beff`, `eneff`, `Dip`, or total energy.

3. **P2 — Independent `NA>1` physics remains deferred to WP6.**  
   The basis-two tests in WP4b are FFT-plumbing tests; the plan explicitly
   assigns independent two-basis physics oracle coverage to WP6.

4. **P2 — Geometry validation is incomplete in the descriptor.**  
   Macrocell centres are staged but not checked against the expected regular target/source basis layout before construction.

The CUDA convolution, independent oracle, and sanitizer results satisfy the
WP4b standalone acceptance criteria. They do not approve production wiring.
