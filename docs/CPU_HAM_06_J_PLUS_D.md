# CPU-HAM-06 — CPU exchange plus DMI backends

**Date:** 2026-09-02
**Source revision:** `061f078` plus the CPU-HAM-06 changes
**Status:** complete for the retained REDUCED-DIRECT and CONVOLUTION backends.

## Decision

CPU-HAM-05 did not show a robust SPARSE crossover: SPARSE lost on the long-
range Nd workload and its isolated bcc-Fe result had a fit residual larger than
the fitted slope. SPARSE therefore remains scalar-J-only and was not extended
to DMI. REDUCED-DIRECT remains an explicit experimental backend, and
CONVOLUTION remains the long-range periodic backend; both were extended to
the canonical scalar `J+D` operator.

## Frozen production convention

`HamiltonianActions:dzyaloshinskii_moriya_field` defines the active CPU
convention:

\[
  \mathbf B_i^D = \sum_{j\in\mathcal N_D(i)}
      \mathbf D_{ij}\times\mathbf M_j,
  \qquad \mathbf D_{ji}=-\mathbf D_{ij}.
\]

For the reciprocal directed list, the pair energy is obtained from the
complete pair field as

\[
  E_{J+D}=-\frac{\mu_B}{2\,mry}\sum_i
       \mathbf M_i\cdot\mathbf B_i^{J+D}.
\]

The field has no `mub/mry` conversion; the conversion is applied only by the
canonical global energy reduction. This agrees with the closed DMI evidence in
[`RCG-02_DMI_HANDEDNESS_EVIDENCE.md`](RCG-02_DMI_HANDEDNESS_EVIDENCE.md),
including its analytic dimer and sign controls.

## Implementation

- `ReducedStencil` now stores independent compact exchange and oriented DMI
  records. Independent lists are required because the production J and DMI
  neighbor counts need not match.
- Reduced setup validates translated DMI source indices, vectors, and the
  reciprocal/antisymmetric relation before exposing the optimized target path.
- The reduced target kernel evaluates the canonical component form
  `(Dy*Mz-Dz*My, Dz*Mx-Dx*Mz, Dx*My-Dy*Mx)` without atomics or scatter writes.
- CPU convolution persists four spectral basis-pair families: `J(q)`,
  `Dx(q)`, `Dy(q)`, and `Dz(q)`. The same transforms are reused for all spin
  components, and the DMI contribution is evaluated as `D(q) × M(q)`.
- `HamiltonianActions` treats the convolution result as the complete bilinear
  pair field, so the existing field-derived energy path remains canonical.
- Reduced setup was moved after DMI mounting so a requested reduced J+D
  representation is validated from the complete production data.

## Correctness evidence

`tests/hamiltonian/test_cpu_ham06_j_plus_d.f90` exercises:

- a periodic 2D radial Néel skyrmion-like fixture (`n3=1`) and a periodic 3D
  two-basis chiral fixture;
- canonical DIRECT versus REDUCED-DIRECT field and energy parity;
- direct CPU convolution versus the same stencil and production convolution
  dispatch;
- field-derived pair-energy parity;
- global DMI sign reversal, which fails field parity;
- component-transposed DMI operator, which fails field parity;
- an incorrectly antisymmetric reciprocal list, which is rejected during
  stencil construction.

Local FFTW results:

```text
2D reduced field max_abs = 0.0000e+00
2D convolution field max_abs = 5.5511e-16
3D reduced field max_abs = 0.0000e+00
3D convolution field max_abs = 3.3307e-16
CPU-HAM-06 J+D tests passed
```

The existing scalar-J field/energy, sparse, reduced-stencil, and convolution
tests also passed as a focused 5/5 CTest subset. The reduced-stencil test was
also built successfully in a CPU-only `USE_FFTW=OFF` configuration.

The B02 production template is intentionally not an eligible convolution
fixture because it has an open third boundary (`BC P P 0`). The synthetic 2D
periodic fixture supplies the required discriminating 2D check without
changing that admitted workload. The B03 genuinely 3D short-range J+D
template was run through all three production modes; each mode reported the
expected activation/fallback behavior.

## Performance evidence

One local Release production run used B03 (`4096` atoms, six exchange and six
DMI directed neighbors per atom), `Nstep=1000`, and FFTW with one provider
thread. Reported wall time for the measurement phase was:

| Backend | `OMP_NUM_THREADS=1` | `OMP_NUM_THREADS=8` |
|---|---:|---:|
| DIRECT | 0.43 s | 0.44 s |
| REDUCED-DIRECT | 0.59 s | 0.49 s |
| CONVOLUTION | 0.43 s | 0.57 s |

These are single local samples and are not a promotion claim. On this
short-range workload the reduced and FFT representations do not beat DIRECT
robustly. The earlier CPU-HAM-05 Nd measurements remain the evidence for
CONVOLUTION's long-range crossover. No hardware-counter data was available.

No MKL test run was needed: this change uses the already available local FFTW
provider and does not add an MKL-specific implementation. A separate MKL
machine is only needed for a future provider-specific comparison.

## Tensor recommendation

Do not implement full general tensor exchange in this slice. If later
measurements justify it, extend the same persistent operator boundary with a
validated BSR/block representation; retain the scalar `J+D` path as the
specialized low-overhead case.

## Commit

`CPU-HAM-06: extend CPU pair backends to DMI`
