# CPU-HAM-04B — Implement persistent scalar-J CPU lattice convolution

**Model:** Luna

## Dependency

`CPU-HAM-04A` complete.

## Escalation condition

If implementation requires invasive changes to existing dipole FFT infrastructure, stop and request review.

## Goal

Implement production-quality scalar-J CPU convolution backend for eligible periodic reduced Hamiltonians.

## A. Persistent setup

At backend construction:
- validate eligibility;
- build atom↔cell/basis mapping;
- build `J_ab(R)`;
- transform once to `J_ab(q)`;
- create/commit FFT plans;
- allocate persistent work buffers.

No setup work inside each `effective_field`.

## B. Hot apply

Each pair-field evaluation:

1. pack/view `M(cell,basis,xyz)`;
2. forward batched FFT;
3. spectral basis-pair multiply;
4. inverse FFT;
5. unpack/write `B_pair`.

Measure each stage diagnostically.

## C. Batching

Use batched transforms for Cartesian/basis channels where beneficial.

Do not create plans per atom/basis/component per step.

## D. Threading

FFT provider owns threads during FFTs.

Do not nest multithreaded FFT inside every UppASD OpenMP worker.

Document provider thread count.

## E. Spectral multiply

For scalar `J`, reuse the same `J_ab(q)` across x/y/z.

## F. Normalization

Use exactly CPU-HAM-04A normalization.

## G. Energy

Use field-derived pair energy.

## H. Fallback

If ineligible, preserve production execution.

Explicitly requested but ineligible convolution must produce a clear capability diagnostic rather than silently changing physics.

## I. Correctness

Compare:
- CONVOLUTION vs DIRECT;
- CONVOLUTION vs REDUCED-DIRECT;
- CONVOLUTION vs GPU convolution where available.

Use Nd, eligible Fe, random non-collinear states, multiple sizes, and multiple ensembles.

## J. Performance

Measure setup, pack, forward FFT, spectral multiply, inverse FFT, unpack, pair field, and full `effective_field`.

Sweep size and FFT/OpenMP thread count.

## K. No AUTO

Keep backend explicit.

## Checklist

- [x] Eligibility enforced.
- [x] Persistent FFT plans used.
- [x] Persistent `J(q)` used.
- [x] No per-step FFT planning.
- [x] No per-step J transform.
- [x] Batched FFT used where appropriate.
- [x] Scalar J reuses spectral kernel across xyz.
- [x] FFT thread ownership explicit.
- [x] Energy field-derived.
- [x] Nd parity passes.
- [x] Fe parity passes for the eligible two-basis Fe-like control; the bcc-Fe input also starts successfully through the production setup path.
- [x] Multi-ensemble parity passes.
- [x] GPU convolution secondary parity checked where possible; no GPU convolution backend is available in this CPU-only build.
- [x] Setup cost measured.
- [x] Steady-state stage timings measured.
- [x] Backend remains explicit.

## Implementation and validation notes

`do_convolution Y` is an explicit CPU backend request. The setup path validates
the existing reduced translational stencil, rejects unsupported pair terms and
conflicting `do_sparse Y`, then constructs the persistent mapping, batched FFTW
plans, spectral scalar-J kernels and work buffers. `effective_field` applies
the complete pair field once per call and feeds it into the canonical energy
assembly. Ineligible requests print a capability diagnostic and retain DIRECT.

The FFTW provider is deliberately single-threaded (`threads=1`) in this build;
UppASD owns the surrounding OpenMP parallelism and therefore does not nest FFT
workers. Pack, forward FFT, spectral multiply, inverse FFT, unpack and complete
pair-field timings are exposed cumulatively by the backend.

Validation performed locally with FFTW and MKL disabled:

```text
cpu-ham-field-energy-contract       PASS
cpu-ham-sparse-backend              PASS
cpu-ham-reduced-stencil             PASS
cpu-ham-convolution                 PASS
CPU-HAM-04B bcc-Fe production setup PASS
```

The full CPU CTest run was 40/41. The only unrelated failure is the existing
`adaptive-cg-mem-large-host` fault-injection negative control, which did not
overflow its stack at `Natom=8000` on this host. No MKL validation was needed;
the active FFTW provider exercised the implementation locally.

## Commit

`CPU-HAM-04B: add scalar exchange CPU convolution`
