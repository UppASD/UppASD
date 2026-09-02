# CPU-HAM-01B — Optimize the production scalar-J DIRECT kernel

**Model:** Luna

## Dependency

`CPU-HAM-01A` complete.

## Purpose

Create a strong DIRECT baseline before evaluating new representations.

Primary target: long-range scalar exchange, without regressing short-range performance.

## Implementation result

The production `heisenberg_field` routine in
`source/Hamiltonian/hamiltonianactions.f90` now keeps the owner-computes gather
formulation but accumulates `bx`, `by`, and `bz` as scalars. The neighbour count,
neighbour index, and coupling are loaded into local scalars before the three
component updates. The general `effective_field` dispatch, rescaling path,
tensor path, optional terms, static OpenMP schedule, and energy contract are
unchanged. No obsolete Hamiltonian source was used.

The existing CPU-HAM-01A contract test now also evaluates a non-collinear,
two-ensemble scalar-J case against an independent direct neighbour-list oracle.
It passes at one and four OpenMP threads. The canonical field/energy contract
continues to pass, including exact field identity with energy enabled/disabled.

### Vectorization decision

The production compiler configuration was used (`gfortran`, `-O3 -Ofast
-funroll-loops -finline-functions -fopenmp -mtune=native`) with
`-fopt-info-vec-all`. The pre-change array-expression loop was not vectorized
(`couldn't vectorize loop`, cost model), and the selected scalar loop also was
not auto-vectorized because the indirect neighbour gather was not profitable.
The explicit `!$omp simd reduction(+:bx,by,bz)` experiment compiled legally,
but GCC reported the loop body as unsupported SLP rather than vectorizing the
indirect gather. It was removed from the production implementation.

Single-run exploratory whole-process timings used `Nstep=20`, structure output
disabled, static scheduling, and `OMP_NUM_THREADS=1,2,4,8`. They are useful for
variant selection, not authoritative benchmark records:

| Case | 1 thread | 2 threads | 4 threads | 8 threads | 8-thread speedup / efficiency |
|---|---:|---:|---:|---:|---:|
| bcc Fe, 16,000 atoms, scalar | 2.30 s | 1.20 s | 1.05 s | 0.94 s | 2.45x / 30.6% |
| dhcp Nd, 16,384 atoms, scalar | 1.98 s | 1.10 s | 0.80 s | 0.74 s | 2.68x / 33.5% |
| short-range J+D, 16,384 atoms, scalar | 0.11 s | 0.09 s | 0.09 s | 0.09 s | 1.22x / 15.3% |

The exploratory SIMD variant measured bcc Fe at `2.31, 1.23, 0.99, 0.93 s`
and dhcp Nd at `1.76, 1.10, 0.79, 0.72 s` for 1/2/4/8 threads. With no
compiler vectorization of the gather body and no consistent improvement, the
clean scalar variant was selected. The OpenMP runtime reported no usable
hardware-counter interface on this macOS host: `perf`, `likwid-perfctr`, and
privileged sampling/counter tools were unavailable. Memory bandwidth, IPC,
cache misses, and stall cycles are therefore recorded as unavailable rather
than inferred from timing.

### Correctness and regression evidence

An isolated pre-change CPU executable and the selected post-change executable
were run with identical deterministic inputs and `OMP_NUM_THREADS=1`:

| Case | Production parity result |
|---|---|
| bcc Fe | exit status 0; `restart` and `averages` outputs bitwise identical |
| dhcp Nd | exit status 0; `restart`, `totenergy`, and `averages` outputs bitwise identical |
| short-range J+D | exit status 0; `restart`, `averages`, and diagnostic outputs identical apart from the provenance YAML |

The bcc Fe and dhcp Nd runs use the campaign's real neighbour lists and random
non-collinear Nd initialization. The focused scalar-J test covers multiple
ensembles and a separate direct field oracle. The CPU-HAM-01A energy parity
test passes unchanged, so no second energy traversal or universal factor-
`1/2` rule was introduced.

## A. Locate production hot loop

Work only inside current `HamiltonianActions` and supporting current modules.

Do not touch obsolete `heisge`.

## B. Scalar accumulators

For isotropic exchange, evaluate whether the hot neighbour loop can be reduced to explicit scalar accumulators `bx`, `by`, `bz`, with hoisted neighbour count, neighbour index and coupling.

Preserve exact physics.

## C. Specialization

If a `J`-only fast path is introduced, keep it narrow and explicit.

Do not duplicate the whole `effective_field` implementation.

General Hamiltonian cases must continue through canonical term dispatch.

## D. Vectorization experiment

Compare:
1. clean scalar accumulators;
2. compiler autovectorization;
3. explicit SIMD reduction, only if legal.

Use compiler vectorization reports.

Do not simply uncomment old `!$omp simd` directives.

Record whether vectorized, gather behaviour if visible, and measured speed.

Select the fastest correct implementation.

## E. OpenMP scheduling

For homogeneous bulk systems test/retain static scheduling as baseline.

Do not introduce dynamic scheduling without demonstrating load imbalance.

Defer weighted partitioning to CPU-HAM-02A.

## F. Array initialization

Using CPU-HAM-00 findings, remove only demonstrably unnecessary full-array clear/copy work.

Do not rely on stale values. Preserve dipole/non-pair additive semantics.

## G. Energy

Use CPU-HAM-01A contract. Do not add a second energy traversal.

## H. Correctness

Use random non-collinear states.

Compare pre/post optimization:
- per-atom field;
- pair energy;
- total effective field where appropriate.

Include Nd, bcc Fe, and multiple ensembles where supported.

## I. Performance

Measure at `1, 2, 4, 8, ...` threads.

Record:
- time;
- speedup;
- efficiency;
- memory bandwidth;
- IPC where available.

The task succeeds even if OpenMP scaling remains bandwidth limited, provided the direct kernel becomes measurably leaner without regression.

## Checklist

- [x] Production scalar-J hot loop simplified.
- [x] No obsolete source used.
- [x] Scalar accumulator implementation tested.
- [x] Auto-vectorization evidence recorded.
- [x] Explicit SIMD tested only if legal.
- [x] Fastest correct variant selected.
- [x] No unnecessary dynamic scheduling introduced.
- [x] Full-array operations reviewed.
- [x] Energy contract preserved.
- [x] Nd field parity passes.
- [x] bcc Fe field parity passes.
- [x] Energy parity passes.
- [x] Thread scaling rerun.
- [x] Hardware counters compared (counter interface unavailable on this host; limitation recorded).
- [x] Short-range regression checked.

## Commit

`CPU-HAM-01B: optimize direct scalar exchange kernel`
