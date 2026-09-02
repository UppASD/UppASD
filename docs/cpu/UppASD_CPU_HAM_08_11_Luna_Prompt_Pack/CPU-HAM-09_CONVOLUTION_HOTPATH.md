# CPU-HAM-09 — Optimize the CPU Convolution Hot Path

**Model:** Luna

## Dependency

CPU-HAM-08 complete.

## Purpose

Optimize the already-correct CPU lattice-convolution backend without changing its physics or eligibility.

This is the primary performance task of the second HAM phase.

## Scope

Initial mandatory optimization scope:

- scalar J;
- existing J+D convolution;
- CPU convolution module and narrow supporting code.

Do not change:
- reduced-stencil physics;
- DMI sign convention;
- FFT normalization;
- backend eligibility;
- public Hamiltonian semantics.

## A. Establish baseline

Before changing code, record on the available machine:

- scalar-J long-range Nd convolution;
- medium-range eligible Fe case;
- one short-range case where convolution loses;
- J+D 2D/3D correctness/performance case.

Record:
- total convolution apply time;
- pack;
- forward FFT;
- spectral multiply;
- inverse FFT;
- unpack;
- thread count/provider;
- raw samples.

Use at least several repetitions for hot-path microtimings.

## B. Remove provably redundant clears

Audit work-buffer initialization.

Potential examples:
- `real_work = 0` immediately before every packed element is overwritten;
- `field_spectral = 0` immediately before every output element is assigned.

Do not remove a clear unless complete overwrite is proven for:
- all basis indices;
- all ensembles;
- all Cartesian components;
- all spectral/grid entries.

Add a test-only poison/NaN mode where useful:
- initialize the buffer to poison;
- run apply;
- verify no poison reaches output.

This proves complete overwrite.

## C. Optimize packing

Current pack/unpack work may be serial.

Investigate:
- OpenMP parallel loops;
- collapse over ensemble/basis/cell where safe;
- contiguous copy opportunities;
- avoiding repeated index reconstruction.

Do not parallelize tiny loops blindly.

Measure whether pack/unpack matters at realistic Nd sizes.

## D. Optimize unpacking

Apply the same discipline as packing.

If output layout allows direct write without an intermediate copy, evaluate it.

Do not sacrifice code clarity or introduce aliasing assumptions without evidence.

## E. Joint Cartesian spectral multiply for J+D

Current spectral work may reload the same:

`J, Dx, Dy, Dz, mx, my, mz`

for each output Cartesian component.

Restructure the innermost computation so for each output basis `a`, spectral point `q`, and input basis `b`:

1. load:
   - `J_ab(q)`
   - `Dx_ab(q)`
   - `Dy_ab(q)`
   - `Dz_ab(q)`
   - `mx_b(q)`
   - `my_b(q)`
   - `mz_b(q)`

2. update all:
   - `Bx`
   - `By`
   - `Bz`

in one pass.

Preserve the exact HAM-06 DMI cross-product convention.

This should reduce repeated spectral-memory traffic.

## F. Scalar-J specialization

For pure scalar J:
- bypass all DMI spectral arrays/branches;
- use the minimal spectral multiply;
- reuse one `J_ab(q)` for x/y/z.

Avoid runtime conditionals in the innermost spectral loop if a clean specialized path is possible.

Do not duplicate the entire convolution implementation.

Share setup/mapping/FFT lifecycle.

## G. Loop order and cache locality

Benchmark legal loop nest choices.

Candidates may include different orderings of:
- ensemble;
- output basis;
- spectral q;
- input basis.

Use measured cache/vectorization evidence.

Do not assume the current loop order is optimal.

## H. Vectorization

Use compiler reports.

Test whether the spectral multiply vectorizes over `q`.

If explicit SIMD is legal and beneficial, add it narrowly.

Do not force vectorization of loops with adverse gather/scatter behavior.

## I. FFT batching

Audit current FFT batching.

Ensure transforms are batched over appropriate:
- Cartesian components;
- basis channels;
- ensembles;

without creating excessive descriptor count.

Do not create/destroy plans per timestep.

## J. FFT threading

Benchmark provider thread counts:

`1, 2, 4, 8, ...` physical cores as appropriate.

Avoid nested oversubscription.

Explicitly record:
- outer UppASD OpenMP threads;
- FFT provider threads.

If transforms are called outside a parallel region, let the provider own those threads.

## K. FFTW versus oneMKL preparation

Do not yet make oneMKL mandatory.

If current abstraction permits both providers:
- preserve both;
- ensure HAM-10 can compare them.

If provider abstraction is incomplete, perform only the minimal cleanup necessary to permit HAM-10 provider benchmarking.

Do not redesign existing dipole FFT code.

## L. Correctness

After each substantive optimization run:
- DIRECT vs CONVOLUTION field parity;
- field-derived energy parity;
- scalar J;
- J+D;
- multi-basis;
- multiple ensembles where supported;
- 2D/3D periodic cases.

Required DMI negative controls from HAM-06 must remain discriminating.

## M. Performance acceptance

Report speedups relative to pre-HAM-09 convolution, not only relative to DIRECT.

For each optimization identify:
- gain;
- neutral;
- regression.

Do not retain complicated transformations that are neutral within timing noise.

## N. Deliverable

Create:

`docs/CPU_HAM_09_CONVOLUTION_OPTIMIZATION.md`

with:
- baseline;
- each tested optimization;
- retained/rejected changes;
- stage timing;
- overall speedup;
- correctness evidence.

## Checklist

- [x] Baseline frozen.
- [x] Redundant clears audited.
- [x] Complete-overwrite poison tests added where needed.
- [x] Proven redundant clears removed.
- [x] Packing benchmarked.
- [x] Packing parallelized only if beneficial.
- [x] Unpacking benchmarked.
- [x] Unpacking parallelized only if beneficial.
- [x] J+D Cartesian outputs updated jointly.
- [x] Scalar-J specialized spectral path implemented if beneficial.
- [x] Loop order benchmarked.
- [x] Vectorization report inspected.
- [x] FFT batching audited.
- [x] FFT thread scaling measured.
- [x] No per-step FFT planning introduced.
- [x] Scalar-J parity passes.
- [x] J+D parity passes.
- [x] Multi-basis parity passes.
- [x] Multi-ensemble parity passes.
- [x] DMI fault-injection tests remain discriminating.
- [x] Overall convolution speedup reported.
- [x] Unhelpful complexity removed/reverted.

## Commit

`CPU-HAM-09: optimize CPU convolution hot path`
