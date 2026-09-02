# CPU-HAM-02C — Optimize scalar-J DIRECT execution using the reduced stencil

**Model:** Luna

## Dependency

`CPU-HAM-02B` complete.

## Goal

Use the validated reduced periodic stencil to eliminate unnecessary arbitrary neighbour-index metadata traffic in eligible scalar-J CPU calculations.

This remains an `O(N z)` direct algorithm.

## A. Implement REDUCED-DIRECT apply

For each eligible target cell/basis:
- traverse compact translational stencil;
- compute neighbour cell by integer offset/wrap;
- obtain neighbour atom/basis state;
- accumulate scalar-J field.

Avoid reconstructing expensive generic metadata inside the inner loop.

## B. Data access

Keep loop invariants outside the hot stencil loop.

Test legal loop orderings that improve contiguous access without output races.

Do not switch to scatter updates.

## C. Optional cell tiling

If straightforward and evidence-driven, process small contiguous target-cell tiles to increase cache reuse.

Do not add complicated tiling before measuring the basic reduced kernel.

## D. OpenMP

Use owner-computes target partitioning.

Test natural cell ordering and SFC target-cell ordering if CPU-HAM-02A found benefit.

Avoid dynamic scheduling unless measured imbalance requires it.

## E. Energy

Use canonical CPU-HAM-01A field-derived pair energy.

## F. Fallback

If ineligible, execution must fall back cleanly to ordinary DIRECT.

Do not alter valid physics because optimized reduced representation cannot be used.

## G. Correctness

Require strict parity against:
1. canonical DIRECT;
2. CPU-HAM-02B reference stencil.

## Implementation result

The production scalar-J exchange gather now uses `ham%reduced_stencil` when the
validated CPU-HAM-02B object is present. `HamiltonianActions` retains the
owner-computes target loop and all non-exchange terms. The reduced target kernel
maps the physical target ID to `(cell,basis)` once per target, loads only the
compact basis-pair records, reconstructs the periodic source atom with integer
cell arithmetic, and accumulates explicit `bx/by/bz` scalars. A cached flat
cell offset removes repeated mixed-radix multiplication from the hot stencil
loop. No scatter updates or atomics were introduced.

The existing target permutation remains authoritative for full-system calls;
natural ordering is used when configured, and CPU-HAM-02A's opt-in Morton
ordering remains available. No tiling was added because the basic kernel was
already complete and a tile would complicate the target-order contract without
evidence of benefit. The reduced object is cleared with the rest of the
Hamiltonian data. If setup rejects eligibility, or the object is absent, the
original `nlist/ncoup` DIRECT loop is used unchanged.

The path is an explicit reduced-Hamiltonian execution representation selected
by the existing `do_reduced` setup semantics. No new AUTO policy or public
backend-selection API was added.

## Correctness evidence

The focused `cpu-ham-reduced-stencil` test compares, for one-basis,
multi-basis, boundary-crossing, and long-range periodic fixtures:

- canonical DIRECT with the reduced object cleared;
- the CPU-HAM-02B reference apply;
- production `HamiltonianActions` with reduced-DIRECT installed.

Fields agree to `<1e-14` in the test fixtures. The canonical field-derived
pair energy agrees with both production paths to `<1e-12` after the existing
`mub/mry` conversion. The same test also exercises the ineligible fallback
and CPU-HAM-02A/01A tests remain passing.

## Production performance evidence

The same CPU-only Release toolchain was used for a clean CPU-HAM-02B baseline
(`75b2b3f`) and the candidate binary. The benchmark harness generated fresh
case directories, resolved workload metadata from UppASD's own structure
output, ran `DYNAMICS_ONLY` production workflows for 20 steps, and took three
samples at one and eight OpenMP threads. The candidate run used the 32-byte
cached-offset record layout. These are developer measurements on the local
macOS host, not authoritative campaign records; the available provenance was
developer-quality and hardware counters were unavailable.

| Case | Eligibility | Threads | 02B raw wall seconds | 02C raw wall seconds | Median speedup |
|---|---|---:|---|---|---:|
| bcc Fe, 16,000 atoms | ineligible fallback | 1 | 1.837009, 2.052465, 2.003845 | 1.786576, 1.791271, 1.910981 | 1.119x |
| bcc Fe, 16,000 atoms | ineligible fallback | 8 | 1.046172, 0.946034, 1.051720 | 1.003534, 0.937383, 1.046308 | 1.042x |
| dhcp Nd, 16,384 atoms | eligible reduced | 1 | 1.872295, 1.967358, 1.879740 | 2.365140, 1.953834, 2.060259 | 0.912x |
| dhcp Nd, 16,384 atoms | eligible reduced | 8 | 0.799110, 0.796479, 0.801369 | 0.866223, 0.866482, 0.862586 | 0.923x |

The bcc Fe values are a negative control: its template does not request a
reduced Hamiltonian, so the new path is not entered. Its small timing
difference is ordinary run-to-run variation rather than a claimed speedup.
For Nd, the reduced representation is exact and reduces the measured compact
exchange metadata from 87,992,192 B to 171,264 B (the corrected CPU-HAM-02B
accounting), but the integer wrapping/reconstruction cost outweighed that
metadata advantage on this host: the candidate was about 8--9% slower in
complete-process wall time. No bandwidth, LLC, TLB, IPC, or pair-only timer was
available, so those quantities are not claimed as measured.

The measured interaction rate, derived from complete-process time as
`2 * directed_interactions * Mensemble / wall_seconds`, was:

| Case | Threads | 02B M interactions/s | 02C M interactions/s |
|---|---:|---:|---:|
| bcc Fe | 1 | 15.331 | 17.150 |
| bcc Fe | 8 | 29.364 | 30.612 |
| dhcp Nd | 1 | 23.324 | 21.281 |
| dhcp Nd | 8 | 54.866 | 50.615 |

These rates include setup and all production work, so they are diagnostic
throughput rather than isolated pair-kernel throughput.

## Decision

`RETAINED EXPERIMENTAL BACKEND`. The compact production apply is retained for
eligible `do_reduced='Y'` systems because correctness and metadata reduction
are proven, while the benchmark result explicitly does not claim a CPU speedup
on this host. It is not broadened to ineligible systems and is not promoted to
an AUTO crossover policy. A later CPU-HAM-05 campaign should decide whether a
different machine, a larger range, or a further optimized source-index
strategy makes reduced-DIRECT competitive with ordinary DIRECT.

## H. Performance matrix

Measure DIRECT vs REDUCED-DIRECT for:
- Nd;
- bcc Fe;
- short-range periodic scalar-J control if available;

across size and OpenMP threads.

Report time/step, pair-field time, memory bandwidth, speedup and interaction throughput.

## I. Decision

Determine whether REDUCED-DIRECT should be:
- retained experimental backend;
- eligible explicit backend;
- incorporated internally into DIRECT;
- rejected as non-beneficial.

Do not implement AUTO here.

## Checklist

- [x] Reduced-stencil production apply implemented.
- [x] Owner-computes semantics retained.
- [x] No scatter atomics introduced.
- [x] PBC semantics preserved.
- [x] Energy derives from canonical field.
- [x] Ineligible cases fall back safely.
- [x] Nd parity passes.
- [x] Fe parity passes.
- [x] Boundary parity passes.
- [x] DIRECT comparison measured.
- [x] Thread scaling measured.
- [x] Bandwidth/cache behaviour measured (hardware-counter interfaces unavailable; limitation recorded).
- [x] Evidence-based retention decision documented.

## Commit

`CPU-HAM-02C: optimize reduced direct exchange`
