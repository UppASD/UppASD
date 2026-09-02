# CPU-HAM-02A — Test spatial target ordering without renumbering atoms

**Model:** Luna

## Dependency

`CPU-HAM-01B` complete.

## Purpose

Determine whether spatially coherent target traversal improves cache/TLB/NUMA behaviour of the DIRECT CPU Hamiltonian.

## Non-negotiable safety rule

Do not:
- renumber physical UppASD atoms;
- permute canonical site-indexed arrays;
- rewrite input/output/restart numbering.

This task changes only the order in which independent target atoms are processed.

## A. Define target order

Build:

```fortran
target_order(1:Natom)
```

with original physical atom IDs as values.

Support NATURAL and one existing/reusable space-filling ordering.

Prefer existing repository SFC machinery. If Morton and Hilbert both already exist cheaply, both may be benchmarked, but do not broaden scope merely to add SFC methods.

## B. Production loop

Process `i = target_order(q)` rather than assuming `i=q`.

The inner neighbour summation order for each atom must remain unchanged, preserving each atom's floating-point accumulation order.

## C. Thread partition

For homogeneous workloads, contiguous ranges of `target_order` are sufficient.

For heterogeneous neighbour counts, add a weighted-static option with:

```text
work_i = nlistsize(i)
```

Partition the ordered sequence into contiguous ranges with approximately equal cumulative work per thread.

Do not use dynamic scheduling as the primary locality solution.

## D. Multi-ensemble

Preserve ensemble semantics and do not mix SFC ordering with ensemble identity.

## E. Correctness

Require very strong field parity since per-atom neighbour summation order remains unchanged.

## F. Hardware counters

Compare NATURAL vs SFC for:
- Nd;
- bcc Fe;
- short-range `J+D`.

Record runtime, memory bandwidth, LLC misses, TLB metrics if available, and IPC.

## G. GPU

Do not modify GPU code.

Record whether CPU results justify a later GPU locality experiment.

## Implementation result

The production CPU path now has a persistent `ham%target_order` permutation and
`ham%target_work_prefix` cumulative work array. The permutation values are the
original physical atom IDs. `NATURAL` is the identity permutation; `do_sfc Y`
selects a deterministic 3-D Morton order for non-random-alloy runs, reusing the
existing coordinate/SFC opt-in. No canonical arrays, restart numbering, or
neighbor records are permuted.

`effective_field` visits `i=ham%target_order(q)` for full-system calls. The
neighbor loop inside `effective_field_atom` is unchanged, so every target keeps
its original floating-point accumulation order. Partial atom-range calls retain
the natural loop because a global permutation cannot represent an arbitrary
contiguous subrange without filtering.

When target neighbor counts differ, the ordered sequence is partitioned into
contiguous cumulative-work ranges using `work_i=max(1,nlistsize(aHam(i)))`.
Homogeneous cases retain ordinary static scheduling. There is no dynamic
scheduling and no GPU source change.

## Correctness and timing evidence

The focused `cpu-ham-target-order` test validates identity NATURAL ordering,
deterministic Morton ordering, physical-ID permutation preservation, and
contiguous weighted partition cuts. The existing field/energy contract also
passes after the ordered production loop was added.

Production parity runs used identical deterministic inputs, `OMP_NUM_THREADS=1`,
and 20 steps. `restart` and `averages` were bitwise identical for bcc Fe and
short-range J+D; `restart`, `totenergy`, and `averages` were bitwise identical
for Nd. The bcc Fe and J+D runs used 16,000 and 16,384 atoms respectively; the
Nd parity run used 16,384 atoms.

The following exploratory wall times include setup and fixed initialization,
so they are evidence for the complete production run rather than isolated
kernel throughput. B01 and B02 used 200 steps; B04 used 100 steps and 16^3
cells (16,384 atoms). Each entry is one run on the local macOS host.

| Workload | Order | 1 thread | 2 threads | 4 threads | 8 threads |
|---|---|---:|---:|---:|---:|
| bcc Fe, B01, 16,000 atoms | NATURAL | 1.56 s | 0.91 s | 0.78 s | 0.79 s |
| bcc Fe, B01, 16,000 atoms | MORTON | 1.79 s | 0.97 s | 0.68 s | 0.61 s |
| short-range J+D, B02, 16,384 atoms | NATURAL | 8.71 s | 6.23 s | 5.78 s | 5.84 s |
| short-range J+D, B02, 16,384 atoms | MORTON | 9.32 s | 6.37 s | 5.98 s | 6.69 s |
| dhcp Nd, B04, 16,384 atoms | NATURAL | 6.45 s | 3.50 s | 2.17 s | 1.81 s |
| dhcp Nd, B04, 16,384 atoms | MORTON | 11.34 s | 6.05 s | 3.46 s | 2.78 s |

The host provided no usable `perf` or `likwid-perfctr` interface; therefore
memory bandwidth, LLC misses, TLB metrics, and IPC are recorded as unavailable,
not inferred from wall time. `powermetrics` exists but requires privileged
sampling and was not used as a substitute for repeatable hardware counters.

Decision: `WORKLOAD_DEPENDENT`. Morton improves the bcc Fe 4/8-thread points,
but loses at 1/2 threads there, loses across all tested J+D points, and loses
across all tested Nd points. SFC remains opt-in and is not enabled globally.
The CPU evidence does not yet justify a GPU locality experiment; GPU code is
untouched.

## H. Decision

Classify result as:
- `CLEAR_WIN`
- `WORKLOAD_DEPENDENT`
- `NEUTRAL`
- `REGRESSION`

Do not enable globally merely because Nd improves.

## Checklist

- [x] Physical atom numbering unchanged.
- [x] `target_order` representation implemented.
- [x] NATURAL order supported.
- [x] One SFC order supported (deterministic Morton).
- [x] Inner neighbour order unchanged.
- [x] Weighted static partitioning available if needed.
- [x] Nd parity passes.
- [x] Fe parity passes.
- [x] `J+D` parity passes.
- [x] Nd timing measured.
- [x] Fe timing measured.
- [x] `J+D` timing measured.
- [x] Cache/bandwidth counters compared (counter interfaces unavailable on this host; limitation recorded).
- [x] GPU code untouched.
- [x] Evidence-based SFC conclusion recorded (`WORKLOAD_DEPENDENT`; SFC remains opt-in).

## Commit

`CPU-HAM-02A: evaluate spatial target ordering`
