# UppASD CPU Pair-Hamiltonian Optimization and Backend Campaign

## Repository
UppASD

## Working branch
`gpu_hip_cu_ab_cg`

## Primary production source
`source/Hamiltonian/hamiltonianactions.f90` and directly related current Hamiltonian data/setup modules.

## Mandatory source-scope rule

`ApplyHamiltonian` / `applyhamiltonian.f90` and `heisge()` are obsolete, historical implementations.

They are **not production code for this campaign**.

Do not:
- optimize them;
- reconnect them;
- use them as the production oracle;
- copy their implementation into new kernels;
- use them to justify production semantics.

`HamiltonianActions` is the production CPU Hamiltonian implementation and must remain the authority for current CPU Hamiltonian semantics.

## 1. Purpose

Improve CPU performance and OpenMP scaling of the production UppASD Hamiltonian while preserving exact Hamiltonian conventions.

Investigate four execution representations:

1. **DIRECT** — current neighbour-list evaluation, optimized carefully.
2. **REDUCED-DIRECT** — translationally reduced regular-stencil representation for eligible periodic reduced Hamiltonians.
3. **SPARSE** — persistent sparse-linear-algebra representation, initially scalar `J`.
4. **CONVOLUTION** — persistent lattice-FFT convolution for eligible periodic reduced systems, initially scalar `J`.

The objective is **not** to force one backend to win everywhere.

Expected qualitative regimes:
- short-range / irregular → DIRECT likely strongest;
- periodic reduced / short-medium range → REDUCED-DIRECT may be strongest;
- irregular medium/long range → SPARSE may be competitive;
- periodic reduced / very long range → CONVOLUTION may be strongest.

Only measurements may determine the actual crossover.

## 2. Physics authority

`HamiltonianActions` owns the physics semantics.

Execution backends are alternative ways of applying the **same Hamiltonian**.

Backends must not independently redefine:
- exchange signs;
- DMI conventions;
- tensor conventions;
- moment factors;
- pair counting;
- energy factors.

## 3. Field and energy contract

Preserve the energy capability of `effective_field()` and related production field routines.

Long-term architectural objective:

> one canonical implementation of each Hamiltonian term → field → energy derived from the same term implementation

Do not create a separate independent global-energy implementation merely to make the field kernel faster.

However, energy calculation must be optional at execution time.

Conceptually:

```fortran
call effective_field(..., measure_energy=.false.)
! canonical field calculation only

call effective_field(..., measure_energy=.true.)
! same canonical field calculation + corresponding energy extraction
```

For ordinary LLG predictor/corrector steps that do not request energy, no unnecessary energy reduction should execute.

For bilinear reciprocal pair interactions,

\[
\mathbf B_i^{\rm pair}=\sum_j \mathsf K_{ij}\mathbf m_j
\]

and the canonical global pair energy may be obtained as

\[
E_{\rm pair}=-\frac12\sum_i\mathbf m_i\cdot\mathbf B_i^{\rm pair}
\]

provided the interaction term obeys the corresponding Hamiltonian reciprocity contract.

This relation covers the eligible bilinear pair operator as a whole, including exchange/DMI/tensor components under their proper conventions.

Do **not** blindly apply the factor `1/2` to:
- Zeeman;
- arbitrary onsite terms;
- higher-order interactions;
- nonlinear/multispin terms.

Those terms retain term-specific canonical energy expressions derived from the same `HamiltonianActions` implementation.

## 4. MC exception

Single-spin-flip Monte Carlo and similar local-update algorithms may retain specialized local `ΔE` kernels when recomputing the global Hamiltonian would be computationally inappropriate.

Such special kernels must be validated against the canonical global Hamiltonian.

Required conceptual oracle:
1. calculate canonical total energy `E_before`;
2. perform one test spin change;
3. calculate canonical total energy `E_after`;
4. compare `E_after-E_before` against specialized MC `ΔE`.

This campaign does not require rewriting MC unless an actual convention defect is discovered.

## 5. DIRECT backend principles

The production direct kernel is an owner-computes gather formulation:
- each target atom `i`;
- traverses neighbours `j`;
- reads `m_j`;
- accumulates `B_i`.

This has a major CPU advantage: each OpenMP worker writes only its own target field.

Do not replace this casually with a unique-pair scatter algorithm requiring atomics, locks, colouring, or large thread-private field buffers unless measured evidence demonstrates a net benefit.

For long-range Nd, use the working hypothesis:

> memory/cache/gather limited

rather than:

> primarily OpenMP scheduling limited

until hardware-counter evidence proves otherwise.

## 6. Locality policy

Space-filling/locality work is allowed, but begin with non-invasive forms.

- **Level 1, allowed early:** change target traversal order only. Physical atom indices and canonical arrays remain unchanged.
- **Level 2, later if justified:** reorder complete neighbour records consistently within each target.
- **Level 3, later:** computationally permuted/mirrored representation.
- **Level 4:** globally renumber UppASD atoms — **out of scope**.

No global atom renumbering is allowed.

Low-risk SFC experiments must preserve:
- physical atom identity;
- neighbour-list semantics;
- input/output ordering;
- restart semantics;
- chemistry/disorder mappings;
- all current functionality.

GPU SFC/locality changes are deferred until CPU evidence establishes that spatial ordering is beneficial.

## 7. Reduced-stencil principle

For an eligible reduced periodic Hamiltonian, replace the hot-loop view

> target atom → arbitrary neighbour atom index

with a translational stencil such as

> source basis, target basis, lattice-cell displacement, coupling.

Conceptually:

```text
(b_out, b_in, dR, K)
```

rather than storing/reading a full arbitrary `nlist` target for every equivalent cell.

The objective is to reduce:
- neighbour-index memory traffic;
- irregular indexing;
- metadata footprint;

and expose more regular spin-access patterns.

The reduced-stencil backend remains `O(N z)` and does **not** approximate or truncate the Hamiltonian.

Eligibility must be established explicitly.

## 8. Sparse principle

The historical `do_sparse` machinery must be treated as archaeology until current production connectivity and correctness are demonstrated.

Do not assume historical sparse support is production-ready.

Initial sparse target: scalar isotropic `J`.

Represent

\[
B_x=JM_x,\quad B_y=JM_y,\quad B_z=JM_z.
\]

Prefer one persistent sparse matrix applied to three RHS vectors where the library supports this efficiently.

Do not initially build a redundant `3N × 3N` scalar-exchange matrix.

For later `J+D`/tensor support, BSR or equivalent block-operator representation may become appropriate.

Sparse handles/inspector state must be persistent across timesteps. No sparse structure construction is allowed in the hot timestep.

## 9. Convolution principle

Initial CPU lattice-convolution scope:
- fully periodic;
- regular replicated lattice;
- reduced Hamiltonian;
- translationally invariant pair interactions;
- scalar isotropic `J`.

Do not broaden scope before scalar `J` is proven useful.

For scalar `J` and `N_A` basis atoms, use approximately `N_A^2` spectral `J_ab(q)` kernels, not a generic `9N_A^2` tensor representation.

The three Cartesian spin components reuse the same scalar `J_ab(q)`.

FFT plans/descriptors, spectral interaction kernels, mappings and work buffers must be persistent.

No plan creation, allocation, kernel transform or descriptor setup inside a production timestep.

Initial CPU FFT provider may use the project's existing MKL infrastructure.

The Hamiltonian design should avoid spreading raw MKL DFTI calls throughout `HamiltonianActions` so an FFTW provider remains possible later.

## 10. Oracle policy

The current canonical DIRECT `HamiltonianActions` implementation is the primary physics oracle before optimization.

Every new representation must compare against it on identical states.

Required field metric:

```text
max_abs_field_difference
```

and appropriate relative metrics.

For eligible bilinear pair terms compare

\[
E_{\rm pair}=-\frac12\sum \mathbf m\cdot\mathbf B_{\rm pair}
\]

across backends.

Fixtures should include:
- random/non-collinear states;
- multi-basis periodic cells;
- long-range scalar `J`;
- multiple ensembles where supported;
- later `J+D`;
- deliberately ineligible convolution cases.

Where applicable, the validated GPU convolution implementation may serve as an independent secondary oracle.

## 11. Benchmark policy

Use the general UppASD production benchmark framework for final decisions.

Representative workloads:
- bcc Fe — medium interaction range;
- 2D skyrmion — short `J+D`;
- 3D skyrmion/chiral model — short `J+D`;
- dhcp Nd — very long-range interactions.

Do not optimize only for Nd.

For low-level CPU profiling collect hardware counters where available:
- memory bandwidth;
- IPC;
- LLC misses;
- cache-miss/stall metrics;
- vectorization information.

The objective is to distinguish:
- OpenMP/fork-join limitation;
- load imbalance;
- memory-bandwidth saturation;
- cache/TLB limitation;
- instruction/vectorization limitation.

## 12. Task sequence

1. `CPU-HAM-00` — Production CPU profiling and call-path baseline
2. `CPU-HAM-01A` — Canonical field/energy execution contract
3. `CPU-HAM-01B` — Lean scalar-J DIRECT kernel
4. `CPU-HAM-02A` — Low-risk target-order/SFC locality experiment
5. `CPU-HAM-02B` — Reduced-stencil representation and oracle
6. `CPU-HAM-02C` — Optimized reduced-stencil DIRECT kernel
7. `CPU-HAM-03A` — Sparse backend archaeology and design decision
8. `CPU-HAM-03B` — Persistent scalar-J sparse backend
9. `CPU-HAM-04A` — CPU lattice-convolution design + spectral oracle
10. `CPU-HAM-04B` — Persistent scalar-J CPU convolution backend
11. `CPU-HAM-05` — Backend crossover campaign
12. `CPU-HAM-06` — Extend worthwhile global backends to `J+D`
13. `CPU-HAM-07` — Backend policy and optional AUTO selection

## 13. Parallel execution

Do not run implementation tasks that modify `HamiltonianActions` concurrently.

Sequential execution is preferred.

`CPU-HAM-03A` is mainly archaeology and could theoretically overlap with locality work in a separate worktree, but sequential execution is safer.

## 14. Model policy

Tasks are sliced for Luna.

If a task encounters:
- ambiguous physics convention;
- invasive public API redesign;
- broad FFT abstraction affecting existing dipole functionality;
- sparse-library ABI/build-system redesign;
- need to change canonical atom numbering;
- unexplained DIRECT/backend parity failure;

**STOP and report the issue.**

Do not resolve such ambiguity by broad refactoring. Escalate that subproblem for higher-capability review.

## 15. General completion rules

For every task:
- keep changes focused;
- avoid unrelated cleanup;
- preserve current code style;
- run relevant correctness tests;
- add discriminating negative controls where requested;
- benchmark before/after when performance is part of acceptance;
- record raw measurements;
- tick checklist items only after evidence exists;
- make one focused commit with the supplied one-line commit message.

At completion report:
1. files changed;
2. implementation/design result;
3. correctness tests;
4. performance evidence;
5. negative controls;
6. remaining risks;
7. whether the task is honestly closable.

## 16. CPU-HAM-07 reconciliation

CPU-HAM-05/06 evidence retains DIRECT as the safe general backend, SPARSE as
an explicit portable scalar-J option without a robust measured crossover, and
CONVOLUTION for fully periodic reduced scalar-J/J+D systems where setup and
steady-state measurements justify it. REDUCED-DIRECT remains internal to the
DIRECT family and experimental; it is not exposed as another production
choice.

CPU-HAM-07 adds the explicit `cpu_ham_backend` selection (`direct`, `sparse`,
or `convolution`), startup requested/resolved/reason diagnostics, capability
rejection for unsupported explicit requests, and compatibility handling for
the historical `do_sparse`/`do_convolution` flags. AUTO remains disabled
because the measured crossover is machine/provider dependent and no robust
portable rule was established. Unsupported physics is never silently
approximated or redirected by the production selector.

The CPU policy was validated locally with FFTW and without MKL. An MKL run is
not required for this task; provider-specific MKL benchmarking remains open.
