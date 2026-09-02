# CPU-HAM-08 — Policy and correctness closeout

**Status:** complete
**Date:** 2026-09-02
**Commit:** `CPU-HAM-08: close CPU Hamiltonian backend policy gaps`

## Audit of the live implementation

The production authority is `source/Hamiltonian/hamiltonianactions.f90`.
`source/Hamiltonian/hamiltonianinit.f90` builds a reduced scalar-J/DMI stencil
for eligible periodic `do_reduced=Y` systems. Before this closeout,
`effective_field` consumed that allocated stencil from the ordinary DIRECT
atom path. Sparse setup also required `size(nlistsize) >= Natom`, which
excluded valid reduced multi-basis maps even though the canonical DIRECT loop
uses `aHam(i)` to select the basis row. Target order was always allocated, so
the old full-range loop could use an indirection even for natural order.

The audit also found that SPARSE and CONVOLUTION only report themselves
applicable for the full atom/ensemble range. Partial requests then entered the
per-atom path without a documented policy. The CPU-HAM-05 driver computed its
`Spin-steps/s` field as `1 / steady_step_seconds`; the shared generic
`benchmarks/analysis/throughput.py` helper already used `Natom / t_step`.

## Implemented policy

### DIRECT and reduced data

`backend=direct` now always uses the canonical neighbour-list implementation.
`do_reduced=Y` may still build and retain the reduced representation for
CONVOLUTION and correctness fixtures, but it cannot promote ordinary DIRECT
execution. REDUCED-DIRECT is available only through the internal
`set_reduced_direct_testing` test hook and is not exposed by the production
backend selector or CPU-HAM-05 campaign driver.

The reduced-stencil tests and the J+D tests continue to compare the internal
oracle against canonical DIRECT. A production B04 dhcp-Nd run with the
template's `do_reduced=Y` printed both:

```text
Validated reduced scalar-J stencil available
CPU Hamiltonian backend: requested=direct resolved=direct ...
```

and completed successfully. The 16³, eight-thread fit from
`build/bin/sd.f95.cuda` was `0.08333038 s/step` in this local run. This
re-check demonstrates that the reduced representation is present without
changing the selected DIRECT execution.

### SPARSE multi-basis

The persistent CSR builder now validates physical-atom coverage, basis-row
maps, coupling dimensions, neighbour counts, and physical neighbour indices.
It accepts reduced multi-basis scalar-J data when those checks pass. The
capability is promoted to `SUPPORTED` in the HAM-07 table.

`tests/hamiltonian/test_sparse_backend.f90` adds a 12-atom, three-basis,
three-ensemble fixture with deterministic seeded random non-collinear moments
and 24 directed entries. SPARSE matches canonical DIRECT for total/internal/
external fields and field-derived energy to `1e-13`.

### Partial-range semantics

Policy 2 is adopted: SPARSE and CONVOLUTION remain full-range operators. An
explicit partial-range request intentionally falls back to canonical DIRECT,
emits a one-time diagnostic identifying the requested range and backend, and
does not alter field or energy results. The sparse and production convolution
tests exercise this behavior and compare the selected range against DIRECT.

### Natural target order

The full-range ordered path is now enabled only when `target_order_sfc` is
true. Natural order uses the simple collapsed atom/ensemble loop; SFC and
weighted traversal remain available. The target-order test benchmarks 40
repetitions of a 4096-atom, four-ensemble scalar-J fixture. This build
measured:

```text
native = 9.0000e-03 s   ordered = 9.0000e-03 s   ratio = 1.00
```

The four-ensemble traversal is therefore checked, and no measurable natural
ordering regression remains to justify the extra indirection.

### Throughput reporting

The CPU-HAM-05 driver now reports `atom_steps_per_second` as
`Natom / Tsteady`. Its interaction field is explicitly named
`directed_interactions_processed_per_asd_step_million_per_second` and is
defined as:

```text
2 * Ndirected * Nensemble / Tsteady
```

The factor two represents the two Hamiltonian field evaluations in a complete
ASD step. Existing JSON values were relabelled/recomputed, and the Markdown
tables now use the corrected labels and values. REDUCED-DIRECT rows remain in
that historical HAM-05 artifact only; the current campaign driver runs
DIRECT, SPARSE, and CONVOLUTION production backends.

## Regression evidence

Focused CPU backend execution passed:

```text
build/bin/cpu_ham07_backend_policy_tests
build/bin/sparse_backend_tests
build/bin/reduced_stencil_tests
build/bin/cpu_ham06_j_plus_d_tests
build/bin/cpu_convolution_tests
build/bin/target_order_tests
build/bin/field_energy_contract_tests
```

This covers scalar-J DIRECT, the reduced oracle, sparse scalar-J, convolution
scalar-J, J+D convolution, 2D DMI, 3D J+D, multiple ensembles, multi-basis
sparse parity, partial-range fallback, target traversal, and explicit backend
selection/rejection behavior.

The complete repository CTest invocation was also attempted. The CPU-HAM
backend tests passed; four unrelated pre-existing integration/reference tests
failed in the shared parallel run (`regression-test`, `asd-tests`,
`cuda-tests`, and `adaptive-cg-production-e2e`). Their failures were missing
or nondeterministic measurement/coarse-dipole outputs, not CPU-HAM-08 failures.

## Reconciled capability table

| Backend | scalar J | DMI | tensor exchange | periodic | nonperiodic | reduced Hamiltonian | non-reduced Hamiltonian | disorder | multi-basis | multiple ensembles |
|---|---|---|---|---|---|---|---|---|---|---|
| DIRECT | SUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED |
| SPARSE | SUPPORTED | UNSUPPORTED | UNSUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED | UNSUPPORTED | SUPPORTED | SUPPORTED |
| CONVOLUTION | SUPPORTED | SUPPORTED | UNSUPPORTED | SUPPORTED* | UNSUPPORTED | SUPPORTED | UNSUPPORTED | UNSUPPORTED | SUPPORTED | SUPPORTED |

`*` CONVOLUTION means a fully periodic translational reduced Hamiltonian.
There are no remaining `NOT_VALIDATED` entries in the production table.

## Checklist

- [x] Live backend routing audited.
- [x] DIRECT no longer silently means REDUCED-DIRECT.
- [x] Reduced stencil retained for legitimate uses.
- [x] Nd DIRECT regression removed/rechecked.
- [x] Sparse multi-basis fixture added.
- [x] Sparse multi-basis promoted to `SUPPORTED`.
- [x] Partial-range production callers audited.
- [x] Partial-range backend policy documented.
- [x] Partial-range tests added.
- [x] Natural target-order overhead benchmarked.
- [x] Natural-order fast path restored.
- [x] Multi-ensemble target traversal checked.
- [x] Spin-throughput formula fixed.
- [x] Interaction-throughput definition fixed.
- [x] HAM-07 capability table reconciled.
- [x] Focused CPU backend regression suite passes.
