# CPU-HAM-08 — Policy and Correctness Closeout

**Model:** Luna

## Dependency

HAM-00 through HAM-07 complete and merged.

## Purpose

Resolve the remaining inconsistencies found by the adversarial post-HAM audit before performing further optimization.

This task should be surgical.

Do not introduce new Hamiltonian algorithms.

## A. Re-audit current production state

Before editing, verify from current source:

1. when the reduced stencil is built;
2. when DIRECT actually calls reduced-stencil apply;
3. sparse eligibility checks;
4. convolution eligibility checks;
5. partial-range `effective_field` behavior;
6. target-order setup and natural-order execution;
7. benchmark throughput calculations.

Document the observed current behavior in:

`docs/CPU_HAM_08_POLICY_CLOSEOUT.md`

Do not trust stale HAM-05/HAM-07 prose over live code.

## B. Decouple DIRECT from automatic REDUCED-DIRECT

HAM-05 concluded REDUCED-DIRECT was not justified as a general automatic replacement.

Enforce that conclusion.

Requirements:

- `backend=direct` must mean canonical neighbour-list DIRECT unless an explicitly documented internal policy has measured justification.
- `do_reduced=Y` must not silently imply REDUCED-DIRECT execution.
- The reduced stencil may still be built when required for convolution or tests.
- Do not remove the reduced-stencil representation; it remains valuable infrastructure and convolution input.

Choose the least invasive design.

Acceptable outcomes include:
- REDUCED-DIRECT becomes an explicit internal/testing mode;
- REDUCED-DIRECT remains available only through a dedicated developer option;
- REDUCED-DIRECT is retained as infrastructure but no longer used in ordinary DIRECT execution.

Do not expose unnecessary new user-facing options unless needed.

### Correctness

DIRECT and REDUCED-DIRECT parity must remain tested.

### Performance

Re-run long-range Nd DIRECT after the policy fix and confirm that `backend=direct` no longer inherits the known REDUCED-DIRECT regression.

## C. Sparse multi-basis capability

Current policy marks multi-basis sparse support as not fully validated.

Resolve this.

Preferred path:

1. construct a small periodic `N_A > 1` scalar-J fixture;
2. compare SPARSE against canonical DIRECT on random non-collinear states;
3. compare field and field-derived pair energy;
4. test multiple ensembles if supported.

If parity is strong and no structural limitation exists:
- promote sparse multi-basis to `SUPPORTED`.

Otherwise:
- add an explicit sparse capability rejection for unsupported multi-basis use.

Do not leave `NOT_VALIDATED` as silently user-selectable production behavior.

## D. Partial-range backend semantics

Audit every production caller of partial-range `effective_field` requests.

Determine whether such calls can occur with explicit SPARSE or CONVOLUTION selection.

Choose one documented policy:

### Policy 1
Global backends explicitly support only full-range field evaluation.
Partial-range requests under explicit sparse/convolution selection reject clearly.

### Policy 2
Partial-range requests deliberately fall back to DIRECT.
If so:
- document this as part of backend semantics;
- emit an appropriate diagnostic at developer/debug verbosity;
- ensure the fallback never changes physics.

Do not leave the behavior accidental.

Add tests for whichever policy is selected.

## E. Natural target-order overhead

Audit the SFC target-order integration.

If target ordering is NATURAL and no weighted/SFC traversal is requested:
- benchmark the current ordered path against the original/simple native atom/ensemble traversal;
- examine multiple ensembles explicitly.

If the extra target-order indirection or changed loop structure causes measurable regression:
- bypass target-order machinery for natural order.

Keep SFC/weighted ordering available when explicitly selected or beneficial.

Do not remove locality support merely for code simplification.

## F. Benchmark throughput bug

Audit HAM-05 and the shared benchmark analysis code.

The quantity previously labelled `Spin-steps/s` was found to be approximately:

`1 / steady_step_time`

rather than:

`Natom / steady_step_time`.

Correct:
- the calculation;
- the label;
- existing generated evidence where practical;
- any shared benchmark helper using the same incorrect expression.

Also define interaction throughput precisely.

If a full ASD timestep contains two Hamiltonian field evaluations, state whether the metric is:

`directed_interactions_per_field_second`

or:

`directed_interactions_processed_per_ASD_step_second`

and include the appropriate factor.

Do not use ambiguous names.

## G. Capability table reconciliation

Update the HAM-07 backend capability table from live source and new tests.

Use only:
- `SUPPORTED`
- `UNSUPPORTED`
- `NOT_VALIDATED`

But a `NOT_VALIDATED` capability must not be silently accepted as unrestricted production behavior.

## H. Regression suite

At minimum run:
- scalar J DIRECT;
- REDUCED-DIRECT oracle;
- SPARSE scalar J;
- CONVOLUTION scalar J;
- J+D convolution;
- 2D DMI fixture;
- 3D J+D fixture;
- multi-ensemble cases;
- backend explicit-selection rejection tests.

## I. Completion criteria

HAM-08 is complete only when:
- backend names match actual execution;
- no known slower backend is silently substituted for DIRECT;
- sparse multi-basis state is resolved;
- partial-range behavior is intentional;
- natural-order overhead is measured/resolved;
- throughput reporting is correct.

## Checklist

- [x] Live backend routing audited.
- [x] DIRECT no longer silently means REDUCED-DIRECT.
- [x] Reduced stencil retained for legitimate uses.
- [x] Nd DIRECT regression removed/rechecked.
- [x] Sparse multi-basis fixture added.
- [x] Sparse multi-basis promoted or rejected explicitly.
- [x] Partial-range production callers audited.
- [x] Partial-range backend policy documented.
- [x] Partial-range tests added.
- [x] Natural target-order overhead benchmarked.
- [x] Natural-order fast path restored if beneficial.
- [x] Multi-ensemble target traversal checked.
- [x] Spin-throughput formula fixed.
- [x] Interaction-throughput definition fixed.
- [x] HAM-07 capability table reconciled.
- [x] Full backend regression suite passes.

## Commit

`CPU-HAM-08: close CPU Hamiltonian backend policy gaps`
