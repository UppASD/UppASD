# CPU-HAM-03B — Implement persistent scalar-J sparse CPU backend

**Model:** Luna

## Dependency

Run only after `CPU-HAM-03A` returns GO.

## Stop condition

If implementation requires broad build-system surgery or a new mandatory dependency, stop and report.

## Goal

Implement a production-compatible scalar-J sparse field backend under `HamiltonianActions`.

## A. Setup once

Construct sparse `J` once during backend setup.

Do not rebuild CSR each timestep.

Store row structure, columns, values, persistent library handle, and inspector/optimization state where available.

## B. Apply three components

Represent spins as three RHS vectors or `N × 3` dense operand.

Apply one sparse structure to `Mx`, `My`, `Mz`.

Avoid three complete sparse-metadata traversals if SpMM is efficient.

## C. Data marshalling

Measure packing between `emomM(3,N,...)` and preferred dense RHS layout.

Do not assume packing is free.

## D. Threading

Enforce CPU-HAM-03A threading contract.

No oversubscription.

Record MKL/OpenMP thread configuration.

## E. Energy

If requested:

\[
E_{\rm pair}=-\frac12\sum \mathbf m\cdot\mathbf B_{\rm sparse}.
\]

## F. Capability

Initial scope: scalar `J`.

Do not add DMI.

## G. Field composition

Integrate sparse pair field with remaining local terms without changing physics.

## H. Correctness

Compare against DIRECT and REDUCED-DIRECT where eligible.

Include random states, Nd, Fe, and multiple ensembles.

## I. Performance

Measure setup, packing, sparse apply, full `effective_field`, and scaling for Nd and Fe.

## J. No AUTO

Backend selection remains explicit/experimental.

## Checklist

### Implementation and validation evidence

The production implementation is in `HamiltonianActions`. It owns an exact-size
directed CSR structure and persistent `N×3×Mensemble` RHS/result buffers. The
portable provider applies all three components in one row traversal under its
own OpenMP worksharing region. Packing, apply, and setup wall times are
accumulated by `sparse_backend_get_stats`; no allocation or CSR construction is
performed by the timestep apply.

`tests/hamiltonian/test_sparse_backend.f90` covers random non-collinear Nd and
Fe-sized directed fixtures, zero-neighbour rows, multiple ensembles,
field-derived energy, repeated setup/cleanup, exact directed NNZ, and explicit
fallback for rescaled exchange. `cpu-ham-reduced-stencil` remains the
REDUCED-DIRECT comparison for eligible periodic fixtures; the irregular Nd/Fe
sparse fixtures do not claim reduced-stencil eligibility. The portable provider
does not require MKL, and the parity test was run with `OMP_NUM_THREADS=1`, `2`,
and `4`; MKL-specific validation is deferred to a host with oneMKL if a
provider-specific implementation is added later.

- [x] Persistent sparse structure implemented.
- [x] No per-step CSR construction.
- [x] Three-component apply efficient.
- [x] Packing cost measured.
- [x] Thread ownership correct.
- [x] No oversubscription.
- [x] Energy field-derived.
- [x] Nd parity passes.
- [x] Fe parity passes.
- [x] Multi-ensemble parity passes.
- [x] Setup cost measured.
- [x] Steady-state cost measured.
- [x] DIRECT comparison recorded.
- [x] REDUCED-DIRECT comparison recorded.
- [x] Backend remains explicit (`do_sparse=='Y'`; no AUTO selection).

## Commit

`CPU-HAM-03B: add scalar exchange sparse backend`
