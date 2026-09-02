# CPU-HAM-03A — Audit historical sparse Hamiltonian support

**Model:** Luna

## Dependency

`CPU-HAM-02C` complete.

## Purpose

Determine exactly what remains of historical `do_sparse` support before creating or reviving a production sparse backend.

No production implementation in this task.

## A. Search current source

Trace:
- `do_sparse` declaration;
- input handling;
- setup handling;
- historical source references;
- build-system references;
- MKL sparse calls;
- documentation/changelog statements;
- tests/examples.

Classify sparse functionality:
- `LIVE`
- `PARTIALLY_CONNECTED`
- `DEAD_BUT_RECOVERABLE`
- `HISTORICAL_ONLY`

## B. Version history

Use git history where helpful.

Do not restore old source blindly.

Document:
- representation of `J`;
- representation of DMI;
- MKL API used;
- API obsolescence;
- threading model;
- setup lifecycle;
- field/energy semantics.

## C. Current architecture fit

Design minimal scalar-J sparse backend:
- persistent sparse `J`;
- three spin-component RHS;
- sparse matrix × dense 3-column matrix.

Prefer CSR + SpMM where supported.

Do not create a `3N × 3N` matrix for scalar `J`.

## D. Build dependency

Determine:
- whether oneMKL sparse is already optional;
- non-MKL behavior;
- alternative sparse providers.

Do not make MKL mandatory.

## E. Threading contract

Define who owns threads during sparse apply.

Avoid nested pattern where every OpenMP worker calls a multithreaded MKL sparse operation.

## F. Energy

Sparse pair energy must use canonical field-derived semantics.

## G. Deliverable

Create `docs/CPU_HAM_03A_SPARSE_AUDIT.md`.

## H. Stop/Go gate

End with `GO` or `NO_GO` for CPU-HAM-03B.

A NO_GO is acceptable.

## Checklist

- [x] `do_sparse` traced.
- [x] Historical implementation located if possible.
- [x] Current connectivity status classified.
- [x] Historical MKL API assessed.
- [x] Scalar-J sparse representation designed.
- [x] Three-RHS approach assessed.
- [x] Build dependency policy defined.
- [x] Thread ownership defined.
- [x] Energy semantics defined.
- [x] GO/NO_GO recommendation produced.
- [x] No production sparse implementation committed.

## Commit

`CPU-HAM-03A: audit sparse Hamiltonian backend`
