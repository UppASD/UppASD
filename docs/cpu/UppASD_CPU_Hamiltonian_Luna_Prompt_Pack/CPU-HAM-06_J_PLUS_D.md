# CPU-HAM-06 — Extend viable CPU global backends to exchange plus DMI

**Model:** Luna

## Dependency

`CPU-HAM-05` complete.

## GO/NO-GO gate

Only extend backends competitive enough to justify complexity.

DIRECT already supports canonical `J+D`.

Possible extensions:
- REDUCED-DIRECT;
- SPARSE;
- CONVOLUTION.

Do not extend clearly dominated backends.

## Stop condition

Any DMI convention ambiguity must stop the task.

## A. Freeze DMI convention

Document canonical `HamiltonianActions` DMI field and energy convention from production code and established tests.

Do not infer signs from textbook notation.

## B. Reduced stencil

Represent each DMI interaction with the exact oriented displacement/pair convention required to reproduce canonical DMI.

Validate reciprocal/antisymmetric relations.

## C. Sparse

If retained, prefer `3×3` block operator representation for `J+D`.

Use:

\[
\mathsf C(\mathbf D)\mathbf m=\mathbf D\times\mathbf m.
\]

The block matrix must reproduce canonical endpoint convention.

## D. Convolution

For each basis pair build:
- `J(q)`;
- `Dx(q)`;
- `Dy(q)`;
- `Dz(q)`.

Compute:

\[
\mathbf B_J(q)=J(q)\mathbf m(q)
\]

and

\[
\mathbf B_D(q)=\mathbf D(q)\times\mathbf m(q).
\]

Do not jump to generic 9-component tensor unless required.

## E. Energy

Use complete canonical pair field:

\[
E_{\rm pair}=-\frac12\sum \mathbf m\cdot\mathbf B_{\rm pair}.
\]

Validate against direct `HamiltonianActions`.

## F. Required negative controls

At minimum:
1. reverse DMI sign;
2. transpose/reverse antisymmetric operator incorrectly.

Both must make parity fail.

## G. Fixtures

Use:
- 2D skyrmion `J+D`;
- 3D skyrmion/chiral `J+D`;
- discriminating DMI small periodic fixture.

## H. Performance

Benchmark only after parity.

Determine whether backend remains competitive with DMI overhead.

## I. Tensor exchange

Do not implement full general tensor exchange.

Recommend later extension only if justified.

## Checklist

- [x] Canonical DMI convention documented.
- [x] REDUCED-DIRECT DMI implemented and retained.
- [x] SPARSE DMI not retained: CPU-HAM-05 showed no robust crossover.
- [x] CONVOLUTION DMI implemented and retained.
- [x] `J(q)+D(q)` representation validated.
- [x] Field parity passes 2D.
- [x] Field parity passes 3D.
- [x] Energy parity passes.
- [x] Sign-reversal negative control fails.
- [x] Antisymmetry/transposition negative control fails.
- [x] `J+D` performance measured.
- [x] Tensor-extension recommendation recorded.

## Commit

`CPU-HAM-06: extend CPU pair backends to DMI`
