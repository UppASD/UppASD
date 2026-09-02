# CPU-HAM-02B — Build and validate a reduced periodic pair stencil

**Model:** Luna

## Dependency

`CPU-HAM-02A` complete.

## Stop condition

If reduced-Hamiltonian semantics, periodic images, or translation equivalence are ambiguous, stop and report rather than guessing.

## Purpose

For eligible reduced periodic Hamiltonians, represent repeated pair geometry once per basis-pair/lattice displacement rather than relying on a full arbitrary atom-index neighbour representation in the hot loop.

This task builds the representation and oracle. Do not optimize aggressively yet.

## A. Establish eligibility

Trace existing production reduced-Hamiltonian semantics.

Document exact requirements for a translational stencil, expected to include:
- periodic lattice;
- regular replicated cells;
- reduced Hamiltonian;
- equivalent coupling pattern under cell translation.

Do not infer them merely from terminology.

## B. Define canonical stencil record

For scalar `J`, use a compact record equivalent to:

```text
output_basis
input_basis
delta_cell_1
delta_cell_2
delta_cell_3
J
```

If wrapped displacement/image information is needed to distinguish interactions, include it.

## C. Cell/basis mapping

Build deterministic:

```text
physical atom <-> (cell1, cell2, cell3, basis)
```

mapping and prove round-trip correctness.

Do not globally reorder atoms.

## D. Construct stencil from production Hamiltonian

Derive the reduced stencil from already accepted Hamiltonian data.

Do not independently rederive interaction physics.

## E. Periodic wrapping

Implement exact cell-index wrapping consistent with production PBC semantics.

Add boundary-crossing fixtures.

## F. Oracle comparison

Write a host/reference stencil application and compare with canonical DIRECT `HamiltonianActions` on:
- one-basis periodic system;
- multi-basis periodic system;
- long-range `J`;
- boundary-crossing interactions;
- random non-collinear spins.

## G. Negative eligibility tests

Reject or decline:
- nonperiodic case;
- non-reduced/ineligible Hamiltonian;
- translation-invariance violation if safely constructible.

Do not silently approximate.

## H. Memory evidence

Compare metadata footprint of full direct `nlist/ncoup` used by field evaluation against reduced-stencil metadata for Fe and Nd.

Do not claim total Hamiltonian memory savings if original arrays remain required elsewhere.

## I. No production switch yet

Do not make reduced stencil the production default.

## Implementation result

The production reduced-Hamiltonian trace is:

- do_reduced='Y' sets nHam=NA;
- setup_aHam maps each physical atom to its compact basis Hamiltonian entry;
- setup_neighbour_hamiltonian keeps the full physical nlist, while the
  reduced nlistsize/ncoup entries are populated for the basis entries used
  by the field kernel.

For the regular non-random-alloy layout, physical atom a maps to
(cell1,cell2,cell3,basis) using:

    a = basis + NA * (cell1 + N1 * (cell2 + N2 * cell3))

with zero-based cell coordinates and one-based basis/atom IDs. The new
ReducedStencil module derives records from the accepted production nlist,
nlistsize, aHam, and scalar ncoup arrays. Each record contains
output_basis, input_basis, delta_cell(1:3), and J. It does not rederive
exchange physics. Eligible production setup calls materialize this object in
ham%reduced_stencil after the canonical exchange lists are mounted.

Eligibility is explicit and conservative: do_reduced='Y', no random alloy,
nHam=NA, a regular NA*N1*N2*N3 layout, fully periodic P P P boundaries,
and a single scalar-J coupling configuration. Construction validates every
translated physical target against the compact first-cell stencil. The
reference apply uses modulo wrapping and owner-computes target accumulation.
No production dispatch or default was changed.

## Correctness and memory evidence

cpu-ham-reduced-stencil passes deterministic random non-collinear fixtures
for:

- one basis with periodic boundary crossings;
- two basis atoms with distinct basis-pair and cell displacements;
- a seven-record-per-basis long-range periodic case;
- explicit cell/basis round-trip mapping and negative nonperiodic/non-reduced
  eligibility cases;
- a deliberately corrupted translated nlist, which is declined with a
  translation-invariance diagnostic.

The host stencil field matches canonical production HamiltonianActions
DIRECT field evaluation, and the pair energy matches
-1/2 sum(m dot B_pair) after the production mub/mry conversion for all
three positive fixtures. The existing CPU-HAM-01A and CPU-HAM-02A tests also
pass. A real reduced FeCo production startup and 5,000-step run reports
“Validated reduced scalar-J stencil available” and completes successfully.

The metadata probe counts the direct nlist, nlistsize, aHam, and reduced
scalar ncoup allocation shapes against packed 28-byte stencil records:

| Reference | Full direct metadata | Reduced-stencil metadata | Ratio |
|---|---:|---:|---:|
| bcc Fe, 16,000 atoms, 96 neighbours, NA=2 | 6,273,536 B | 5,376 B | 1167.0x |
| dhcp Nd, 16,384 atoms, 1,340 allocated max neighbours, NA=4; 1,338 records | 87,992,192 B | 149,856 B | 587.2x |

These are metadata representation figures only. The current production full
arrays remain allocated elsewhere and other Hamiltonian terms are outside
this scalar-J comparison; total Hamiltonian memory savings are not claimed.

## Checklist

- [x] Reduced-Hamiltonian semantics documented.
- [x] Eligibility rules explicit.
- [x] Cell/basis mapping proven.
- [x] Canonical stencil record defined.
- [x] Stencil derived from production Hamiltonian.
- [x] PBC wrapping validated.
- [x] One-basis field oracle passes.
- [x] Multi-basis field oracle passes.
- [x] Long-range field oracle passes.
- [x] Boundary-crossing fixture passes.
- [x] Nonperiodic case rejected/declined.
- [x] Ineligible Hamiltonian rejected/declined.
- [x] Translation-invariance assumptions tested.
- [x] Metadata footprint measured.
- [x] No production default changed.

## Commit

`CPU-HAM-02B: add reduced periodic pair stencil`
