# CPU-HAM-07 — Define production CPU pair-backend policy

**Model:** Luna

## Dependencies

`CPU-HAM-05` complete. `CPU-HAM-06` complete if pursued.

## Purpose

Turn measured backend capabilities and crossovers into a conservative production selection policy.

Do not invent heuristics from intuition.

## Implementation result

The production selector is the `cpu_ham_backend` input keyword. It accepts
`direct`, `sparse`, and `convolution` (case-insensitively). The shorter
`ham_backend` spelling is accepted as an input alias. The historical
`do_sparse` and `do_convolution` flags remain compatibility aliases when the
canonical request is `direct`; selecting both, or mixing either alias with a
different canonical backend, is rejected.

Production setup now resolves and reports the requested and resolved backend
before the simulation starts. DIRECT is the safe general path. An explicit
SPARSE or CONVOLUTION request must complete its capability checks and
persistent setup; an ineligible request terminates with a capability
diagnostic instead of silently falling back to DIRECT. The backend state is
fixed for the run, so changing legacy flags during a timestep cannot change
the execution representation.

AUTO is intentionally disabled. CPU-HAM-05 found a robust convolution win for
the measured long-range dhcp Nd workload, but no portable rule spanning Nd,
bcc Fe, and the short-range control. The crossover is also provider- and
architecture-dependent. No MKL run is required for this policy: the existing
CPU convolution implementation was validated through FFTW, while MKL remains
an unvalidated provider on this machine.

### Capability table

`SUPPORTED` means the production selector admits the capability and the
backend implementation has a defined path. `UNSUPPORTED` is rejected rather
than approximated. `NOT_VALIDATED` means the implementation boundary is
present but the campaign has not yet supplied an end-to-end fixture for that
combination.

| Backend | scalar J | DMI | tensor exchange | periodic | nonperiodic | reduced Hamiltonian | non-reduced Hamiltonian | disorder | multi-basis | multiple ensembles |
|---|---|---|---|---|---|---|---|---|---|---|
| DIRECT | SUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED |
| SPARSE | SUPPORTED | UNSUPPORTED | UNSUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED | SUPPORTED | UNSUPPORTED | NOT_VALIDATED | SUPPORTED |
| CONVOLUTION | SUPPORTED | SUPPORTED | UNSUPPORTED | SUPPORTED* | UNSUPPORTED | SUPPORTED | UNSUPPORTED | UNSUPPORTED | SUPPORTED | SUPPORTED |

`*` CONVOLUTION means fully periodic, regular replicated cells with a
translationally reduced Hamiltonian. It is currently provided through FFTW;
an MKL CPU convolution provider is not validated here. SPARSE is the
portable persistent scalar-J CSR implementation; CPU-HAM-06 deliberately did
not retain a sparse DMI extension.

## A. Capability table

For each backend document:
- scalar `J`;
- DMI;
- tensor exchange;
- periodic;
- nonperiodic;
- reduced Hamiltonian;
- non-reduced Hamiltonian;
- disorder;
- multi-basis;
- multiple ensembles.

Use:
- `SUPPORTED`
- `UNSUPPORTED`
- `NOT_VALIDATED`

No ambiguous entries.

## B. User-selection API

Define explicit backend choices following UppASD style, conceptually:
- `direct`;
- `sparse`;
- `convolution`.

REDUCED-DIRECT may remain internal to DIRECT if evidence supports it.

Do not expose unnecessary implementation detail.

## C. DIRECT default safety

If explicitly requested backend is unsupported, reject clearly.

AUTO may fall back to DIRECT internally.

Never silently approximate unsupported physics.

## D. AUTO mode

Implement AUTO only if crossover evidence is robust.

AUTO must first check correctness eligibility, then performance.

Potential metadata:
- `Natom`;
- `Ndirected`;
- mean neighbours;
- basis size;
- periodic/reduced state;
- Hamiltonian terms.

Do not use simplistic `z` threshold without evidence.

## E. Machine dependence

Recognize crossover is CPU-architecture dependent.

If no portable rule exists, leave AUTO disabled and retain explicit selection.

That is acceptable.

## F. Startup diagnostic

Print requested backend, resolved backend, reason and eligibility information.

## G. Tests

For every backend:
- valid explicit selection works;
- invalid explicit selection rejects;
- AUTO never selects ineligible backend;
- field/energy parity holds.

## H. Documentation

Document likely regimes without universal guarantees:
- DIRECT → general/default;
- SPARSE → optional irregular large-pair workload where benchmarked;
- CONVOLUTION → periodic reduced translational Hamiltonians, especially long range.

## I. Final report

Update master blueprint with retained backends, rejected experiments, measured crossover conclusions and remaining work.

## Checklist

- [x] Capability table complete.
- [x] Explicit backend selection defined.
- [x] DIRECT remains safe general path.
- [x] Unsupported explicit selection rejects.
- [x] No backend approximates unsupported physics.
- [x] AUTO implemented only if evidence warrants; AUTO remains disabled.
- [x] AUTO eligibility precedes performance heuristic; no heuristic is enabled.
- [x] Machine dependence addressed.
- [x] Startup diagnostics clear.
- [x] Field parity tests pass.
- [x] Energy parity tests pass.
- [x] User documentation updated.
- [x] Master blueprint reconciled.

### Validation evidence

`tests/hamiltonian/test_cpu_ham07_backend_policy.f90` covers the public
vocabulary, case normalization, invalid names, AUTO rejection, legacy alias
resolution, and conflicting requests. The existing scalar-J SPARSE parity
test and periodic scalar-J/J+D CONVOLUTION parity tests now exercise the
production selector. DIRECT remains covered by the canonical field/energy
contract and all backend parity tests compare fields and field-derived energy
against the DIRECT oracle.

The local FFTW build passed the CPU-HAM-07 policy test, SPARSE parity,
CONVOLUTION parity, J+D parity, and the field/energy contract. No MKL-specific
validation was needed.

## Commit

`CPU-HAM-07: define CPU Hamiltonian backend policy`
