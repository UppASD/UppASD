# CPU-HAM-07 — Production CPU pair-backend policy

## Selection

Set the production CPU pair backend with:

```text
cpu_ham_backend direct
cpu_ham_backend sparse
cpu_ham_backend convolution
```

`direct` is the default and general path. `sparse` selects the persistent
scalar-J CSR backend. `convolution` selects the persistent periodic reduced
scalar-J/DMI FFT backend when its capability checks pass. The legacy
`do_sparse` and `do_convolution` flags remain accepted as compatibility
aliases, but conflicting requests are rejected.

The startup diagnostic reports the requested backend, resolved backend,
reason, and eligibility result. Explicit unsupported requests terminate with a
clear diagnostic; no backend approximates unsupported physics. `auto` is
rejected because the CPU-HAM-05 measurements do not establish a portable
crossover rule across machines and workloads.

## Capability summary

| Backend | Supported scope | Rejected scope |
|---|---|---|
| DIRECT | General current `HamiltonianActions` physics, periodic/nonperiodic, reduced/non-reduced, disorder, multi-basis, ensembles | None within the canonical production scope |
| SPARSE | Static scalar isotropic J, periodic/nonperiodic, reduced/non-reduced, ensembles | DMI, tensor exchange, disorder/LSF; multi-basis remains not validated end-to-end |
| CONVOLUTION | Fully periodic reduced translational scalar J or validated J+D, multi-basis, ensembles | Nonperiodic/non-reduced, tensor exchange, disorder/LSF, unsupported onsite/pair extensions |

The measured campaign retained DIRECT for general and short-range workloads,
SPARSE as optional, and CONVOLUTION for the measured long-range dhcp Nd
regime. These are likely regimes, not universal performance guarantees.

## Validation

The policy resolver test covers valid names, invalid names, AUTO rejection,
legacy aliases, and conflicts. Existing SPARSE and CONVOLUTION parity tests
exercise the production selector and compare fields and field-derived energy
against canonical DIRECT. Local validation used FFTW with MKL disabled; no
MKL-specific run was required.
