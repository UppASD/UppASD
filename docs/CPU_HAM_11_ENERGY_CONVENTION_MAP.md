# CPU-HAM-11 — Energy convention map

**Date:** 2026-09-02
**Status:** complete
**Authority:** `source/Hamiltonian/hamiltonianactions.f90`

## Global energy call graph

`Energy::calc_energy` is the production measurement/output entry point. For
non-LSF measurements it now requests `term_fields` from
`HamiltonianActions::effective_field`; the same call performs the field
assembly and the backend-specific pair application. `DIRECT`, `SPARSE`, and
eligible `CONVOLUTION` therefore enter the same term-energy reduction
interface. LSF measurement remains delegated to `LSF::totalenergy_LSF`, which
has moment-size-dependent semantics.

The optional scalar-energy interface of `effective_field` remains backward
compatible for diagnostics and replica drivers. Its field-only path is
unchanged and still returns zero energy when `measure_energy=.false.`.

The remaining global-looking routes are `SpinIce::calculate_spinice_energy`
(a loop-algorithm reference over argument-owned spin-ice data),
`LSF::totalenergy_LSF`, and the induced-moment routines. They are included in
the exception inventory below rather than silently routed through the
atomistic state-owned helper.

| energy term | current routine | canonical field routine | same physics implementation? | can derive/reuse? | performance concern | action |
|---|---|---|---|---|---|---|
| isotropic exchange | `Energy::calc_energy` | `heisenberg_field` / `heisenberg_rescaling_field` | Yes for non-LSF measurement | Yes | One target traversal on requested measurement; sparse can supply the pair field | Route through `term_fields`; remove the old global exchange loop |
| DMI | `Energy::calc_energy` | `dzyaloshinskii_moriya_field` | Yes | Yes | Convolution returns combined J+D, so term separation needs one canonical DMI evaluation | Split combined J+D into exchange and DMI term fields |
| symmetric anisotropic exchange | `Energy::calc_energy` | `symmetric_anisotropic_field` | Yes | Yes | Only evaluated when enabled and requested | Preserve `ene_sa`; map it to the canonical term field |
| tensor exchange | `Energy::calc_energy` | `tensor_field` | Yes | Yes | Tensor traversal is measurement-only | Preserve `ene_pair` and the tensor term slot |
| pseudo-dipolar exchange | `Energy::calc_energy` | `pseudo_dipolar_field` | Yes | Yes | Term-specific neighbour traversal | Preserve `ene_pd`; use its field projection |
| BIQDM | `Energy::calc_energy` | `dzyaloshinskii_moriya_bq_field` | Yes | Yes | Higher-order field remains term-specific | Preserve `ene_bqdm`; retain its 1/2 relationship |
| biquadratic exchange | `Energy::calc_energy` | `biquadratic_field` | Yes | Yes | Higher-order local-neighbour work | Preserve `ene_bq`; retain its 1/4 relationship |
| four-spin ring | `Energy::calc_energy` | `ring_field` | Yes | Yes | Multispin field is not a pair traversal | Preserve `ene_ring`; retain its 1/4 relationship |
| scalar chirality | `Energy::calc_energy` | `chirality_field` | Yes | Yes | Multispin field is term-specific | Preserve `ene_chir`; retain its 1/2 relationship |
| dipole | `Energy::calc_energy` | `DipoleManager::dipole_field_calculation` | Yes for field/energy result | Yes | Dipole manager owns its method and macrocell treatment | Expose the incoming dipole field as a term; keep manager accumulation authoritative |
| onsite anisotropy | `Energy::calc_energy` | `uniaxial_anisotropy_field` / `cubic_anisotropy_field` | Yes after explicit expression review | Yes, via `canonical_onsite_energy` | No generic field prefactor is assumed | Preserve `ene_ani`; use the explicit uniaxial/cubic polynomial |
| static and time-dependent Zeeman | `Energy::calc_energy` | external-field assembly in `HamiltonianActions` | Yes | Yes | O(N) requested measurement work | Use `-m·B_ext` with factor 1 |
| spin-ice loop reference total | `SpinIce::calculate_spinice_energy` | no state-owned `HamiltonianActions` equivalent | No | No | Argument-owned loop/vertex contract | Retain as a specialized MC global reference |
| LSF energy | `LSF::totalenergy_LSF` | LSF field/interpolation routines | Specialized | No generic atomistic reduction | Moment-size interpolation and LSF tables | Retain as an explicit LSF exception |

## Canonical relationships

For a reciprocal bilinear pair field, the global energy is

\[
 E_\mathrm{pair}=-\frac12\sum_i \mathbf m_i\cdot\mathbf B_i^\mathrm{pair}.
\]

This applies separately to isotropic exchange, DMI, SA, PD, tensor exchange,
BIQDM, dipole, and chirality where the production field is the derivative of
the corresponding global term. BQ and ring retain their established 1/4
field projections. The Zeeman term is

\[
 E_Z=-\sum_i \mathbf m_i\cdot\mathbf B_i^\mathrm{ext}.
\]

Onsite anisotropy is not inferred from homogeneity. `canonical_onsite_energy`
uses the production expression directly. For a primary uniaxial axis with
`c=m·e`, that expression is

\[
 E_\mathrm{ani}=K_1c^2+2K_2c^2-K_2c^4,
\]

with the corresponding cubic and optional secondary-axis expressions kept
explicit in the helper. This preserves the production field convention while
removing the former MC polynomial drift.

## Result decomposition and measurement cadence

`term_fields(3,HAM_TERM_COUNT,Natom,Mensemble)` is an optional output of the
full `effective_field` interface. It carries the minimum decomposition needed
to preserve the existing `ene_t` and `localenergy` outputs:

| slot | term |
|---:|---|
| 1 | scalar exchange |
| 2 | DMI |
| 3 | symmetric anisotropic exchange |
| 4 | pseudo-dipolar |
| 5 | BIQDM |
| 6 | biquadratic |
| 7 | ring |
| 8 | onsite anisotropy |
| 9 | Zeeman |
| 10 | dipole |
| 11 | chirality |
| 12 | tensor exchange |

The global energy routine allocates this decomposition only on a non-LSF
measurement call. Normal LLG field calls continue to pass
`measure_energy=.false.` and use the lean field-only path. No term result is
required on the hot dynamics path.

## MC exception matrix

The specialized `montecarlo_common::calculate_energy` routine remains a
single-flip transition/acceptance kernel. It is allowed to use local outgoing
neighbour work instead of a global traversal, but its supported terms are
checked against canonical global before/after energies in
`tests/hamiltonian/test_cpu_ham11_energy.f90`.

| local path | status | reason |
|---|---|---|
| exchange | parity fixture passes | Local outgoing contribution equals the reciprocal global delta |
| DMI | parity fixture passes | Frozen HAM-06 handedness and directed-list convention are preserved |
| onsite anisotropy | parity fixture passes | MC now calls `canonical_onsite_energy` |
| SA, PD, BIQDM, BQ, ring, chirality | retained specialized path | Local updates remain term-specific; broad MC rewrite is out of scope |
| dipole/macrocell | retained specialized path | Macrocell trial state and manager layout require a local algorithm |
| LSF/induced moments | retained specialized path | Moment-size/interpolation state cannot use the atomistic global helper |
| heat-bath `calculate_efield` | retained local field path | Transition field construction is not a global energy implementation |

The parity test includes controlled negative calculations: negating the
exchange delta and doubling the DMI delta both disagree with the canonical
global delta. These are test-side mutations; no production mutation remains.

## Backend parity and performance review

The same `term_fields` request is exercised for the retained backends:

- `cpu_ham11_energy_tests`: DIRECT canonical global/MC parity;
- `sparse_backend_tests`: scalar-J term fields match DIRECT for Nd, Fe, and
  the multi-basis fixture;
- `cpu_convolution_tests`: convolution term fields match DIRECT on the
  production periodic scalar-J fixture.

The HAM-09/HAM-10 representative Nd, Fe, and J+D measurements remain valid:
energy decomposition is measurement-gated, while the normal LLG path is not
given an additional traversal or allocation. The new parity fixture and
focused backend checks complete in the existing Release build without a
separate global pair-energy pass.

## Explicit exceptions

`LSF::totalenergy_LSF`, `LSF::calculate_energy_wLSF`, induced-moment energy,
and lattice Hamiltonian energy belong to distinct state/physics contracts and
are not silently redirected to the atomistic helper. Historical
`applyhamiltonian.f90` / `heisge()` code is not a production authority and was
not reconnected.
