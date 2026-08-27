# Benchmark cases

Immutable case definitions: one directory per case family, each holding a case
manifest and the input templates it owns.

Planned families (blueprint §3):

| Case | Workload | Status |
| --- | --- | --- |
| [`B01_bccFe`](B01_bccFe/README.md) | bcc Fe, medium-range interactions — central reference case | Admitted (WP-08A) |
| [`B02_skyrmion2D`](B02_skyrmion2D/README.md) | 2D skyrmion, short-range J+D | Admitted (WP-08B) |
| [`B03_skyrmion3D`](B03_skyrmion3D/README.md) | 3D skyrmion/chiral magnet, short-range J+D | Admitted (WP-08C) |
| [`B04_dhcpNd`](B04_dhcpNd/README.md) | dhcp Nd, very long-range interactions | Admitted (WP-08D) |
| [`B05_dipoleFFT`](B05_dipoleFFT/README.md) | open-boundary dipole/FFT workload | Admitted (WP-08E) |

Rules that already apply:

* **Templates are immutable inputs.** A benchmark run never rewrites anything
  here; it copies what it needs into its own generated work directory
  (`benchmarks/harness/cases.py::generate_run_directory`), and it never
  executes inside a tracked template directory.
* **No generated output in this directory.** Results and work directories live
  under gitignored paths — see the result-data policy in [`../README.md`](../README.md).
* **No physics simplification for convenience.** Interaction shells are not
  truncated, cutoffs are not reduced, weak interactions are not discarded and
  the basis is not simplified in order to make benchmarking easier.

All five planned case families (blueprint §3) are now admitted.
`B01_bccFe` was admitted by **WP-08A**; see
[`B01_bccFe/README.md`](B01_bccFe/README.md) for its template audit,
size-ladder rationale and scaling/sanity validation. `B02_skyrmion2D` was
admitted by **WP-08B**; see
[`B02_skyrmion2D/README.md`](B02_skyrmion2D/README.md) for its template
audit (including a rejected sibling template) and two real backend-
compatibility findings (a GPU-unsupported `skyno` mode, and a cross-binary
chiral-domain degeneracy in the maintainer's own field-free initial-phase
search). `B03_skyrmion3D` was admitted by **WP-08C**; unlike `B01`/`B02`,
no existing maintainer example deck is a genuinely bulk (three-
dimensionally periodic) short-range J+D system, so this case's template is
a maintainer-approved *derived* construction rather than a verbatim copy
-- see [`B03_skyrmion3D/README.md`](B03_skyrmion3D/README.md) for the full
audit, the exact derivation from `B02_skyrmion2D`'s template, and why this
case's simplified initial condition makes it a *stronger* CPU/GPU
bit-exactness discriminator than `B02_skyrmion2D`'s. `B04_dhcpNd` was
admitted by **WP-08D**; unlike `B01`-`B03`, its ~1338
directed-interactions/atom exchange cutoff is large enough to alias onto
itself through the periodic boundary below a real, measured minimum
supercell size, and large enough that CPU (not only GPU) memory
availability becomes a genuine per-size classification question -- see
[`B04_dhcpNd/README.md`](B04_dhcpNd/README.md) for the aliasing-floor
measurement, the new `harness/host_memory.py` CPU memory classifier this
admission added (mirroring WP-06's GPU-side `gpu_memory.py`), and why this
case's size ladder deliberately stops short of the ~4M-atom tier every
other admitted case reaches. `B05_dipoleFFT` was admitted by **WP-08E**;
unlike every other case, its dominant per-step cost is a global FFT-based
dipole calculation rather than a neighbour-list sum -- no maintainer
`examples/` deck exercises this, so its template is derived (with two
disclosed, minimal changes) from this project's own internal dipole-
validation base fixtures; see
[`B05_dipoleFFT/README.md`](B05_dipoleFFT/README.md) for that audit, a
confirmed-broken legacy CPU FFT-dipole implementation this case
deliberately avoids using, a real "silently drops the dipole term with no
error" backend-override trap, and why this case's size ladder is
thin-film/racetrack-shaped rather than the isotropic cubes every other
case uses.

## Case manifest (WP-02)

Each case is one directory holding a `case.yaml` (or `.json`) manifest,
validated against
[`../schema/case_manifest.v1.schema.json`](../schema/case_manifest.v1.schema.json),
plus the `template_directory` it points to. See
[`benchmarks/harness/cases.py`](../harness/cases.py) for the loader and
[`INFRA_test_only/case.yaml`](INFRA_test_only/case.yaml) for a complete
worked example.

| Field | Meaning |
| --- | --- |
| `id`, `description` | Case identity; `id` matches TERMINOLOGY.md `case_id` and the manifest's directory name. |
| `infrastructure_test_only` | `true` only for the one synthetic fixture case (see below). Omitted/`false` for every real case. |
| `template_directory` | Path to the immutable input template, relative to the manifest. Never written to. |
| `workload_class` | Case-level interaction-range tag: `short_range_neighbor`, `medium_range_neighbor`, `long_range_neighbor`, `fft_dipole`, `neighbor_plus_dipole`. Maps onto the per-record `workload_class` enum (`*_neighbor` → `NEIGHBOR_LIST`, `fft_dipole` → `FFT_DIPOLE`, `neighbor_plus_dipole` → `NEIGHBOR_LIST_PLUS_FFT_DIPOLE`). |
| `dimensionality` | `3D` (scales `[Nx, Ny, Nz]`), `2D` (scales `[Nx, Ny]`, fixes `thickness_cells`), or `case_specific` (an explicit legal replication tuple whose meaning is the case's own). |
| `scalable_dimensions` | Names of the tuple `sizes[].replication` varies, in order. |
| `variants` | Alternative physics/runtime configurations of the case (TERMINOLOGY.md VARIANT), each an `id`, `description` and an `overrides` map. Every override key must also appear in `allowed_input_overrides`. |
| `sizes` | The case's legal size ladder (TERMINOLOGY.md SIZE): `id` plus a `replication` tuple. Not auto-generated — a case is not forced onto a power-of-four ladder. |
| `expected_physics` | Documentation (`lattice`, `hamiltonian_terms`, `notes`) auditors can check overrides against. Never machine-enforced beyond the allow-list. |
| `allowed_input_overrides` | The closed list of `inpsd.dat` keywords this case permits changing. Must be a subset of the harness-wide `cases.GLOBALLY_SAFE_OVERRIDE_KEYS` and must include `ncell` (every size is applied through it). |
| `workload_metadata_method` | Name of the parser in [`../harness/workload_metadata.py`](../harness/workload_metadata.py) used to determine `natom`/`directed_interactions`/`mean_neighbors`/`max_neighbors` or `natom`/`fft_grid`/`fft_grid_padded`/`fft_grid_points`. |

### Allow-listed input overrides

`cases.GLOBALLY_SAFE_OVERRIDE_KEYS` is the harness-wide ceiling on what any
case may ever declare overridable: `ncell` (supercell replication), `Nstep`
(step count), `tseed` (RNG seed), `temp` (temperature, only when a case
variant defines it), `avrg_step`/`cumu_step`/`tottraj_step`/`ene_step`
(measurement cadence), `do_prnstruct` (requests UppASD's own
`struct.<simid>.out` diagnostic dump, which the neighbour-workload metadata
parser reads), `gpu_mode` (CPU/GPU backend dispatch — see
[`B01_bccFe/README.md`](B01_bccFe/README.md#backend-dispatch-gpu_mode) for
why this is a per-run override rather than template content), and `skyno`
(skyrmion-number diagnostic method — see
[`B02_skyrmion2D/README.md`](B02_skyrmion2D/README.md#backend-dispatch-gpu_mode-skyno)
for why the GPU backend needs this overridden away from the maintainer's
own `T` (triangulation) value, which it does not implement), and
`do_dip`/`gpu_dipole_mode` (CPU/GPU dipole Hamiltonian on/off — see
[`B05_dipoleFFT/README.md`](B05_dipoleFFT/README.md#a-template-audit) for
why a GPU dipole sample needs both set together with `gpu_mode`, and what
happens if one is forgotten). `temp`/`do_dip`/`gpu_dipole_mode` are the
only entries here that change genuine Hamiltonian content rather than pure
backend dispatch or diagnostics — each is safe-listed only because a case's
own `variants` define the legal values, the same "a case variant defines
it" justification `temp` already established. None of these changes any
lattice or moment parameter. A case's
`allowed_input_overrides` only ever narrows this set. Anything outside it —
exchange/DMI parameters, interaction cutoffs, lattice structure
(`cell`/`BC`/`Sym`), magnetic moments, or any other Hamiltonian term — is
rejected outright, at manifest-load time for variant overrides and at
override-build time for any other requested override; it is never silently
dropped and never silently applied.

### The infrastructure fixture

`INFRA_test_only/` is the one case with `infrastructure_test_only: true`: a
tiny synthetic 2-atom-basis template that exists solely to exercise this
machinery (manifest parsing, template copying, size generation, input
overrides, fingerprints — see `benchmarks/tests/test_case_manifest.py`). Its
physics is not validated and not meaningful.
`cases.filter_cases_for_campaign(..., authoritative=True)` and
`cases.require_not_infrastructure_only` are how a later campaign layer (WP-09)
must refuse it when generating an authoritative report.
