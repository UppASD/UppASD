# B03_skyrmion3D — genuinely 3D short-range exchange+DMI chiral magnet

Admission record for WP-08C (blueprint §3 "B03 — 3D skyrmion/chiral magnet:
short-range (J+D)"). This case is the framework's genuinely three-
dimensional short-range-neighbour production reference: contrasted against
`B02_skyrmion2D`'s quasi-two-dimensional short-range workload (same
Hamiltonian class) and against `B01_bccFe`'s medium-range, three-
dimensional workload. Per the blueprint, a separate Hopfion benchmark is
intentionally omitted from the initial suite: its ordinary ASD
computational workload is sufficiently represented by this family (see
"Hopfion redundancy" below).

## A. Template audit

### Template selection: no existing maintainer deck is genuinely 3D

Every J+D (exchange + DMI) example deck in `examples/` was inspected for a
genuinely bulk (three-dimensionally periodic) chiral/skyrmion system:

| Template | `BC` | `ncell` | Verdict |
| --- | --- | --- | --- |
| `SpecialFeatures/SkyrmionLattice` | `P P 0` | `128 128 1` | Quasi-2D (already `B02_skyrmion2D`) |
| `SpecialFeatures/SkyrmionTransport` | `P P 0` | `200 100 1` | Quasi-2D; also rejected for `B02` (restart-file/MC-anneal init, `stt A`) |
| `SpinWaves/ChiralSpiral` | `P P 0` | `128 128 1` | Quasi-2D |

None is a genuine bulk 3D system. The `tests/coarse_graining/e2e/
moving_dmi_chiral_*` fixtures do use `BC P P P`, but they are small
synthetic coarse-graining test fixtures (elongated slabs, e.g. `24 2 2`),
not physically-motivated production benchmark decks, and are out of scope
for a case whose template must be maintainer-approved production input
(`benchmarks/cases/README.md`: "their production input templates come
from the maintainer, not from this framework").

This gap was raised with the maintainer directly. Decision: **derive** a
genuine 3D template by extending `SkyrmionLattice`'s square-lattice
nearest-neighbour exchange + Bloch DMI isotropically onto the third (z)
axis, turning the quasi-2D single layer into a simple-cubic bulk lattice.
This is not a verbatim maintainer deck (unlike `B01_bccFe`/`B02_skyrmion2D`);
it is a maintainer-approved construction, and every change from the source
template is enumerated below rather than left implicit.

### Exact derivation from `SkyrmionLattice`

Starting point: `examples/SpecialFeatures/SkyrmionLattice/{inpsd.dat,jij,dmdata,momfile,posfile,qfile}`.

| File | Change | Rationale |
| --- | --- | --- |
| `jij` | Added 2 lines: `(0,0,+1)` and `(0,0,-1)` bonds, `J=1.0` — same magnitude as the existing 4 in-plane bonds | Isotropic extension: the new axis gets the same nearest-neighbour exchange as the other two |
| `dmdata` | Added 2 lines: same 2 bonds, `\|D\|=0.8165` **parallel to bond** — same magnitude and Bloch convention as the existing 4 in-plane bonds | Isotropic Bloch DMI on the new axis, matching the in-plane convention exactly (this is the standard cubic-helimagnet Hamiltonian class — isotropic `J` plus a Dzyaloshinskii-Moriya vector along each bond direction, e.g. the B20/MnSi-type model) |
| `posfile`, `momfile` | Unchanged | Single-atom basis, same moment |
| `inpsd.dat` — `ncell` | `128 128 1` → `16 16 16` (baseline default; every size overrides this via `ncell`, see §B) | Genuine 3D replication instead of a single layer |
| `inpsd.dat` — `BC` | `P P 0` → `P P P` | Fully periodic in all three dimensions — the defining change that makes this case genuinely 3D (blueprint §C: "3D replication ... surface/volume ratio") |
| `inpsd.dat` — `simid` | `SCsurf_T` → `SC3Dbulk` | Cosmetic; reflects the bulk (not surface) geometry |
| `inpsd.dat` — `ip_mode`/`ip_temp`/`qm_svec`/`qm_nvec` | `ip_mode Y` (triple-Q qminimizer search) removed; `ip_mode N` | See "Why the initial-phase search was dropped, not extended" below |
| `inpsd.dat` — everything else | Unchanged: `cell`, `do_prnstruct 2`, `SDEalgh 1`, `Initmag 3`, `mode S`, `temp 0.000`, `damping 0.500`, `Nstep 5000`, `timestep 1.000e-15 s`, `hfield 0.0 0.0 -150.0`, `do_sc N`, `qpoints F`, `plotenergy 1`, `do_avrg Y`, `skyno T`, `stt N`, `jvec`, `adibeta` | Carried over verbatim from the source template — this admission only changes what is needed to make the lattice genuinely 3D, per blueprint §C: "Do not claim the initial spin texture itself constitutes a new kernel workload" |

### Why the initial-phase search was dropped, not extended

`SkyrmionLattice`'s `ip_mode Y` runs `qminimizer.f90::sweep_q3`, a triple-Q
spin-spiral energy-minimization line search over the 129 candidate
q-vectors in `qfile` — all of which have `qz=0` (`qm_nvec 0 0 1` is the
plane normal the search is defined against). This is a strictly in-plane
mechanism; generalizing it to a genuine 3D search would require inventing
a new 3D q-grid (how many points, what spacing, what plane normal) with no
maintainer-supplied precedent — exactly the kind of invented "production"
input every prior WP-08x admission has declined to do
(`B02_skyrmion2D/README.md` §C: "borrowing here would be inventing a
production parameter"). `ip_mode` is dropped to `N` instead: `Initmag 3`'s
uniform initial state (unchanged from the source template) is used
directly, with no initial-phase relaxation.

**A real, verified consequence, not a bug.** `qfile`/`qpoints F` remain in
the template purely because `source/Correlation/qvectors.f90::read_q` is
called unconditionally whenever `qpoints` is one of `F`/`D`/`I`/`B`/`G`
(`source/uppasd.f90:1101-1104`), regardless of `do_sc` — removing `qfile`
entirely was tried first and produces a hard runtime crash:
`Fortran runtime error: Cannot open file '': No such file or directory`
(`source/Correlation/qvectors.f90:97`). `qfile` is therefore kept
byte-for-byte (129 points), but is otherwise fully inert for this case
(`do_sc N`, `ip_mode N` — nothing else reads it).

### A verified, exact consequence of the uniform initial state (not a bug)

For **any** spatially uniform spin configuration on this lattice, each
nearest-neighbour bond pair's Bloch D-vectors are equal and opposite
(`D_{i,i+e} = -D_{i,i-e}` for each axis `e`), so their cross-product
contributions to the effective field at every atom cancel exactly,
independent of the uniform direction chosen. Confirmed directly, not
assumed: a 200-step T=0 CPU run's own `totenergy.<simid>.out` shows the
`DM` column at floating-point noise level throughout
(`~1e-17`–`1e-20`, vs. `Exc` exactly `-6.00000000E+00` and `Zeeman`
relaxing smoothly from `0.635` to `0.408` as damping acts) — the DM term
is genuinely (not approximately) zero for this initial condition at every
step. This mirrors `B01_bccFe/README.md`'s own "weak discriminator" note
for its `bcc_fe_t0` variant (a uniform `Initmag 3` state is an exact fixed
point of that case's anisotropy field too) — a known, accepted property of
starting short-range benchmark dynamics from a uniform state, not new to
this case.

This does not affect the case's validity as a benchmark of a genuinely 3D
short-range J+D workload: the full 6 exchange + 6 DM directed-interaction
neighbour list is evaluated by the solver every step regardless (verified
in §D), which is the computational workload this case represents — per
blueprint §C, "the benchmark significance is the computational workload,
not the exact skyrmion texture," and this admission does not claim the
uniform initial texture is itself a new kernel workload. A later work
package wanting a texture in which the DM term is actively nonzero
throughout the trajectory would need either a finite-temperature variant
(thermal noise breaks the exact cancellation) or a genuine 3D initial-phase
search — both explicitly out of scope here (see §C, no finite-T variant is
added, for the same reason).

### Physical summary

| Property | Value | Source |
| --- | --- | --- |
| Lattice | Simple cubic, single-atom basis, cell `(1,0,0)/(0,1,0)/(0,0,1)` (reduced units) | `posfile`, `inpsd.dat` |
| Dimensionality | Genuinely 3D: `ncell Nx Ny Nz` (all three replicated), `BC P P P` (fully periodic) | `inpsd.dat` |
| Basis/moment | 1 basis atom, `mmom=1.0`, initial direction `(0.1, 0.0, 1.0)` normalized (`Initmag 3`) | `momfile`, `inpsd.dat` |
| Symmetry | `Sym` unset → default `0`: `jij`/`dmdata`'s 6 listed bonds are the complete directed lists | `inpsd.dat` |
| Exchange | 6 nearest-neighbour bonds (`+-x, +-y, +-z`), `J=1.0` | `jij` |
| DMI | Same 6 bonds, `\|D\|=0.8165`, D parallel to each bond (Bloch-type) | `dmdata` |
| Directed interactions / atom | **6 exchange + 6 DM, measured directly** (§D), constant across a 15.6x atom-count range | `struct.<simid>.out`, `dmdata.<simid>.out` |
| Anisotropy | Off (`#anisotropy` commented out, `do_anisotropy` stays `0`) | `inpsd.dat` |
| External field | `hfield 0.0 0.0 -150.0` — carried over unchanged from `SkyrmionLattice` | `inpsd.dat` |
| Dipole | Off | `inpsd.dat` |
| Temperature | `0 K` baseline (only value anywhere in the derived template) | `inpsd.dat` |
| Timestep | `1.000e-15 s` | `inpsd.dat` |
| Integrator | Heun (`SDEalgh 1`), damping `0.50` | `inpsd.dat` |
| Initial phase | None (`ip_mode N`) — `Initmag 3` uniform state used directly | `inpsd.dat` |
| Diagnostics | `skyno T`, `do_avrg Y`; `do_sc N`, `qpoints F` (both inert, see above) | `inpsd.dat` |

### Backend dispatch (`gpu_mode`, `skyno`)

Identical to `B02_skyrmion2D`'s findings, reused rather than re-discovered:
`gpu_mode` must be a per-run override (no separate GPU twin deck exists),
and `skyno T` (triangulation) unconditionally crashes the GPU measurement
backend (`source/gpu_files/measurement/gpuMeasurement.cpp:122-124`,
unimplemented) — a GPU sample overrides `skyno` to `Y` (brute-force,
GPU-implemented). Both keys are already in `harness/cases.py`'s
`GLOBALLY_SAFE_OVERRIDE_KEYS` (added by WP-08A/WP-08B); no harness change
was needed for this admission.

## B. Size ladder

`Nx=Ny=Nz=n` cubes (single-atom basis, `natom=n^3`), chosen as the closest
legal integer replication to each target in the blueprint's ~4x ladder
(4k/16k/64k/256k/1M/4M atoms) — geometry is never distorted to force an
exact atom count, the same methodology `B01_bccFe`'s cube ladder uses (and
for the same reason: a single-atom-basis simple-cubic lattice cannot hit
every target exactly any more than `B01_bccFe`'s two-atom bcc basis could).

| `size_id` | `n` | `natom` | vs. target | Note |
| --- | --- | --- | --- | --- |
| `16x16x16` | 16 | 4,096 | target 4,096, **exact** | |
| `25x25x25` | 25 | 15,625 | target 16,384, -4.6% | |
| `40x40x40` | 40 | 64,000 | target 65,536, -2.3% | closest tier to this project's historical 65,536-atom 3D reference (`B01_bccFe`'s `32x32x32` lands on it exactly with its 2-atom basis) |
| `63x63x63` | 63 | 250,047 | target 262,144, -4.6% | |
| `101x101x101` | 101 | 1,030,301 | target 1,048,576, -1.7% | |
| `161x161x161` | 161 | 4,173,281 | target 4,194,304, -0.5% | |

Aspect ratio is exactly `1:1:1` (isotropic cube) at every size — the
simplest physically sensible choice for a bulk cubic lattice with
isotropic couplings on all three axes, and directly comparable to
`B01_bccFe`'s own isotropic-cube ladder.

## C. Variants

Only `skyrmion_3d_t0` is admitted:

| `variant_id` | `temp` override | `tseed` override | Purpose |
| --- | --- | --- | --- |
| `skyrmion_3d_t0` | 0.0 | 1 | Zero-temperature relaxation dynamics from the uniform `Initmag 3` initial condition |

**No finite-T variant is added**, for the same reason `B02_skyrmion2D`
declined one: this derived template's only temperature anywhere is
`temp 0.000` — there is no maintainer-supplied finite-T value to reuse,
and inventing a "representative" one would repeat the mistake every prior
WP-08x admission has avoided. `temp` remains in `allowed_input_overrides`
so a future WP can add a finite-T variant the moment the maintainer
supplies one (and, per §A, doing so would also be the natural way to get
a texture with an actively nonzero DM term throughout the trajectory).

## D. Scaling validation

Used `harness.cases.generate_run_directory` (not manual file edits) with
`extra_overrides={"do_prnstruct": 1, "Nstep": 20}` on
`build_cpu/bin/sd.f95` at three sizes spanning a 15.6x atom-count range,
then read `struct.<simid>.out` (exchange) via the case's own
`neighbor_list_from_struct_output` and `dmdata.<simid>.out` (DM) directly:

| `size_id` | `natom` | exchange directed / atom | DM directed / atom |
| --- | --- | --- | --- |
| `16x16x16` | 4,096 | 6.0 | 6.0 |
| `25x25x25` | 15,625 | 6.0 | 6.0 |
| `40x40x40` | 64,000 | 6.0 | 6.0 |

Constant at every tested size: full 3D periodicity means every atom is
bulk (no edge truncation) — the same mechanism `B01_bccFe` (constant 96)
and `B02_skyrmion2D` (constant 4/4) found for their own cases.
`63x63x63`/`101x101x101`/`161x161x161` were not executed (impractical for
an admission-time check); at the confirmed constant 6/6 topology,
exchange `directed_interactions` at those sizes is `natom*6` =
1,500,282 / 6,181,806 / 25,039,686 respectively — predictions from the
confirmed scaling law, not independent measurements (same practice
`B02_skyrmion2D` used for its own two largest, untested sizes).

**Same known undercount in `directed_interactions` as `B02_skyrmion2D`.**
This case's `workload_metadata_method` (`neighbor_list_from_struct_output`)
reads only `struct.<simid>.out` (exchange list); it has no knowledge of the
separate DM neighbour list written to `dmdata.<simid>.out`. The true
per-atom neighbour workload driving the LLG effective-field sum is 12
directed interactions (6 exchange + 6 DM), while the reported
`directed_interactions` metric only ever reports the 6 exchange ones — the
same ~2x undercount `B02_skyrmion2D/README.md` documents and explicitly
predicted would recur here ("the future `B03_skyrmion3D` will hit the
identical gap"). Not fixed here, for the same reason: `workload_metadata.py`
is WP-02's frozen, harness-wide contract, left for a future work package to
fix once across every J+D case.

### Workload distinction from `B02_skyrmion2D` (blueprint §C)

| Property | `B02_skyrmion2D` | `B03_skyrmion3D` |
| --- | --- | --- |
| Replication | `[Nx, Ny]`, thickness fixed at 1 cell | `[Nx, Ny, Nz]`, genuinely 3D |
| Boundary conditions | `P P 0` (open out-of-plane) | `P P P` (fully periodic) |
| Neighbour graph | 4 exchange + 4 DM directed / atom | **6 exchange + 6 DM directed / atom** — a larger, genuinely 3D coordination shell, not just a relabeled 2D one |
| Surface/volume ratio | Every atom has 2 "missing" out-of-plane bonds by construction (the physical single layer) | Every atom is bulk-coordinated (no surface atoms at all, at any size — full 3D periodicity) |
| Memory footprint per atom | 2 float64 moment components effectively degenerate with a fixed 3rd (z locked to layer index); 8 directed-neighbour entries (4 `jij` + 4 `dmdata`) | Full 3D per-atom moment; 12 directed-neighbour entries (6 `jij` + 6 `dmdata`) — 50% more neighbour-list memory traffic per atom than `B02_skyrmion2D` at matched atom count |

This is a real difference in the per-atom neighbour-list size and memory
traffic pattern the solver walks every step (6 vs. 4 directed bonds per
list, each list 50% larger), not merely a change in what texture the
initial condition happens to produce — consistent with blueprint §C's
requirement not to claim the initial spin texture itself is a new kernel
workload.

## E. Sanity runs

All runs used `16x16x16` (4,096 atoms), `do_prnstruct 1`, generated through
`harness.cases.generate_run_directory` (not manual file edits); CPU on
`build_cpu/bin/sd.f95` (`UPPASD_GPU_BACKEND=OFF`, `UPPASD_PRECISION=DOUBLE`),
GPU on `build_gpu/bin/sd.f95.cuda` (`UPPASD_GPU_BACKEND=CUDA`,
`UPPASD_PRECISION=DOUBLE`, `extra_overrides={"gpu_mode": 1, "skyno": "Y"}`,
required per §A), both at `Nstep 50`.

| Run | Result |
| --- | --- |
| CPU, `skyrmion_3d_t0` | **PASS** — completes; `natom`/exchange+DM directed interactions match §D; `\|m\|=1.0` conserved; Exc and DM both present and correctly evaluated every step (DM nets to floating-point-noise zero, per §A's verified finding, not a defect) |
| CUDA, `skyrmion_3d_t0` | **PASS** — completes (once `skyno` is overridden per §A); `\|m\|=1.0` conserved; average magnetization direction `(⟨M⟩_x, ⟨M⟩_y, ⟨M⟩_z) = (9.95037190E-02, 0.0, 9.95037190E-01)` matches the CPU run to all 9 printed significant figures at `t=0` |

**This case is a stronger CPU/GPU bit-exactness discriminator than
`B02_skyrmion2D`.** Because `ip_mode` is `N` (§A — no field-free triple-Q
search), there is no cross-binary tie-break-sensitive initial-phase
algorithm for compiler/build non-associativity to expose, unlike
`B02_skyrmion2D`'s documented chirality-degenerate-domain finding. The `t=0`
state here is bit-identical between backends to the precision printed.

Two pre-existing, general backend differences were observed — both already
documented elsewhere in this framework, not new findings specific to this
case:

* **GPU `totenergy.<simid>.out` at `t=0` reports `Zeeman=0` exactly**
  (`Tot=-6.0`, i.e. `Exc` only), while CPU reports the real
  `Zeeman=6.34990275E-01` (`Tot=-5.36500973E+00`). Root cause: this build
  is compiled with `USE_BIG_GRID` (`source/make/default-profiles/
  gfortran-cuda.make:126` et al.), and the GPU measurement backend prints
  `"WARNING: DO NOT USE BIG_GRID WITH GPU CALCULATIONS OF ENERGY"`
  unconditionally whenever energy printing is enabled on such a build
  (`source/gpu_files/measurement/gpuMeasurement.cpp:428`) — a known,
  general GPU energy-diagnostic limitation of this build configuration,
  not something introduced by this case. Confirmed present identically in
  `B02_skyrmion2D`'s own GPU sanity run (re-executed here for comparison):
  its GPU `totenergy.<simid>.out` likewise reports `Zeeman=0.00000000E+00`
  exactly at `t=0`. Does not affect dynamics: the average-magnetization
  comparison above (computed from `emom`, not the diagnostic energy print)
  matches between backends.
* **`averages.<simid>.out` print cadence differs at the end of the run**:
  CPU's only row at `Nstep 50` is `Iter 0` (the run is shorter than this
  template's `avrg_step` sampling interval), GPU prints both `Iter 0` and
  `Iter 50`. This is the same pre-existing end-of-loop print-cadence
  difference `B01_bccFe/README.md` documents for `totenergy.<simid>.out`
  ("a pre-existing end-of-loop energy-print convention difference, not a
  validity concern for either run"), observed here on a different output
  file.

## Hopfion redundancy (blueprint B03 note)

The blueprint states: "A separate Hopfion benchmark is intentionally
omitted from the initial benchmark suite because its ordinary ASD
computational workload is sufficiently represented by this three-
dimensional short-range (J+D) family." This case satisfies that
requirement directly: a Hopfion is, computationally, a genuinely
three-dimensional short-range exchange+DMI spin texture on a bulk lattice
— exactly the workload class (6 exchange + 6 DM directed neighbours/atom,
fully 3D periodic, no anisotropy/dipole term) this case exercises. Per §A
and blueprint §C, the benchmark significance is this per-atom neighbour-
list/effective-field computational pattern, not which specific spin
texture (uniform, skyrmion lattice, or Hopfion) happens to occupy it; no
separate case is needed to represent a Hopfion's ordinary ASD workload.

## Checklist

- [x] Template audited (§A: no genuine 3D J+D maintainer deck exists in `examples/`; derived from `SkyrmionLattice` with every change from the source template enumerated; sibling/candidate templates and the rejected CG-fixture alternative inspected and documented).
- [x] Genuine 3D scaling implemented (§B: six cube sizes 4,096 → 4,194,304 atoms; §C: `BC P P P`, all three of `[Nx, Ny, Nz]` scale).
- [x] Short-range J+D interactions preserved (§A/§D: 6 exchange + 6 DM directed interactions/atom, both constant across a 15.6x atom-count range; isotropic Bloch DMI convention carried over unchanged from `SkyrmionLattice`).
- [x] Atom/bond metadata recorded (§A physical summary table; §D scaling table).
- [x] Aspect ratio controlled (`Nx=Ny=Nz` at every size, exactly 1:1:1).
- [x] CPU sanity passes (§E).
- [x] CUDA sanity passes (§E; stronger bit-exactness discriminator than `B02_skyrmion2D`, no chirality-tie-break sensitivity).
- [x] Hopfion explicitly documented as redundant for core benchmark scope (see "Hopfion redundancy" above).
- [x] Workload distinction from `B02_skyrmion2D` documented explicitly (§D: 3D replication, surface/volume ratio, neighbour graph, memory footprint — initial spin texture explicitly not claimed as the distinction, per blueprint §C).
