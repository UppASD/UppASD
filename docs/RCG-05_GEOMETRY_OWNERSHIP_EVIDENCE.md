# RCG-05 geometry/ownership evidence

**Status: RCG-05A only. Evidence-contract and defect-reproduction slice —
no tracked fixture, generator, comparator, or production fix was added or
changed.** This document is created fresh; no prior RCG-05 evidence
document existed.

**Base commit:** `59e0464b7eab78f9a1a5e4573ac32fbfa506f7cc` ("RCG-04: record
Human closure decision (Anders Bergman, 2026-08-10)"), the exact commit
named by `docs/RCG-05_GEOMETRY_OWNERSHIP_PROMPT_PACK.md` §4 as the required
starting point. `git status --short` at session start showed no staged or
modified tracked files, only pre-existing untracked build directories and
unrelated `ASD_GUI/` generated UI files, none of which are touched by this
slice.

`docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md` (Task RCG-05
entry and the RCG-04 follow-up tasks subsection),
`docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md`, and
`docs/RCG-04_MOVING_E2E_EVIDENCE.md` §17 (RCG-04I) were read in full before
this reproduction was performed. RCG-04 is closed (Human decision,
2026-08-10); its own evidence explicitly stopped short of RCG-05's scope
(`docs/RCG-04_MOVING_E2E_EVIDENCE.md` §17.1).

---

## 1. Scope of this slice

RCG-05A is an evidence-contract and defect-reproduction slice. It:

- reproduces, by direct execution of the real `sd.f95`/`sd.f95.cuda`
  binaries against a disposable, hand-built anisotropic-block-shape input,
  the buffer-width scalarization defect the RCG-05 prompt pack's §1
  inventory identified from source reading alone;
- independently re-confirms, by reading the current source fresh in this
  session, the `dilateAdaptiveState` read/write hazard the same inventory
  flagged;
- inventories the existing reusable evidence infrastructure (workspace
  hash-verification pattern, directed-bond neighbour list, static-topology
  oracle) that RCG-05B/C must build on rather than duplicate;
- defines the canonical ownership-map evidence record schema RCG-05B
  onward must produce.

It does **not** implement a tracked geometry generator (RCG-05B), a
tracked CPU/GPU comparator (RCG-05C), or any production fix (RCG-05D/E).
The reproduction fixture used below is intentionally disposable: it lives
only in this session's scratch directory, was never copied into the
tracked tree, and is not required for any later slice to re-derive its own
tracked fixtures independently.

---

## 2. Buffer-width scalarization: reproduced through the real executable

### 2.1 Disposable fixture

A hand-built `inpsd.dat` plus a custom single-shell `jfile` (kept only in
a scratch directory, never added to the tracked tree), reusing the shared
two-atom `posfile`/`momfile` basis from `tests/coarse_graining/e2e/`:

```text
simid rcg05arepro
alat 1.0
ncell 12 8 8
BC P P P
cell 1.0 0.0 0.0
     0.0 1.0 0.0
     0.0 0.0 1.0
Sym 1
posfile posfile
momfile momfile
exchange jfile
SDEalgh 1
Initmag 1
mseed 92731
Mensemble 1
ip_mode N
mode S
temp 0.0
damping 0.05
Nstep 1
do_avrg N
do_cumu N
do_tottraj N
plotenergy 0
gpu_mode 1
do_gpu Y
do_gpu_llg Y
do_gpu_measurements N
do_gpu_correlations N
do_adaptive_cg Y
block_size_x 1
block_size_y 2
block_size_z 4
cg_operator PROJECTED
cg_mask_mode STATIC
cg_static_mask_file mask.dat
cg_buffer_blocks 0
cg_diagnostics 3
```

`jfile` (single shell, deliberately anisotropic only in the sense that
matters here — the *block shape*, not the bond direction; the CPU formula
maps an isotropic physical interaction radius into per-axis block counts
using the block-vector metric, so any nonzero radius exposes the defect
once `block_size_x/y/z` differ):

```text
1 1 4.0 0.0 0.0 0.01
```

`mask.dat`:

```text
1 FINE
```

`block_size_x/y/z = 1/2/4` gives `block_grid = 12 4 2`
(`ncell 12 8 8` divides exactly, satisfying
`adaptivecgproduction.f90:740-742`'s requirement). By the CPU formula
(`statichybridoperator.f90:180-187`,
`buffer_width_blocks(axis) = ceiling(max(0, fractional_radius(axis) - 64*eps)) + safety_dilation_blocks`,
where `fractional_radius(axis) = radius * |inverse_block_vectors(axis,:)|`
and, for this diagonal/orthogonal block metric,
`inverse_block_vectors(axis,:)` reduces to `1/block_size_axis`;
`safety_dilation_blocks = 0` here since `cg_buffer_blocks 0`), with
`radius = 4.0` (the single shell's displacement norm) the `64*eps` guard
is negligible at this scale:

- x: `4.0 / 1 = 4.0` → `ceil = 4`
- y: `4.0 / 2 = 2.0` → `ceil = 2`
- z: `4.0 / 4 = 1.0` → `ceil = 1`

predicting a genuinely anisotropic `buffer_width_blocks = (4, 2, 1)`, three
distinct values, one per axis.

### 2.2 Temporary diagnostic prints (added, exercised, then reverted)

Two one-line prints were added for this reproduction only and removed
before any commit (confirmed by `git diff --stat` showing no residual
change to either file — see §6):

1. In `source/CoarseGraining/adaptivecgproduction.f90`, inside
   `print_resolved_configuration`, immediately after the existing
   `block_size=`/`block_grid=` diagnostic line: printed
   `adaptive_cg_state%runtime%hybrid%buffer_width_blocks` (the CPU
   three-component vector) and `adaptive_cg_state%gpu_buffer_dilation`
   (the already-scalarized value staged for the GPU descriptor, computed
   at `adaptivecgproduction.f90:615-616` via `int(maxval(...))`).
2. In `source/gpu_files/gpuSimulation.cpp`, immediately after
   `adaptiveSelectorPolicy.bufferDilationBlocks = *FortranData::adaptive_buffer_dilation;`:
   printed `adaptiveSelectorPolicy.bufferDilationBlocks`, the exact
   C++ field the CUDA/HIP `dilateAdaptiveState` kernel launch consumes
   (`gpuAdaptiveRuntime.cpp:328`, `policy.bufferDilationBlocks`).

### 2.3 Fresh builds and runs (see §6 for full commands/environment)

CPU run (`do_gpu N`; fresh out-of-tree `build_rcg05a_cpu`, GNU Fortran
13.3.0, Release, fp64):

```text
AdaptiveCG: block_size=1 2 4  block_grid=12 4 2
RCG05A_TEMP: cpu_buffer_width_blocks=4 2 1  gpu_buffer_dilation_scalar=0
```

(`gpu_buffer_dilation_scalar=0` here only because `gpu_requested` is false
when `do_gpu N`; `allocate_gpu_staging` — and therefore the collapse — is
skipped entirely on this path, confirming the collapse is GPU-staging-only
code, not a universal CPU-side computation.)

GPU run (`do_gpu Y do_gpu_llg Y`, same fixture; fresh out-of-tree
`build_rcg05a_cuda`, `UPPASD_GPU_BACKEND=CUDA`, NVIDIA CUDA 13.3.73,
`CMAKE_CUDA_ARCHITECTURES=native` against an NVIDIA RTX A4000, compute
capability 8.6):

```text
AdaptiveCG: block_size=1 2 4  block_grid=12 4 2
RCG05A_TEMP: cpu_buffer_width_blocks=4 2 1  gpu_buffer_dilation_scalar=4
...
RCG05A_TEMP: gpu_policy.bufferDilationBlocks=4 (scalar)
```

This is the defect, reproduced through the real executable rather than
cited from source reading: the CPU side computes and holds the genuinely
anisotropic `(4, 2, 1)`; the exact same run, on the exact same fixture,
stages `4` (the `maxval`) into the GPU descriptor field consumed by the
dilation kernel — confirmed twice, once at the Fortran staging site and
independently again at the C++ consumption site, so the value observed is
not an artifact of one print but the same scalar traced end to end through
the real GPU path. Per §3.4 of the prompt pack, the reproduction shows
*identity* of what CPU computed being discarded, not merely a difference
in some derived count.

Both runs completed successfully (`AdaptiveCG: GPU lifecycle complete`),
so this is not a setup-rejection artifact — the anisotropic geometry is
accepted and dispatched by production code today, with the GPU side
silently using an isotropic value.

---

## 3. Dilation kernel read/write pattern: re-confirmed from current source

Independently re-read in this session, not assumed from the prompt pack's
own inventory:

`dilateAdaptiveState` (`source/gpu_files/gpuAdaptiveRuntime.cpp:322-350`)
is launched with one GPU thread per block (`target = adaptiveThreadIndex()`).
Each thread that finds `runtime.pendingState[target] == coarseState`
scans a periodic box of neighbouring blocks and, for each neighbour
`source`, reads `runtime.pendingState[source]`; if it finds
`fineState`, it writes `runtime.pendingState[target] = bufferState` and
returns. Both the read (`runtime.pendingState[source]`) and the write
(`runtime.pendingState[target]`) target the same array,
`runtime.pendingState`, within the same kernel launch, across threads
running concurrently with no documented ordering guarantee between them —
confirmed present, unchanged, in the current source at the cited lines.
`proposeAdaptiveState` (`gpuAdaptiveRuntime.cpp:302-320`), the kernel that
populates `pendingState` from `blockState`, is a separate, distinct kernel
launch that (given normal CUDA/HIP stream semantics) completes and
synchronizes before `dilateAdaptiveState` begins, so the specific hazard is
whether one `dilateAdaptiveState` thread's write of `bufferState` into
`pendingState[target]` can be observed as a spurious `fineState`-adjacent
signal by another thread's read of that same slot converted from
`coarseState`, within the same kernel launch — exactly the "no double
buffer and no invariant documented in source" characterization the prompt
pack's §1 item 2 gives. `bufferState`, not `fineState`, is written, so a
naive worst case is a false transitive dilation (a block seeing an
already-diluted neighbour and diluting again) rather than the *wrong*
class being propagated; this is a plausibility observation from reading
the state values, not proof of correctness or of a benign bound — RCG-05E
is the slice authorized to determine actual safety with sanitizer
evidence, and this session performs no such run.

No production or CPU-analogue equivalent of this same-launch self-referential
read/write exists: `rebuild_static_hybrid_ownership`
(`source/CoarseGraining/statichybridoperator.f90:198-256`) computes each
block's `atomistic_block` state from the immutable `fine_mask` input and
the geometry-only `buffer_width_blocks`, never reading another block's
freshly-written output within the same rebuild — confirmed by re-reading
lines 232-248 in this session (reproduced in full at
`docs/RCG-05_GEOMETRY_OWNERSHIP_PROMPT_PACK.md` §1 item 2 and verified
unchanged here).

---

## 4. Existing reusable infrastructure inventoried (not duplicated)

Read in full this session before any new evidence was written, per the
prompt pack §3.1's explicit instruction not to write a second ownership-map
dumper or duplicate mechanism:

- **`tests/coarse_graining/static_topology_oracle.py`** (217 lines):
  independently reproduces the CPU ownership algorithm (`buffer_width_blocks`
  formula, periodic Chebyshev-box dilation, block numbering) from first
  principles, but is explicitly restricted at lines 23-28
  (`UnsupportedTopologyError`/`NotImplementedError`) to
  `block_grid_y == block_grid_z == 1`. Confirmed unchanged at this commit.
  This is RCG-05B's generalization target, not something to reimplement.
- **`tests/coarse_graining/torque_oracle.py`** (686 lines):
  `build_geometric_bonds` (line 213) is the directed-bond neighbour-list
  builder RCG-05C/F must reuse for cross-interface bond coverage, rather
  than re-deriving bond geometry independently.
- **`tests/coarse_graining/coarse_torque_oracle.py`** (209 lines):
  block-level averaging/stiffness helpers for the coarse side; reusable by
  RCG-05F's transition-invariant checks.
- **`tests/coarse_graining/run_moving_backend_parity.py`**
  (`prepare_workspace`, lines 215-249; `sha256_file`, lines 211-212): the
  per-backend isolated-workspace pattern (copy the full `e2e` tree, verify
  every non-GPU-appended file is byte-identical by SHA-256 before any run)
  that RCG-05C's own CPU/GPU comparator workspace must reuse rather than
  reinvent.
- **`tests/coarse_graining/trajectory_evidence.py`**: transition-history
  parsing (`parse_mag_conf_text` family) RCG-05F should reuse for
  resolution-transition ownership-invariant checks, per the prompt pack's
  explicit pointer to RCG-04G's infrastructure.
- **`tests/coarse_graining/moving_state_generator.py`**'s `Geometry`/
  manifest/provenance dataclass pattern (`_manifest`, `manifest_json`,
  `content_hash`): the convention RCG-05B's new skew/unequal-width
  generators must follow for deterministic, tracked provenance.

No new dumper, comparator, or oracle was written in this slice; the above
is a reading-based inventory only, matching the RCG-05A prompt's explicit
scope limit.

---

## 5. Canonical ownership-map evidence record (schema for RCG-05B onward)

Later slices must produce a machine-readable evidence record per compared
geometry (JSON/YAML sidecar or an equivalently structured, parseable stdout
block) containing at minimum:

1. **Provenance:** source commit, `git status --short`, generator
   name/version and parameters (RCG-05B), input content hash (jfile/
   posfile/momfile/mask, or generator-produced equivalents),
   configure/build/run commands, compiler, backend, device, precision —
   mirroring `docs/RCG-04_MOVING_E2E_EVIDENCE.md` §5 item 1 exactly, since
   this is the same provenance contract applied to a different observable.
2. **Geometry identity:** full cell vectors (not just an
   orthogonal/skew flag — the actual 3x3 matrix, since RCG-05B must cover
   non-orthogonal cells), `ncell`/`na`, `block_size_x/y/z`, the derived
   `block_grid` (nb1, nb2, nb3), and periodicity (`BC` per axis).
3. **Per-block ownership record:** for every block id (1-based, matching
   `block_id = 1 + bx + grid_x*(by + grid_y*bz)`), its
   `block_grid_coordinate` triple and ownership class
   (`FINE`/`BUFFER`/`COARSE`, matching
   `STATIC_HYBRID_FINE`/`BUFFER`/`COARSE`), plus — where meaningful — the
   per-atom `atomistic_atom` flag and owning block id for every atom.
4. **Dipole-all-grid ownership tag:** an explicit, separate boolean/enum
   field (not inferred from the fine/buffer/coarse class) recording that a
   given atom's uniform FFT dipole contribution is owned independently of
   its mask class, per `coarsetensoroperator.f90:345-346` and
   `statichybridoperator.f90:18/115-118` — so a comparator can assert
   "included exactly once" without conflating dipole ownership with the
   mask-derived classes.
5. **Source tag:** which backend produced this specific record instance
   (`CPU`/`CUDA`/`HIP`), with that backend's own build identity (`git
   describe`, compiler, GPU device/architecture, precision) — since
   RCG-05C's whole purpose is comparing two *different* records of the
   same geometry, one per backend.
6. **Buffer-width identity used to produce this record:** the CPU
   three-component `buffer_width_blocks(3)` and, on a GPU-sourced record,
   the scalar `bufferDilationBlocks` actually staged — so that any
   observed mismatch is directly attributable to the exact scalarization
   defect this slice reproduced, both before RCG-05D's fix (expected
   mismatch) and after (expected match).
7. **Comparator verdict fields (RCG-05C onward):** per-block/per-atom
   match or mismatch (index-set identity, never a bare count), periodic-
   wrapping coverage confirmed per axis, and cross-interface bond coverage
   confirmed via the reused `torque_oracle.build_geometric_bonds` list —
   per the prompt pack §3.3's explicit rule against aggregate-only
   ownership comparison.
8. **Regeneration provenance:** enough of the above (generator name/
   version/parameters plus content hash) that the exact same fixture can
   be regenerated byte-identically, matching
   `docs/RCG-04_MOVING_E2E_EVIDENCE.md` §5's existing regeneration
   requirement.

This schema is a specification for RCG-05B/C to implement, not new code
introduced here — consistent with how `docs/RCG-04_MOVING_E2E_EVIDENCE.md`
§5 was written by RCG-04A for RCG-04B/C to implement.

---

## 6. Fresh build/test evidence (this slice)

**Environment:** GNU Fortran 13.3.0, GNU C/C++ 12.4.0, CMake 3.28.3,
NVIDIA CUDA 13.3.73 (`nvcc`), NVIDIA RTX A4000 (compute capability 8.6,
`CMAKE_CUDA_ARCHITECTURES=native`), Release build type, fp64 (default
precision). No HIP toolchain is present on this host (`hipcc` not found),
matching RCG-04-FU1's still-open deferral.

**Fresh out-of-tree CPU configure/build (with the temporary diagnostic
prints from §2.2 present):**

```text
$ cmake -DCMAKE_BUILD_TYPE=Release -DUPPASD_GPU_BACKEND=OFF \
    -S . -B build_rcg05a_cpu
...
-- Git tag found: VERSION="v6.0.2-459-g59e0464b".
-- Output binary:  .../build_rcg05a_cpu/bin/sd.f95
-- Configuring done
$ cmake --build build_rcg05a_cpu -j2
... exit 0
```

**Disposable-fixture run (CPU, `do_gpu N`):** see §2.3 for full annotated
output; `AdaptiveCG: GPU lifecycle complete` not reached (no GPU path
requested), `RCG05A_TEMP: cpu_buffer_width_blocks=4 2 1
gpu_buffer_dilation_scalar=0` recorded.

**Fresh out-of-tree CUDA configure/build (with the temporary diagnostic
prints present):**

```text
$ cmake -DCMAKE_BUILD_TYPE=Release -DUPPASD_GPU_BACKEND=CUDA \
    -S . -B build_rcg05a_cuda
...
-- The CUDA compiler identification is NVIDIA 13.3.73
-- CMAKE_CUDA_ARCHITECTURES: native
-- Output binary:  .../build_rcg05a_cuda/bin/sd.f95.cuda
-- Configuring done
$ cmake --build build_rcg05a_cuda -j2
... exit 0
```

**Disposable-fixture run (CUDA, `do_gpu Y do_gpu_llg Y`):** see §2.3;
`RCG05A_TEMP: cpu_buffer_width_blocks=4 2 1  gpu_buffer_dilation_scalar=4`
and `RCG05A_TEMP: gpu_policy.bufferDilationBlocks=4 (scalar)` both
recorded from the real GPU path; run completed
(`AdaptiveCG: GPU lifecycle complete`).

**Diagnostic prints reverted; confirmed clean:**

```text
$ git diff --stat -- source/CoarseGraining/adaptivecgproduction.f90 \
    source/gpu_files/gpuSimulation.cpp
(no output)
```

**Baseline regression evidence, against the reverted source** (full
`make -j2` rebuild picked up the revert incrementally in both trees before
these runs):

```text
$ ctest --test-dir build_rcg05a_cpu -L cg13-cpu --output-on-failure
...
100% tests passed, 0 tests failed out of 22
```

```text
$ ctest --test-dir build_rcg05a_cuda -L cg13-cuda --output-on-failure
...
100% tests passed, 0 tests failed out of 21
(includes adaptive-cg-moving-backend-parity ... Passed 76.98 sec)
```

**Worktree check after the runs:**

```text
$ git status --short --porcelain=v1 | grep -v '^??'
(no output after restoring test-run byproducts, see below)
```

`adaptive-cg-production-e2e` (run as part of both `cg13-cpu`/`cg13-cuda`
labels) runs the real binary against tracked
`examples/AdaptiveCoarseGraining/*` inputs, which receive their own
provenance-stamped `uppasd.*.yaml` output (`date`/`git_revision` fields
only) as a side effect, exactly as `docs/RCG-04_MOVING_E2E_EVIDENCE.md` §8
already documented for RCG-04A; these were restored with `git checkout --`
after the runs, confirmed by the empty `git status --short` above other
than pre-existing untracked build directories and unrelated `ASD_GUI/`
files. The disposable reproduction fixture (`inpsd.dat`, `jfile`,
`mask.dat`, copied `posfile`/`momfile`) lived only in this session's
scratch directory, never inside the tracked tree, and no tracked or
untracked repository file was created, modified, or deleted by this slice
other than this new evidence document.

---

## 7. RCG-05A checklist

- [x] Buffer-width scalarization is reproduced through the real executable, not only cited from source reading.
- [x] A deliberately anisotropic block shape is used, with the per-axis CPU widths shown to differ.
- [x] The GPU descriptor's collapse to a single scalar is shown from actual GPU-path output.
- [x] The dilation kernel's read/write pattern is independently re-confirmed by reading the current source (not assumed from this document).
- [x] The canonical ownership-map evidence schema is documented.
- [x] Existing reusable infrastructure (`static_topology_oracle.py`, `torque_oracle.py`/`coarse_torque_oracle.py`, `run_moving_backend_parity.py`'s hash-verification pattern) is inventoried for reuse, not duplicated.
- [x] No permanent fixture, generator, or production fix is introduced in this slice.
- [x] Fresh configure/build/test commands and results are recorded.
- [x] Unrelated worktree changes remain untouched and unstaged.

## Open items (carried forward, not blocking RCG-05B)

- The dilation kernel's read/write hazard is re-confirmed present by
  source reading only; whether it is an actual race (vs. benign under
  CUDA/HIP's specific scheduling) is explicitly **not** determined by this
  slice — that determination, with sanitizer evidence, is RCG-05E's own
  exit criterion, per the prompt pack's governing rule against clearing a
  race "based on source reading alone."
- No skew (non-orthogonal) cell was exercised in this slice's disposable
  reproduction; the fixture used is orthogonal-only. RCG-05B is
  responsible for both the unequal-width orthogonal and skew-cell cases.
- HIP was not exercised (no toolchain on this host); recorded as a
  deferral, not a pass, matching RCG-04-FU1.

---

## 8. RCG-05B: skew and unequal-width geometry generators

**Base commit:** `d190a169f6052b77c73184e259e156620f989a55` ("RCG-05A:
reproduce buffer-width scalarization and define geometry evidence
contract"), the accepted RCG-05A commit. `git status --short` at session
start showed no modified tracked files.

### 8.1 Inspection of existing conventions (before implementation)

Read before writing any code, per the prompt pack's explicit instruction
not to write a second ownership-map/geometry mechanism:

- `tests/coarse_graining/moving_state_generator.py` in full — the
  `Geometry`/manifest/provenance dataclass pattern RCG-04B established.
  `Geometry.cartesian()` was hard-restricted to the identity cell; nothing
  else in that module assumes identity beyond that one method.
- `tests/coarse_graining/static_topology_oracle.py` in full — the exact 1D
  restriction (`block_grid_y == block_grid_z == 1`, `UnsupportedTopologyError`)
  and its algorithm, confirmed still present and unchanged at this commit.
- `source/CoarseGraining/blocktopology.f90:build_block_topology` (lines
  134-437) — the authoritative 3D block-numbering algorithm
  (`regular_spatial_block_id`,
  `block_id = 1 + coord(1) + grid(1)*(coord(2)+grid(2)*coord(3))`),
  atom-to-block floor division on every axis independently, and
  `block_vectors(:,axis) = block_shape(axis)*cell_vectors(:,axis)` for a
  general (including skew) cell — this is the metric the buffer-width
  formula must actually invert, not a per-axis shortcut.
- `source/CoarseGraining/statichybridoperator.f90:180-187` and its
  `inverse3` (lines 414-430) — the exact per-axis buffer-width formula
  (`radius * |inverse3(block_vectors)(axis,:)|`) already partially
  transcribed by RCG-05A's inventory and independently re-confirmed here by
  direct reading, term-for-term, before writing `_inverse3` in Python.
- `source/Input/inputdata.f90:335` (`posfiletype = 'C'` default) and
  `source/Input/inputhandler_ext.f90:1079-1080` (`posfiletype=='C'` ->
  `r_red=r_tmp`, no cell-vector conversion of the basis) — confirms a
  fixture's basis offsets are valid unchanged under any `cell_vectors`,
  skew or not, so no new basis/posfile generator was needed; the existing
  shared `tests/coarse_graining/e2e/posfile` basis is reused directly.
- `source/Hamiltonian/neighbourmap.f90:266-290` — production's own
  periodic-image convention for matching a declared jfile bond to a real
  neighbour atom (`dfrac = dcart @ inverse(cell_matrix)`, `itrans =
  nint(dfrac)`, `rescart = dcart - itrans @ cell_matrix`) — the convention
  `_lattice_minimum_image` reproduces exactly, term-for-term, rather than
  inventing an equivalent one.
- `tests/coarse_graining/torque_oracle.py` in full — `build_geometric_bonds`
  and `_minimum_image` are calibrated against real production exchange-energy
  output (module docstring), but restricted to a plain per-axis orthogonal
  box wrap; not modified (regression risk on RCG-04's calibration), reused
  unchanged for every orthogonal-cell case, and used as the regression
  target for the new general bond builder (§8.3).
- `tests/coarse_graining/run_moving_static_mixed.py`,
  `run_moving_dmi_chiral.py`, `run_moving_all_coarse.py` — grepped for
  every real call site of `compute_expected_topology`/`interface_bond_count`
  before changing either function's signature. Every call site passes
  `block_size_y=2, block_size_z=2` (matching `GEOMETRY.n2=n3=2` exactly, the
  degenerate 1D case) and no `cell_vectors`, and depends on the
  `.nblocks_x` attribute name directly (`run_moving_static_mixed.py:197`).
  This fixed the compatibility contract the generalization below had to
  preserve exactly: same signature plus one new optional `cell_vectors`
  keyword (default identity), same `.nblocks_x` attribute name (now holding
  the total block count, numerically identical to the old per-axis count
  whenever `grid_y == grid_z == 1`).
- `tests/coarse_graining/e2e/*/inpsd.dat` — grepped for every distinct
  `(ncell, block_size_x/y/z)` combination any tracked fixture actually uses.
  Every one uses `ncell _ 2 2` / `block_size_y=block_size_z=2`; the distinct
  `(ncell_x, block_size_x)` pairs are `(10,1) (24,1) (24,2) (24,4) (24,8)
  (6,1) (6,2)`, plus one deliberate negative control
  (`invalid_partial_block`, `ncell 5`/`block_size_x 2`, non-divisible).

### 8.2 What was implemented

**`tests/coarse_graining/moving_state_generator.py`** (extended, not
replaced): `Geometry` gained an optional `cell_vectors` field (default the
identity cell every prior caller implicitly used), and `cartesian()` was
generalized to the full `i1*C1 + i2*C2 + i3*C3 + Bas(i0)` formula
(`source/System/geometry.f90:445`), reducing exactly to the old
`(i1+bx, i2+by, i3+bz)` for the default. No RCG-04B generator function
(`conical_spiral_state` etc.) passes a non-default `cell_vectors`, so none
of RCG-04's own fixtures are affected — checked by re-running
`test_moving_state_generator.py`'s full 40-test suite unmodified (§8.4).

**`tests/coarse_graining/static_topology_oracle.py`** (generalized):

- `_inverse3` — a general 3x3 matrix inverse, transcribed term-for-term
  from Fortran's `inverse3` (three independent copies in
  `statichybridoperator.f90`/`coarsetensoroperator.f90`/
  `multichannelcoarsetensoroperator.f90`, confirmed identical to each
  other before transcribing one), raising a new `DegenerateGeometryError`
  for a singular (zero-determinant) matrix rather than propagating a
  divide-by-zero.
- `buffer_width_blocks(shells, cell_vectors, block_shape, cg_buffer_blocks)`
  — the general per-axis dilation width for an arbitrary (including skew)
  metric. `buffer_width_blocks_x` (every existing caller's entry point) is
  now defined *in terms of* this general function on an identity cell with
  `block_shape=(block_size_x, 1, 1)` — a regression check baked into the
  implementation itself, not merely asserted by a separate test.
- `compute_expected_topology` — same name, same required-argument
  signature, plus one new optional `cell_vectors` keyword (default
  identity). Internals rewritten to a genuine 3D block grid: `_block_id`/
  `_block_coord` implement `regular_spatial_block_id`'s exact numbering and
  its inverse; the divisibility check, the FINE-seed periodic dilation, the
  atom-to-block floor division, and `distance_from_boundary_blocks` (now a
  true 3D Chebyshev/L-infinity distance in grid-index space, reducing to
  the old scalar distance when only one axis is nontrivial) all generalize
  to all three axes. `UnsupportedTopologyError` is no longer raised (kept,
  empty, only so `sto.UnsupportedTopologyError` remains importable).
  `ExpectedTopology.nblocks_x` keeps its exact name (every real caller
  depends on it) but documents that it now holds the *total* block count;
  `.nblocks`/`.buffer_width_x` properties were added as clearer aliases
  without removing the old names.
- `_lattice_minimum_image` — the general periodic-image reduction
  (fractional coordinates via `_inverse3`, round-to-nearest integer
  translation, reduce), matching `neighbourmap.f90`'s own convention
  exactly. `build_geometric_bonds_general`/`interface_bond_count_general`
  generalize `torque_oracle.build_geometric_bonds`/`interface_bond_count`
  for a non-orthogonal cell, wrapping over the *simulation box*
  (`geometry.n1/n2/n3` replicas of `cell_vectors`) exactly as
  `torque_oracle._minimum_image`'s calibrated `box=(n1*alat,n2*alat,
  n3*alat)` convention does — an implementation bug (wrapping over the bare
  unit cell instead) was caught by the regression test in §8.3 before this
  document was written; see that section for what the bug looked like and
  how it was found and fixed. `torque_oracle.py` itself was not modified.
- `unequal_width_orthogonal_fixture`/`skew_cell_fixture` — the new
  geometry generators, matching the `moving_state_generator._manifest`
  provenance-record shape (reimplemented locally, not imported, only
  because this module's own `GENERATOR_VERSION` must not be shadowed by
  RCG-04B's; `content_hash` itself is imported and reused directly). Each
  generator is *self-validating*: it calls `compute_expected_topology` on
  its own chosen parameters and raises `ValueError` unless every one of
  fine/buffer/coarse is nonempty, rather than asserting the host is "large
  enough" without checking. `skew_cell_fixture` additionally rejects a
  `cell_vectors` that is not actually non-orthogonal (every pairwise dot
  product ~0), so it cannot silently degrade to an orthogonal fixture that
  merely calls itself "skew". Both reuse the existing shared 2-atom
  `e2e/posfile` basis (§8.1) rather than generating a new one.

Per the prompt's explicit scope limit, neither generator's output was
written into a tracked `e2e/` fixture directory, registered with CTest, or
run through the real executable in this slice — RCG-05B produces
generator *functions*, self-validated by this module's own oracle, exactly
as RCG-04B produced moving-state generator functions without materializing
fixture directories itself (RCG-04B's one exception, `initmag_restart_atomistic`,
was triggered by an unplanned production capability-gap fix, not a general
requirement; no such gap was found in this slice).

### 8.3 A bug found and fixed by this slice's own regression testing

While writing `GeneralBondBuilderRegressionTests` (§8.4), the first version
of `build_geometric_bonds_general` reduced every displacement using the
*bare unit cell* (`cell_vectors` as passed in), not the simulation box. On
`GEOMETRY` (`na=2, n1=24, n2=2, n3=2`, identity cell) this produced 4480
matched bond-endpoint pairs where `torque_oracle.build_geometric_bonds`
(the real, production-calibrated function) finds only 176 — every
single-unit-cell-step image was being treated as a distinct "nearby" atom
instead of being folded back into the finite periodic box. The fix scales
`cell_vectors` by `(geometry.n1, geometry.n2, geometry.n3)` before calling
`_lattice_minimum_image`, exactly matching `torque_oracle._minimum_image`'s
own `box=(n1*alat, n2*alat, n3*alat)` convention. After the fix,
`build_geometric_bonds_general(GEOMETRY, SHELLS, _IDENTITY_CELL) ==
torque_oracle.build_geometric_bonds(GEOMETRY, SHELLS)` exactly (dict
equality, not just equal counts) — this is now a permanent regression test
(`test_matches_orthogonal_build_geometric_bonds_exactly`), and is stronger
evidence than a small synthetic case since `GEOMETRY`/`SHELLS` here are the
exact fixture `torque_oracle.py`'s own docstring calibrated against real
production exchange-energy output.

### 8.4 Tests (`tests/coarse_graining/test_static_topology_oracle.py`)

Every one of the original 19 tests (`BufferWidthTests`,
`ExpectedTopologyTests`, `MaskFileTests`, `InterfaceBondCountTests`) still
passes **unmodified**, except the one test that specifically exercised the
now-removed 1D restriction
(`test_rejects_block_grid_not_spanning_yz`), which is replaced by
`test_genuine_3d_block_grid_is_now_supported` (same call, now asserting the
real, hand-derived 2D partition it produces: `block_grid=(24,2,1)`,
`width=(2,2,1)`, block 2/3/25 all `BUFFER`, hand-derived from
`block_vectors=diag(1,1,2)` before running any code) and one new negative
control for the axis-y divisibility check. 26 new tests were added:

- `GeneralBufferWidthTests` — the orthogonal anisotropic case reproduces
  RCG-05A's own real-hardware result (`block_size=(1,2,4)`, radius 4.0 ->
  `(4,2,1)`, the exact numbers CPU and CUDA printed in
  `docs/RCG-05_GEOMETRY_OWNERSHIP_EVIDENCE.md` SS2.3); the skew case
  (`cell=((1,0,0),(0.5,1,0),(0,0,1))`, radius=|C2|=1.118034) uses widths
  independently computed by a standalone script *before* this module's
  generalization was written (`(1,2,2)->(2,1,1)`, `(2,2,2)->(1,1,1)`);
  `buffer_width_blocks_x`'s equivalence to the general formula is asserted
  directly; malformed/singular-cell inputs are rejected.
- `Full3DOwnershipTests` — a genuine `grid=(5,5,5)` case (`na=1`, single
  distance-1.0 shell, `block_size=(1,1,1)` -> `width=(1,1,1)`), hand-derived
  (not printed from the module under test): a 3x3x3=27-block box (1 FINE +
  26 BUFFER) out of 125, 98 COARSE; six specific block states checked by
  coordinate, plus the atom-to-block floor-division and the 3D Chebyshev
  `distance_from_boundary_blocks` generalization.
- `SkewCellOwnershipTests` — the same skew cell at `block_size=(1,2,2)`,
  `ncell=(8,8,8)`: `block_counts()=={"fine":1,"interface":44,"coarse":83}`,
  matching a standalone independent computation performed before writing
  this test.
- `LatticeMinimumImageTests` — two hand-derived skew reductions
  (`1*C2 + residual`, `2*C1 - C3 + residual`, each reducing to exactly the
  residual), an identity-cell sanity check, and a singular-cell rejection.
- `GeneralBondBuilderRegressionTests` — the bug-and-fix described in §8.3,
  now a permanent exact-equality regression against
  `torque_oracle.build_geometric_bonds`/`interface_bond_count`.
- `RealFixtureGeometryRegressionTests` — every real `(ncell_x,
  block_size_x)` combination from §8.1's grep, run through
  `compute_expected_topology` and checked for exact equality (block grid,
  width, full `block_state_by_id` dict) against
  `_reference_1d_topology`, an independent transcription of the *original*
  (removed) 1D algorithm written directly into the test file rather than
  calling anything from `static_topology_oracle` — an oracle-of-the-oracle.
  Also confirms this file's `JFILE_TEXT` constant is byte-identical to the
  tracked `e2e/jfile` every one of these real fixtures actually reads.
- `BlockGeometryGeneratorTests` — both generators produce all-nonempty
  partitions matching their independently-computed widths; deterministic
  regeneration (byte-identical `jfile_text`/`mask_text`, identical
  `content_sha256`, byte-identical `manifest_json`); a changed parameter
  changes the hash; malformed inputs (isotropic block size for the
  "unequal-width" fixture, non-divisible `ncell`, an orthogonal
  `cell_vectors` for the "skew" fixture, a host too small to leave every
  state nonempty) all raise `ValueError` with a specific message.

Total: 45 tests (19 original, unmodified except the one repurposed
restriction test and one added negative control, plus 26 new), all
passing, run both directly (`python3 test_static_topology_oracle.py -v`,
from both `tests/coarse_graining/` and a build directory, after a
CTest-only path bug was found and fixed -- see §8.5) and via the existing
CTest target `coarse-graining-static-topology-oracle`.

### 8.5 A second bug found by fresh out-of-tree CTest evidence

`RealFixtureGeometryRegressionTests` opened `e2e/jfile` with a bare
relative path, which only resolves when the interpreter's current working
directory happens to be `tests/coarse_graining/` -- true when the test is
run directly from that directory, false under `ctest`, whose working
directory is the build tree. This was caught by running the actual
CTest target from a fresh build (§8.6), not just direct `python3 ...`
invocation from the source directory, exactly the reason the governing
rules require CTest-based evidence rather than only direct script runs.
Fixed by resolving the path relative to `Path(__file__)` (`E2E_ROOT =
Path(__file__).with_name("e2e")`), the same convention
`run_moving_static_mixed.py`'s `ROOT = Path(__file__).with_name("e2e")`
already uses -- reused, not invented. Re-verified passing both from
`tests/coarse_graining/` and from the build directory directly (§8.6).

### 8.6 Fresh build/test evidence (this slice)

**Environment:** GNU Fortran 13.3.0, GNU C/C++ 12.4.0, CMake 3.28.3, CPU
backend, fp64 (default precision), Release build type. No CUDA/HIP
evidence was gathered in this slice (not required: RCG-05B produces no
GPU-path fixture or claim, per the prompt's explicit "do not yet run these
fixtures through the real executable's GPU path" scope limit).

```text
$ python3 tests/coarse_graining/test_static_topology_oracle.py -v
...
Ran 45 tests in 2.7s
OK

$ python3 tests/coarse_graining/test_moving_state_generator.py -v
...
Ran 40 tests in 0.03s
OK

$ python3 tests/coarse_graining/test_torque_oracle.py -v
...
Ran 29 tests in 0.11s
OK

$ python3 tests/coarse_graining/audit_fixture_dependencies.py
adaptive-CG fixture dependency audit: PASS (58 fixture directories, 118 input paths)
```

(58 directories / 118 input paths -- identical to the count before this
slice, confirming no new tracked fixture or dangling reference was
introduced, consistent with RCG-05B producing generator functions only.)

**Fresh out-of-tree configure/build:**

```text
$ cmake -DCMAKE_BUILD_TYPE=Release -DUPPASD_GPU_BACKEND=OFF \
    -S . -B build_rcg05b_cpu
...
-- Git tag found: VERSION="v6.0.2-460-gd190-dirty".   # dirty solely from
                                                        # this slice's own
                                                        # uncommitted
                                                        # changes
-- Configuring done
$ cmake --build build_rcg05b_cpu -j2
... exit 0
```

**Test run (`cg13-cpu` label, full production regression):**

```text
$ ctest --test-dir build_rcg05b_cpu -L cg13-cpu --output-on-failure
...
100% tests passed, 0 tests failed out of 22
 3: coarse-graining-block-topology ................... Passed
 ...
17: coarse-graining-static-topology-oracle ........... Passed   # this
                                                                  # slice's
                                                                  # 45 tests
18: adaptive-cg-production-e2e ....................... Passed
19: adaptive-cg-setup-rejection-matrix ............... Passed
20: adaptive-cg-moving-off-fine ...................... Passed
21: adaptive-cg-moving-all-coarse .................... Passed
22: adaptive-cg-moving-static-mixed .................. Passed
23: adaptive-cg-moving-adaptive-wall ................. Passed
24: adaptive-cg-moving-dmi-chiral .................... Passed
```

Tests 21/22/24 (`adaptive-cg-moving-all-coarse`/`-static-mixed`/
`-dmi-chiral`) are the real production harness scripts identified in §8.1
as live callers of `compute_expected_topology`/`interface_bond_count`;
their unmodified pass, through the real `sd.f95` binary, is the decisive
regression evidence that the generalization did not change behaviour for
any existing RCG-04 caller (a stronger check than re-running this file's
own unit tests alone, since those callers run the actual executable and
compare against its actual diagnostic output, not just this module's
self-consistency).

**Worktree check after the run:**

```text
$ git status --short --porcelain=v1 | grep -v '^??'
 M examples/AdaptiveCoarseGraining/adaptive/uppasd.adaptive.yaml         # test byproduct, restored below
 M examples/AdaptiveCoarseGraining/initial_phase_texture/uppasd.adaptive.yaml
 M examples/AdaptiveCoarseGraining/static_mixed/uppasd.adaptive.yaml
 M tests/coarse_graining/moving_state_generator.py
 M tests/coarse_graining/static_topology_oracle.py
 M tests/coarse_graining/test_static_topology_oracle.py
$ git checkout -- examples/AdaptiveCoarseGraining/adaptive/uppasd.adaptive.yaml \
    examples/AdaptiveCoarseGraining/initial_phase_texture/uppasd.adaptive.yaml \
    examples/AdaptiveCoarseGraining/static_mixed/uppasd.adaptive.yaml
$ git status --short --porcelain=v1 | grep -v '^??'
 M tests/coarse_graining/moving_state_generator.py
 M tests/coarse_graining/static_topology_oracle.py
 M tests/coarse_graining/test_static_topology_oracle.py
```

The three `uppasd.*.yaml` diffs were only `date`/`git_revision` provenance
stamps (identical pattern to RCG-05A's own evidence §6), confirmed before
restoring. No other tracked or untracked file was created, modified, or
deleted by this slice other than this evidence document and the three
files listed above.

### 8.7 RCG-05B checklist

- [x] Existing RCG-04B generator conventions were inspected and reused, not duplicated. (SS8.1, SS8.2)
- [x] A deterministic unequal-block-width orthogonal generator is implemented. (SS8.2, `unequal_width_orthogonal_fixture`)
- [x] A deterministic skew-cell (non-orthogonal lattice vector) generator is implemented. (SS8.2, `skew_cell_fixture`)
- [x] Atom ordering, basis/material identity, and moment magnitudes are preserved. (SS8.1: reuses the existing shared `e2e/posfile` basis unchanged; no momfile/state claim is made by this slice)
- [x] `static_topology_oracle.py`'s block-grid restriction is removed and replaced with a general 3D formula. (SS8.2)
- [x] The generalized oracle is regression-tested against every existing RCG-04 fixture geometry. (SS8.4, `RealFixtureGeometryRegressionTests`; SS8.6, real production harnesses)
- [x] The generalized oracle is tested against at least one hand-verified anisotropic/skew case. (SS8.4, `Full3DOwnershipTests`, `SkewCellOwnershipTests`, `LatticeMinimumImageTests`)
- [x] Generator parameters and provenance are tracked in a manifest. (SS8.2, `_manifest`/`content_hash`)
- [x] Malformed or degenerate generator requests fail clearly. (SS8.4, `BlockGeometryGeneratorTests` malformed-input cases)
- [x] Tracked-fixture/package audit covers the new generators and required inputs. (SS8.6: `audit_fixture_dependencies.py` unaffected, 58/118 unchanged, consistent with no new tracked fixture directory)
- [x] No CPU/GPU equivalence or production-fix claim is made in this slice. (SS8.2 explicit scope note)
- [x] Unrelated worktree changes remain untouched and unstaged. (SS8.6)

---

Shall I create the focused RCG-05B commit with the one-line message
`RCG-05B: add skew and unequal-width geometry generators`?
