# RCG-03 polarization gate and anisotropy capability audit — evidence

**Status: RCG-03 is closed (2026-08-08, Human decision: Anders Bergman).**
Sections 1–14 below (through "Open limitations") are the original
2026-08-06 exploratory session, performed while RCG-02 was still open, and
are retained unchanged as the historical record of that work. The
"RCG-03 closure (2026-08-08)" section at the end of this document records
the rebase onto the now-closed RCG-02, the three previously-open checklist
items being delivered, and the Human closure decision including the
HIP/independent-review deferral. Read that final section for the current
status; do not treat the paragraph below (written 2026-08-06) as current.

**Original 2026-08-06 status (superseded, kept for history):** exploratory
work performed while RCG-02 remains open (see
`docs/RCG-02_DMI_HANDEDNESS_EVIDENCE.md`). Per
`docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md` section 8/10, a
later task may be implemented and tested ahead of its dependency closing,
but it is not accepted until the dependency closes and this work is
rebased onto that accepted commit and rerun. **Eleven of the fourteen
RCG-03 checklist items in the remediation blueprint are ticked, each on
its own individually demonstrated evidence; the remaining three are left
open (see that document's evidence paragraph for exactly which and why).
Ticking those items is not a claim that RCG-03 as a whole is closed** —
closure still requires RCG-02, independent review, HIP evidence, and the
production `ANI-NONUNIFORM-REJECT` fixture.

**Base commit:** `352d695ee16e160ba6a675661972e3432ac383a0` ("Fix atomistic
DMI handedness"), `git status --short` clean for tracked source before this
session's edits.

**Environment (CPU, Patches 1–2):** GNU Fortran 13.3.0, GNU C/C++ 12.4.0,
CMake 3.28.3, CPU backend, fp64.

**Environment (CUDA GPU port, added 2026-08-06):** same host toolchain plus
CUDA 13.3.73 (`nvcc`), NVIDIA driver 610.57.04, two NVIDIA RTX A4000 GPUs
(`nvidia-smi`), configured with `-DUPPASD_GPU_BACKEND=CUDA`, fp64. HIP was
explicitly not touched this session (deferred by instruction) — no HIP
build, kernel-name check, or device evidence exists for the polarization
gate.

---

## Patch 1: polarization safety gate (F-14)

### Mechanism

`source/CoarseGraining/blockselector.f90` already provides a non-overridable
interlock: `advance_selector_state` checks `hard_atomistic_mask(block)`
before any dwell/hysteresis-gated refine/coarsen decision and forces
`RESOLUTION_ATOMISTIC` with reason `'hard-atomistic-exclusion'`. This patch
makes `adaptive_cg_state%hard_fine_mask` dynamic: on every synchronization
step, `update_adaptive_mask` (`source/CoarseGraining/adaptivecgproduction.f90`)
now calls the existing `restrict_channel_moments` (unchanged) into persistent
per-block scratch arrays, then a new pure function
`evaluate_polarization_gate` (`source/CoarseGraining/adaptivehybridsolver.f90`)
marks a block unsafe whenever any channel/ensemble has an undefined
direction (zero/near-zero resultant — reusing `restrict_channel_moments`'s
own "never normalize noise" convention) or a resultant/moment-sum ratio
below `cg_polarization_threshold`. `hard_fine_mask` is recomputed each step
as `static_hard_fine_mask .or. polarization_unsafe_block`, preserving any
pre-existing static-mask-file exclusion.

### New input keyword

`cg_polarization_threshold` (real, default `0.9`, validated to lie in
`(0,1]`), following the existing `cg_refine_threshold` pattern exactly:
`source/Input/inputdatatype.f90` (`adaptive_cg_config_t`),
`source/Input/inputhandler.f90` (keyword parsing), `validate_configuration`
in `adaptivecgproduction.f90` (bounds check).

### Fixtures and results

`POL-THRESHOLD` and `POL-CANCELLATION` are
`tests/coarse_graining/test_polarization_gate.f90`, calling
`evaluate_polarization_gate` directly:

- ratios strictly below, exactly at, strictly above, and roundoff-scale
  (`threshold - 1e-12`) below the threshold each select the documented
  `<` boundary;
- an exactly-cancelled channel resultant (`direction_defined=.false.`)
  is always unsafe;
- invalid-argument branches (not-ready topology, mismatched shapes,
  threshold outside `(0,1]`) are rejected.

**Negative control:** the threshold comparison was temporarily replaced
with an always-false condition (`ratio < 0.0_dblprec`) and rebuilt; the two
threshold-sensitive assertions failed as expected
(`ratio strictly below threshold is unsafe`,
`roundoff-scale below-threshold ratio is unsafe`). The source was restored
and the rebuilt suite passes again.

**Production regression (`POL-CANCELLATION`, operator layer):** a new
fixture, `tests/coarse_graining/e2e/polarization_gate_cpu`, uses the same
geometry as `parity_adaptive_cpu` with the misalignment selector
deliberately neutralized (`cg_refine_threshold 2.00`,
`cg_coarsen_threshold 1.99`, saturating the criterion so it always requests
coarsening) but a genuinely incoherent (`Initmag 1`, random) initial state.
`run_production_e2e.py` asserts `accepted_transitions == 0` and the final
`resolution_counts coarse=0`, i.e. every block remains atomistic. The
existing `parity_adaptive_cpu` fixture (identical thresholds, but an
aligned `momfile` initial state) does coarsen
(`accepted_transitions > 0`, already asserted) — same selector
configuration, different initial polarization, opposite outcome, isolating
the polarization gate as the differentiating mechanism.

### Fixture collateral: `parity_adaptive_cpu` / `parity_adaptive_gpu`

These two fixtures previously used `Initmag 1` (random) with the same
neutralized selector thresholds, asserting `accepted_transitions > 0`. That
assertion started failing under the new gate — correctly, since a random
initial state genuinely has low per-block polarization for this geometry
(4 atoms/block). Random noise coarsening into a macrospin is exactly what
this task prevents. Both fixtures were changed to `Initmag 3` (the same
aligned `momfile` already used by the working `adaptive_mixed` fixture) so
the parity/transition-mechanics assertions they exist for are exercised
under a physically legitimate precondition. `parity_adaptive_gpu` was
changed identically to `parity_adaptive_cpu` so a future CPU/GPU parity run
still compares matched inputs; it was not executed in this CPU-only
session.

### Commands and results

```text
cmake -S . -B <build> -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build <build> --target sd.f95 <all coarse_graining test targets> -j4
ctest --test-dir <build> -L coarse-graining
# 12/13 passed. The one failure (adaptive-cg-setup-rejection-matrix) is the
# pre-existing "mode:" rejection-matrix gap, reproduced identically on the
# unmodified base commit in a separate clean build — confirmed unrelated to
# this patch (see "Pre-existing failure" below).
```

### GPU (CUDA) port

Investigation found the GPU backend is not a thin compute-offload shim:
`GpuAdaptiveRuntime` (`source/gpu_files/gpuAdaptiveRuntime.cpp`) is a fully
autonomous adaptive-CG engine that makes its own per-block resolution-state
decisions in device kernels, driven from its own native time-step loop
(`gpuSDSimulation.cpp` → `GpuSimulation::advanceAdaptiveStep`). It never
read the Fortran `hard_fine_mask`, so the CPU gate above had **zero effect**
on GPU runs until this port. The anisotropy-uniformity and exchange-tensor-
symmetry checks (Patch 2) needed no GPU work: they run once at Fortran-side
setup, before GPU dispatch begins.

The existing kernel `proposeAdaptiveState` already accepted a `hardMask`
parameter with exactly the CPU `hard_atomistic_mask` bypass semantics, and
`restrictAdaptiveMoments` already computed the channel-resolved resultant/
moment-sum every step — both were already invoked, just never connected to
a gate. This port added:

- a new kernel `evaluateAdaptivePolarizationGate` (one thread per block,
  looping internally over channel/ensemble to avoid any multi-writer race)
  and its host method `GpuAdaptiveRuntime::evaluatePolarizationGate`,
  mirroring the Fortran `evaluate_polarization_gate` contract exactly,
  including the same `64*epsilon*scale` "never normalize noise" convention
  (factored into a shared `epsilonDevice()` helper also now used by
  `restrictAdaptiveMoments`, removing a previously duplicated `#ifdef`);
- a new device buffer (`polarizationUnsafeBlockMask_`) and its lifecycle
  (allocate/release/`estimateBytes`/`refreshDeviceDescriptors`), and a new
  `GpuAdaptivePhaseMetrics::polarizationMilliseconds` field for honest phase
  accounting;
- wiring in `GpuSimulation::advanceAdaptiveStep`: the new gate call now
  runs between the existing `restrictMoments`/`evaluateSelectorScores` and
  `proposeSelectorState`, whose previously-unused `hardAtomisticBlockMask`
  argument now receives the computed mask;
- a `cg_polarization_threshold` bridge extension through five files
  (`source/chelper.f90`, `source/gpu_files/fortranData.hpp`/`.cpp`,
  `source/gpu_files/gpuSimulation.hpp`/`.cpp`) — no new keyword, this only
  makes the CPU keyword also reach the GPU path.

**Device-level test (`tests/coarse_graining/test_gpu_adaptive_runtime.cpp`,
`testPolarizationGate`):** builds one block with two exactly-aligned atom
directions (high polarization, expected safe) and one block with two
exactly-opposing atom directions (undefined direction, expected unsafe),
calls `restrictMoments` then `evaluatePolarizationGate(0.9)`, downloads
`polarizationUnsafeBlockMask()`, and asserts both outcomes.

**Negative control (device level):** the kernel's unsafe-flagging condition
was temporarily replaced with `if(false) unsafe = 1;`, rebuilt, and rerun on
the RTX A4000: the test failed with `exactly-cancelled block was not
flagged unsafe`, confirming the assertion is defect-sensitive rather than
vacuously true. Reverted and confirmed the suite passes again.

**Production e2e (`tests/coarse_graining/e2e/polarization_gate_gpu`):**
same construction as the CPU `polarization_gate_cpu` fixture (neutralized
selector thresholds, random `Initmag 1`) with `do_gpu Y`/`gpu_mode 1`
added, run through the real `sd.f95.cuda` binary via
`run_production_e2e.py --gpu-binary`. Asserts `accepted_transitions == 0`
and the final GPU-printed `resolution_counts coarse=0` (the GPU engine
prints its own `Gpu: AdaptiveCG resolution_counts ...` line — the CPU
Fortran summary printer explicitly skips this output for GPU-dispatched
runs, so the harness pattern match target differs between the CPU and GPU
fixtures even though the assertion intent is identical).

**Negative control (production e2e):** with the same kernel fault injected,
`adaptive-cg-production-e2e` failed at
`assert metric(polarization_gate_gpu.stdout, "accepted_transitions") == 0`
— i.e. the real GPU production run actually coarsened the low-polarization
block when the gate was disabled. Reverted and confirmed the full suite
passes again.

**Commands and results:**

```text
nvidia-smi   # confirms 2x NVIDIA RTX A4000, driver 610.57.04
nvcc --version   # CUDA 13.3.73
cmake -S . -B <build> -DUPPASD_GPU_BACKEND=CUDA -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build <build> --target sd.f95.cuda gpu_adaptive_runtime_tests gpu_dmi_dimer_tests <CPU coarse_graining targets> -j8
ctest --test-dir <build> -L coarse-graining --exclude-regex 'adaptive-cg-setup-rejection-matrix|dipole'
# 14/14 passed: all CPU coarse-graining targets, coarse-graining-gpu-dmi-dimer,
# coarse-graining-gpu-adaptive-runtime (including testPolarizationGate),
# and adaptive-cg-production-e2e (including the new GPU fixture, run with
# --gpu-binary automatically supplied by CMake when USE_CUDA is set) all
# pass on real hardware. CPU/GPU parity for the pre-existing
# parity_adaptive_cpu/gpu pair (accepted_transitions equality, energy,
# field-checksum, and restart-state comparisons) also passed unchanged.
```

`docs/CG-10_GPU_ADAPTIVE_KERNELS.md` was checked for a policy-scalar
inventory table to update; it does not enumerate individual bridge scalars
in a way that needed a matching entry, so it was left unchanged.

---

## Patch 2: anisotropy capability audit + exchange tensor symmetry (F-06)

### Non-uniform anisotropy rejection

`build_production_anisotropy` (`adaptivecgproduction.f90`) samples
anisotropy from one central translated copy per basis index and broadcasts
it to every block. A new check, inserted between the existing per-atom
staging loop and that central-cell read, compares every atom sharing a
basis index (`anumb(atom)`) against the first such atom's staged
`taniso`/`eaniso`/`kaniso` (and `_diff` variants); any divergence is
rejected via the existing `reject(...)` helper (matching the neighboring
`random_anisotropy` capability rejection) rather than silently sampled.
Uniform anisotropy (the only case reachable through the ordinary,
non-ralloy, non-random kfile path, where `setup_anisotropies` sets every
atom's values purely as a function of `anumb(i)`) is unaffected.

**Tolerance bug found and fixed during verification:** the first version of
this check used `max(1.0_dblprec, ...)` as the comparison scale for the
Joule-scale `k1_j`/`k2_j` constants (typically `1e-21`–`1e-24` J). The `1.0`
floor completely swamped any real divergence, making the check vacuous — a
controlled fault injection (see below) initially passed through silently
and only surfaced as a wrong downstream energy value, not a rejection. The
scale was corrected to `max(tiny(1.0_dblprec), ...)`, matching the existing
single-channel pattern in `coarsetensoroperator.f90:203-211`.

**Negative control:** `ham%kaniso(1,2)` was temporarily perturbed
(`+1.0_dblprec`, a physically enormous but structurally realistic stand-in
for atom-indexed divergence such as a cluster-embedding override) right
before the per-atom staging loop, and the tree rebuilt. The corrected check
correctly rejects with `AdaptiveCG setup rejected: anisotropy: basis 2 is
not cell-periodic; atom 4 differs from atom 2`. Reverted and confirmed the
existing `anisotropy_uniform_fine`/`anisotropy_uniform_coarse` fixtures
(genuinely per-basis-uniform `kfile_cg_x`) still pass unchanged
(`ANI-UNIFORM-TRANSLATED` regression).

**`ANI-NONUNIFORM-REJECT` fixture — not delivered this session.** The
concrete atom-indexed divergence source identified during planning
(`do_cluster='Y'` overriding specific atom indices in
`source/Clusters/clusters.f90`) was traced further: `setup_clus_geometry`
increases `Natom` beyond `NA*N1*N2*N3` whenever `clus_expand > 0`, which the
pre-existing `validate_configuration` geometry check
(`Natom /= NA*N1*N2*N3`) already rejects before `build_production_anisotropy`
runs — for that common case, `do_cluster` is transitively unreachable and
would not exercise this new check at all. Whether a `clus_expand == 0`
cluster configuration (atoms replacing rather than adding to host sites)
can still reach this code path was not verified in this session, and the
`Clusters` subsystem has no existing example or test fixture in this
repository to adapt from. Rather than ship an unverified or possibly-dead
fixture, this negative control rests on the fault-injection evidence above,
which is one of the three forms of evidence the remediation blueprint's
section 2.3 accepts. A production `ANI-NONUNIFORM-REJECT` fixture (via a
verified `do_cluster` configuration or another genuine trigger) remains
open.

### Exchange tensor symmetry assertion

`coarsetensoroperator.f90` already asserted Cartesian symmetry of the
single-channel exchange stiffness tensor.
`multichannelcoarsetensoroperator.f90` asserted symmetry of the scalar
`local_exchange` matrix but never checked its sibling per-channel-pair
spatial tensor `exchange_stiffness(3,3,2,2)`, which the same mixed-derivative
gradient discretization relies on. A matching assertion was added,
looped over the four channel pairs, using the same `tiny()`-floored
relative tolerance as the single-channel reference (avoiding the same
scale-floor mistake described above).

**Negative control:** `tests/coarse_graining/test_multichannel_coarse_tensor_operator.f90`
gained `test_asymmetric_exchange_stiffness_rejected`, which perturbs
`material%exchange_stiffness(1,2,1,1)` after extraction (a state
unreachable through the ordinary atomistic extraction path, which always
yields a symmetric displacement-outer-product sum) and asserts setup is
rejected with a diagnostic containing `'Cartesian-symmetric'`.

### Commands and results

```text
cmake --build <build> --target sd.f95 multichannel_coarse_tensor_operator_tests -j4
ctest --test-dir <build> -L coarse-graining
# 12/13 passed; same pre-existing rejection-matrix gap, unrelated.
```

---

## Pre-existing failure (not caused by this session)

`adaptive-cg-setup-rejection-matrix` fails with
`AssertionError: mode: expected setup-time rejection (mode:)` on this
branch. This was reproduced identically after `git stash`-ing every change
made in this session and rebuilding in a second, independent out-of-tree
directory from the unmodified `352d695e` commit — confirming it predates
and is unrelated to RCG-03. It is the same class of pre-existing gap the
RCG-01 evidence already flagged ("the broader setup-rejection label retains
pre-existing failures outside RCG-01").

---

## Open limitations

- **HIP is untouched and untested**, by explicit instruction this session
  ("save HIP tests for later"). The `hipLaunchKernelGGL` branches for the
  new kernel/host method were written symmetrically with the CUDA branches
  (following the file's existing dual-path convention) but never compiled
  or run under `-DUPPASD_GPU_BACKEND=HIP`. This is a real gap, not merely an
  untested formality — the HIP path should be built and run before any HIP
  claim is made.
- **`ANI-NONUNIFORM-REJECT` has no production fixture** (see above);
  evidence rests on source-level fault injection only.
- **Independent review has not occurred.** Per the remediation blueprint's
  delegation guide, the implementer of a capability-safety change should
  not be its sole reviewer; that separation has not happened here.
- **RCG-02 is still open**, so per the blueprint's own sequencing rules this
  entire document records exploratory, not accepted, evidence. Before
  acceptance, this work must be rebased onto the commit that closes RCG-02
  and every fixture above rerun.

---

## RCG-03 closure (2026-08-08)

**RCG-03 is closed (2026-08-08, Human decision: Anders Bergman).** This
section covers the remaining work after RCG-02 closed: the mandatory
rebase/rerun, the three checklist items left open above, and the closure
decision itself, including HIP and independent review.

### Rebase and rerun

RCG-02 closed at commit `fae4c413` ("RCG-02: close, with HIP and
independent review deferred by Human decision"). The two RCG-03 patches
above (`20ef6f4d`, `3733b6ce`) are chronological ancestors of RCG-02's
closing commits, not the reverse, so no git rebase was needed — RCG-02's
DMI/Monte Carlo fixes already sit on top of RCG-03's polarization/
anisotropy work in a single line of history. HEAD at the start of this
session was `d8b8c5ab` ("Updating LSF moment file reading"), with a clean
`git status --short` for tracked source.

Per the blueprint's explicit instruction not to assume nothing regressed,
both suites were rerun from fresh out-of-tree builds on `d8b8c5ab` before
any new edits:

```text
CPU:  cmake -S . -B <build> -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
      cmake --build <build> -j2
      ctest --test-dir <build> -L cg13-cpu
      -> 12/12 passed, including adaptive-cg-production-e2e and
         adaptive-cg-setup-rejection-matrix (see note below).

CUDA: cmake -S . -B <build> -DUPPASD_GPU_BACKEND=CUDA -DBUILD_TESTING=ON \
        -DCMAKE_BUILD_TYPE=Release
      cmake --build <build> -j2
      ctest --test-dir <build> -L cg13-cuda
      -> 15/15 passed, including the GPU production e2e and GPU
         adaptive-runtime fixtures.
```

Environment: GNU Fortran 13.3.0, GNU C/C++ 12.4.0, CMake 3.28.3; CUDA
13.3.73 (`nvcc`), NVIDIA driver 610.57.04, two NVIDIA RTX A4000 GPUs
(`nvidia-smi`), `-DUPPASD_GPU_BACKEND=CUDA`, fp64.

**`adaptive-cg-setup-rejection-matrix` no longer has the pre-existing
`mode:` failure** that both the RCG-01 and 2026-08-06 RCG-03 evidence
flagged as unrelated and pre-existing: it passed cleanly in this rerun,
before any RCG-03-closure edits were made. This session did not
investigate why it now passes (out of scope), but it is recorded here
since it changes the honest baseline this closure work rebases onto.

**Note on tracked-file cleanliness:** running `adaptive-cg-production-e2e`
executes the real binary directly against `examples/AdaptiveCoarseGraining/*`
(tracked input directories that also receive the run's own provenance-
stamped `uppasd.*.yaml` output in place), which modifies those three
tracked files as a side effect of testing, not of any source edit. They
were restored with `git checkout` after every test run in this session so
`git describe --tags --dirty` reports a clean, non-`-dirty` tree between
edits.

### Item 1 — already-coarse unsafe blocks refine at an accepted sync point

**Finding before implementing a fixture:** production reconstruction
(`reconstruct_coarse_atoms` in `adaptivecgproduction.f90`) runs every step
for every currently-coarse block, before the polarization gate is ever
evaluated that same step, and it always rebuilds the dormant atoms to
exactly match the coarse channel's own moment/direction (both
`RECONSTRUCTION_ALIGNED` and `RECONSTRUCTION_CONSTRAINED_CONE` solve for a
resultant equal to `channel_moment_sum_mub * coarse_direction`, i.e. ratio
`1.0` to `1e-12`, by construction — see `reconstruct_block_constrained_cone`'s
bisection search on `cone_longitudinal` against `requested_norm`). This
means a block that has been coarse for at least one production step can
never be observed polarization-unsafe through the ordinary
`adaptive_cg_cpu_step` loop: reconstruction actively "launders" its
polarization back to `1.0` before the gate ever sees it, every step. A
production-executable fixture for this specific transition (coarse →
unsafe → forced atomistic, purely from step-loop dynamics) is therefore
not achievable honestly; building one would require either disabling
reconstruction (not a real production configuration) or waiting on an
event the architecture is specifically designed to prevent.

**Fixture delivered (operator layer):**
`tests/coarse_graining/test_adaptive_hybrid_solver.f90`,
`test_polarization_forces_refine_of_coarse_block`, matching the blueprint's
fixture-layer taxonomy (section 6.1, layer 2: "operator fixtures ... with
exact or high-accuracy references"), which does not go through
`reconstruct_coarse_atoms` at all:

1. A new single-block fixture helper, `make_single_block_fixture` (width
   `=n`, unlike the file's existing `make_chain_fixture`'s width `=1`), so
   a block can contain more than one atom per channel and therefore have a
   genuinely variable polarization ratio.
2. `setup_adaptive_hybrid_runtime` is confirmed to leave the block
   `RESOLUTION_COARSE` by construction (`fine_seed` all false) — "already
   coarse before any evolution," matching the checklist wording directly.
3. The block's dormant atoms are hand-set to a genuine, honest low-
   polarization state (three of four members aligned, one flipped),
   giving an exact `0.5` resultant/moment-sum ratio.
4. `restrict_channel_moments` + `evaluate_polarization_gate` (threshold
   `0.9`) are called directly on this state and confirmed to measure
   exactly `0.5` and flag the block unsafe.
5. `apply_adaptive_transitions` is called with that mask/ratio at
   `ADAPTIVE_STAGE_PREDICTOR`: rejected (`ADAPTIVE_HYBRID_INVALID_STAGE`),
   state unchanged — "not mid-integrator-stage."
6. The same call at `ADAPTIVE_STAGE_COMPLETE_STEP`: accepted, block
   transitions `RESOLUTION_COARSE -> RESOLUTION_ATOMISTIC`, and the logged
   transition event carries `selector_reason == 'polarization-unsafe'` and
   `polarization_ratio == 0.5` (both new fields; see Item 2).

**Negative control:** the `polarization_caused` branch in
`advance_selector_state` (`blockselector.f90`) was temporarily forced to
`.false.` (`if (.false. .and. present(...))`), rebuilt, and rerun: only the
new "own reason and ratio" assertion failed
(`FAIL: the forced transition is logged with its own reason and the
measured polarization ratio`); every other assertion, including the
state-transition and stage-gating ones, still passed. Reverted and
confirmed `adaptive_hybrid_solver_tests`, `polarization_gate_tests`, and
`block_selector_tests` all pass again.

### Item 2 — selector diagnostics report polarization and the reason

**CPU (Fortran):**

- `evaluate_polarization_gate` (`adaptivehybridsolver.f90`) gained an
  optional `block_ratio(:)` output: the worst (minimum) resultant/moment-
  sum ratio observed across every channel/ensemble at each block,
  regardless of whether the block is flagged unsafe, with the documented
  `0.0` floor for an undefined (near-zero, never-normalized) direction.
- `advance_selector_state` and `apply_adaptive_transitions`
  (`blockselector.f90`, `adaptivehybridsolver.f90`) gained optional
  `polarization_unsafe_mask`/`polarization_ratio` parameters. When
  supplied, a hard-forced transition is logged with reason
  `'polarization-unsafe'` instead of the generic `'hard-atomistic-
  exclusion'` whenever polarization was (at least partly) the cause;
  `'hard-atomistic-exclusion'` is kept, unchanged, for a static-mask-only
  exclusion. All new parameters are optional, so every pre-existing call
  site (both production and the other unit tests) is unaffected.
- `selector_transition_event_type` and `adaptive_transition_event_type`
  gained a `polarization_ratio` field (sentinel `-1.0`, never a valid
  ratio, when not supplied), threaded through
  `append_transition`/`append_adaptive_event`.
- Production wiring: `adaptive_cg_state` gained a persistent
  `polarization_ratio_block(:)` array, filled by the same
  `evaluate_polarization_gate` call `update_adaptive_mask` already made,
  and passed into `apply_adaptive_transitions` alongside the existing
  `polarization_unsafe_block` mask. `print_transition_events` now prints a
  `polarization_ratio=<value>` field on every `AdaptiveCG: transition ...`
  line.

**GPU (CUDA; HIP deferred, see below):**

- `evaluateAdaptivePolarizationGate` (`gpuAdaptiveRuntime.cpp`) gained the
  same worst-ratio tracking, writing to a new
  `GpuAdaptiveDeviceRuntime::polarizationRatio` device buffer
  (`polarizationRatioBlock_`, allocated/wired/freed alongside the existing
  `polarizationUnsafeBlockMask_`). The loop can no longer short-circuit on
  the first unsafe hit the way the mask-only version could, since every
  channel/ensemble must be visited to find the true minimum.
- `estimateBytes` was updated to include the new buffer; omitting this
  caused `gpu_adaptive_runtime_tests` to abort with "CG-10 scratch or
  construction storage bypassed memory preflight" until fixed.
- `GpuAdaptiveDiagnosticSnapshot` gained a `polarizationRatio` vector,
  downloaded in `diagnosticSnapshot()`; `gpuSimulation.cpp`'s
  `release()` diagnostic block prints a new
  `Gpu: AdaptiveCG polarization_ratio values=...` line alongside the
  existing `final_state` line.
- Unlike the Fortran side, the GPU's `hardAtomisticBlockMask` (consumed by
  `proposeAdaptiveState`) is *only* ever populated from
  `polarizationUnsafeBlockMask()` (`gpuSimulation.cpp`) — there is no
  separate static-mask concept ported to the GPU path yet — so a distinct
  GPU-side reason enum would have exactly one possible value today. A code
  comment records this explicitly rather than adding a redundant enum; the
  device-reachable gap this checklist item is actually about (no numeric
  ratio surfaced) is closed by the buffer/snapshot/print above.
- `test_gpu_adaptive_runtime.cpp`'s `testPolarizationGate` gained a new
  `polarizationRatioBlock()` accessor call, asserting the aligned block's
  ratio is `1.0` and the exactly-cancelled block's ratio is `0.0`
  (matching the CPU floor convention).

**Commands and results:**

```text
cmake --build <cpu build> -j2 --target sd.f95 adaptive_hybrid_solver_tests \
  polarization_gate_tests block_selector_tests
ctest --test-dir <cpu build> -L cg13-cpu           # 12/12
cmake --build <cuda build> -j2 --target gpu_adaptive_runtime_tests sd.f95.cuda
ctest --test-dir <cuda build> -L cg13-cuda         # 15/15
```

HIP was not touched, per the same explicit deferral as the 2026-08-06
session (see "HIP and independent review" below): the `hipLaunchKernelGGL`
launch call for `evaluateAdaptivePolarizationGate` is unchanged (the kernel
body it launches is shared source, compiled by either backend, so the
change automatically applies to a future HIP build), but it has not been
built or run under `-DUPPASD_GPU_BACKEND=HIP`.

### Item 3 — ANI-NONUNIFORM-REJECT production fixture

**New fixture:** `tests/coarse_graining/e2e/ani_nonuniform_reject_cluster/`
(full provenance in that directory's `README.md`), wired into
`tests/coarse_graining/run_setup_rejection_matrix.py` as a 31st rejection
case and into `tests/coarse_graining/fixture_dependencies.py`'s
`all_e2e_cases()`/`audit_fixture_dependencies.py` tracked-fixture audit.

**Mechanism (reachable, unlike `do_cluster`'s original candidate path):**
a `do_cluster` embedding of exactly one atom, placed at the same fractional
coordinate as an existing host lattice site (`ncell_clus 1 1 1`,
`posfile_clus` at `(0,0,0)`), so `clus_expand` (`clus_geometry.f90`) stays
`0` and `Natom` is unchanged. The pre-existing `Natom /= NA*N1*N2*N3`
geometry rejection in `validate_configuration` — which the 2026-08-06
session found transitively blocks the `clus_expand>0` case — therefore
never fires. `anisotropy_clus` then gives that one embedded atom (basis 1)
a divergent `k1` (`-0.005`) from the host's uniform kfile value for basis 1
(`-0.002`, from `../kfile_cg_x`, the same file `anisotropy_uniform_fine`
already uses). Running the ordinary `sd.f95` executable against this input
produces:

```text
AdaptiveCG setup rejected: anisotropy: basis 1 is not cell-periodic; atom 3 differs from atom 1
```

— `build_production_anisotropy`'s cell-periodicity check
(`adaptivecgproduction.f90`), reached through the ordinary executable and
an ordinary capability, exactly as the checklist requires.

**Two pre-existing, unrelated setup-order obstacles were found and worked
around from the fixture's own input files, with no source changes:**

1. `Landeg_glob` defaults to `2.0`, but the ordinary host moment path
   halves it (`Landeg(i)=Landeg_ch(...)*0.5`,
   `source/System/magnetizationinit.f90`) while `cluster_creation`
   (`source/Clusters/clusters.f90`) copies the cluster's `Landeg_ch_clus`
   into `Landeg(iatom)` directly, with no such factor. With the default
   `set_landeg 0`, the embedded atom would silently end up with `Landeg=2.0`
   against every host atom's `1.0`, tripping the unrelated
   `Landeg/do_site_damping: ... requires uniform gamma and damping`
   capability guard before the anisotropy check is ever reached. The
   fixture sets `set_landeg 1` (making the momfiles' gyromagnetic column
   authoritative for both host and cluster) and gives `momfile_clus` a
   `1.0` column so the embedded atom's `Landeg` matches the host exactly.
2. `cluster_creation`'s exchange-overwrite loop unconditionally consumes
   `ham_clus%nlistsize_clus`/`ncoup_clus` regardless of whether the
   cluster's anisotropy is even in play, so `exchange_clus` (`jfile_clus`)
   is required input; the fixture supplies one self-referencing bond with
   zero coupling so the host's real exchange is left untouched.

Both are genuine, narrow inconsistencies in the pre-existing `do_cluster`
embedding path, unrelated to RCG-03's capability-safety scope. Neither
required a source change — a candidate fix to the apparently-missing
`ham_clus%taniso_clus`/`kaniso_clus`/`eaniso_clus` population was drafted
and then reverted after finding that `hamiltonianinit.f90:299-313` already
calls the shared `setup_anisotropies` routine for the cluster path
correctly; that candidate fix was an unnecessary duplicate based on an
incomplete search and never appears in the final diff.

**Negative control:** replacing `kfile_clus`'s `k1=-0.005` with the host's
own `k1=-0.002` (anisotropy-matched embedding) was run manually and
completes normally (`AdaptiveCG: capability accepted`, `returncode=0`)
instead of rejecting, confirming the rejection is produced by the genuine
per-atom divergence and not merely by `do_cluster`'s presence.

**Commands and results:**

```text
python3 tests/coarse_graining/run_setup_rejection_matrix.py --binary <sd.f95>
# CG-13 setup-rejection matrix passed (31 cases)
python3 tests/coarse_graining/audit_fixture_dependencies.py
# adaptive-CG fixture dependency audit: PASS (38 fixture directories, 60 input paths)
```

### HIP and independent review

Per Human decision (Anders Bergman, 2026-08-08), RCG-03 closes with the
same two deferrals RCG-02 closed with, for the same reasons:

- **HIP execution evidence** is deferred: no HIP toolchain or device
  exists in any environment used so far. The `hipLaunchKernelGGL` launch
  call for the polarization-gate kernel, and the (unchanged, shared)
  kernel body it launches, remain untested under
  `-DUPPASD_GPU_BACKEND=HIP`.
- **A separate independent Opus/Terra or Sol adversarial physics/
  implementation review**, distinct from Human approval, is deferred to a
  later stage of the remediation program.

Neither deferral reflects a physics disagreement or an unresolved
correctness question — every fixture this document and the remediation
blueprint's RCG-03 checklist require passes now on CPU and CUDA.

### Closure summary

All fourteen RCG-03 checklist items in the remediation blueprint are now
ticked on delivered, tested evidence. Exit evidence:

- `POL-THRESHOLD`, `POL-CANCELLATION` — passing since 2026-08-06, rerun
  clean in this session's fresh builds.
- `ANI-UNIFORM-TRANSLATED` — passing since 2026-08-06 (Fortran-side
  setup, no GPU port needed), rerun clean.
- `ANI-NONUNIFORM-REJECT` — new production fixture delivered above.
- Already-coarse-unsafe-refine and polarization-reason/ratio diagnostics —
  new operator fixture and CPU+CUDA diagnostic wiring delivered above.

`ctest -L cg13-cpu` (12/12) and `ctest -L cg13-cuda` (15/15) both pass on
fresh out-of-tree builds of the final commit; `run_setup_rejection_matrix.py`
passes 31/31 cases; `audit_fixture_dependencies.py` passes.
