# RCG-03 polarization gate and anisotropy capability audit — exploratory evidence

**Status:** exploratory work performed while RCG-02 remains open (see
`docs/RCG-02_DMI_HANDEDNESS_EVIDENCE.md`). Per
`docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md` section 8/10, a
later task may be implemented and tested ahead of its dependency closing,
but it is not accepted until the dependency closes and this work is
rebased onto that accepted commit and rerun. **No RCG-03 checklist box in
the remediation blueprint is ticked by this document.**

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
