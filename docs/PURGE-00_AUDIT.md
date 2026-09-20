# PURGE-00 — Freeze baseline and establish purge ownership map

> Historical disposition note: this audit records the pre-purge state at
> `438d14d7`. For the final disposition, see
> [`ADAPTIVE_CG_PURGE.md`](ADAPTIVE_CG_PURGE.md).

Date: 2026-09-20
Branch: `gpu_hip_cu_ab_dip`
HEAD: `438d14d721a2aac82745413000a84501f50b0f8e` (`438d14d7 Added type projection for S(r)`)

This is an evidence-only audit. It freezes the current baseline and identifies the adaptive exchange/coarse-graining implementation that may be removed in the later purge. No production source, build file, test, input, reference output, GPU file, or legacy file was changed for this audit. The only intended repository change is this document.

## A. Branch and worktree evidence

The required branch and commit were verified:

```text
git branch --show-current
gpu_hip_cu_ab_dip

git rev-parse HEAD
438d14d721a2aac82745413000a84501f50b0f8e

git log -1 --oneline
438d14d7 Added type projection for S(r)
```

`gpu_hip_cu_ab_cg` points to the same commit as the current branch. Its merge-base with `gpu_hip_cu_ab_dip` is also `438d14d721a2aac82745413000a84501f50b0f8e`, and the committed branch-to-branch diff is empty. The CG branch reference was not modified.

The worktree was already substantially dirty before this audit: `git status --short` reported 20 modified tracked paths and 41 untracked paths. The existing changes include example inputs, ordinary CPU-Hamiltonian sources/tests, build artifacts, GUI artifacts, and additional files. They were preserved and not used as a reason to reset, clean, or overwrite the worktree. The build evidence below therefore distinguishes a clean `HEAD` snapshot from the existing dirty checkout.

## B. Baseline build and test evidence

### Clean `HEAD` build

To avoid mixing pre-existing worktree changes into the baseline, `HEAD` was exported with `git archive` to `/private/tmp/purge00-source` and configured in `/private/tmp/purge00-build`.

Configuration:

```text
cmake -S /private/tmp/purge00-source -B /private/tmp/purge00-build \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBUILD_TESTING=ON -DBUILD_PYTHON=OFF \
  -DUPPASD_GPU_BACKEND=OFF -DUPPASD_PRECISION=DOUBLE \
  -DRUN_SLD_TESTS=ON -DRUN_GNEB_TESTS=ON
```

Result: configure succeeded. The environment provided GNU Fortran 16.2.0, AppleClang, OpenMP, FFTW, BLAS/LAPACK, and a CPU-only DOUBLE configuration. CUDA/HIP were not built or run because no GPU toolchain/hardware was required for this CPU baseline.

Build:

```text
cmake --build /private/tmp/purge00-build -j2
```

Result: succeeded at 100%.

The existing checkout was also rebuilt with `cmake --build build -j2`; it succeeded at 100%, but that build is not the clean baseline because it incorporates the pre-existing worktree changes.

### CTest inventory

The clean `HEAD` inventory contains 46 registered tests:

```text
regression-test
asd-tests
sld-tests
gneb-tests
coarse-graining-block-topology
coarse-graining-stiffness-material
coarse-graining-dmi-dimer-energy
cpu-ham-field-energy-contract
cpu-ham-energy-convention
cpu-ham-backend-policy
cpu-ham-sparse-backend
cpu-ham-target-order
cpu-ham-reduced-stencil
coarse-graining-adaptive-hamiltonian-contract
coarse-graining-tensor-operator
coarse-graining-dispersion
coarse-graining-multichannel-tensor-operator
coarse-graining-smooth-projected-operator
coarse-graining-static-hybrid-operator
coarse-graining-static-hybrid-dilation-scaling
coarse-graining-adaptive-setup-scaling
coarse-graining-block-selector
coarse-graining-adaptive-hybrid-solver
coarse-graining-polarization-gate
coarse-graining-moving-state-generator
coarse-graining-trajectory-evidence
coarse-graining-torque-oracle
coarse-graining-coarse-torque-oracle
coarse-graining-static-topology-oracle
coarse-graining-ownership-map-comparator
adaptive-cg-mem-large-host
adaptive-cg-reconstruction-rng-spatial-stats
adaptive-cg-production-e2e
adaptive-cg-setup-rejection-matrix
adaptive-cg-moving-off-fine
adaptive-cg-moving-all-coarse
adaptive-cg-moving-static-mixed
adaptive-cg-moving-adaptive-wall
adaptive-cg-moving-dmi-chiral
adaptive-cg-timing-reconciliation
adaptive-cg-fixture-dependencies
adaptive-cg-transition-ownership-invariants
dipole-open-fft-oracle
dipole-open-host-builder
dipole-open-host-goldens
dipole-ewald-host-builder
```

The clean baseline non-adaptive selection was:

```text
ctest --test-dir /private/tmp/purge00-build --output-on-failure \
  -R '^(regression-test|asd-tests|sld-tests|gneb-tests|cpu-ham-|dipole-open-|dipole-ewald-)'
```

Results:

| Result | Tests | Interpretation |
|---|---|---|
| PASS | `cpu-ham-field-energy-contract`, `cpu-ham-energy-convention`, `cpu-ham-backend-policy`, `cpu-ham-sparse-backend`, `cpu-ham-target-order`, `cpu-ham-reduced-stencil`, `dipole-open-fft-oracle`, `dipole-open-host-builder`, `dipole-open-host-goldens`, `dipole-ewald-host-builder` | Clean CPU Hamiltonian and surviving open/FFT/Ewald dipole host baseline passes. |
| FAIL — environment | `regression-test`, `asd-tests`, `sld-tests`, `gneb-tests` | All fail before the physics checks because `tests/bergtest.py` cannot import the unavailable Python module `yaml` (`ModuleNotFoundError`). |
| NOT RUN | CUDA/HIP tests | GPU toolchains/hardware were not part of this CPU baseline. |
| NOT RUN | Adaptive/coarse-graining tests | Intentionally not used as non-adaptive acceptance evidence; their later disposition is recorded below. |

There were no silent omissions in the selected clean CPU run. The dirty checkout’s `ctest -N` inventory additionally contains two current-worktree CPU-Hamiltonian tests (`cpu-ham-convolution`, `cpu-ham-j-plus-d`) and GPU/dipole additions; those were not treated as clean-`HEAD` evidence.

## C. Deterministic baseline manifest

The following clean-`HEAD` executables were run directly. Their stdout was captured and hashed as an exact, deterministic oracle:

| Test executable | Exact stdout | SHA-256 of stdout |
|---|---|---|
| `/private/tmp/purge00-build/bin/dmi_dimer_energy_tests` | `DMI-DIMER-ENERGY passed` | `13d26598b2c29b9809250943583218e4dc48c9bcda2ffcf2a4d2b74a1d0fc80b` |
| `/private/tmp/purge00-build/bin/block_topology_tests` | `block topology tests passed` | `24bc68909ac0b8170d1196a3961fd1a7457e11627b9140ecef24b9663ac3772f` |
| `/private/tmp/purge00-build/bin/cpu_ham11_energy_tests` | `CPU-HAM-11 energy convention parity passed` | `b33b000f51467d181b624951528466bbcacf67bb27313fc69096eaffc6c32646` |

These three checks are deliberately retained as purge guardrails: the DMI dimer test protects independent DMI correctness, the block-topology test protects surviving macrocell/coarse-block infrastructure, and CPU-HAM-11 protects the independent CPU Hamiltonian energy convention.

The existing Heisenberg-chain workflow produces timestamped and stochastic output files; its current acceptance is tolerance-based rather than exact byte/hash equality. It was not promoted into this exact manifest.

## D. Ownership map

Classification vocabulary:

- `KEEP`: surviving ordinary ASD/LLG, GPU, FFT/Ewald dipole, CPU-Hamiltonian, macrocell, or coarse-block infrastructure.
- `REMOVE`: adaptive exchange/coarse-graining implementation with no non-CG production consumer found.
- `REHOME`: test or documentation content that is currently grouped under coarse-graining but protects an independent surviving feature.
- `INSPECT DURING PURGE`: mixed files where only an adaptive slice may be removed.
- `ARCHIVE`: historical adaptive-CG documentation retained outside the active implementation path.

### Coarse-graining source ownership

| Component/file | Current role | Non-CG consumer/evidence | Classification | Planned stage | Evidence |
|---|---|---|---|---|---|
| `source/CoarseGraining/blocktopology.f90` | Regular spatial block IDs and FFT-grid channel mapping | `macrocells.f90` and `hamiltonianmacroblocks.f90` use `BlockTopology`; macrocell paths are used by ordinary drivers, Hamiltonian, dipole, SLD/GNEB, MC, and related production code | KEEP | Preserve in the coarse-block/macrocell pass | Direct `use BlockTopology` consumers found outside adaptive operators |
| `source/CoarseGraining/macrocells.f90` | Macrocell layout and ordinary coarse-block/macrocell support | `uppasd`, `sd_driver`, `mc_driver`, `sld_driver`, GNEB, MS/PT/WL/SX, Hamiltonian and dipole code | KEEP | Preserve | Broad production fan-out; this is not adaptive atom ownership |
| `source/CoarseGraining/hamiltonianmacroblocks.f90` | Hamiltonian macroblock layout and block-level Hamiltonian support | `source/Hamiltonian/hamiltonianinit.f90` calls `build_macroblock_layout` | KEEP | Preserve | Direct non-CG production caller |
| `source/CoarseGraining/blockselector.f90` | Adaptive active-block selection | Adaptive modules/tests only | REMOVE | Remove after adaptive callers are removed | No surviving non-CG consumer found |
| `coarsetensoroperator.f90` | Adaptive projected tensor operator | Adaptive modules/tests only | REMOVE | Remove with adaptive operator stack | No surviving non-CG consumer found |
| `multichannelcoarsetensoroperator.f90` | Adaptive multi-channel projected operator | Adaptive modules/tests only | REMOVE | Remove with adaptive operator stack | No surviving non-CG consumer found |
| `smoothprojectedoperator.f90` | Adaptive smooth/projected operator | Adaptive modules/tests only | REMOVE | Remove with adaptive operator stack | No surviving non-CG consumer found |
| `statichybridoperator.f90` | Adaptive static hybrid operator | Adaptive modules/tests only | REMOVE | Remove with adaptive operator stack | No surviving non-CG consumer found |
| `adaptivehybridsolver.f90` | Adaptive hybrid solver and transitions | Adaptive production/tests only | REMOVE | Remove with adaptive runtime | No surviving non-CG consumer found |
| `adaptivecgproduction.f90` | Adaptive CG setup, step, diagnostics, and teardown | `uppasd`, `sd_driver`, `chelper`, GPU bridge, and adaptive tests | REMOVE | Remove after mixed production hooks are excised | Dedicated adaptive owner; all callers are purge targets |
| `source/CoarseGraining/CMakeLists.txt` | Registers all coarse-graining objects | Build-system-only dependency | INSPECT DURING PURGE | Keep the three surviving objects; remove the seven adaptive objects after callers are gone | Current list contains both surviving macrocell infrastructure and adaptive operators |

The seven adaptive objects planned for removal from the coarse-graining build lists are `blockselector`, `coarsetensoroperator`, `multichannelcoarsetensoroperator`, `smoothprojectedoperator`, `statichybridoperator`, `adaptivehybridsolver`, and `adaptivecgproduction`. `blocktopology`, `macrocells`, and `hamiltonianmacroblocks` remain.

### Stiffness ownership split

`source/SpinWaves/stiffness.f90` is mixed and must not be deleted wholesale.

| Slice | Current role | Evidence | Classification | Planned stage |
|---|---|---|---|---|
| Legacy stiffness arrays, `do_stiffness`, `stiffness_wrapper`, `init_stiffness` | Ordinary legacy stiffness and spin-wave support | `source/uppasd.f90` uses `Stiffness`; `source/SpinWaves/prn_micromagnetic.f90` provides the legacy routines; `source/Input/inputhandler.f90` parses `do_stiffness` | KEEP | Preserve | Independent production consumers exist |
| Continuum-material/CG extraction types and procedures: `coarse_lattice_input_type`, `coarse_lattice_sums_type`, `coarse_material_metadata_type`, `coarse_material_diagnostics_type`, `coarse_material_type`, `calculate_coarse_lattice_sums`, `fit_coarse_material`, `extract_coarse_material`, `extract_coarse_material_from_uppasd`, `validate_coarse_material_small_q`, `validate_coarse_material_two_sublattice_modes`, `coarse_material_runtime_status`, and associated `COARSE_MATERIAL_*`/`COARSE_RUNTIME_*` constants | Adaptive CG material extraction and validation | Consumers are adaptive tensor/operator modules and adaptive tests; no non-adaptive production consumer found | REMOVE slice | Inspect and remove only this CG slice during the stiffness pass | Keep the legacy stiffness portion intact |

### Production hook map

| File | Adaptive hook | Surviving neighboring role | Classification / planned stage |
|---|---|---|---|
| `source/uppasd.f90` | Adaptive setup, preflight, handoff, runtime gating, cleanup, and summary around the `adaptive_cg` state | Ordinary initialization, macrocell setup, dipole, and the main ASD/LLG workflow | INSPECT DURING PURGE; remove only adaptive blocks |
| `source/sd_driver.f90` | `use AdaptiveCGProduction`, `adaptive_cg_is_enabled()`, and `adaptive_cg_cpu_step` | Ordinary SD integration and macrocell-aware execution | INSPECT DURING PURGE; remove adaptive branch only |
| `source/chelper.f90` | Adaptive topology/kernel C ABI bridge (`FortranData_setAdaptiveTopology`, `clearAdaptiveTopology`, `setAdaptiveKernels`) | `setMacrocell` and `setPmeMacrocell` bridge surviving macrocell/dipole infrastructure | INSPECT DURING PURGE; retain ordinary bridge |
| `source/Input/inputdatatype.f90` | Adaptive configuration/state type and defaults | Ordinary block-size, macrocell, dipole, and GPU input state | INSPECT DURING PURGE; remove adaptive fields only |
| `source/Input/inputdata.f90` | Adaptive state reset/default handling | Ordinary input state handling | INSPECT DURING PURGE; remove adaptive reset only |
| `source/Input/inputhandler.f90` | Adaptive keyword parser cases | Ordinary macrocell/dipole/GPU keywords | INSPECT DURING PURGE; remove exact adaptive cases only |
| `source/measurement` and `source/Measurement` | No adaptive-CG production hook found | Ordinary measurement; unrelated adaptive-timestep wording is not adaptive CG | KEEP |
| `source/gpu_files/gpuSimulation.cpp/.hpp` | Adaptive runtime state, setup, step, masks, reconstruction, and diagnostics | Ordinary GPU simulation, measurement, integrator, and dipole flow | INSPECT DURING PURGE |
| `source/gpu_files/gpuSDSimulation.cpp` | Adaptive-enabled branch and `gpuAdaptiveMomentUpdater` integration | Ordinary GPU SD path and integrator | INSPECT DURING PURGE |
| `source/gpu_files/gpuHamiltonianCalculations.cpp/.hpp` | `UpdateAdaptiveMacroMoments`, adaptive FFT-dipole methods/state | Ordinary Hamiltonian, `GpuDipoleConvolution`, and FFT/Ewald dipole code | INSPECT DURING PURGE; keep ordinary backend/dipole APIs |
| `source/gpu_files/fortranData.cpp/.hpp` | Adaptive state and Fortran data exchange | Ordinary Fortran/GPU data exchange | INSPECT DURING PURGE |
| `source/gpu_files/fort_helper.cpp` | `gpusim_updateadaptivemask_` and related adaptive bridge | Ordinary GPU wrappers and Fortran ABI | INSPECT DURING PURGE |
| `source/gpu_files/measurement` | No adaptive-specific path found | Ordinary GPU measurement | KEEP |

### GPU helper spillover

| File | Adaptive-only or mixed symbol family | Non-CG evidence | Classification |
|---|---|---|---|
| `gpuAdaptiveRuntime.cpp/.hpp` | Dedicated adaptive runtime | `gpuSimulation` and adaptive tests only | REMOVE |
| `gpuAdaptiveMomentUpdater.cpp/.hpp` | Dedicated adaptive moment updater | Adaptive GPU path/tests only | REMOVE |
| `gpuAdaptiveReconstructionRng.hpp` | Adaptive reconstruction RNG | Adaptive reconstruction/tests only | REMOVE |
| `gpuAtomicDouble.hpp` | Generic-looking atomic double helper | Only adaptive runtime/energy tests use it; no surviving production caller found | REMOVE |
| `gpuDepondtIntegrator.cpp/.hpp` | Active-list overloads and `CommitActivePredictor`, `RestoreActiveInitial`, `CommitActiveCorrector` | Full-range `evolveFirst/evolveSecond` are used by ordinary GPU simulation/SD | INSPECT DURING PURGE; keep full-range integrator, remove active subset slice |
| `gpuThermfield.cpp/.hpp` | Active-list randomization overload | Full `randomize(mmom)` is used by ordinary Depondt integration | INSPECT DURING PURGE; keep ordinary overload |
| `gpuParallelizationHelper.hpp/.tpp` | `active_atom_kernel`, `active_atom_site_kernel`, `gpuActiveAtomCall`, `gpuActiveAtomSiteCall` | Full-range helpers have ordinary consumers | INSPECT DURING PURGE; remove active subset helpers only |
| `gpuCommon.hpp` | `AddAtoms`, `AvgAtoms`, `AddToAtoms` active-subset helpers | Ordinary `Add`, `AddTo`, `Avg`, `ScalarMult`, `Inv` remain in use | INSPECT DURING PURGE; preserve shared ordinary helpers |
| `gpuDipoleConvolution.cpp/.hpp` | Adaptive callers may use `devicePaddedField()`; “active grid” descriptors also occur in ordinary macrocell/dipole projection | Ordinary open FFT/Ewald, `addFieldsToAtoms`, `accumulateEnergy`, and coarse-block/macrocell projection are surviving production roles | KEEP file; inspect/remove only an adaptive-only method after caller proof |
| `gpu_wrappers.h` | Shared GPU abstraction | Used by ordinary and adaptive GPU code | KEEP; remove only a proven adaptive bridge |

The GPU CMake inventory confirms that `gpuAdaptiveMomentUpdater.cpp`, `gpuAdaptiveRuntime.cpp`, `gpuAdaptiveRuntime.hpp`, and `gpuAdaptiveMomentUpdater.hpp` are dedicated entries. Mixed GPU files must be edited surgically later; the presence of the word “adaptive” or “active” alone is not sufficient evidence for file removal.

## E. Exact adaptive input-keyword inventory

The following 18 parser keywords are the adaptive-CG inventory from `source/Input/inputhandler.f90` and the corresponding `adaptive_cg_config_t` state:

```text
do_adaptive_cg
cg_operator
cg_mask_mode
cg_selector
cg_refine_threshold
cg_coarsen_threshold
cg_polarization_threshold
cg_update_interval
cg_minimum_dwell_updates
cg_buffer_blocks
cg_channel_mode
cg_channel_file
cg_reconstruction
cg_cone_angle
cg_static_mask_file
cg_energy_jump_gate
cg_energy_jump_limit
cg_diagnostics
```

These keywords and their adaptive state/default/reset handling are purge targets. The following similarly named or neighboring inputs are not adaptive-only and must remain unless a later dependency audit proves otherwise:

```text
block_size_x  block_size_y  block_size_z  block_size
do_macro_cells  prn_dip_subset  dip_file
GPU dipole / FFT / Ewald keywords
```

In particular, macrocell/block geometry and dipole controls belong to the surviving macrocell, coarse-block, CPU-Hamiltonian, and FFT/Ewald infrastructure.

## F. Test disposition

The disposition is based on dependency and protected behavior, not on directory or test-name matching alone.

### KEEP

- `regression-test`, `asd-tests`, `sld-tests`, `gneb-tests`
- CPU-Hamiltonian tests: `cpu-ham-field-energy-contract`, `cpu-ham-energy-convention`, `cpu-ham-backend-policy`, `cpu-ham-sparse-backend`, `cpu-ham-target-order`, `cpu-ham-reduced-stencil`, plus the current-worktree `cpu-ham-convolution` and `cpu-ham-j-plus-d`
- Surviving dipole tests: `dipole-open-fft-oracle`, `dipole-open-host-builder`, `dipole-open-host-goldens`, `dipole-ewald-host-builder`, and conditional GPU/open-FFT tests such as `dipole-gpu-fft-convolution`, `dipole-open-fft-layout`, `dipole-gpu-wp5-e2e`, and `dipole-open-fft-coarse-e2e`

### REHOME

- `coarse-graining-block-topology` — move/classify with surviving `BlockTopology`, macrocell, and dipole/coarse-block infrastructure.
- `coarse-graining-dmi-dimer-energy` — move/classify with independent DMI correctness.
- Conditional `coarse-graining-gpu-dmi-dimer` — move/classify with independent GPU DMI correctness.

Moving adaptive DMI/chiral tests remain adaptive-path tests and are not protected merely because they contain “DMI”.

### REMOVE

The following tests exercise adaptive operators, adaptive transitions, adaptive ownership, or adaptive GPU runtime and have no independent non-CG contract:

```text
coarse-graining-stiffness-material
coarse-graining-adaptive-hamiltonian-contract
coarse-graining-tensor-operator
coarse-graining-dispersion
coarse-graining-multichannel-tensor-operator
coarse-graining-smooth-projected-operator
coarse-graining-static-hybrid-operator
coarse-graining-static-hybrid-dilation-scaling
coarse-graining-adaptive-setup-scaling
coarse-graining-block-selector
coarse-graining-adaptive-hybrid-solver
coarse-graining-polarization-gate
coarse-graining-moving-state-generator
coarse-graining-trajectory-evidence
coarse-graining-torque-oracle
coarse-graining-coarse-torque-oracle
coarse-graining-static-topology-oracle
coarse-graining-ownership-map-comparator
adaptive-cg-mem-large-host
adaptive-cg-reconstruction-rng-spatial-stats
adaptive-cg-production-e2e
adaptive-cg-setup-rejection-matrix
adaptive-cg-moving-off-fine
adaptive-cg-moving-all-coarse
adaptive-cg-moving-static-mixed
adaptive-cg-moving-adaptive-wall
adaptive-cg-moving-dmi-chiral
adaptive-cg-timing-reconciliation
adaptive-cg-fixture-dependencies
adaptive-cg-transition-ownership-invariants
```

Conditional adaptive GPU tests/benchmarks with the same dependency are also purge targets: `adaptive-cg-energy-fp32-accum`, `adaptive-cg-energy-hierarchical-precision`, `coarse-graining-gpu-adaptive-runtime`, `gpu-depondt-active-atoms`, `gpu-adaptive-moment-updater`, `adaptive-cg-dilation-sanitizer`, `coarse-graining-adaptive-asd-parity`, `gpu_adaptive_runtime_benchmark`, `gpu_inactive_runtime_overhead_benchmark`, and `gpu_thermfield_rng_benchmark`. The independent `dipole_gpu_fft_benchmark` remains.

## G. Build-system and documentation inventory

### Build systems

The following build files contain the adaptive/coarse-graining registration that will require a later surgical edit:

- `source/CoarseGraining/CMakeLists.txt`: remove the seven adaptive objects listed in section D; keep `blocktopology.f90`, `macrocells.f90`, and `hamiltonianmacroblocks.f90`.
- `source/Makefile.gfortran`, `source/Makefile.legacy`, `source/Makefile.legacy.gfortran`: each lists the same coarse object group and a generic `CoarseGraining/%.o` rule. Remove only the seven adaptive object entries later; keep the three surviving objects and the generic rule if still needed.
- `source/gpu_files/CMakeLists.txt`: remove dedicated adaptive GPU source/header entries later and preserve ordinary GPU, dipole, FFT, and macrocell entries. Mixed source files require code-level surgery.
- Top-level `CMakeLists.txt`: contains adaptive tests and conditional GPU adaptive tests/benchmarks as well as surviving CPU-Hamiltonian and dipole registrations. Edit by target dependency, not by broad text substitution.

No build file was edited in PURGE-00.

### Documentation

Archive as historical adaptive-CG material rather than presenting it as active behavior:

- `docs/ADAPTIVE_COARSE_GRAINING_BASELINE_20260730.md`
- `docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md`
- `docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md`
- Root `docs/CG-*.md`, `docs/CGP-*.md`, and `docs/RCG-*.md`
- `docs/cg14/**`, `docs/rcg08/**`, `docs/rcg09/**`, and `docs/rcg10/**`

Keep active or historical documentation for surviving behavior, including:

- `docs/CPU_HAM_*` and `docs/cpu/**`
- `docs/FFT-dipole_implementation_plan.md`
- `docs/FFT-dipole_status.md`
- `docs/GPU_FFT_DIPOLE_DESIGN.md`
- `docs/WP10_OPEN_FFT_BLUEPRINT.md`
- `docs/WP10_OPEN_FFT_LUNA.md`
- `docs/sol_wp10_7a_open_coarse_projection.md`
- `docs/terra_wp10_7b_open_coarse_integration.md`
- `docs/terra_wp10_8_open_performance.md`
- `docs/luna_acceptance_review.md`, `docs/luna_wp10_6b_open_acceptance.md`, `docs/luna_wp10_7c_open_coarse_acceptance.md`, `docs/luna_wp10_final_acceptance.md`, and `docs/luna_wp10_open_acceptance.md`
- General production/build/feature/benchmark documentation such as `docs/FEATURES.md`, `docs/UppASD_GPU_features.md`, `docs/BUILDING_GPU.md`, and production benchmark material

The word “coarse” is not sufficient to archive a document: the open FFT and coarse-block projection documents describe surviving dipole infrastructure and must remain available.

### GPU newline convention

No line-ending normalization was performed. The current files were measured as follows:

| File | Observed convention |
|---|---|
| `source/gpu_files/gpuSimulation.cpp` | CRLF-only |
| `source/gpu_files/gpuSimulation.hpp` | Mixed: 104 CRLF records and 2 LF-only records |
| `source/gpu_files/gpuSDSimulation.cpp` | CRLF-only |
| `source/gpu_files/gpuSDSimulation.hpp` | Not present in the tree |
| `source/gpu_files/gpu_wrappers.h` | CRLF-only |

Later GPU edits must preserve each file’s existing convention and must not normalize whole files.

## H. Open items for PURGE-01

These are genuine implementation-stage checks, not reasons to broaden the purge boundary:

1. After adaptive callers are removed, re-run the symbol search for `gpuDipoleConvolution::devicePaddedField()` before deleting that method; the current evidence points to adaptive ownership, but the surviving dipole file must remain intact.
2. Rehome the three independent topology/DMI tests into a surviving test grouping and update their CTest labels/paths without changing their protected physics contract.
3. Re-check the split edit in `stiffness.f90` after removing the CG continuum slice; legacy stiffness symbols and `do_stiffness` must still compile and run.
4. Because the starting worktree is dirty, run the post-purge build/test matrix from a fresh commit or clean export. Do not use a worktree cleanup as a shortcut.

## Acceptance record

- Required branch and expected HEAD verified.
- CG branch reference left untouched.
- Existing dirty worktree preserved.
- Clean CPU-only `HEAD` configure and build passed.
- Full clean CTest inventory captured; selected non-adaptive suite classified as PASS/FAIL-environment/NOT RUN.
- Deterministic stdout/hash manifest recorded for DMI, BlockTopology, and CPU-HAM energy guardrails.
- Adaptive ownership map, mixed production hooks, GPU spillover, exact input keywords, build registrations, test disposition, documentation disposition, and newline constraints recorded.
- No adaptive implementation, production path, build registration, test, input, reference physics, or line ending was changed by PURGE-00.
