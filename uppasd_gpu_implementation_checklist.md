# UppASD GPU port — implementation checklist

Companion to `uppasd_gpu_audit.md`. This document turns the audit findings, plus four additional focus areas (memory footprint, CMake modernization, convolution backend, flexible precision), into self-contained tasks suitable for handing to a coding assistant (Claude Code / Opus / Sonnet) one at a time. Cross-references like (C1), (P4), (I2) point to sections of the audit.

---

## Instructions for the implementing assistant

Read this section before every task.

1. **Scope discipline.** Do exactly one task per session/branch. Do not refactor code outside the task's file list, do not "fix" unrelated style issues, and do not change any physics formulas unless the task explicitly says so. If a task turns out to require touching files outside its list, stop and report why instead of expanding scope.
2. **Fortran is the reference implementation.** The CPU code in `source/*.f90` defines correct physics. When in doubt about a sign, prefactor, or convention, find the corresponding Fortran routine (`hamiltonianactions.f90`, `depondt.f90`, `montecarlo.f90`, `prn_averages.f90`, `correlation*.f90`) and match it. Never modify the Fortran physics to match the GPU code.
3. **Build both backends.** Every change must compile under CUDA (`-DUSE_CUDA=ON`) and HIP (`-DUSE_HIP=ON`). If you cannot compile HIP locally, at minimum ensure no CUDA-only API is introduced outside `gpu_wrappers.h` / `real_type.h`.
4. **Validate.** After Phase 0 lands, run the T=0 regression (`tests/gpu_regression`, task V2) after every task and report the diff summary. Tasks that change numerics on purpose (precision modes) state their own tolerances.
5. **Update this file.** Mark the task's checkbox `[x]`, and append a one-line note after the task (`> Done <date>, PR #, notes`). If a task is found to be invalid or already fixed, mark `[x]` with a note explaining why.
6. **Commit style.** One task = one PR/branch named `gpu/<task-id>-short-slug`. Commit messages reference the task ID.
7. **Line references are approximate.** Re-locate by symbol name, not line number.

Priorities: **P0** = correctness/crash, do first. **P1** = high-value. **P2** = valuable, do after P0/P1. **P3** = cleanup/opportunistic.
Effort: S (< half day), M (a day-ish), L (multi-day).

### Checkbox reconciliation, 2026-07-20

The checkboxes had gone stale: ~200 commits had landed against 4 ticked items, which made
finished work look pending (M3, for instance, still read as "blocked on CV4" after CV4 was
already satisfied). Items marked **verified landed 2026-07-20** were reconciled in a sweep.

**What that marker means, and does not mean.** It means the implementing code was located and
inspected for that item. It does **not** mean the item's stated Acceptance criterion was
re-executed — no profiling traces were captured, and no per-item benchmarks were re-run. Treat
those ticks as "the work exists" rather than "the acceptance bar was re-measured." PF3 is the
exception: it was validated end-to-end (see its entry).

**Correction, same day — Phase 4 (PR1–PR7) was over-ticked and has been reverted to `[ ]`.**
The first pass ticked the PR cluster on the strength of `SINGLE_PREC` code existing. That was
the exact error this note warns about. Phase 4 is specified in terms of a three-mode
`UPPASD_PRECISION` (DOUBLE/SINGLE/MIXED) with `storage_t` / `accum_t` / `accum_cplx` typedefs —
**none of those identifiers exist anywhere in the tree.** What landed is a *binary* `SINGLE_PREC`
switch reaching five files (`real_type.h`, `tensor.hpp`, `gpu_wrappers.h`, `gpuFftWrapper.hpp`,
`gpuSimulation.cpp`). Real work, but not the specified design, and MIXED — the mode PR3 exists
to serve — is absent entirely. PR7's tolerances are also a single 1e-8 tier, not per-precision.
Phase 4 is best read as "prototyped, not designed in."

**Correction 2026-07-20 — V3 was over-ticked and has been un-ticked.** `28f4b5a` marked it
"verified landed"; the tree has no `try`/`catch` in `fort_helper.cpp` and no `MemGetInfo`
anywhere under `source/gpu_files/`. Same failure mode that pass had with PR1–PR7. It is a P0
and it **blocks M1**. If you find another item ticked by `28f4b5a` that you are about to build
on, re-verify it against the tree before trusting it — that pass is not reliable.

**Deliberately left unticked despite plausible evidence:**

- **CM1 is now done** (`b426193`) — `UPPASD_GPU_BACKEND`/`UPPASD_PRECISION` exist, the GPU
  `add_compile_definitions` are target-scoped, and `-DDEBUG` injection is gone. See the CM1
  entry for what was verified and the one deliberate deviation (C++20 kept, not 17).
- **CM2 is now done** (`5515dba`) — legacy CUDA discovery, `select_compute_arch` and the
  blanket `CUDA::toolkit` link are gone. See the CM2 entry, including the NVTX dependency the
  blanket link had been hiding.
- **CM4 is now done** (`e32c11b`) — one source list per directory, and the reason the
  duplication existed (directory-scoped source properties) is written up in the CM4 entry,
  along with a latent HIP no-op it exposed.
- **CM5 is now done** (`d5ebd17`) — presets, `docs/BUILDING_GPU.md`, and CI switched over. The
  hand-configured `build_fastcopy`/`build_ptds`/`build_ab` now have preset equivalents.
- **CM3 and CM6 remain not done** (line numbers below predate CM1/CM2/CM4 and have shifted —
  re-locate by symbol): `CMAKE_CXX_COMPILER`/`CMAKE_HIP_COMPILER` set inside
  CMakeLists (:420-421), which CM3 explicitly forbids; hardcoded `/opt/rocm` search paths
  (:449 and `gpu_files/CMakeLists.txt:66`). **CM3 was skipped deliberately: no AMD hardware
  here.** Note that `.github/workflows/gpubuilds.yml` has a `rocm/dev-ubuntu-24.04` job, so
  CM3's acceptance ("configures in the ROCm container") is reachable through CI even though it
  is not reachable locally.
- **PR5** — the fp32 infrastructure exists, but that does not establish that the audit of
  numerically sensitive spots was actually carried out. Its precision-independent piece (the
  versine formula) has since landed in `9b1011c`; the rest is still written against `accum_t`,
  which does not exist until PR1 lands. PR1 is itself no longer blocked now that CM1 is in.
- **F5, F6, PF1, PF2** — MC items. F5's exact anisotropy ΔE and F6's sublattice coloring both
  appear implemented, and `single_spin_mc_uniaxial_easy_axis` passes, but MC is gated out of
  scope and MC semantics were not validated. Confirm before ticking.

**Not tracked as an item but should be:** `sc_dm_uniaxial` segfaults in
`gpusim_initiatematrices` (`sd_driver.f90:1185`) in every configuration — the only failing
regression case. Also `cpuRestMeasurement.cpp:166` reads `fortran_beff.copy_sync(emom)`,
copying `emom` into the `beff` buffer; likely a typo, unfixed.

---

## Phase 0 — Safety net (do these before anything else)

### - [x] V1 — GPU compile-only CI  (P0, M)
**Files:** `.github/workflows/` (new file `gpubuilds.yml`), top-level `CMakeLists.txt` only if a flag is broken.
**Task:** Add a CI workflow that compiles (not runs) the GPU configurations on plain runners:
- Job 1: Ubuntu + NVIDIA CUDA toolkit (apt `nvidia-cuda-toolkit` or NVIDIA's repo), `cmake -DUSE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=80`, full build. nvcc compiles without a GPU present.
- Job 2: container `rocm/dev-ubuntu-24.04` (or current LTS tag), `cmake -DUSE_HIP=ON -DCMAKE_HIP_ARCHITECTURES=gfx90a`, full build.
- Both jobs must also compile `chelper.f90` (they will, since `USE_CUDA/USE_HIP` pulls it in) with gfortran — this is what catches Fortran breakage in GPU-only files.
**Acceptance:** Both jobs green on the branch after C3 is fixed; a deliberately introduced syntax error in `chelper.f90` fails CI.
**Note:** Expect the CUDA job to *fail* on the current branch until C3/CM tasks land; that is the point. Land V1 together with C3.

> Done 2026-07-18, local change: added CUDA 13.3 and ROCm compile-only jobs, both serializing the build because the current Fortran module output directory is shared. CUDA CMake now preserves an explicitly requested architecture for GPU-less runners. Remote CI execution remains to be observed.

### - [x] V2 — Deterministic T=0 regression harness  (P0, L)  — **verified landed 2026-07-20**
**Files:** new directory `tests/gpu_regression/` (inputs + Python driver), optionally a tiny hook in `chelper.f90`/drivers if a "dump moments after K steps" switch doesn't already exist (it does: restart/moment printing via `prn_mag_conf`).
**Task:** Build a pytest (or plain Python) harness that, for a matrix of small systems (~10³–10⁴ atoms), runs **the same input** through the Fortran path (`do_gpu N`) and the GPU path (`do_gpu Y`) and compares final moment configurations and energy outputs.
- Cases: bcc Fe Heisenberg-only; +DM; +uniaxial aniso; +cubic aniso; tensor exchange; reduced Hamiltonian (`do_reduced`); MC mode; convolution backend on/off (`do_gpu_convolution`); M=1 and M=4 ensembles; one case with N not a multiple of 256 (e.g. N=1000) and one with odd 3·N·M — these two specifically exercise audit bugs C1 and C2.
- All runs at `temp 0` / `ip_temp 0` with thermal noise irrelevant, `SDEalgh 1`, short (50–200 steps) → trajectories are deterministic up to floating-point reassociation. Compare with rtol ~1e-10 (fp64) on emom and per-term energies.
- Emit a compact pass/fail table.
**Acceptance:** Harness runs on a GPU machine with one command; currently-failing cases are recorded (expected: energy cases with N%256≠0 fail until C1). This harness is the acceptance gate for most later tasks.

> Implemented 2026-07-18, local change: added a standard-library Python harness and JSON case matrix. It runs CPU and GPU variants in isolated copied fixtures, compares final restarts and all energy rows at fp64 tolerance, and includes scalar, DM/aniso, tensor, reduced, M=4, N=1000, odd-count, convolution, and MC variants. Mark complete after its first run on a host with a visible device.

### - [x] V3 — Memory & error reporting baseline  (P0, S)  — **part 1 (catch + report) landed 2026-07-20; part 2 folded into M1**

**Correction.** This was ticked "verified landed" by `28f4b5a` (the same reconciliation pass
that over-ticked PR1–PR7). It is not implemented. Evidence, re-checked against the tree:
`fort_helper.cpp` contains **no `try`, no `catch`, and no `MemGetInfo`** — the four `extern "C"`
entry points at `fort_helper.cpp:42-56` are bare one-line forwarders. `grep -rn "MemGetInfo"`
over `source/gpu_files/` returns **nothing at all**, so part 2 (the post-init summary line) does
not exist either. `git log -S "catch (const std::exception" -- source/gpu_files/fort_helper.cpp`
returns no commit, so it was never added and later removed — it simply never landed.

The premise of the task is intact and still correct: `ASSERT_GPU` (`base.hpp:32-36`) does
`throw std::runtime_error`, and nothing between there and Fortran catches it.

**This blocks M1**, whose task text opens "V3 delivers the catch + report. Extend with a
pre-allocation estimate." Do V3 first, or fold it into M1 as its first commit.

**Files:** `gpu_files/fort_helper.cpp`, `gpu_files/gpuSimulation.cpp`, `gpu_files/measurement/memoryMeasurement.cpp`.
**Task:** The tree already has `TensorMemoryTracker` (host+device byte counters, `printResults()`, `saveToFile()`), but nothing catches the `std::runtime_error` that `ASSERT_GPU` throws on allocation failure — it propagates through `extern "C"` into Fortran and terminates without diagnostics (this is the "very large systems crash" symptom, see M1). Two small changes:
1. Wrap the bodies of all `extern "C"` entry points in `fort_helper.cpp` (`gpusim_initiateconstants_`, `gpusim_initiatematrices_`, `gpusim_gpurunsimulation_`, `gpusim_release_`) in `try { ... } catch (const std::exception& e)` that prints the exception message, calls `TensorMemoryTracker::printResults()` and the free/total device memory from `cudaMemGetInfo`/`hipMemGetInfo`, then `std::exit(EXIT_FAILURE)` with a clear "GPU out of memory / GPU error" banner.
2. After a successful `initiateMatrices()`, print one summary line: device bytes allocated per `TensorMemoryTracker`, and free/total from `MemGetInfo`.
**Acceptance:** Running a deliberately oversized system produces a readable OOM report with the tracker breakdown instead of a bare `terminate called after throwing...`.

> Done 2026-07-20 (commit `0582ac5`, part 1 of 2). All four `extern "C"` entry points in
> `fort_helper.cpp` now wrap their body in `try { ... } catch (const std::exception& e)` calling a
> single `reportGpuFailureAndExit(where, e)` helper. The helper flushes stdout, prints a banner
> that says **GPU OUT OF MEMORY** when `e.what()` contains "out of memory" and **GPU ERROR**
> otherwise (the exception only carries the driver string, so the substring match is how OOM is
> distinguished — addresses the "don't swallow non-OOM errors" risk), then `TensorMemoryTracker::printResults()`
> for the peak host/device breakdown, then free/total from a new `GPU_MEM_GET_INFO` wrapper
> (`cudaMemGetInfo`/`hipMemGetInfo`, added to `gpu_wrappers.h`), then `std::exit(EXIT_FAILURE)`.
> Drive-by: fixed the botched `"...allocate meFlags.mory: %s"` string at `gpuSimulation.cpp:377`.
>
> **Demonstrated** against a real OOM: a 2 000 000-atom bccFe run (`ncell 100 100 100`, needs
> 1.87 GB device) on a GPU pre-occupied by a 15 000 MiB holder (0.65 GB free). Output — exit
> code 1, no `terminate called...`:
> ```
> ========================================================================
>  GPU OUT OF MEMORY in gpusim_initiatematrices
>    out of memory
> ------------------------------------------------------------------------
> ------------------TENSOR MEMORY CONSUMPTION REPORT----------------------
> CATEGORY                PEAK CONSUMPTION
> Max CPU                 0                        B
> Max GPU                 8                        MB
>  GPU memory: 0.65 GB free / 16.75 GB total
> ========================================================================
> ```
> The banner fires inside `gpusim_initiatematrices`, before any kernel launch, and the exit is
> clean (no deadlock against the measurement thread — that thread never touches the tracker mutex;
> only `tensor.hpp` Allocate/Free on the main thread do — so `printResults()` from the handler
> is safe). CUDA-only path exercised; **HIP path unvalidated (no AMD hardware).**
> Regression suite unchanged: 11 PASS / 1 FAIL on `cuda-debug`.
>
> **Part 2 not done here.** The task's second item — a post-`initiateMatrices` success summary
> line (allocated device bytes + free/total) — is deferred into **M1**, where the `MemGetInfo`
> budget reporting around `initiateMatrices` lives naturally; doing it twice would duplicate the
> reporting. M1 must add it or re-scope it explicitly.
>
> **Correction to the handover's `sc_dm_uniaxial` claim:** the persistently failing regression
> case does **not** segfault in `gpusim_initiatematrices`. Observed under this build, the CPU
> *reference* binary aborts with a Fortran runtime bounds error in `read_exchange`
> (`hamiltonianinit.f90:1326`, "Index '5' of dimension 1 of array 'nlist' above upper bound of
> 4") during Hamiltonian setup — before any GPU work and unrelated to OOM or the GPU boundary.
> V3/M1 cannot address it; it is a separate Fortran/input bug, out of scope for M1.
>
> **Exit-path hardening (2026-07-21, commit `7620a9d`).** The six "Uncaught GPU error" handlers in
> `gpuSDSimulation.cpp` / `gpuMCSimulation.cpp` did `GPU_DEVICE_RESET(); std::exit()` — the reset
> tears down the context, then `std::exit` runs the static `GpuSimulation` destructor →
> `release()` → `cudaFree` on the dead context → `terminate: driver shutting down` → SIGABRT (this
> is the second failure in the large-system crash log). Replaced with flush + `std::_Exit`, which
> skips static destructors and any further CUDA calls — robust against sticky errors too. V3's own
> `reportGpuFailureAndExit` was changed from `std::exit` to `std::_Exit` for the same reason (a
> post-init error caught at the boundary would otherwise hit the same static-dtor `cudaFree`).
> Demonstrated by temporarily reintroducing the `sumOverAtoms` overflow (commit `4d7af10` fixed a
> separate large-system bug: the atom-sum finalize kernel sized dynamic shared memory by the
> atom-block count, overflowing the 48 kB/block limit past ~1.5 M atoms): the run now prints
> "Uncaught GPU error 1: invalid argument" and exits code 1 cleanly, no terminate/SIGABRT.

---

## Phase 1 — Correctness fixes (audit C1–C10)

### - [x] F1 — Fix the Hamiltonian energy reduction  (P0, M) — audit C1  — **verified landed 2026-07-20**
**Files:** `gpu_files/gpuHamiltonianCalculations.cpp` (function `sum_warp_energy` and all `Heisge*::each` energy blocks), possibly `gpu_files/measurement/kernels.hpp/.cpp` (to reuse `block_reduce_sum_1d`).
**Background:** Three bugs: (a) divergent `__syncthreads()` because out-of-range threads (`site >= N` in the tail block) never enter `each()`; (b) `warp_sums[w]` read for warps that never wrote it this launch (stale shared memory); (c) no `__syncthreads()` between the three consecutive `sum_warp_energy` calls (shared-memory race).
**Task (recommended shape):** Replace the in-kernel block reduction + `atomicAdd` with the two-phase pattern already correctly implemented in `measurement/kernels.cpp`: the Heisge kernel writes per-atom energy terms (or per-block partial sums computed with a *correct* reduction) to a small global scratch buffer; a separate finalize kernel reduces per ensemble into `energyM`. If you keep an in-kernel reduction instead, you must: move the bounds check inside `each()` so **all** threads of the block execute the reduction (out-of-range threads contribute 0), zero-initialize `warp_sums` before use with a `__syncthreads()`, and add a trailing `__syncthreads()` at the end of `sum_warp_energy`. Also replace `int mInd = blockIdx.y` with the `ensemble` argument that `each(atom, site, ensemble)` already receives (this simultaneously fixes audit C7 for the energy path).
**Acceptance:** V2 energy cases pass for N=1000 (not a multiple of 256), for all Hamiltonian variants, and repeated runs give bit-identical energies (no race). No performance regression > 5 % on the measure-step (energies are measured rarely; correctness dominates).

### - [x] F2 — RNG length + status checking  (P0, S) — audit C2  — **verified landed 2026-07-20**
**Files:** `gpu_files/cudaThermfield.cu`, `gpu_files/hipThermfield.cpp` (or their unified successor after U1).
**Task:** `curand/hiprandGenerateNormal[Double]` require an even sample count; `field.size() = 3·N·M` can be odd and the status is unchecked. Allocate the field buffer with one element of slack when `size()` is odd (`n_gen = size + (size & 1)`), generate `n_gen` values, and wrap **every** `curand*/hiprand*` call in a status check macro (add `ASSERT_CURAND`/`ASSERT_HIPRAND` next to `ASSERT_GPU` in `base.hpp`, throwing with the status code). Do the same for the generator creation/seed/stream calls in `initiate`.
**Acceptance:** A run with N·M odd (e.g. N=1001, M=1) at finite temperature produces a fluctuating thermal field (compare magnetization variance against the Fortran run statistically); status-check macro verified by forcing an error in a scratch test.

> Test note, 2026-07-18: the former `sc_odd_rng_count` GPU-regression case ran at `temp=0` and changed the fixture geometry to 1001×1×1, so it was removed from the default suite rather than serving as F2 coverage. Replace it with an odd-sized finite-temperature statistical case when GPU-MC/thermal validation is revisited.

### - [x] F3 — `chelper.f90` USE statement  (P0, S) — audit C3
**Files:** `source/chelper.f90` line ~20.
**Task:** `use prn_averages,` has a trailing comma and no only-list — invalid Fortran. Determine what symbols from `prn_averages` are actually referenced in `Chelper` (grep the module; likely the measurement buffer arrays passed in `FortranData_setMeasurables`: `mavg_buff`, `mavg2_buff`, `mavg4_buff`, `mavg_buff_proj`, `indxb_ac`-style arrays — check `prn_averages` module contents) and restore `use prn_averages, only : <list>`. Compile with gfortran under `-DUSE_CUDA=ON` to verify.
**Acceptance:** GPU configuration compiles with gfortran (CI job V1 green).

> Done 2026-07-18, local change: restored an explicit `only` import list for every `prn_averages` symbol used by `Chelper`; compiled `source/chelper.f90` through the CUDA CMake target with CUDA 13.3. CI coverage is provided by V1.

### - [x] F4 — Ensemble index under `USE_BIG_GRID`  (P0, S) — audit C7  — **verified landed 2026-07-20**
**Files:** `gpu_files/gpuHamiltonianCalculations.cpp` (all `blockIdx.y` uses inside functors), `gpu_files/gpuMetropolis.cpp` (verify only — `MCSweep` launches its own grid where `blockIdx.y` *is* the ensemble by construction; leave it but add a comment).
**Task:** Inside `each(atom, site, ensemble)` functors, replace every `int mInd = blockIdx.y;` with the `ensemble` parameter. Grep the whole tree for `blockIdx` inside `each(` methods and fix all instances — functors must never read launch geometry directly.
**Acceptance:** V2 passes with `-DUSE_BIG_GRID=true` forced in a test build.

### - [ ] F5 — MC: exact on-site anisotropy ΔE  (P0, M) — audit C5
**Files:** `gpu_files/gpuMetropolis.cpp` (`MCSweep`), `gpu_files/gpuMCSimulation.cpp` (pass aniso arrays through), `gpu_files/gpuMetropolis.hpp`.
**Background:** The acceptance uses `ΔE ∝ eneff·(S_new−S_old)`, exact for bilinear terms under sublattice decomposition but wrong for on-site anisotropy (energy is quadratic/quartic in the flipped spin).
**Task:** Pass `kaniso/eaniso/taniso/sb` (device pointers) and the aniso flag into `MCSweep`. Compute the exact anisotropy energy difference for old vs new spin in-kernel using the same formulas as Fortran `montecarlo.f90` (uniaxial: `E = k1·(1−(S·e)²) + k2·(1−(S·e)²)²` — verify exact convention against `calculate_energy`/`effective_field` in the Fortran MC path; cubic analogous) and add it to ΔE, while keeping the field-based term for exchange/DM/Zeeman **with the anisotropy contribution removed from `eneff` in MC mode** (cleanest: in MC mode call the Heisge variant *without* aniso to build `eneff`, and handle aniso entirely in `MCSweep`). Alternative minimal fix if that's too invasive: `if (do_aniso) { fprintf(stderr, "GPU MC with anisotropy not yet supported"); exit; }` — acceptable as an interim commit but the full fix is the task.

> Progress, 2026-07-18: strict GPU-MC now uses the Fortran `calculate_energy` uniaxial convention, `K1(S·e)^2 + K2(S·e)^4`, for the exact on-site trial-spin ΔE and obtains its field term without anisotropy. Cubic (`taniso=2`) and combined (`taniso=7`) anisotropy are rejected before GPU-MC initialization. The brute-force GPU-MC path still rejects all anisotropy.

> **Deferred / WIP:** Revisit GPU-MC as a complete unit before declaring this task finished: reconcile GPU anisotropy measurement energies with the Fortran convention, support or explicitly scope the remaining anisotropy modes, add the finite-temperature statistical CPU-vs-GPU regression, and complete the linked F6/PF1/PF2 work.

**Acceptance:** New V2 case: MC at low finite T with uniaxial anisotropy, compare equilibrium ⟨m_z⟩ and energy against Fortran MC statistically (same T, enough steps; 3σ agreement). Note this task naturally merges with performance task PF1 (local ΔE) — if PF1 is done first, F5 becomes "add the on-site terms to the local ΔE", which is the preferred order.

### - [ ] F6 — MC sublattice coloring on union of neighbor lists  (P1, S) — audit C6

> **Deferred / WIP:** Revisit together with the rest of the GPU-MC path (F5/PF1/PF2).
**Files:** `gpu_files/gpuMetropolis.cpp` (`split_lattice` and its callers).
**Task:** Build the independence graph from the union of `nlist` and `dmlist` (and any future lists) so that parallel updates never touch interacting pairs. Combine with PF2 (coloring rewrite) if convenient.
**Acceptance:** A test system with DM neighbors ⊄ exchange neighbors (construct via input with different cutoffs) yields a coloring where no sublattice contains a DM-coupled pair (assert in a debug check).

### - [x] F7 — Upload staging: kill the in-place transpose  (P1, M) — audit C4 (pairs with PR2)  — **verified landed 2026-07-20**
**Files:** `gpu_files/gpuSimulation.cpp` (`copyFromFortran`), `gpu_files/tensor.hpp` (new helper), later `real_type.h`.
**Task:** Add a `Tensor<T,2>::transposed_copy_to(GpuTensor<T,2>&)` (host-side transpose into a scratch `std::vector`/pinned buffer, then upload) or a small device transpose kernel; use it for `ncoup`/`nlist` so the Fortran-owned arrays are **never mutated**. Delete the transpose/upload/transpose-back dance. Design the helper signature so PR2 (precision conversion on upload) can extend it with a source/destination type pair.

> Progress, 2026-07-18: `transposed_copy_to` now stages the transposed rank-2 data in a host vector and uploads it directly; `copyFromFortran` no longer mutates `ncoup` or `nlist`.
**Acceptance:** V2 passes; add an assertion test that `ham%ncoup` on the Fortran side is bit-identical before/after `gpuSim_initiateMatrices` (e.g., checksum printed from Fortran around the call in a debug run).

### - [x] F8 — `fastCopy` event ordering  (P1, S) — audit C8
**Files:** `gpu_files/measurement/fortranMeasurement.cpp` (`copyQueueFast`), same pattern in `correlations/fortranCorrelation.cpp` if present.
**Task:** Record an additional event `d2hDone` on `copyStream` **after** the three pinned D2H copies, and enqueue the queue callback behind it (either `GPU_STREAM_WAIT_EVENT(workStream, d2hDone)` immediately before the callback, or move the callback to `copyStream`). Keep the existing early-release of `workStream` after the D2D event — that part is correct and intentional. While here, migrate `GPU_STREAM_ADD_CALLBACK` to `cudaLaunchHostFunc`/`hipLaunchHostFunc` in `gpu_wrappers.h`.
**Acceptance:** Build with `-DUSE_FAST_COPY`; run a measurement-heavy V2 case under `compute-sanitizer --tool racecheck` (CUDA) with no reports from this path; averages match the non-fastCopy run bitwise.

> Progress, 2026-07-20: validated on top of `8f1b657`. The `d2hDone` choreography itself was
> correct, but three surrounding defects were found and fixed:
> (a) `flushMeasurements`/`flushCorrelations` still synchronized only `workStream`, while C8 had
> moved the consumer callback to `copyStream`, so flushes no longer waited for pending pushes;
> (b) `SDmphase` called `mqueue.finish()` *before* the final fast-copy callback had run, joining
> the consumer thread and destroying its mutex with a push still outstanding;
> (c) `SDmphase` never exported final state to Fortran — it relied on `copyQueueSlow` refreshing
> `fortran_emomM` as a side effect, so with `-DUSE_FAST_COPY` the restart file held a stale
> `printMdStatus` snapshot. `GpuEventPool::Event::active` is now `std::atomic<bool>` (it is
> cleared from the driver callback thread and read from the host thread).
> Validated with a dedicated `-DUSE_FAST_COPY` build (`build_fastcopy`): full GPU regression suite
> matches the default build, and the fast-copy cases are bitwise-identical to CPU across repeats.
> Note `USE_FAST_COPY` is not set by any CMake option, so the fast path is dead code by default.

### - [x] F9 — Move device query out of static init; guard zero-norm rotation; `initiate` idempotence  (P1, S) — audit C9, C10  — **verified landed 2026-07-20**
**Files:** `gpu_files/gridHelper.hpp`, `gpu_files/gpuParallelizationHelper.{hpp,cpp,tpp}`, `gpu_files/gpuDepondtIntegrator.cpp`, `gpu_files/gpuHamiltonianCalculations.cpp`.
**Task:** (a) Move `cudaGetDeviceProperties` from the `GridHelper` constructor into `GpuParallelizationHelper::initiate()` (store `maxGridSize1` there, pass to `GridHelper`), with error checking. (b) In `Rotate::each`, guard `norm == 0` (skip rotation, `mrod = emom`). (c) In `GpuHamiltonianCalculations::initiate`, add a `listsPrepared` flag on the device Hamiltonian struct set by `SetupNeighbourList*`; if `initiate` is called again without a fresh upload (flag still set), either skip the setup kernels or error out loudly — never double-decrement. (d) Remove the unused `delta_t` parameter of `GpuDepondtIntegrator::rotate` or actually use it.
**Acceptance:** Compiles and V2 passes; a scratch test calling `initiate` twice does not corrupt neighbor lists.

---

## Phase 2 — Memory footprint  ★ user emphasis 1

**Context (from the audit of allocations).** Current per-run device footprint, fp64, in units of `NM = N·M`:
`gpuLattice`: beff, b2eff, eneff, emomM, emom, emom2 (6 × 3NM·8 B) + extfield (3NM·8 B, always) + btorque (3NM·8 B if STT) + mmom/mmom0/mmom2/mmomi (4 × NM·8 B).
`GpuDepondtIntegrator`: mrod, blocal, bdup (3 × 3NM·8 B) + thermfield field (3NM·8 B) → integrator alone adds ~4 × 3NM.
Total streaming arrays: **11–12 × 3NM doubles ≈ 288 NM bytes** (+32 NM for the mmom family). For N=10⁷, M=1 that is ~3.2 GB before the Hamiltonian.
`Hamiltonian`: nlist `mnn×N×4 B` + ncoup `mnn×NH×8 B` (+ DM: `mnndm×N×4` + `3·mnndm×NH×8`; + tensor: `9·mnn×NH×8`). For mnn=30, N=10⁷: nlist ≈ 1.2 GB.
`Correlations` (do_sc=Q/Y): `qt` = 3·nq·sc_max_nstep complex (16 B) and `w_block` = 3·(blocks)·nq·nw complex — multi-GB for dense q-meshes/long windows.
`Convolution`: spin_fft + field_fft = 2 × cells·3·NA·M·16 B (C2C double), kernel_fft = cells·9·NA²·16 B, **plus** cuFFT/hipFFT internal workspace of comparable size, **plus** a second plan (`kernel_plan`) kept alive forever.
`MC`: d_state = max_spins·M × sizeof(curandStateXORWOW≈48 B).
The crash mode for big systems is M1 below, not a leak.

### - [x] M1 — Graceful OOM + upfront budget check  (P0, S/M)  — **landed 2026-07-20**
**Files:** `gpu_files/fort_helper.cpp`, `gpu_files/gpuSimulation.cpp`, plus the per-component
estimators in `gpu_files/measurement/memoryMeasurement.{h,cpp}`, `gpu_files/gpuThermfield.{hpp,cpp}`,
`gpu_files/gpuDepondtIntegrator.{hpp,cpp}`, `gpu_files/gpuLatticeConvolutionHamiltonian.{hpp,cpp}`,
`gpu_files/correlations/gpuCorrelations.{hpp,cpp}`.
**Task:** V3 delivers the catch + report. Extend with a **pre-allocation estimate**: before `initiateMatrices` allocates anything, compute the projected device bytes from `N, M, NH, mnn, mnndm`, flags (jtensor/dm/aniso/stt/do_ene), convolution descriptor, and correlation settings (`nq, nw, sc_max_nstep`), compare against `MemGetInfo` free bytes, and if projected > ~90 % of free, abort with a table of the top consumers and actionable hints ("disable eneff via do_ene 0 saves X GB; single precision saves Y GB; reduce sc_max_nstep..."). Keep the estimator as one function so it stays in sync with allocations (unit-test it against `TensorMemoryTracker` totals after a real init: agreement within 5 %).
**Acceptance:** Oversized run aborts before any kernel launch with the budget table; estimator matches tracker within 5 % on three V2 cases with different flags.

> Done 2026-07-20 (commits `c0cb704`, `ba75419`, `0adf7ab`; V3's part-2 summary line folded in).
>
> **Design — per-component, not one function.** The task says "keep the estimator as one function",
> but the allocations live in five places, so a single function would have to reach into all of
> them and would rot exactly like the CMake source list did (CM4). Instead each owner exposes a
> `static estimateBytes(...)` next to its `Allocate` calls — `GpuThermfield`, `GpuDepondtIntegrator`
> (which sums its own thermfield), `GpuLatticeConvolutionHamiltonian`, `GpuCorrelations` — and the
> `initiateMatrices` allocations (Hamiltonian/Lattice/Energies) are covered by file-local helpers
> in `gpuSimulation.cpp`. `computeAndCheckDeviceBudget()` sums them in **one place**. Every estimator
> uses `sizeof(real)` / `sizeof(unsigned int)` / `sizeof(GpuFftComplex)` directly, so the byte math
> is identical to what `TensorMemoryTracker` records and auto-adapts to single vs double precision.
>
> **The 5% self-check runs on every run, not once by hand** (Risk 1). `release()` prints
> `estimate vs measured peak_device() (delta%)`. Added `TensorMemoryTracker::peak_device()` for it.
> Measured pairs on large systems (small regression cases are dominated by fixed overhead and are
> noisy for a ratio, so these use scaled inputs), all four **+0.0%** to displayed precision:
>
> | Case | Distinct flags | Atoms | Estimate | Peak |
> |---|---|---|---|---|
> | scalar SD | ncoup, none | 2 000 000 | 1.872 GB | 1.872 GB |
> | tensor | `do_jtensor` | 1 920 000 | 1.229 GB | 1.229 GB |
> | multi-ensemble | `Mensemble 4` | 500 094 × 4 | 0.948 GB | 0.948 GB |
> | energy | `do_ene 1` | 2 000 000 | 1.872 GB | 1.872 GB |
>
> **Budget abort demonstrated.** With a helper process holding all but 0.698 GB of a 16 GB card, a
> 2 M-atom run aborts *before any device allocation* (no `TENSOR MEMORY` report is printed) with a
> sorted top-consumers table (Hamiltonian 1.264 GB / Lattice 352 MB / integrator 256 MB / ...) and
> hints (single precision, fewer ensembles, smaller system). Threshold is estimate > 90 % of
> `MemGetInfo` free; the banner states the estimate is approximate and excludes FFT workspace.
>
> **Escape hatch:** `UPPASD_GPU_SKIP_BUDGET=1` bypasses the abort (projected line is tagged
> `[OVER BUDGET - bypassed]`); the run then falls through to V3's graceful OOM banner rather than a
> raw crash. Documented in the abort message itself (Risk 3).
>
> **Two gates had to mirror runtime behaviour, or the estimate over-counts:**
> - *Convolution:* counted only when it will actually engage — `willUseConvolution()` replicates
>   `canUseLatticeConvolution`'s parameter checks (`do_gpu_convolution`, `N == N1·N2·N3·NA`,
>   `NH == NA`). Without this, a `do_gpu_convolution Y` run that falls back to neighbour lists (the
>   common `NH != NA` case) was estimated +41 % high.
> - *Correlations:* the raw `Flags.do_gpu_correlations` is unreliable — line 54 does
>   `static_cast<bool>(*FortranData::do_gpu_correlations)`, and `'N'` (0x4E) casts to **true**. The
>   estimate is gated on the real trigger `*FortranData::do_gpu_correlations == 'Y'` (what
>   `CorrelationFactory` uses). That phantom `coord` (3·N reals) was the entire residual error;
>   removing it took every case from ~+2.6 % to +0.0 %. **The underlying char→bool bug at
>   `gpuSimulation.cpp:54` is left as-is — out of M1's scope, flagged here.**
>
> **Not addressed / observations:**
> - `sc_dm_uniaxial` is unchanged: its failure is a Fortran `read_exchange` bounds error in the CPU
>   reference binary before any GPU init (see V3 note), which the budget check never reaches.
> - M2 is ticked "verified landed" but `initiateMatrices` still allocates `eneff` **and** `extfield`
>   unconditionally in the tree — M2's conditional-skip/alias is not active. The estimator counts
>   all seven 3NM arrays accordingly, which is why agreement is exact; if M2 is genuinely completed
>   later, `latticeBytes`/`hamiltonianBytes` must drop those terms to stay at 0 %.
> - HIP path unvalidated (no AMD hardware); the `GPU_MEM_GET_INFO` wrapper and estimators compile
>   for HIP but were only exercised under CUDA.
> - Preserved **11 PASS / 1 FAIL** on `cuda-debug`, `cuda-fastcopy`, and `cuda-ptds`.

### - [x] M2 — Eliminate always-allocated arrays that are conditionally needed  (P1, M)  — **partial landed 2026-07-21 (eneff + mmom0/mmom2); extfield + blocal deferred**
**Files:** `gpu_files/gpuSimulation.cpp` (`initiateMatrices`), `gpu_files/gpuHamiltonianCalculations.cpp`, `gpu_files/gpuStructures.hpp`.
**Task:**
1. **`eneff` (3NM):** only meaningful when energies are measured (`do_ene>0` or GPU-MC). When not needed, don't allocate; make the Heisge kernels write `eneff` only under the `Measure` template flag (after U2) — interim: point `gpuLattice.eneff` at `gpuLattice.beff` (alias) when `do_ene==0 && !MC`, since the kernels currently write identical values in the non-measure path for the non-aniso variants — **verify per-variant before aliasing** (the aniso variants write different `*_en` prefactors; for those, allocation stays whenever aniso+MC/energy).
2. **`extfield` (3NM):** detect the all-zero external field at upload (scan the Fortran array once on host); if zero, skip allocation and pass `nullptr`/a flag so kernels skip the loads (compile-time after U2, runtime flag before).
3. **`blocal` (3NM):** removed by kernel fusion — defer to PF5, but note the dependency here.
4. **`mmom0`, `mmom2`:** only used by `mompar!=0` paths and MC `moms`; allocate conditionally.
**Acceptance:** V2 passes with all flag combinations; tracker report for a plain Heisenberg SD run (no field, no energy) shows ≥ 3 × 3NM·8 B reduction vs baseline.

> **History.** Ticked "verified landed" by `28f4b5a` with nothing in the tree (third over-tick
> from that pass, after V3 and PR1–PR7); un-ticked 2026-07-20 with the analysis below; items 1 and
> 4 then implemented 2026-07-21 by first plumbing the run mode into init. Items 2 (`extfield`) and
> 3 (`blocal`) remain deferred — see the end.
>
> **Why it was blocked as specced (the analysis that shaped the fix):**
> - `gpuSim_initiateMatrices()` is called from **both** `sd_driver.f90` (1087, 1185) and
>   `mc_driver.f90` (578, 662), and the run mode (SD vs MC) is only passed later, as `whichsim`,
>   to `gpuRunSimulation`. So `initiateMatrices` **cannot tell SD from MC**, and the MC kernels
>   (`gpuMetropolis*::moms`, `MCSweep`) write `eneff`, `mmom0`, `mmom2` unconditionally. Skipping
>   any of them at allocation time would corrupt/crash an MC run. This blocks items 1 and 4.
> - `eneff` is in fact read **only** by the MC kernels (`gpuMetropolis.cpp:74-76`,
>   `gpuMetropolis_bruteforce.cpp:67-69`) — never by energy measurement, which builds `energyM`
>   directly in the `Measure` branch of `Heisge`. So the task's "`do_ene`" condition is not even
>   the right trigger; `eneff` is dead for *all* SD runs and needed only by MC. The honest fix is
>   to move `eneff` allocation into the MC path (or gate it behind a mode flag / U2 template),
>   **not** to alias it — aliasing to `beff` is unsafe anyway because the aniso `Heisge` writes
>   `eneff != beff` (`gpuHamiltonianCalculations.cpp:342-344`, `ax_en != ax`) and
>   `gpuMultiscaleField` *adds* to both (double-add on aliased memory).
> - `extfield` (item 2) is the only mode-independent item, but the field-present path has **zero
>   regression coverage** (no GPU case sets a nonzero field; `SCsurf`'s `hfield` is `0 0 -0.000`
>   and that case fails on the CPU side anyway). Shipping a kernel guard that skips the `ext_f`
>   load in `Heisge`, `HeisgeJijElement`, and the convolution `apply` kernel without a with-field
>   reference case would add an untested hot-path branch.
>
> **What landed (2026-07-21, commit `0e6f1d2`) — items 1 and 4, via run-mode plumbing:**
> - **Mode plumbing.** `gpusim_initiatematrices_` now takes `int *is_mc`;
>   `initiateMatrices(int is_mc)` stores `runIsMC`. `sd_driver.f90` passes `0`, `mc_driver.f90`
>   passes `1`. Safe per-phase because each phase does its own
>   `initiateConstants → initiateMatrices → run → release()` cycle (verified in both drivers), so
>   there is no phase-ordering trap — an SD phase and an MC phase in the same run allocate
>   independently. The Fortran↔C binding is an implicit interface (no `bind(c)` block), so a missed
>   call site degrades to a garbage `is_mc` — usually nonzero → allocate-everything (the safe
>   fallback), except MC callers which must pass 1; all four sites updated.
> - **`eneff` (item 1)** — aliased to `beff` (shallow `GpuTensor` copy, shares the device buffer)
>   when `!runIsMC && do_aniso == 0`; allocated normally otherwise. Aliasing is provably safe there:
>   `Heisge` writes `eneff == beff` when `!HasAniso`, `HeisgeJijElement` and the convolution kernel
>   always write them equal, `eneff` is never *read* in SD, and `gpuMultiscaleField` (the double-add
>   risk) is a no-op stub today (`canUseMultiscaleDipole` returns false — noted as a dependency).
>   `release()` skips the free for the aliased view; the tracker never double-counts.
> - **`mmom0`/`mmom2` (item 4)** — allocated only when `runIsMC || mompar != 0`; the four
>   copy sites (`copyFromFortran`/`copyToFortran`) and the `gpuHasNoData` checks are guarded with
>   `!empty()` so a plain SD run neither allocates nor touches Fortran's `mmom0`/`mmom2`.
> - **Estimator kept in lockstep** (same commit): `latticeBytes(F, P, is_mc)` counts `eneff` only
>   for `is_mc || aniso` and `mmom0/mmom2` only for `is_mc || mompar`, matching the allocations.
>
> **Evidence (measured, not asserted):**
> - Plain SD, 2 000 000 atoms: peak device **1.872 → 1.792 GB**, an 80 MB drop =
>   `eneff` 48 MB + `mmom0` 16 MB + `mmom2` 16 MB. M1 self-check **+0.0 %** (estimate 1.792 =
>   peak 1.792). The MC branch's pre-alloc estimate is 1.872 GB (eneff+mmom counted), confirming the
>   `is_mc` gate.
> - Regression **11 PASS / 1 FAIL bitwise** (`abs=0.00e+00`) — SD cases exercise the alias/skip
>   paths, MC (`bcc_mc`, `single_spin_mc`, the latter MC+aniso) exercise full allocation; a large MC
>   run reached the measurement phase with no crash. compute-sanitizer: **0 errors** across
>   memcheck/racecheck/synccheck/initcheck on the aliasing case (`bcc_scalar`).
>
> **Still deferred, with reasons:**
> - **`extfield` (item 2)** — mode-independent but needs a guard in three hot kernels (`Heisge`,
>   `HeisgeJijElement`, convolution `apply`) *and* a new nonzero-field GPU regression case (no
>   existing case sets a field, so the field-present branch is untested). Left out to avoid an
>   untested hot-path change.
> - **`blocal` (item 3)** — deferred to PF5 (kernel fusion) as the task states. So this delivers
>   `eneff (3NM) + mmom0 + mmom2 (2×NM)`, i.e. ~1.67 × 3NM, short of the acceptance's "≥ 3 × 3NM"
>   which assumed `blocal` and `extfield` too. Acceptance partially met; recorded honestly rather
>   than ticked clean.

### - [ ] M3 — Convolution mode: drop the sparse Hamiltonian from device  (P1, M; depends on **CV2**)
**Files:** `gpu_files/gpuSimulation.cpp`, `gpu_files/gpuHamiltonianCalculations.cpp`.
**Dependency corrected 2026-07-20:** this item read "depends on CV4", but CV4 is the stream/plan
association audit and does not affect which steps the convolution backend serves. The blocker is
**CV2 — Energy from the convolved field**, which is still open. M3 is therefore still blocked;
ticking CV4 does not release it.
**Task:** Once CV2 makes the convolution backend serve `measure=true` steps too, SD runs on the convolution backend no longer need `nlist/ncoup/dmvect/dmlist` resident (they're still needed for MC and as fallback). Add a mode where, if the convolution kernel builds successfully, the sparse arrays are freed after kernel construction (they are consumed during `buildIsotropicDmKernel`/`buildTensorKernel`). Keep a flag so a fallback to sparse triggers a clear error instead of touching freed memory.
**Acceptance:** Large Bravais SD run shows nlist/ncoup absent in the tracker; V2 convolution cases still pass including energy measurement.

### - [x] M4 — Correlations: remove the `sc_max_nstep` dimension from device (base S(q,t))  (P1, L)  — **landed 2026-07-21; new host transform bit-identical to prior device result (T=0 bccFe_sqw, all 40k sqw values, max Δ = 0)**
**Files:** `gpu_files/correlations/{correlation_types.h, gpuCorrelations.{hpp,cpp}}`.
**Task:** `qt(3, nq, sc_max_nstep)` stores the full time series on device only so the t→ω transform can run at flush. Implement **chunked streaming with the transform finished on host** (option (b) below); option (a) was considered and rejected — see note.

**Chunked streaming (implemented for the base `sc` path):**
- Device `sc.qt` is a rolling `(3, nq, K)` chunk (`K = min(SC_CHUNK, sc_max_nstep)`, `SC_CHUNK = 64`). The dynamic sampler writes slot `t_cur % K`; each slot is zeroed before reuse (the FinalSum kernels accumulate with `+=`).
- Every `K` samples the full chunk is copied `D2H` into the host `m_kt(3, nq, sc_max_nstep)` buffer at the current time offset (column-major → a chunk of consecutive time slices is contiguous, single memcpy). The final partial chunk is flushed at `flushCorrelations`.
- The windowed t→ω DFT is done **on host** (`transform_kt_to_kw_host`, an exact port of `GPUSwSum`/`corr_kernel_time`: `m_kw(c,q,ω) = Σ_t m_kt(c,q,t)·exp(+i·t·dt(t)·w(ω))·win(t+1, sc_max_nstep)`, OpenMP over ω,q). This also matches Fortran better than the old device tree-reduction did, since both are sequential sums.
- Device `sc.qw`, `sc.w_block`, and the device `dt`/`w` scratch are **dropped** for the base path. Base dynamic-correlation device footprint drops from `3·nq·(sc_max_nstep + nw + blocks·nw)` complex to `3·nq·K`.

**Note — why not (a) incremental device DFT:** it keeps `qw(3, nq, nw)` on device, which for our regime (`nw ≈ 1–2·sc_max_nstep`) is no smaller than `qt`, so ~no win; and it *discards* `m_kt`, which is a required output (`sqt` file) for `do_sc T`/`Y`. The maths is exact either way (the transform is a linear windowed DFT, no origin-averaging), so (a) loses no statistics — it simply doesn't fit the memory regime or the output contract.

**Still open (not done here):** the projected paths (`sc_proj`/`sc_projch`, `qt(3, nproj, nq, sc_max_nstep)`) still transform on device and keep the full series; apply the same chunked-streaming treatment there. Also `w_block` right-sizing from actual launch geometry is now moot for the base path (buffer removed) but still applies to the projected paths.
**Acceptance:** Dynamic-correlation case (`do_gpu_correlations Y`, `do_sc Q`/`Y`) matches the pre-change device result and Fortran `m_kw` within tolerance; tracker shows the sc_max_nstep-proportional base allocation gone.

### - [x] M5 — FFT workspace hygiene  (P2, S)  — **verified landed 2026-07-20**
**Files:** `gpu_files/gpuLatticeConvolutionHamiltonian.cpp`.
**Task:** (1) Destroy `kernel_plan` immediately after the kernel FFT is done in `build*Kernel` — it is never needed again and holds a workspace. (2) Investigate `cufftSetAutoAllocation(plan, 0)` + `cufftSetWorkArea` sharing one work area between forward/backward use of the single `plan` (they already share the plan — verify no extra workspace duplication) and, if kernel_plan must persist for rebuilds, share the work area with `plan`. (3) After CV3 (R2C), sizes drop further.
**Acceptance:** Tracker/MemGetInfo shows workspace reduction on a convolution run; V2 convolution cases pass.

### - [x] M6 — 32-bit index overflow audit  (P2, S)  — **verified landed 2026-07-20**
**Files:** `gpu_files/gpuParallelizationHelper.{hpp,tpp}`, `gpu_files/gridHelper.hpp`, `gpu_files/gpuLatticeConvolutionHamiltonian.cpp`.
**Task:** `op.NM3 = N*M*3` (unsigned int) overflows at N·M > 1.43e9; `X*Y` in `dim2d` big-grid path similarly; `pack_spins_to_fft` total and cuFFT `int` dims cap the convolution grid. Convert the helper structs' size fields and index math to `size_t`/`unsigned long long` where kernels permit (index arithmetic in kernels can stay 32-bit when safe — gate with a checked cast at launch: if the product exceeds `UINT_MAX`, error out with a clear message rather than wrapping). Document the hard limits that remain (cuFFT int dims).
**Acceptance:** A mock launch with N·M·3 > 2³² triggers the checked-cast error, not silent wraparound (unit test the check function on host).

### - [ ] M7 — Precision-driven halving  (P1, tracked in Phase 4)
Single/mixed storage (PR1–PR6) halves every streaming array and the FFT buffers. Listed here for completeness of the memory story; implemented in Phase 4.

---

## Phase 3 — CMake modernization  ★ user emphasis 2

Goal: `cmake --preset <name> && cmake --build` works out of the box on (a) a laptop with a consumer NVIDIA GPU, (b) a laptop/desktop with an RDNA/CDNA AMD GPU, (c) LUMI/Dardel-class HPC with Cray PrgEnv + ROCm, (d) NVIDIA HPC with modules — with no hand-editing of CMakeLists. Backends and precision are cache options; nothing is hardcoded.

### - [x] CM1 — Option surface and target hygiene  (P1, M)
**Done** (`b426193`). `UPPASD_GPU_BACKEND` (OFF|CUDA|HIP) and `UPPASD_PRECISION`
(DOUBLE|SINGLE|MIXED) are cache STRINGs resolved before the binary is named; `USE_CUDA`/
`USE_HIP` survive as deprecated aliases that emit `message(DEPRECATION ...)` and are derived
back from the resolved backend, so `gpu_files/*/CMakeLists.txt` and CI keep working until
CM2–CM4 convert them. Disagreeing spellings (`UPPASD_GPU_BACKEND=CUDA` + `USE_HIP=ON`) are a
hard error rather than a silent pick. `UPPASD_PRECISION=MIXED` is a `FATAL_ERROR` — the mode is
specified against `storage_t`/`accum_t`, which PR1 has not created yet; SINGLE maps to the
existing `SINGLE_PREC` macro. `USE_FAST_COPY` promoted to a real option (it had rotted
undetected precisely because it was not one); verified end-to-end, not just on the compile line.
`CUDA_V`/`HIP_V`/`ON_LUMI` moved from `add_compile_definitions` to
`target_compile_definitions(${UppASD_LIB} PRIVATE ...)`; the five remaining global ones
(`USE_OVF`, `VSL`, `USE_MKL_FFT`, `USE_FFTW`) are non-GPU and belong to CM6. `-DDEBUG` removed
from the CUDA, HIP and CXX flag strings — its only consumers are under `gpu_files/old/`, which
is not built.

**Deviation from the task text:** it asks for `cxx_std_17`. The repo already sets
`CMAKE_CXX_STANDARD`/`CMAKE_CUDA_STANDARD` to **20** with `_REQUIRED TRUE`; downgrading would be
a regression, so 20 was kept and `CMAKE_HIP_STANDARD`/`_REQUIRED` added alongside, which were
missing.

**Verified:** `-DUPPASD_GPU_BACKEND=CUDA`, `=OFF`, and legacy `-DUSE_CUDA=ON` all configure
clean; `=OPENCL`, `MIXED`, and the CUDA/HIP conflict all fail with the intended message;
`CUDA_V`+`SINGLE_PREC` confirmed present in the `CUDA_DEFINES`, `CXX_DEFINES` *and*
`Fortran_DEFINES` lines of `asdlib`'s `flags.make`. `build_ab`, `build_fastcopy`, `build_ptds`
all rebuild and still give 11 PASS / 1 FAIL (`sc_dm_uniaxial`, pre-existing), as does a fresh
tree configured with `-DUSE_FAST_COPY=ON`. **HIP was not built — no AMD hardware here.**

**Files:** top-level `CMakeLists.txt`, `source/CMakeLists.txt` (if present), `source/gpu_files/CMakeLists.txt`, `source/gpu_files/{measurement,correlations}/CMakeLists.txt`.
**Task:** Replace `USE_CUDA`/`USE_HIP` booleans with a single `UPPASD_GPU_BACKEND` cache STRING = `OFF|CUDA|HIP` (keep the old booleans as deprecated aliases mapping onto it, with a warning). Add `UPPASD_PRECISION` = `DOUBLE|SINGLE|MIXED` (default DOUBLE; consumed by Phase 4). Convert all `add_compile_definitions` / global `CMAKE_CXX_FLAGS` string surgery for GPU defines into `target_compile_definitions(${UppASD_LIB} PRIVATE CUDA_V …)` and `target_compile_features(${UppASD_LIB} PRIVATE cxx_std_17)` (and `cuda_std_17` / `hip_std_17` via `CMAKE_CUDA_STANDARD`/`CMAKE_HIP_STANDARD` on the target). Remove `-DDEBUG` injection (nothing consumes it); ensure `NDEBUG` follows build type as normal so device asserts vanish in Release (audit P7).
**Acceptance:** `cmake -DUPPASD_GPU_BACKEND=CUDA` and `=HIP` and `=OFF` all configure+build; `grep -r add_compile_definitions` in GPU sections returns nothing global.

### - [x] CM2 — CUDA path  (P1, M)
**Done** (`5515dba`). Legacy `find_package(CUDA)`, the `CMAKE_VERSION VERSION_LESS 3.17` fork
and the whole `FindCUDA/select_compute_arch` detection block are gone; only
`enable_language(CUDA)` + `find_package(CUDAToolkit REQUIRED)` remain. GPU builds now hard-fail
below CMake 3.24 (what `native` needs); CPU-only builds keep the 3.13 project minimum.
`CMAKE_CUDA_ARCHITECTURES` defaults to `native` when unset. The manual
`-I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}` splicing into `CMAKE_CUDA_FLAGS`/`CMAKE_CXX_FLAGS`
is replaced by `target_include_directories`; the CCCL probes were kept (they solve a real
layout problem) but now use `${UppASD_LIB}` rather than a hardcoded `asdlib`.

**Two findings worth carrying forward:**

- **NVTX was a hidden dependency of the blanket `CUDA::toolkit` link.** `gpuEventPool.cpp` and
  `measurement/fortranMeasurement.cpp` include `nvToolsExtCuda.h` and call `nvtxRangePush` /
  `nvtxNameOsThread`. Naively linking only `cudart/curand/cufft` as the task text says would
  have dropped it. Handled: `CUDA::nvtx3` if present, else `CUDA::nvToolsExt` (renamed across
  toolkit versions). On CUDA 13.3 here `CUDA::nvtx3` is an INTERFACE (header-only) target, so
  nothing appears on the link line — verified by probe, not assumed.
- **The `-D__${CMAKE_Fortran_COMPILER_ID}__` flag splicing was dead** and was deleted: zero
  occurrences of `__GNU__`/`__Intel__` anywhere under `source/`, including `gpu_files/old/`.

**On the broken `MATCHES "^[__86+PTX]"` regex (audit B2):** confirmed it does not exist in this
tree — the previous handover's doubt was right. The block that was deleted used
`CUDA_DETECT_INSTALLED_GPUS` with a hardcoded fallback of `80`. Do not claim that regex as
fixed.

**Verified:** `native` configures, builds, and bakes `sm_86` — matching the actual card (RTX
A4000, cc 8.6), checked with `cuobjdump --list-elf`. Note the task's acceptance text ("not 86")
assumed a different machine; 86 *is* the correct answer here, and the old hardcoded fallback
was 80. CI-style `-DCMAKE_CUDA_ARCHITECTURES=80` also configures and builds, baking `sm_80`, so
V1 keeps working — and `.github/workflows/gpubuilds.yml` already passes an explicit arch, so
the `native` default never reaches the GPU-less runner. `build_ab`, `build_fastcopy`,
`build_ptds` and a fresh `native` tree all give 11 PASS / 1 FAIL.

**Files:** top-level `CMakeLists.txt`.
**Task:** Require CMake ≥ 3.24 for GPU builds. Use only `enable_language(CUDA)` + `find_package(CUDAToolkit REQUIRED)` (delete the legacy `find_package(CUDA)` fallback and `FindCUDA/select_compute_arch` entirely, including the broken `MATCHES "^[__86+PTX]"` regex — audit B2). Default `CMAKE_CUDA_ARCHITECTURES=native` when not set by the user; document overriding for cross-compilation (`-DCMAKE_CUDA_ARCHITECTURES=80;90`). Link explicitly: `target_link_libraries(${UppASD_LIB} PRIVATE CUDA::cudart CUDA::curand CUDA::cufft)` — today curand comes in implicitly. Replace the manual `-I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}` flag splicing with the imported targets (they carry includes).
**Acceptance:** Clean build on a machine with only the toolkit (no GPU) using `-DCMAKE_CUDA_ARCHITECTURES=80`, and on a GPU machine with `native`; `nvcc` line in verbose build shows the arch actually detected, not 86.

### - [ ] CM3 — HIP path  (P1, M)
**Files:** top-level `CMakeLists.txt`, `source/gpu_files/CMakeLists.txt`.
**Task:** Use `enable_language(HIP)` with `CMAKE_HIP_ARCHITECTURES` (default: let CMake/hipcc detect; presets pin gfx90a/gfx942 for LUMI-class, gfx11xx for RDNA laptops). **Do not set `CMAKE_CXX_COMPILER`/`CMAKE_HIP_COMPILER` inside CMakeLists** — remove those lines; the compiler comes from the environment/toolchain (`CC=hipcc` or CMake's ROCm detection). Replace the hardcoded `find_library(HIPFFT_LIBRARY PATHS /opt/rocm/lib …)` with `find_package(hipfft CONFIG REQUIRED)` and `find_package(hiprand CONFIG REQUIRED)` (+ `find_package(hip CONFIG REQUIRED)`), linking `hip::hipfft hip::hiprand hip::host`. ROCm config packages are found via `CMAKE_PREFIX_PATH=$ROCM_PATH` — document that one env var instead of baking paths.
**Acceptance:** Configure+build inside `rocm/dev-ubuntu` container with only `CMAKE_PREFIX_PATH` set; no absolute `/opt/rocm-*` strings remain in any CMakeLists.

### - [x] CM4 — Single source of truth for GPU sources; per-backend language  (P1, M)
**Done** (`e32c11b`). The top-level `set_source_files_properties(... LANGUAGE CUDA)` blob is
deleted. Each of the three `gpu_files` directories now declares its own `GPU_KERNEL_SOURCES`
and calls a new `uppasd_mark_gpu_sources()` helper, so adding a `.cpp` means editing exactly
one file. All commented-out source entries are gone.

**The non-obvious part — why the duplication existed at all.** Source-file properties are
directory-scoped, and `${UppASD_LIB}` is defined in the *top-level* directory. A plain
`set_source_files_properties()` inside `source/gpu_files/` is therefore **silently ignored** for
that target, which is why the list had to be repeated at top level. The helper uses
`TARGET_DIRECTORY ${UppASD_LIB}` (CMake >= 3.18, and CM2 already requires 3.24 for GPU builds)
to aim the property at the directory that owns the target. **This also means
`correlations/CMakeLists.txt`'s `set_source_files_properties(correlation_kernels.cpp PROPERTIES
LANGUAGE HIP)` was a no-op** — the HIP build was compiling that file as plain C++. Unverifiable
here without AMD hardware, but worth knowing when CM3 is picked up.

**Drift found and removed:** the blob listed `${SRC}/gpu_files/gpuMeasurement.cpp`, which **does
not exist** (the real file is `measurement/gpuMeasurement.cpp`; the `gpu_files/` one is
commented out of the source list). `set_source_files_properties` does not error on missing
files, so this sat there harmlessly and invisibly.

**Verified** by dumping `compile_commands.json` and classifying every `gpu_files` entry by
whether it invokes `nvcc`: exactly 20 sources compile as CUDA and 8 as plain C++ — matching the
old blob's 20 real entries and the 8 files it omitted, so the split is byte-for-byte the
previous behaviour. Confirmed independently that none of those 8 contain `__global__` or
`<<<`. Full rebuild + 11 PASS / 1 FAIL on `build_ab`, `build_fastcopy`, `build_ptds` and a
fresh tree.

**Files:** `source/gpu_files/CMakeLists.txt`, top-level `CMakeLists.txt`.
**Task:** Today the GPU file list exists twice (subdir `target_sources` + a drifted `set_source_files_properties(... LANGUAGE CUDA)` blob at top level) with commented-out entries disagreeing. Restructure: `source/gpu_files/CMakeLists.txt` defines one variable/OBJECT-library `uppasd_gpu` with the complete source list; a single loop sets `LANGUAGE` to `CUDA` or `HIP` on the kernel-containing sources depending on `UPPASD_GPU_BACKEND` (after U1 lands there is exactly one list — until then, include the per-backend thermfield/correlation file selection here and nowhere else). Delete the top-level `set_source_files_properties` block and all commented-out source entries. Headers/`.tpp` need not be listed as sources (list them via `target_sources(... FILE_SET HEADERS)` if IDE visibility is wanted).
**Acceptance:** Adding a new `.cpp` under gpu_files requires editing exactly one file; both backends build.

### - [x] CM5 — Presets  (P1, S)
**Done** (`d5ebd17`). `CMakePresets.json` now provides `cpu`, `cpu-debug`, `cuda`,
`cuda-debug`, `cuda-single`, `cuda-fastcopy`, `cuda-ptds`, `hip-cdna`, `hip-rdna` and `lumi`,
all inheriting a hidden `base` and setting cache variables only — no paths. Matching
`buildPresets` mean `cmake --build --preset <name>` works. `docs/BUILDING_GPU.md` covers the
four target platforms, the `CMAKE_PREFIX_PATH=$ROCM_PATH` note and the LUMI module list.
`.github/workflows/gpubuilds.yml` (V1) now configures via `--preset cuda` / `--preset hip-cdna`.

**The hand-configured build dirs collapse into presets**, which was the real measure of done:
`build_fastcopy` → `cuda-fastcopy`, `build_ptds` → `cuda-ptds`, `build_ab` → `cuda-debug`. The
old `debug`/`release`/`debugCuda`/`releaseCuda` presets are kept as deprecated aliases pointing
at their original `binaryDir`s, so existing trees and habits do not break.

**Verified:** `cpu`, `cuda`, `cuda-fastcopy`, `cuda-ptds` each configure *and build* from
scratch; the three GPU ones give 11 PASS / 1 FAIL. The flags are confirmed to actually reach
the compiler rather than being decorative — `-DUSE_FAST_COPY` appears in `cuda-fastcopy`'s
`flags.make`, `--default-stream per-thread` in `cuda-ptds`', and neither in plain `cuda`'s. The
exact CI invocation (`cmake --preset cuda -DCMAKE_CUDA_ARCHITECTURES=80 -DBUILD_TESTING=OFF
...`) was run locally and builds, confirming a command-line `-D` still overrides a preset cache
variable. `hip-cdna` and `hip-rdna` resolve the backend correctly and then fail at
`find_package(HIP)` — the expected outcome with no ROCm present, not a preset error. **Neither
HIP preset has been built.** Note the LUMI module versions in the doc are unverified; the doc
says so.

**Files:** new `CMakePresets.json` at repo root, `docs/` snippet.
**Task:** Provide configure presets: `cpu`, `cuda` (native arch, Release), `cuda-debug`, `hip-cdna` (gfx90a;gfx942), `hip-rdna` (gfx1100;gfx1101;gfx1102), `lumi` (hip-cdna + hints matching PrgEnv-cray/amd: `CMAKE_Fortran_COMPILER=ftn` etc. — keep minimal, rely on modules), plus precision variants via `UPPASD_PRECISION`. Each preset sets only cache variables, never paths. Add a short `docs/BUILDING_GPU.md`: three commands per platform, the LUMI module list, and the `CMAKE_PREFIX_PATH=$ROCM_PATH` note.
**Acceptance:** `cmake --preset cuda && cmake --build --preset cuda` works on a stock CUDA box; `cmake --preset hip-cdna` configures in the ROCm container; CI (V1) switched to use the presets.

### - [ ] CM6 — BLAS/LAPACK and "local library shenanigans"  (P2, S)
**Files:** top-level `CMakeLists.txt`.
**Task:** Keep `find_package(BLAS REQUIRED)`/`find_package(LAPACK REQUIRED)` but respect `BLA_VENDOR` (document `-DBLA_VENDOR=Intel10_64lp` for MKL, `=OpenBLAS`, Cray's libsci needs nothing under PrgEnv). Remove the dead commented MKL block. Ensure no `link_directories` remain. If site-local hacks were needed historically, capture each as either a preset cache variable or a documented environment expectation — the CMakeLists itself stays generic.
**Acceptance:** Builds against OpenBLAS (Ubuntu), MKL (`BLA_VENDOR`), and cray-libsci (LUMI preset) without edits.

---

## Phase 4 — Flexible precision: FP64 / FP32 / MIXED  ★ user emphasis 4

**Design decision (adopt this):** precision is a **build-time** choice (`UPPASD_PRECISION=DOUBLE|SINGLE|MIXED`), producing distinctly named binaries (`sd`, `sd_sp`, `sd_mp`). Runtime switching would require compiling every kernel twice and doubles binary/testing surface for little benefit on HPC where builds are per-machine anyway. The Fortran side stays `dblprec` (double) in **all** modes; conversion happens at the upload/download boundary. This simultaneously resolves audit I2 (the current SINGLE_PREC build reinterprets doubles as floats across the bridge — i.e. today SINGLE_PREC is broken, not merely inaccurate).

Semantics of the three modes:
- **DOUBLE:** storage fp64, accumulation fp64 (today's behavior).
- **SINGLE:** storage fp32, accumulation fp32. Fastest; for exploratory dynamics, large-N equilibration, consumer GPUs (where fp64 is 1/32–1/64 rate).
- **MIXED:** storage fp32 (all 3NM streaming arrays, couplings, FFT buffers), accumulation fp64 (all reductions: energies, magnetization sums, Binder cumulants, correlation accumulators), plus fp64 for the handful of scalar-sensitive spots listed in PR5. Recommended default for production on consumer/RDNA hardware — memory and bandwidth of fp32 with fp64 statistics.

### - [ ] PR1 — Type foundation  (P1, M) — **partial: binary `SINGLE_PREC` only, not the specified three-mode design (see Phase 4 note)**
**Files:** `gpu_files/real_type.h`, `gpu_files/base.hpp`, CMake (CM1's `UPPASD_PRECISION`).
**Task:** Introduce two typedefs controlled by the CMake option: `using storage_t = float|double;` and `using accum_t = float|double;` per the table above, and complex counterparts `storage_cplx`, `accum_cplx`. Keep `using real = storage_t;` so the existing code compiles unchanged; new/edited code uses the explicit names. Map the FFT type macros (`GpuFftComplex`, `GPUFFT_EXEC_C2C` → C2C vs Z2Z; after CV3, R2C vs D2Z) off `storage_t`. Delete the `ON_LUMI` thrust fork opportunistically if U1 has landed.
**Acceptance:** All three modes compile for both backends (numbers may be wrong until PR2 — say so in the PR).

### - [ ] PR2 — Boundary conversion layer  (P0-within-phase, M; extends F7) — **partial: binary `SINGLE_PREC` only, not the specified three-mode design (see Phase 4 note)**
**Files:** `gpu_files/tensor.hpp`, `gpu_files/gpuSimulation.cpp`.
**Task:** Generalize F7's staging helper to `upload_convert(const Tensor<double,dim>& src, GpuTensor<storage_t,dim>& dst, bool transpose=false)` and `download_convert(...)` (convert in a pinned staging buffer sized max over uses, reused across calls; async on `copyStream` where callers permit). Route **every** `copyFromFortran`/`copyToFortran` transfer through it. In DOUBLE mode it degenerates to the current memcpy path (zero overhead — specialize). Host-side `cpuLattice`/`cpuHamiltonian` tensors keep aliasing the Fortran doubles; only device tensors change type.
**Acceptance:** DOUBLE mode bit-identical to pre-change (V2). SINGLE mode: V2 with relaxed tolerances (see PR7) passes; no `reinterpret`-style type punning remains (`grep` for `real*` casts on `FortranData` pointers).

### - [ ] PR3 — fp64 accumulation everywhere it matters (MIXED)  (P1, M) — **partial: binary `SINGLE_PREC` only, not the specified three-mode design (see Phase 4 note)**
**Files:** `gpu_files/gpuHamiltonianCalculations.cpp` (energy path — coordinate with F1/U2), `gpu_files/measurement/kernels.cpp` (all `block_reduce_sum_1d` users), `gpu_files/correlations/*` (S(q) accumulators `q/qt/qw` and block scratch), `gpu_files/measurement/gpuMeasurement.cpp` buffers that hold running sums.
**Task:** Change reduction accumulator types and the persistent accumulation buffers (`energyM`, `emomMEnsembleSums`, correlation `q/qt/qw`) from `real` to `accum_t` / `accum_cplx`. Per-element loads stay `storage_t`; the promotion happens at the `+=`. Binder cumulant fourth powers especially need fp64 (catastrophic cancellation in `⟨m⁴⟩−3⟨m²⟩²`-type combinations at fp32).
**Acceptance:** MIXED-mode long run (10⁵ steps, finite T): energy drift and Binder cumulant match DOUBLE within statistical error, where a pure-SINGLE control visibly drifts. Add this comparison as a (manual/nightly, not per-PR) test script.

### - [ ] PR4 — RNG per precision  (P1, S) — **partial: binary `SINGLE_PREC` only, not the specified three-mode design (see Phase 4 note)**
**Files:** thermfield (post-U1), `gpu_wrappers.h`, `gpu_files/gpuMetropolis.cpp`.
**Task:** Under SINGLE/MIXED generate `curandGenerateNormal` (float) into the fp32 field; keep the even-length fix from F2. In `MCSweep`, replace the hardwired `GPU_NORMAL_DOUBLE`/`GPU_RAND_UNIFORM_DOUBLE` with precision-dispatched macros (`GPU_NORMAL_REAL`, `GPU_RAND_UNIFORM_REAL`). MC acceptance comparison (`exp(βΔE) < u`): compute `βΔE` and the `exp` in `accum_t` (double under MIXED) — it's one scalar per attempt, cost-free.
**Acceptance:** SINGLE/MIXED MC reproduces Fortran ⟨m(T)⟩ curve on the standard test system across 5 temperatures within statistics.

### - [ ] PR5 — Numerically sensitive spots under fp32  (P1, M) — **partial: (a)'s versine fix landed; (a) rest, (b), (c) blocked on PR1**

**(a) versine, done** (`9b1011c`). `1 − cos θ` is now `2·sin²(θ/2)` via a `versine()` device
helper, applied at both rotation sites in `gpuDepondtIntegrator.cpp` — the standalone `rotate`
kernel and the fused `EvolveFirst` predictor. Measured: **no change at all** in double
precision — `kagome_tensor` still reports exactly `abs=1.00e-11, rel=3.62e-09` and every other
comparison exactly `0.00e+00`, across `build_ab`/`build_fastcopy`/`build_ptds`. No tolerance was
touched. So this buys nothing today and is purely insurance for SINGLE/MIXED; do not cite it as
a measured accuracy gain. The rest of (a) needs `accum_t`, i.e. PR1, i.e. CM1's
`UPPASD_PRECISION` — now landed, so PR1 is unblocked. (b) and (c) still wait on PR1.

**Files:** `gpu_files/gpuDepondtIntegrator.cpp`, `gpu_files/gpuMomentUpdater.cpp`.
**Task:** (a) Rotation: compute the angle `norm·Δt·γ/(1+α²)` and `sincos` in `accum_t`, cast the rotation matrix to storage — the angle is tiny and fp32 `sincos` of a tiny angle loses the `1−cos` digits (`u = 1−cos` should be computed as `2·sin²(θ/2)` regardless of mode — small formula improvement, do it for all precisions). (b) Spin-length drift: under SINGLE, renormalize `emom` every `K` steps (new input/def, default K=100 for SINGLE, off otherwise) — one cheap element kernel; document. (c) `dp_factor`/`sigma` prefactors: compute on host in double, cast once.
**Acceptance:** SINGLE 10⁶-step conservative test (λ=0, T=0): |m| deviation from 1 stays < 1e-6 with renormalization on; energy conservation comparable to Fortran fp64 within fp32 expectations (document measured drift in the PR).

### - [ ] PR6 — FFT precision  (P1, S; with CV3) — **partial: binary `SINGLE_PREC` only, not the specified three-mode design (see Phase 4 note)**
**Files:** `gpu_files/gpuFftWrapper.hpp`, `gpu_files/gpuLatticeConvolutionHamiltonian.cpp`.
**Task:** Select `CUFFT_C2C/R2C` vs `Z2Z/D2Z` (and hipFFT equivalents) from `storage_t`. The spectral multiply runs in storage precision; kernel construction (done once) always in fp64 on host/device then cast — the kernel entries are sums of couplings and deserve full precision.
**Acceptance:** Convolution V2 case passes in all three modes at mode-appropriate tolerance.

### - [ ] PR7 — Tolerance-tiered validation  (P1, S) — **partial: binary `SINGLE_PREC` only, not the specified three-mode design (see Phase 4 note)**
**Files:** `tests/gpu_regression/`.
**Task:** Parameterize V2 tolerances by precision: DOUBLE rtol 1e-10; MIXED rtol 1e-5 on trajectories / 1e-7 on averaged energies; SINGLE rtol 1e-4 / 1e-5, plus the statistical thermal tests from PR3/PR4. CI (compile-only) builds all three modes; the numeric suite runs on demand/GPU runner.
**Acceptance:** One command runs the full matrix {DOUBLE, MIXED, SINGLE} × {CUDA, HIP if available} and prints a table.

---

## Phase 5 — Convolution backend  ★ user emphasis 3

**Verdict: keep it, invest in it.** The FFT lattice convolution is the right long-term architecture for Bravais(-with-basis) systems: O(cells·log) per field evaluation independent of coordination number, no neighbor-list memory (the biggest single array for high-mnn systems — see M3), naturally batched over ensembles, and — as suspected — it is exactly the infrastructure needed for FFT dipole–dipole (the demagnetizing field is a convolution with the dipole tensor; this is how every micromagnetics code does it). The `NA×NA` block-kernel structure (`kernel_fft` sized `cells·9·NA²`) is already the correct generalization for multi-atom bases. Two things need fixing before it can be trusted broadly (CV1) and made cheap (CV3); then the dipole extension (CV6) is a natural, well-scoped project.

### - [x] CV1 — Boundary-condition guard (correctness)  (P0, S)  — **verified landed 2026-07-20**
**Files:** `gpu_files/gpuHamiltonianCalculations.cpp` (`canUseLatticeConvolution`), `gpu_files/gpuLatticeConvolutionHamiltonian.cpp`.
**Background:** The kernel-construction code respects BCs via `wrapped_coord_diff` (min-image only for periodic dims), but the *application* is a cyclic FFT convolution — intrinsically periodic in all three dims. With any BC = open (`'0'`), boundary sites receive spurious wrapped contributions from the far side. `canUseLatticeConvolution` never checks the BCs.
**Task:** Add `if (BC1!='P' || BC2!='P' || BC3!='P') return false;` (with a printed note that the convolution backend requires full PBC for now). Zero-padding support for open dims is CV6's problem (dipole needs it anyway); don't implement it here.
**Acceptance:** Open-BC input falls back to the sparse path with a message; a new V2 case (open BC film) passes because it takes the sparse path; full-PBC convolution cases unchanged.

### - [ ] CV2 — Energy from the convolved field  (P1, M)
**Files:** `gpu_files/gpuLatticeConvolutionHamiltonian.cpp`, `gpu_files/gpuHamiltonianCalculations.cpp` (`heisge` dispatch).
**Background:** `heisge(measure=true)` currently bypasses the convolution backend entirely, so measuring steps pay the sparse cost and the sparse arrays must stay resident (blocks M3).
**Task:** After `unpack_fft_field`, when `measure`, launch an energy kernel computing per-ensemble `E_exch = −½ Σ_i S_i·B^conv_i` (the convolution field is the bilinear field, so the ½ double-counting factor applies as in the sparse `HeisgeJij` path — verify against `energyM` column conventions), `E_ext = −Σ S·B_ext` (sign per Fortran), plus the on-site anisotropy energy evaluated directly (it is on-site — cheap; note the current convolution kernel folds uniaxial aniso into the q=0 kernel entry via `add_uniaxial_anisotropy_kernel`, which is only valid as a *field*; the energy prefactors differ — evaluate the aniso energy in this kernel from `kaniso/eaniso/taniso` directly and make sure the field-side folding and the energy don't double count; simplest correct arrangement: remove aniso from the convolution kernel and apply aniso field+energy in the unpack/energy kernel on-site). DM energy analogous via the antisymmetric part already in the kernel. Reuse the F1 reduction machinery for the sums.
**Acceptance:** V2 convolution+energy cases match the sparse-path energies to rtol 1e-10 (fp64) for J, J+DM, J+aniso.

### - [x] CV3 — R2C/C2R transforms  (P1, M)  — **verified landed 2026-07-20**
**Files:** `gpu_files/gpuLatticeConvolutionHamiltonian.cpp`, `gpu_files/gpuFftWrapper.hpp`.
**Task:** The spin field and the output field are real; switch spin/field transforms to R2C forward / C2R inverse (`CUFFT_D2Z/Z2D` resp. `R2C/C2R` per precision). Spectral buffers shrink to `(n1/2+1)·n2·n3` complex per component; the kernel FFT likewise (kernel is real in real space). The spectral multiply indexing changes to the half-spectrum layout. Mind hipFFT parity. Keeps ~45 % of the FFT memory and time.
**Acceptance:** Bitwise-tolerant (1e-12 fp64) agreement with the C2C result on V2 convolution cases; tracker shows the buffer reduction.

### - [x] CV4 — Stream + plan association audit  (P2, S)  — **verified landed 2026-07-20**
**Files:** same.
**Task:** ~~Ensure `GPUFFT_SET_STREAM(plan, workStream)` is called~~ — **verified done 2026-07-20**: `gpuLatticeConvolutionHamiltonian.cpp:409-417` sets the stream on `forward_plan`, `backward_plan` and `kernel_plan`. The convolution cases are also bitwise-identical under `--default-stream per-thread` (see PF3), which they would not be if any plan were still on the default stream. Remaining: destroy `kernel_plan` after build (also listed as M5-1 — do once).
**Acceptance:** `nsys`/`rocprof` trace shows FFTs on the work stream.

### - [ ] CV5 — Crossover benchmark + auto-selection groundwork  (P2, S)
**Files:** `tests/gpu_regression/bench_conv.py` (new), docs.
**Task:** Benchmark sparse vs convolution field evaluation across N (10³…10⁷) and mnn (6…50+) on one CUDA and one HIP machine; record the crossover. Print the measured per-step timings at startup when both backends are available (groundwork for later auto-selection; do **not** auto-select yet — keep the input flag authoritative).
**Acceptance:** A table in `docs/` with crossover guidance.

### - [ ] CV6 — FFT dipole–dipole backend  (P1, L — flagship follow-on)
**Files:** new `gpu_files/gpuDipoleConvolution.{hpp,cpp}` (or extend the existing class), `gpu_files/gpuHamiltonianBackend.hpp` (a `FftDipole` slot exists in spirit via `MultiscaleDipole` — add a proper enum value), `chelper.f90`/`fortranData.*` (pass `do_dip` mode + prefactors), dispatch in `heisge`.
**Task (staged):**
1. **Kernel construction:** dipole tensor `D_ab(r) = (3 r_a r_b − r² δ_ab)/r⁵` summed over the lattice into the same `cells·9·NA²` block-kernel storage used by exchange, with the correct `μ0/4π`-equivalent prefactor in UppASD's field units (reference: Fortran `dipole*` / `DipoleManager` modules — match `do_dip=1` brute-force results). Build once on device (a kernel over cells×NA², each thread summing its displacement — for PBC use Ewald-style or direct real-space sum with a documented cutoff/convergence parameter; simplest correct v1: real-space sum over `n_rep` periodic images with convergence check against the Fortran brute-force energy).
2. **Open boundaries:** zero-pad each open dimension to `2·n_i` (standard micromagnetics trick makes the cyclic convolution exact for open BC). This also unlocks lifting CV1's restriction for the exchange backend later. Memory cost up to 8× on the padded spectral buffers — surface in the M1 budget estimator.
3. **Application:** the dipole field is *additive* — accumulate its spectral multiply into the same `field_fft` as exchange when both are convolution-served (one inverse FFT total), or run standalone added onto `beff` when exchange is sparse.
4. **Validation:** V2 case vs Fortran `do_dip=1` (brute force) on a small system, field rtol 1e-8 (fp64) and energy agreement; a thin-film case checking the demag-driven easy-plane behavior qualitatively.
**Acceptance:** All four stages; documented input flag (`do_dip 4` or similar, coordinated with Fortran input parsing); budget estimator updated.

---

## Phase 6 — Performance (audit P1–P8, minus items absorbed above)

### - [ ] PF1 — MC: local ΔE in `MCSweep`, drop per-sublattice heisge  (P1, L) — audit P1, pairs with F5

> **Deferred / WIP:** Revisit together with F5/F6/PF2 after the remaining GPU-MC correctness work is scoped.
**Files:** `gpu_files/gpuMetropolis.cpp`, `gpu_files/gpuMCSimulation.cpp`.
**Task:** Rewrite `MCSweep` to gather the trial spin's neighbors from `nlist/ncoup` (+ `dmlist/dmvect`) directly and compute ΔE = ΔE_exch + ΔE_DM + ΔE_Zeeman + ΔE_aniso(exact, from F5) in-kernel. Remove the `hamCalc.heisge` call from the per-sublattice loop; call `heisge` once per MC step only when a measurement needs `beff`/energies that step. The neighbor gather uses the transposed layout (`pos[site + i·N]`, `coup[rsite + i·NH]`) for coalescing, same as the field kernels. Keep the old path behind the existing `bruteforce` toggle until validated, then delete.
**Acceptance:** V2 MC cases statistically match Fortran across the T sweep; wall-clock per MC step improves by ≈ `num_subL`× on the benchmark system (report numbers).

### - [ ] PF2 — O(N·mnn) coloring, silent  (P1, S) — audit P5, with F6

> **Deferred / WIP:** Revisit together with F5/F6/PF1.
**Files:** `gpu_files/gpuMetropolis.cpp` (`split_lattice` + helpers `fillneighbours_in_play` etc. — delete them).
**Task:** Replace with greedy graph coloring: iterate sites, mark colors of already-colored neighbors (union list per F6) in a small bitset/array, assign the lowest free color; then bucket site indices by color into `subIdx`. Remove all `printf`s of the table and the O(N²) scans.
**Acceptance:** Same or fewer colors than the old code on V2 systems (assert independence property in a debug check); startup time for a 10⁶-site system < 1 s.

### - [x] PF3 — Stream discipline sweep  (P1, M) — audit P2 — **done 2026-07-20**
**Files:** `gpu_files/measurement/fortranMeasurement.cpp`, `gpu_files/measurement/cpuRestMeasurement.cpp`, `gpu_files/correlations/fortranCorrelation.cpp`, `gpu_files/correlations/gpuCorrelations_cucomplex.cu`, `CMakeLists.txt`.

Most of the original task was already complete: all ~30 launches in `gpuCorrelations.cpp` and all of `measurement/*` already pass `workStream`; the remaining `GPU_STREAM_SYNC(workStream)` calls are correctly scoped to host-consumption boundaries; `StopwatchDeviceSync::sync()` is inert (`ASYNC_STOPWATCH`); `gpuDepondtIntegrator.cpp` `GPU_DEVICE_SYNCHRONIZE` is an init-time allocation check.

**Verification method:** built with `-DCMAKE_CUDA_FLAGS="--default-stream per-thread"` (all GPU `.cpp` compile as CUDA, so the flag covers them). This removes the implicit edge between the legacy default stream and the blocking `workStream`, exposing any latent stream-0 dependency.

**Latent bug found and fixed:** the `copyQueueSlow` D2H copies in `fortranMeasurement`, `cpuRestMeasurement` and `fortranCorrelation` ran on the default stream and depended entirely on legacy default-stream semantics to wait for `workStream`. Under PTDS, 7 of 12 regression cases produced wrong energies (e.g. `feco_unreduced` rel 7.9e-01) from moments the integrator had not finished writing. The `*_gpu_measurements` cases passed throughout because `saveToFile` syncs `workStream` explicitly. Note: the old claim in this item that "copy_sync already synchronizes" was **false** — `copy_sync` is an ordinary `cudaMemcpy` and carries no ordering of its own. Fixed by an explicit `GPU_STREAM_SYNC(parallel.getWorkStream())` before each slow-path copy block.

**`USE_CUCOMPLEX_CORRELATIONS` resolved as unsupported:** `gpuCorrelations_cucomplex.cu` launches everything on the default stream, and the source-selection block in `correlations/CMakeLists.txt` that once chose it is commented out — the file is attached to no target and is never compiled, so the option silently built the Thrust path instead. It is now a configure-time `FATAL_ERROR`, with revival requirements documented in the file header.

**Acceptance:** all 12 cases bitwise-identical to CPU across `build_ab`, `build_fastcopy` and the PTDS build (only the pre-existing `sc_dm_uniaxial` crash in `gpusim_initiatematrices`, `sd_driver.f90:1185`, still fails, identically in all three); `compute-sanitizer` racecheck/memcheck/synccheck/initcheck clean on both `build_ab` and `build_fastcopy`.

**Out of scope, still open:** `gpuMetropolis.cpp` / `gpuMetropolis_bruteforce.cpp` launch on the default stream and `gpuMCSimulation` has four `GPU_DEVICE_SYNCHRONIZE` calls (C5/C6/P1/P5). The MC regression cases pass under PTDS only because those device syncs are heavy enough to mask the missing ordering — this is worth revisiting with MC work, not before. `nsys` step-time profiling of a correlation-heavy case was not done (no correlation case exists in `tests/gpu_regression`).

### - [x] PF4 — `j_tensor` device layout transpose  (P2, M) — audit P4  — **verified landed 2026-07-20**
**Files:** `gpu_files/gpuSimulation.cpp` (upload via F7/PR2 staging), `gpu_files/gpuHamiltonianCalculations.cpp` (`HeisgeJijTensor*` indexing, `SetupNeighbourListExchangeTensor`), convolution `buildTensorKernel` reader.
**Task:** Store the device tensor as `[site + NH·(component + 9·neigh)]` (site fastest) and update all readers. Upload-time transform via the staging layer.
**Acceptance:** V2 tensor cases pass; tensor-mode field kernel bandwidth (nsys/rocprof) improves materially (report before/after).

### - [x] PF5 — Integrator kernel fusion  (P2, M) — audit P6, enables M2-3  — **verified landed 2026-07-20**
**Files:** `gpu_files/gpuDepondtIntegrator.cpp`, `gpu_files/gpuCommon.hpp`.
**Task:** Fuse per half-step: one kernel doing local field = beff+btherm (+btorque), damping cross-product, normalization, rotation-matrix build, rotation, and (first half) `emomM = mrod·mmom` — all per-atom register work. `blocal` disappears as an array; `bdup` is still needed across the two half-steps (b2eff avg). Keep the old kernels behind a define for A/B validation for one release.
**Acceptance:** V2 bitwise (fp64, same operation order where feasible — otherwise 1e-12 rtol); integrator section of the stopwatch report ~2× faster; `blocal` gone from the tracker.

### - [x] PF6 — Misc: `__restrict__` uniformity, N%32 note, unbuffered stdout  (P3, S)  — **verified landed 2026-07-20**
**Files:** Heisge functors (after U2 this is one place), `gpuHamiltonianCalculations.cpp`, all `std::setbuf(stdout, nullptr)` sites.
**Task:** Uniform `const __restrict__` pointers in functors; fix the "multiple of 32" message to be wave-size aware; remove `setbuf(…, nullptr)` from production phases (keep behind a debug define).

---

## Phase 7 — Unification & structural refactors (audit A2, A3)

### - [x] U1 — Single-source CUDA/HIP: thermfield + correlation kernels  (P1, L) — audit A2  — **verified landed 2026-07-20**
**Files:** `cudaThermfield.cu` + `hipThermfield.cpp` → `gpuThermfield.cpp`; `correlations/correlation_kernels.cpp` + `.hip.cpp` → one file; `gpu_wrappers.h`, `real_type.h`, both `CMakeLists.txt` (with CM4).
**Task:** Add to `gpu_wrappers.h`: `WAVE_SIZE` compile-time constant (`__AMDGCN_WAVEFRONT_SIZE` on AMD HIP, else 32), RNG generator-API macros (`GPU_RAND_CREATE_GENERATOR`, `GPU_RAND_GENERATE_NORMAL[_REAL]`, `GPU_RAND_SET_STREAM`, …). Merge the file pairs: mechanical diff resolution using the existing `GPU_CREAL/GPU_CIMAG/MAKE_GPU_COMPLEX` macros, `cg::tiled_partition<WAVE_SIZE>`, and the raw-`blockIdx` idioms (pick the raw form — it's what the HIP file uses and it's fine). Pick **one** complex representation for HIP (recommend `thrust::complex` everywhere, i.e., generalize the `ON_LUMI` path and delete the flag — rocThrust ships with ROCm) unless benchmarking shows a reason not to. Delete the twin files.
**Acceptance:** Both backends build from the unified sources; V2 passes; `diff`-able twins no longer exist; `ON_LUMI` removed from CMake and code.

### - [x] U2 — Templated Heisge kernel  (P1, L) — audit A3 + P3  — **verified landed 2026-07-20**
**Files:** `gpu_files/gpuHamiltonianCalculations.cpp/.hpp`.
**Task:** Replace the six `HeisgeJij*` functors with one `template<bool HasDM, bool HasAniso, bool HasTensor, bool Measure> class Heisge : ParallelizationHelper::AtomSiteEnsemble` using `if constexpr`; a 16-way (practically 12-way — tensor excludes scalar-J DM path per current semantics, mirror existing dispatch) dispatch table in `heisge()` selects the instantiation. Fold in: F1's fixed energy path (only compiled when `Measure`), F4's ensemble usage, PF6's `__restrict__`, and resolve the "To Anders: sign????" prefactor questions against the Fortran reference (`effective_field`/`calculate_energy` in `hamiltonianactions*.f90`) — document the resolved conventions in a header comment block and **get sign-off from Anders on the uniaxial `*_en` sign and the cubic ½ and ⅓ energy factors before merging**.
**Acceptance:** V2 passes for every Hamiltonian permutation including energies; non-measure step time unchanged or better; measure-step register count reported (nvcc `-Xptxas -v`) lower than baseline for the plain-J non-measure instantiation; total file shrinks by ~800 lines.

### - [x] U3 — `bind(C)` the whole bridge  (P2, M) — audit I1  — **verified landed 2026-07-20**
**Files:** `source/chelper.f90`, `gpu_files/fortranData.cpp`, `gpu_files/fort_helper.cpp`, `gpu_files/c_helper.h`.
**Task:** Give every `FortranData_set*`, `gpusim_*`, and the C→Fortran callbacks (`fortran_measure_*`, `fortran_calc_simulation_status_variables`, …) explicit `bind(C)` interfaces with `iso_c_binding` types; C symbols lose the trailing underscore (`fortrandata_setconstants` etc. — update both sides); character flags become `character(kind=c_char)` scalars or `integer(c_int)` enums (prefer enums for new code); settle `integer(c_int)` vs the current `unsigned int*` mismatches on `int`. The `fortran_print_measurables` path passing `allocatable` dummies must be audited: if it is ever invoked with a C-constructed actual argument, convert to `bind(C)` + explicit shape/`CFI_cdesc_t`, otherwise document that it is Fortran-called only.
**Acceptance:** Builds with gfortran and (if available) another compiler (nvfortran or Cray ftn via CI/LUMI) without symbol errors; V2 passes; no `extern "C" void *_(...)` hand-mangled symbols remain.

### - [ ] U4 — Measurement thread hardening or retirement  (P2, M) — audit A4
**Files:** `gpu_files/measurement/measurementQueue.cpp`, `fortranMeasurement.cpp`.
**Task:** Minimum hardening: the empty-measurement path must snapshot too or the worker must never receive live `FortranData::` pointers (pass `nullptr`-means-skip instead of fallback-to-live); replace pthreads with `std::thread/std::mutex/std::condition_variable`; document (in code) which Fortran routines the worker may call and assert the main thread doesn't call them concurrently (a simple atomic "fortran_busy" flag with an abort in debug builds). If by this point GPU measurement coverage is complete for all observables used in practice, propose retirement of the thread in a short design note instead.

---

## Phase 8 — Hygiene sweep  (P3, S each — batch into one or two PRs)

- [ ] H1 — Delete `gpu_files/old/`; delete dead code blocks (`null_energy`, TEMPORARY PRINTING sections, disabled `cpuRestMeasurement` wiring, commented twin of `fortran_print_measurables` in `chelper.f90`, dead locals `mnn/l/NH` in both SD phases, the stray `real cv{};` inside the SDmphase loop).
- [ ] H2 — Remove `#pragma once` from `.cpp` files (`fortranData.cpp`, `gpuMetropolis.cpp`, `gpuCorrelations.cpp`, `gpuParallelizationHelper.cpp`).
- [ ] H3 — Normalize line endings to LF (`.gitattributes` with `* text=auto`, one renormalize commit).
- [ ] H4 — `enum class EnergyColumn { Exchange=0, Aniso=1, DM=2, Tensor=3, External=4, Total=5 }` replacing the magic indices into `energyM`.
- [ ] H5 — Rename `redNeibourCount` → `redNeighbourList` (it holds `aHam`, a site→reduced-site map, not a count — fix the semantics in the name), fix "Depontd" comments.
- [ ] H6 — `GridHelper::dim*` error paths: return false and propagate to a caller-side check instead of `exit()`; resolve the "AS:" comments.
- [ ] H7 — Add `clang-format` config matching the dominant existing style and format `gpu_files/` once (after U1/U2 to avoid churn).

---

## Suggested execution order

```
V1 + F3  →  V2, V3  →  F1, F2, F4  →  M1  →  CM1–CM5 (unblocks everyone's builds)
→ F7+PR2, PR1  →  F5+PF1, PF2, F6  →  U1  →  U2 (+F1 energy integration, PF6)
→ PR3–PR7  →  CV1, CV2, CV3, M3  →  PF3, PF5, M2, M4  →  CV6 (dipole)
→ U3, U4, CM6, M5, M6, PF4, CV4, CV5  →  H1–H7
```

Rationale: safety net first; the three P0 correctness bugs and the OOM report make current results trustworthy and crashes diagnosable; CMake early because every subsequent task builds on it; the staging layer (F7/PR2) is shared infrastructure for precision, memory, and layout work; unification (U1/U2) before the wide precision sweep so type changes happen in one copy of each kernel; convolution improvements before the dipole flagship; hygiene last to avoid rebase pain.
