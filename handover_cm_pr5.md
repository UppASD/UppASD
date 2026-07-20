# Handover — CM1–CM6 (CMake modernization) and PR5

Written 2026-07-20 at commit `28f4b5a`, branch `gpu_hip_cu_ab`.

Status was surveyed with real evidence, not inferred from commit titles. **CM1–CM5 are
confirmed not done.** CM6 is partially done. PR5 is not done and is *blocked*. Details and
line references below — re-locate by symbol, line numbers drift.

---

## Repo conventions you must follow

### Commit messages: subject line only

**All 23 commits preceding this handover have an empty body.** The convention is a single
subject line, `area: imperative summary`, roughly <= 72 chars:

```
gpu: size Binder cumulant partials by full launch grid
cmake: make MKL selection compiler-safe
tests: add compute-sanitizer runner for GPU measurement paths
```

Prefixes in use: `gpu:`, `cmake:`, `tests:`/`test:`, `docs:`, `fortran:`.

Two commits from the preceding session (`0c4179d`, `28f4b5a`) shipped 15–20 line prose
bodies. That was a departure from the convention, not a new standard — **do not copy those
two.** Put the reasoning in the checklist entry or in this kind of handover file, not in the
commit body.

The harness requires a `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>` trailer, which
the existing history does not have. Keep the trailer (it is required) but keep the prose body
empty — trailer only.

### Other conventions

- Mixed CRLF/LF exists in the tree (`gpuSDSimulation.cpp` notably). Editing tools may
  normalize whole regions and produce a huge spurious diff. Run `git diff --stat` after every
  edit; keep added lines LF.
- `tests/*/inpsd.dat` and `tests/gpu_regression/cases.json` carry **pre-existing uncommitted
  edits that are not yours**. Never stage them. Stage only files you intend to change.
- `git diff --check` reports a pre-existing trailing-whitespace warning in
  `tests/bccFe_cuda/inpsd.dat`. Scope the check: `git diff --check -- source CMakeLists.txt`.

---

## Validation commands

```
cmake --build build_ab        --parallel 8      # default CUDA build     -> bin/sd.cuda
cmake --build build_fastcopy  --parallel 8      # USE_FAST_COPY          -> bin/sd.f95.cuda
cmake --build build_ptds      --parallel 8      # per-thread default str -> bin/sd.f95.cuda

python3 tests/gpu_regression/run.py --binary build_ab/bin/sd.cuda
python3 tests/gpu_regression/run.py --binary build_fastcopy/bin/sd.f95.cuda
python3 tests/gpu_regression/run.py --binary build_ptds/bin/sd.f95.cuda
python3 tests/gpu_regression/sanitize.py                                  # fast-copy, 4 tools
python3 tests/gpu_regression/sanitize.py --binary build_ab/bin/sd.cuda
```

Expected: **11 PASS, 1 FAIL** in every configuration. The single failure is `sc_dm_uniaxial`,
which segfaults in `gpusim_initiatematrices` (`sd_driver.f90:1185`) in *all* configurations and
predates this work. Confirm it fails identically; do not treat it as a regression.

`USE_FAST_COPY` is **not** a CMake option — `build_fastcopy` was hand-configured:

```
cmake -S . -B build_fastcopy -DCMAKE_BUILD_TYPE=DEBUG -DUSE_CUDA=ON -DUSE_MKL=ON \
  -DCMAKE_CUDA_FLAGS=-DUSE_FAST_COPY -DCMAKE_CXX_FLAGS=-DUSE_FAST_COPY
```

`build_ptds` was configured with `-DCMAKE_CUDA_FLAGS="--default-stream per-thread"`. It is the
regression net for stream-ordering bugs: it removes the implicit edge between the legacy
default stream and the blocking `workStream`. **CM work that touches CUDA/HIP flags must keep
this build working**, or the net is lost. All GPU `.cpp` compile as CUDA
(`CMakeLists.txt:348-356`), which is why a `CMAKE_CUDA_FLAGS` change reaches them.

---

## Part 1 — CM1–CM6

Phase 3 goal (from the checklist): `cmake --preset <name> && cmake --build` works out of the
box on a consumer NVIDIA laptop, an RDNA/CDNA AMD laptop, LUMI/Dardel-class HPC with Cray
PrgEnv + ROCm, and NVIDIA HPC with modules — no hand-editing.

Current reality: `CMakeLists.txt.local` exists in the tree and `build_ab`, `build_fastcopy`,
`build_deb`, `build_gpu_deb`, `build_fastcopy`, `build_ptds` are all hand-configured. That is
the problem Phase 3 exists to fix, and it is a good measure of done: when CM5 lands, those
hand-rolled invocations should collapse into presets.

### Evidence gathered (all confirmed present, 2026-07-20)

| Item | Finding | Location |
|---|---|---|
| CM1 | No `UPPASD_GPU_BACKEND`, no `UPPASD_PRECISION` | absent |
| CM1 | 8 global `add_compile_definitions` | `CMakeLists.txt` |
| CM1 | `-DDEBUG` injected into CUDA, HIP and CXX flags | `:368`, `:429`, `:430` |
| CM2 | Legacy `find_package(CUDA)` still present | `:352` |
| CM2 | `FindCUDA/select_compute_arch` still used | `:384` |
| CM2 | Blanket `CUDA::toolkit` link, not `cudart/curand/cufft` | `:398` |
| CM2 | Manual `-I` flag splicing | `:371`, `:400` |
| CM3 | Sets `CMAKE_HIP_COMPILER` / `CMAKE_CXX_COMPILER` (CM3 forbids) | `:420-421` |
| CM3 | Hardcoded `/opt/rocm` search paths | `:449`, `gpu_files/CMakeLists.txt:66` |
| CM3 | `find_library(hiprand)` rather than `find_package(... CONFIG)` | `:446-449` |
| CM4 | Source list duplicated: subdir `target_sources` + top-level blob | `:348-356` |
| CM4 | Commented-out disagreeing entries | `correlations/CMakeLists.txt:1-11` |
| CM5 | Presets are `debug/release/debugCuda/releaseCuda` | `CMakePresets.json` |
| CM5 | `docs/BUILDING_GPU.md` missing | absent |
| CM6 | `USE_MKL` option exists; some hygiene landed | `:195-196` |

**One thing CM2 asks for that I could not confirm:** the checklist cites a broken
`MATCHES "^[__86+PTX]"` regex (audit B2). I did not find that string — the arch block at
`:383-395` uses `CUDA_DETECT_INSTALLED_GPUS` with a hardcoded fallback of `80`. Either the
regex was already removed or the checklist's description is stale. Verify before writing a
changelog entry claiming you fixed it.

### Ordering

`CM1 -> CM2 -> CM3 -> CM4 -> CM5`, with CM6 independent.

CM1 first: it defines `UPPASD_GPU_BACKEND` and `UPPASD_PRECISION`, which CM2/CM3 branch on and
CM5's presets set. **CM1's `UPPASD_PRECISION` is also the unblocker for all of Phase 4** — see
Part 2. CM4 wants U1 landed for a single source list; U1 *is* landed (correlation and
thermfield kernels are unified), so CM4 is unblocked.

### Risks specific to this work

1. **Do not break the three existing builds.** `build_ab`, `build_fastcopy`, `build_ptds` are
   the only regression net. If CM1 changes how GPU defines are injected, `USE_FAST_COPY`
   (passed via `CMAKE_CXX_FLAGS`/`CMAKE_CUDA_FLAGS`) and `--default-stream per-thread` must
   both still work. Consider promoting `USE_FAST_COPY` to a real CMake option as part of CM1 —
   it rotted undetected precisely because it was not one. That is arguably in CM1's scope
   ("option surface hygiene") and would prevent a recurrence.
2. **Keeping `USE_CUDA`/`USE_HIP` as deprecated aliases is load-bearing**, not politeness.
   `CMakeLists.txt.local`, existing build dirs, and CI all pass the old booleans.
3. **No AMD hardware is available here.** CM3 cannot be validated beyond "configures in a ROCm
   container". Do not claim CM3 acceptance without a real HIP build; say what you tested.
4. CM2's `CMAKE_CUDA_ARCHITECTURES=native` default will change what CI builds. V1 (compile-only
   CI) is ticked and must keep passing.

### Suggested commits

```
cmake: introduce UPPASD_GPU_BACKEND and UPPASD_PRECISION
cmake: modernize CUDA discovery and linking
cmake: remove hardcoded ROCm paths and compiler assignment
cmake: single source of truth for GPU sources
cmake: add per-platform configure presets
cmake: respect BLA_VENDOR for BLAS and LAPACK
```

---

## Part 2 — PR5, and the Phase 4 correction that precedes it

### Read this first: PR1–PR7 were over-ticked and have been reverted

An earlier pass in this session ticked PR1–PR4, PR6 and PR7 on the strength of `SINGLE_PREC`
code existing. **That was wrong and has been corrected** — those items are `[ ]` again, marked
partial.

The whole of Phase 4 is specified against a three-mode `UPPASD_PRECISION` (DOUBLE/SINGLE/MIXED)
with `storage_t`, `accum_t`, `accum_cplx` typedefs. **None of those identifiers exist anywhere
in `source/`.** What exists is a binary `SINGLE_PREC` macro reaching five files: `real_type.h`,
`tensor.hpp`, `gpu_wrappers.h`, `gpuFftWrapper.hpp`, `gpuSimulation.cpp`. Genuine work
(`122f73d`, `18c3cf5`, `d8e1361`), but not the specified design. MIXED — the mode PR3 exists to
serve — is absent. PR7's tolerances are a single 1e-8 tier, not per-precision.

Treat Phase 4 as **prototyped, not designed in**.

### PR5 is blocked

PR5's task is written in terms of `accum_t`:

- (a) Rotation: compute the angle `norm·Δt·γ/(1+α²)` and `sincos` **in `accum_t`**, cast the
  rotation matrix to storage.
- (b) Spin-length drift: under SINGLE, renormalize `emom` every K steps (new input, default
  K=100 for SINGLE, off otherwise).
- (c) `dp_factor`/`sigma` prefactors: compute on host in double, cast once.

`accum_t` does not exist, so (a) cannot be written as specified. **PR5 depends on PR1, which
depends on CM1's `UPPASD_PRECISION`.** The real chain is:

```
CM1 -> PR1 -> PR2/PR3/PR4/PR6 -> PR5 -> PR7
```

This is why doing CM1 first matters beyond CMake tidiness.

### One piece of PR5 is not blocked — do it now

PR5(a) contains a precision-independent improvement, called out in the checklist itself:

> `u = 1−cos` should be computed as `2·sin²(θ/2)` regardless of mode — small formula
> improvement, do it for all precisions

`1−cos θ` at small θ loses digits to cancellation in **any** precision; `2 sin²(θ/2)` does not.
The rotation is in `gpuDepondtIntegrator.cpp`. This is a self-contained, independently testable
change that does not wait on PR1.

Expect it to *change results at the last bits* — it is a genuine accuracy improvement, so
bitwise-identical output is the wrong acceptance bar. Check the regression harness tolerances
(`rtol`/`atol` 1e-8 in `cases.json`) still hold, and state the observed delta. Do not silently
loosen a tolerance to make it pass; if it does not fit, say so and bring it back for a
decision.

Note `kagome_tensor` already reports `abs=1.00e-11, rel=3.62e-09` rather than exact zero, so
that case is the sensitive one to watch.

Suggested: `gpu: compute versine stably in Depondt rotation`

### PR5(b) and (c)

Both need PR1. (c) is nearly precision-independent (compute prefactors on host in double, cast
once) and could be done early, but its benefit only materializes under SINGLE, so there is
little reason to split it out.

---

## Definition of done for this handover

- CM1–CM5 land, and the hand-configured build dirs collapse into presets.
- `build_ptds` (or its preset successor) still passes 11/12 — the stream net survives.
- `docs/BUILDING_GPU.md` exists.
- PR5(a)'s versine fix lands independently, with measured deltas reported.
- PR5(b)/(c) and PR1–PR4/PR6/PR7 are explicitly deferred behind CM1, not silently ticked.
- Checklist updated. If an item is found already done, tick it **with the evidence**, not with
  a general impression — that is the failure mode that produced the correction above.
