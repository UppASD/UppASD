# RCG-09B reduction evidence

**Task:** RCG-09B, replace serialized adaptive FP64 energy reductions.
**Implementation date:** 2026-08-14.

## Current inventory

The live CUDA/HIP source was checked after the RCG-09A refactor. The affected
energy terms and kernels are:

| term | live producer | contribution granularity |
| ---: | --- | --- |
| 0 | `evaluateAdaptiveAtomisticBonds` | one bond/ensemble thread |
| 1 | `evaluateAdaptiveAtomisticOnsite` | one active-atom/ensemble thread |
| 2, 3 | `evaluateAdaptiveCoarseTensor` | one active-block/ensemble thread, with the existing tensor loops |
| 4, 5 | `finalizeAdaptiveCoarseLocal` | one active-block/ensemble thread |
| 6 | `addAdaptiveDipole` or `addAdaptiveBasisResolvedDipole` | one block/ensemble thread |

There are no adaptive production calls to `atomicAddEnergyTerm` after this
change. `source/gpu_files/gpuAtomicDouble.hpp` remains for the historical
RCG-06B microbenchmark and evidence-only controls; it is not used by the
adaptive production kernels.

## Reduction and numerical contract

Each energy-producing launch writes one `double` partial per CUDA/HIP block to
a term-major compact array. Threads contribute to a shared-memory FP64 binary
tree. Out-of-range threads contribute zero, so every launched block writes a
defined partial. A final one-block kernel sums each term's partials in
ascending launch-block order and writes the existing eight FP64 `energyTerms`
slots. The existing fixed term order then computes slot 7.

This preserves FP64 contribution evaluation/storage and the unconditional
energy availability contract. The summation order changes from nondeterministic
atomic-arrival order to a defined block-local tree followed by ascending block
order. Results are therefore reproducible for a fixed launch geometry and
backend/compiler; they are not asserted bitwise identical to the former CAS
arrival order. The accepted energy comparison must allow only the resulting
floating-point reduction-order difference. Field updates remain unchanged
apart from the removal of energy-side CAS traffic.

The scratch allocation is included in `estimateBytes()` and is allocated with
the other CG-10 kernel storage. CUDA and HIP use the same launch geometry and
shared-memory reduction through the shared launch macro.

## Evidence status

Completed locally:

- live inventory verified at eight affected producer groups;
- all eight call-site groups converted, including both atomistic and coarse
  paths and both dipole variants;
- no adaptive production energy CAS call remains;
- CUDA Release/fp64 build completed in `build_rcg09a_cuda_fp64`;
- CUDA Release/fp32-storage/fp64-energy compile completed in
  `/tmp/rcg09b_cuda_single`.

The implementation host has no usable device for independent local execution,
but the CUDA handoff run supplied standalone correctness and diagnostic timing
evidence. The benchmark's `--parity-only` target passed the RCG-09A Hamiltonian
oracle and adaptive-ASD cases, including exchange, DMI, anisotropy, combined,
and finite-temperature cases. Field differences were at roundoff level and
named energy differences were approximately 0--3e-16 against a 1e-12
tolerance. This does not make the full production e2e suite pass.

The original selected ctest command reported one failure,
`adaptive-cg-production-e2e`, before reaching its GPU cases. Its CPU chiral
fixture `dmi_anisotropy_mixed` used `ncell 10 2 2` and listed both `+y` and
`-y` DMI shells. With two periodic cells in `y`, those shells aliased to
canonical pair 1/21; the two loaded vectors were the same rather than
reciprocal. The RCG-09A folded-pair validator correctly rejected this
ambiguous representation. Both CPU and GPU fixture inputs were corrected to
`ncell 10 4 2`: this preserves the block-size divisibility contract and
removes the periodic-image alias. A fresh local CPU production-e2e rerun then
passed all cases. The subsequent CUDA device-side run passed the corrected
GPU ctest: the production e2e, GPU adaptive-runtime, and adaptive-ASD parity
targets all passed (3/3).

The first corrected CUDA ctest then reached the GPU parity cases but exposed a
separate diagnostic-buffer defect: the static CPU/GPU `coarse_sum` comparison
was `535.7009401755746` versus `0.0`. `evaluateHybrid()` writes the assembled
coarse field to `coarseField_`; `diagnosticSnapshot()` had incorrectly read the
zeroed `predictorCoarseField_` integration scratch buffer. The diagnostic now
reads `coarseField_`. This changes only the reported checksum, not field
assembly, integration, or energy reduction. The CUDA build was rerun after
the fix, followed by a passing device-side ctest run.

The next CUDA rerun passed the coarse-field checksum assertion but failed the
following static CPU/GPU restart-state trajectory assertion. The maximum
observed component difference was `6.291162e-3` in the coarse reconstructed
block, against the existing `8.0e-5` budget; adaptive-ASD parity and the GPU
runtime test still passed. The tolerance was not widened. Source inspection
identified the remaining state-handoff defect: the active-subset Depondt
corrector stores accepted fine-atom vectors in `emom2`, while the adaptive
driver reconstructed/published `emom` and then used that stale pre-corrector
buffer for the next step. The driver now reconstructs into `emom2` and swaps
the accepted buffer into `emom` before selector work or the next field
evaluation. The subsequent CUDA production-e2e run passed the existing energy,
field, restart-state, and all-fine parity checks; the `8.0e-5` restart budget
was retained.

The corrected production e2e harness and direct scaling benchmark were then
run on a device accepted by the clean gate: an NVIDIA RTX A4000, driver
610.57.04, 0% utilization, 47 C, SM clock 210/2100 MHz, and raw clock-event
mask `0x0000000000000001`. The idle-only bit `0x1` is NVIDIA's ordinary GPU
idle event; the harness retains the raw mask and rejects actionable power,
thermal, application-clock, and hardware-slowdown bits. The full direct
scaling artifact, including raw samples, medians/MADs, active fractions,
phase timings, and metadata, is
`docs/rcg09/rcg09b_cuda_clean_scaling.txt`.

The production comparison measures only the all-fine fixture
(`active_atoms=16384`, `active_blocks=0`, `interface_atoms=0`,
`adaptive_device_bytes=41861688`):

| texture | atomistic median/MAD (us) | adaptive median/MAD (us) | ratio atomistic/adaptive | crossover |
| --- | ---: | ---: | ---: | --- |
| spiral | 142.733 / 10.987 | 1337.040 / 151.960 | 0.106753 | not observed |
| uniform | 176.164 / 4.758 | 1273.256 / 220.171 | 0.138357 | not observed |

Raw production samples and full metadata are in the operator artifact
`/tmp/rcg09b_cuda_clean.json`; the repository's direct-scaling artifact retains
the corresponding sweep samples. These values show that this all-fine
adaptive production path is approximately 7.23--9.37x slower than the ordinary
atomistic path on the measured RTX A4000. The direct sweep observed no
production crossover at any size or active fraction. Its five-point log-log
OLS exponents are 0.320712 for the production baseline, 0.185922 for
all-fine adaptive, and 0.200597 for zero-fine adaptive; endpoint exponents are
0.333815, 0.202523, and 0.234338 respectively. Coarse and interface phases
remain dominant as the active fraction decreases, rather than disappearing.

Standalone CUDA parity, corrected CPU production-e2e correctness, the
device-side CUDA production-e2e correctness gate, and clean-device CUDA
scaling/fraction measurements are demonstrated. Production crossover and HIP
execution remain open; no speedup claim is restored.

## Negative/control evidence

The accepted pathological comparison remains in
`docs/RCG-09_PRODUCTION_PERFORMANCE_EVIDENCE.md` §3.3 and the raw samples in
`docs/rcg09/scaling_with_energy_atomic.txt` and
`docs/rcg09/scaling_without_energy_atomic.txt`. It restored only the old bond
single-address accumulation behavior, showed the 242x/1,762x/10,489x ratios,
and was reverted before production code. No pathological implementation is
selectable in this tree. This historical control is evidence for causation,
not a correctness comparison for the new reduction.

## GPU handoff commands

On a clean, non-throttled CUDA device, rebuild and run the unchanged RCG-09A
correctness/parity targets:

```bash
cmake --build build_rcg09a_cuda_fp64 -j2
ctest --test-dir build_rcg09a_cuda_fp64 --output-on-failure \
  -R 'coarse-graining-(gpu-dmi-dimer|gpu-adaptive-runtime)|adaptive-cg-(production-e2e|moving-backend-parity|timing-reconciliation|ownership-map-comparator|dipole-ownership-check|adaptive-asd-parity)'
```

Capture device state immediately before the measurement, then run the
production sweep. The harness refuses dirty devices unless explicitly told
otherwise; `--allow-dirty-device` is not acceptance evidence.

```bash
nvidia-smi
nvidia-smi --query-gpu=index,name,driver_version,utilization.gpu,memory.used,memory.total,temperature.gpu,clocks.current.sm,clocks.max.sm,clocks_event_reasons.active,persistence_mode,compute_mode --format=csv
python3 tests/coarse_graining/run_rcg09_perf_e2e.py \
  build_rcg09a_cuda_fp64/bin/sd.f95.cuda \
  --json /tmp/rcg09b_perf.json
```

For the accepted direct scaling points and raw output, run:

```bash
for blocks in 64 256 1024 4096 16384; do
  build_rcg09a_cuda_fp64/bin/gpu_adaptive_runtime_benchmark \
    --blocks "$blocks" --atoms-per-block 4 --warmup 1 \
    --iterations 3 --repetitions 3 --texture spiral \
    | tee -a /tmp/rcg09b_scaling.txt
done
```

The clean CUDA artifact now provides the requested raw samples, median/MAD,
phase times, active fractions, and device clock/utilization/thermal metadata.
Report the measured values only: every adaptive point was slower than the
production baseline, and no crossover was observed. The previously observed
roughly 6x atomistic architectural opportunity remains the reporting ceiling
unless a separate task removes additional O(N) work and measures it.
