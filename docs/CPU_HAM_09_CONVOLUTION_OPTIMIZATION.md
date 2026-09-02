# CPU-HAM-09 — CPU convolution hot-path optimization

**Date:** 2026-09-02
**Status:** complete
**Baseline source:** pre-HAM-09 `290519b`
**Provider:** FFTW 3.3.10, optional `libfftw3_threads`
**Compiler/build:** GNU Fortran 13.3, Release, OpenMP, double precision

## Scope and method

The production source remains `source/Hamiltonian/cpuconvolution.f90`. The
implementation preserves the periodic reduced-stencil mapping, FFT
normalization, HAM-06 DMI convention, backend eligibility, and field-derived
energy path.

The new `cpu_convolution_benchmark` target runs five repetitions after two
warm-up applies and reports raw wall and stage samples. The stage counters now
use `omp_get_wtime()` in OpenMP builds; serial builds retain the `cpu_time()`
fallback. The parent-revision comparison used the same benchmark driver and
linked the pre-HAM-09 `CPUConvolution` object against the current dependency
objects, so the hot-path comparison does not include unrelated source drift.

The host reports eight physical cores, one socket, two hardware threads per
core, and one NUMA node; the active container exposes two CPUs via `nproc`.
The headline timing table therefore uses `OMP_NUM_THREADS=2`, with
`UPPASD_FFT_THREADS=1`. Eight-thread results were also sampled for continuity
with the earlier campaign but are not treated as an uncontended physical-core
claim.

## Baseline and retained changes

| Item | Result | Evidence |
|---|---|---|
| Baseline freeze | Retained as parent-revision source object | `290519b`, same benchmark driver, five samples |
| Work-buffer clears | Retained optimization | `real_work` pack and `field_spectral` spectral output are complete-overwrite paths; NaN-poison tests pass |
| Packing | Retained parallel loop | OpenMP over ensemble/cell/basis tiles; measurable reduction in the Nd and Fe stage samples |
| Unpacking | Retained parallel loop | OpenMP over ensemble/cell/basis tiles; output is assigned exactly once |
| J+D multiply | Retained joint Cartesian update | One load of J/Dx/Dy/Dz and Mx/My/Mz updates Bx/By/Bz together |
| Scalar-J multiply | Retained specialized path | DMI arrays and component branches are bypassed entirely |
| Loop order | Retained `ensemble → output basis → input basis → q` | Parent-order comparison and compiler report support contiguous q traversal |
| Explicit SIMD | Not added | GNU report shows q-loop vectorization; no adverse gather/scatter SIMD directive was needed |
| FFT batching/plans | Preserved | Three persistent `plan_many` plans remain; no apply-time plan creation or kernel transform |
| FFT threading | Retained optional provider knob | `UPPASD_FFT_THREADS=1,2,4,8` measured; default remains one provider thread |

## Stage timing

These are post-HAM-09 medians in milliseconds from the five raw samples at
`OMP_NUM_THREADS=2`, `UPPASD_FFT_THREADS=1`. The benchmark fixtures use the
production neighbour-count classes (`z=1338`, `96`, and `6`) and exercise both
2D and 3D J+D layouts. The stage totals are diagnostic; the raw wall sample is
the timing used for speedup comparisons.

| Case | Natom/grid | Pack | Forward FFT | Spectral | Inverse FFT | Unpack | Apply stages | Parity max |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Nd scalar-J | 62500 / 25³ | 0.159 | 0.934 | 0.309 | 1.067 | 0.152 | 2.598 | 1.5e-11 |
| Fe scalar-J | 16000 / 20³ | 0.065 | 0.244 | 0.063 | 0.267 | 0.069 | 0.704 | 5.0e-14 |
| Short scalar-J | 32768 / 32³ | 0.075 | 0.172 | 0.052 | 0.184 | 0.069 | 0.550 | 5.6e-17 |
| 2D J+D | 1024 / 32²×1 | 0.005 | 0.011 | 0.005 | 0.044 | 0.006 | 0.072 | 2.8e-17 |
| 3D J+D, 2 basis | 8192 / 16³ | 0.057 | 0.123 | 0.121 | 0.134 | 0.056 | 0.496 | 1.4e-16 |

The Nd and Fe microbenchmark parity figures include long neighbour-list
accumulations and are therefore less stringent than the focused field parity
tests. The production-style acceptance tests remain at approximately
`1e-12`–`1e-15` as appropriate for their fixture sizes.

## Overall hot-path speedup

The table compares median raw wall samples from the parent hot path with the
optimized path. A value above one is a speedup. These are convolution apply
microbenchmarks, not complete ASD timestep claims.

| Case | Parent median (ms) | HAM-09 median (ms) | Speedup |
|---|---:|---:|---:|
| Nd scalar-J | 3.044 | 2.599 | 1.17× |
| Fe scalar-J | 0.896 | 0.704 | 1.27× |
| Short scalar-J | 0.769 | 0.551 | 1.40× |
| 2D J+D | 0.131 | 0.072 | 1.82× |
| 3D J+D, 2 basis | 0.635 | 0.496 | 1.28× |

The parent and optimized samples are intentionally reported as local evidence,
not as a portable performance guarantee. For the earlier production
crossover, CPU-HAM-05 measured the 32³ short-range control at `0.001852 s`
for DIRECT versus `0.002984 s` for CONVOLUTION, so the optimized convolution
does not change the backend policy: short-range DIRECT remains preferred.

## FFT provider thread scaling

For the 25³, four-basis scalar-J fixture with outer OpenMP held at one, the
five raw post-HAM-09 wall samples were:

| FFTW provider threads | Raw wall samples (ms) | Median (ms) |
|---:|---|---:|
| 1 | 2.461, 2.376, 2.341, 2.265, 2.211 | 2.341 |
| 2 | 2.217, 2.171, 2.312, 2.350, 2.211 | 2.217 |
| 4 | 2.355, 2.279, 2.207, 2.173, 2.243 | 2.243 |
| 8 | 2.221, 2.171, 2.205, 2.149, 2.129 | 2.171 |

The provider abstraction now reads `UPPASD_FFT_THREADS` once during persistent
plan setup, calls the optional FFTW threaded API, and reports the selected
count. The default remains one to avoid nested oversubscription when UppASD
owns outer OpenMP parallelism. oneMKL remains optional and untouched; HAM-10
can compare providers through the same boundary.

## Correctness evidence

The focused CTest subset passed, 7/7:

```text
cpu-ham-field-energy-contract
cpu-ham-backend-policy
cpu-ham-sparse-backend
cpu-ham-target-order
cpu-ham-reduced-stencil
cpu-ham-convolution
cpu-ham-j-plus-d
```

Coverage includes:

- scalar-J DIRECT/convolution parity and field-derived energy;
- scalar-J multi-basis and three-ensemble parity;
- 2D and 3D J+D parity;
- DMI sign and component-transpose negative controls;
- NaN-poison complete-overwrite checks for scalar-J and J+D;
- persistent plans and repeated applies;
- partial-range fallback behavior;
- benchmark coverage with FFTW provider thread counts 1, 2, 4, and 8.

The compiler report was inspected with GNU `-fopt-info-vec-optimized` and
`-fopt-info-vec-missed`. The spectral q loops are vectorized; missed reports
are associated with OpenMP regions, complex multi-array accesses, or setup and
diagnostic code. No explicit SIMD transformation was retained because the
report did not show a safe, necessary improvement beyond the q-loop result.

## Conclusion

HAM-09 retains the useful changes: redundant clears are gone, packing and
unpacking are parallelized, scalar-J arithmetic is specialized, and J+D
Cartesian outputs are produced in one spectral pass. FFT batching and
persistent planning remain unchanged. The measured improvement is strongest
for basis-rich J+D and medium/short synthetic fixtures; the long-range Nd
fixture improves modestly because FFT execution remains dominant. No physics,
eligibility, normalization, or public backend semantics changed.
