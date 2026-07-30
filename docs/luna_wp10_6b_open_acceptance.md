# WP10.6b Luna acceptance — basis-resolved block-one `NA>1`

Date: 2026-07-28  
Terra commit: `0830897` (`Stabilize OPEN_FFT basis energy derivative E2E`)  
Decision: **GO**

WP10.6b is complete. Terra’s derivative-probe revision clears the prior
numerical gate without changing the independent oracle or threshold. The
independent NA=2 production runtime checks, CUDA multi-basis direct-oracle
suite, full CTest subsets, and sanitizer checks pass. No implementation code
was repaired by Luna; a temporary diagnostic print used to expose the prior
maximum was reverted.

## Scope and authority

The expected values were generated from the independent finite point-dipole
sum in `tests/dipole_validation/open_fft_oracle.py`, using every ordered
target/source basis pair and Cartesian component. No implementation output was
used as expected data.

The checked runtime scope was block one, `NA=2`, fp64, all-open BC, GPU spin
dynamics, and existing basis-resolved macrocell ordering
`macro = basis + NA*cell`. Cases used M=4 nonuniform moments, orthogonal and
skew primitive cells, basis permutations, and centered energy derivatives.

## Independent production results

Actual Fortran `sd.f95.cuda` runs passed against the direct oracle:

| Case | Field maximum (formatted T) | Internal energy maximum | Derivative maximum |
|---|---:|---:|---:|
| Orthogonal 2x1x1, NA=2, M=4 | `4.520126e-38` | `7.444398e-46` | `1.294712e-42` |
| Skew 2x2x1, NA=2, M=4 | `4.951753e-38` | `8.407791e-45` | `5.066219e-42` |

Basis-swapped runs preserved the corresponding field labels and ensemble
energies within the same maxima. The runtime diagnostics reported two basis
channels and four macro channels for the orthogonal 2x1x1 case, confirming the
existing basis-resolved macrocell representation rather than a second
atomistic path.

The accepted NA=1 production matrix remains unchanged after rebuilding from
`0830897`: maximum formatted field error `4.363540e-39` T and maximum relative
formatted energy error `1.163557e-8` in the main matrix.

The independent CUDA multi-basis oracle/plumbing suite also passed its direct
production-builder checks:

```text
multi-basis production L10-skew field_error=6.56858655152258e-47
                           energy_error=7.5265054236196229e-48
multi-basis property ... basis_permutation=6.938893903907228e-18
multi-basis oracle matrix production_field=6.938893903907228e-17
                            production_energy=3.469446951953614e-17
```

That suite exercised every basis/component source impulse on noncubic and
skew grids, M=1/M=4 isolation, reciprocity, basis permutation, and direct
energy. The independent production Fortran runs additionally checked the
input-to-macrocell bridge.

## Revisited CUDA production E2E

After rebuilding `build_gpu_wp9_fp64` from `0830897`, the required CUDA
production layout test passes:

```text
open-production-e2e impulses field_max=3.1263880373444408e-13
                               energy_max=2.7570656868647347e-10
```

The test gate is `field < 2e-11` and `energy < 5e-9`; both pass. Terra changed
only the centered derivative probe from `1e-7` to `1e-5` to avoid fp64
cancellation in the energy difference. The oracle, expected values, and gate
were not weakened.

Commands and outcomes:

```text
cmake --build build_gpu_wp9_fp64 -j2                         PASS
ctest --test-dir build_gpu_wp9_fp64 --output-on-failure \
  -R 'dipole-open-fft-oracle|dipole-open-host-builder|dipole-open-host-goldens|dipole-open-fft-layout|dipole-ewald-host-builder'
                                                               PASS: 5/5
ctest --test-dir build_gpu_wp9_fp64 --output-on-failure \
  -R 'dipole-gpu-fft-convolution|dipole-open-fft-layout|dipole-gpu-wp5-e2e'
                                                               PASS: 3/3
```

The CUDA convolution/multi-basis suite, basis-resolved layout E2E, and
periodic EWALD3D_FFT E2E test all passed.

## Sanitizer and backend

CUDA Compute Sanitizer memcheck on the corrected layout executable reported no
memory errors and the application passed:

```text
ERROR SUMMARY: 0 errors
```

An actual orthogonal NA=2 production `sd.f95.cuda` launch passed CUDA
memcheck:

```text
Gpu: OPEN_FFT geometry staged (... 2 basis channels ...)
Gpu: full-run planned device inventory 27728 B vs TensorMemoryTracker process peak 27728 B
========= ERROR SUMMARY: 0 errors
```

HIP numerical execution was unavailable: no `hipcc`, `rocminfo`, `/dev/kfd`,
or AMD device exists on this runner.

## Ownership and gate

**GO for WP10.6b.** The basis-resolved block-one `NA>1` scope is accepted for
orthogonal and skew cells, M=1/M=4, nonuniform moments, basis permutations,
all source/component impulses, energy, and derivatives.

WP10.7 is now unblocked by WP10.6b. Coarse-block implementation and acceptance
remain separate gates; the already accepted WP10.5 `NA=1` production slice
remains valid.
