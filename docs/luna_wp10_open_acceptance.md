# WP10.5 Luna acceptance — first production `OPEN_FFT` slice

Date: 2026-07-28  
Commit: `0b2ac5619374498e223b422b2138d79bc91d2c11` (`Reconcile OPEN_FFT GPU memory accounting`)  
Decision: **GO**

WP10.5 is complete. The previously blocking full-run memory-accounting defect
was corrected and independently verified. Production `OPEN_FFT` enablement is
approved for the documented block-one, `NA=1`, fp64 GPU-SD slice.
No implementation files were repaired or changed by this acceptance.

## Scope and environment

I read sections 1–7 of `docs/WP10_OPEN_FFT_BLUEPRINT.md`, including the
recorded WP10.0–WP10.4 handoffs, the independent oracle and red tests, and
only the WP10.4 input/dispatch production diff needed to identify surfaces.
The expected-value authority was the independent fp64 direct finite
point-dipole evaluator; CPU `do_dip=1`, `do_dip=3`, periodic Ewald output, and
GPU output were not used as expected data.

Backend/builds exercised:

- Linux, GCC/GFortran 12.4, CUDA toolkit 13.3, CMake CUDA backend.
- `build_gpu_wp9_fp64` was rebuilt successfully at the acceptance commit.
- `build_gpu_wp9_fp32` was rebuilt successfully at the acceptance commit and
  its OPEN_FFT rejection was rechecked.
- External CUDA verification used one available NVIDIA RTX A4000 (two were
  enumerated by `nvidia-smi`); CUDA device memory reported 16.732 GB total.
  CUDA toolkit version was 13.3 and Compute Sanitizer was 2026.2.1.0.
- `rocminfo`, `hipcc`, `/dev/kfd`, and an AMD device were unavailable, so HIP
  numerical execution was not possible.

## Commands and numerical results

Independent and host layer:

```text
python3 tests/dipole_validation/test_open_fft_oracle.py -v                    PASS (9 tests)
ctest --test-dir build_gpu_wp9_fp64 --output-on-failure \
  -R 'dipole-open-fft-oracle|dipole-open-host-builder|dipole-open-host-goldens|dipole-open-fft-layout|dipole-ewald-host-builder'
                                                                              PASS (5/5)
python3 tests/dipole_validation/test_open_kernel_goldens.py \
  --driver build_gpu_wp9_fp64/bin/dipole_open_kernel_host_driver
                                                                              PASS
```

The ten independent golden cases were: 1x1x1 self, axial pair, transverse
pair, nonuniform 2x3x1, G3>G2 (2x1x3), skew cell, nonuniform basis-resolved
case, four distinct ensembles, thin-film in-plane, and thin-film normal.

Observed host maxima:

| Check | Maximum |
|---|---:|
| Host builder versus independent goldens | `5.3290705182007514e-15` |
| Host direct-convolution matrix | `7.8159700933611020e-14` |
| Dimensionless energy versus `-1/2 sum(M dot B)` | `3.5527136788005009e-15` |
| Reciprocity | `7.1054273576010019e-15` |
| Exact two-basis phase error | `0` |
| Layout red-test field matrix | `1.1368683772161603e-13` |
| Layout red-test energy matrix | `4.5796699765787707e-16` |
| Layout finite-difference derivative mismatch | `1.1677344e-9` absolute, below the `5e-9` test gate |

The host/oracle results cover every field component, exact point self, axial
and transverse signs, no opposite-face wraparound, skew geometry, thin-film
orientations, M=4 isolation, changed moments, zero padding, batch order, and
the energy derivative. The CUDA CTest layout seam also passed the active-only
storage, padded-zero, changed-moment, batch, same-kernel energy, and impulse
matrix checks. Its verbose output included:

```text
open-layout finite-difference derivative=-2.0000000011677344 negative_field=-2
open-layout impulses field_max=1.1368683772161603e-13 energy_max=4.5796699765787707e-16
```

Regression layer:

```text
ctest --test-dir build_cpu --output-on-failure -R 'regression-test|asd-tests'
                                                                              PASS (2/2)
python3 tests/dipole_validation/run.py --binary build_cpu/bin/sd.f95        PASS
python3 tests/dipole_validation/test_periodic_ewald_reference.py             PASS
./build_gpu_wp9_fp64/bin/dipole_ewald_host_tests                             PASS
```

The CPU `do_dip=1` comparison passed as an additional formatted-output
regression, with maximum reported field errors from `1.566e-11` to
`4.900e-11` T and energy relative errors up to `1.919e-8`. It is not the
acceptance oracle. A CPU `OFF` one-cell run also completed successfully.

The existing CUDA periodic convolution and EWALD3D_FFT E2E CTest suite passed
in the configured CTest environment, including the recorded production-slice
periodic maxima (`7.006492321624085e-46` T and
`1.0263416486754031e-48` mRy). This is a regression result only; it does not
substitute for a production `OPEN_FFT` run. The repository
`tests/dipole_validation/run_wp5_e2e.py` result is intentionally recorded as
the EWALD3D_FFT regression; its fixture remains periodic.

## Production `OPEN_FFT` result

Valid temporary production inputs used the per-atom restart format, with
`gpu_dipole_mode OPEN_FFT`, `BC 0 0 0`, `do_dip=0`, and unit macrocell blocks.
The production matrix covered 1x1x1, axial and transverse 2x1x1 pairs,
nonuniform 2x3x1, G3>G2 2x1x3, skew 2x2x1, four distinct M=4 ensembles, and
4x3x1 thin films with in-plane and normal magnetization.

Against the independent direct finite point-dipole oracle, all reported field
components and energies passed the formatted-output layer. Maximum observed
errors were:

| Production check | Maximum |
|---|---:|
| Field absolute error (T) | `4.602701e-39` |
| Energy absolute error (mRy/atom) | `7.192777e-41` |
| Energy relative error | `1.923e-8` |
| Main 1/M=1 matrix energy relative error | `1.923e-8` |

The internal CUDA energy trace reported `total == dip` for every tested
ensemble; for example, the distinct-ensemble trace values were
`8.58341964206586393e-32`, `1.93877240866242943e-31`,
`-2.25198810316629452e-31`, and `6.03131197099516712e-31` in the internal
field-energy units. The formatted Fortran layer is lower precision and is
the source of the approximately 2e-8 relative energy ceiling.

The production OFF/ON exchange comparison passed. The ON-minus-OFF field
delta was
`(-2.9145401579330997e-32,-1.9430267760892e-32,-2.4287834601115e-32)` and
`(6.70407928e-32,-6.70407928e-33,3.35203964e-33)` T for the two atoms;
the total-energy delta was `6.35288800735249e-35` mRy/atom, the OPEN_FFT dipole
column was `6.35288801e-35`, and the exchange-energy delta was zero. The
`do_dip=1` conflict was rejected, so the legacy dipole was not applied twice.
The production predictor/corrector and direct operator seams also passed the
changed-moment, zero-padding, and additive scatter checks.

Finite-difference energy/field consistency passed in the internal fp64 layout
test (`-2.0000000011677344` versus `-2`) and in the independent host oracle.
The production formatted layer was checked against the same oracle; no
opposite-face wraparound was observed in the nonuniform and thin-film cases.

Commands for the CUDA production evidence included:

```text
python3 tests/dipole_validation/run_wp5_e2e.py --binary build_gpu_wp9_fp64/bin/sd.f95.cuda
PYTHONPYCACHEPREFIX=/tmp/lunacache python3 /tmp/luna_extra_prod.py
PYTHONPYCACHEPREFIX=/tmp/lunacache python3 /tmp/luna_exchange_open.py
```

The first command is the periodic EWALD3D_FFT regression; the latter two are
temporary independent OPEN_FFT acceptance drivers and were not added to the
production tree.

Terra re-address acceptance was rerun at commit `0b2ac561` with:

```text
cmake --build build_gpu_wp9_fp64 -j2
cmake --build build_gpu_wp9_fp32 -j2
ctest --test-dir build_gpu_wp9_fp64 --output-on-failure \
  -R 'dipole-open-fft-oracle|dipole-open-host-builder|dipole-open-host-goldens|dipole-open-fft-layout|dipole-ewald-host-builder'
ctest --test-dir build_gpu_wp9_fp64 --output-on-failure \
  -R 'dipole-gpu-fft-convolution|dipole-open-fft-layout|dipole-gpu-wp5-e2e'
PYTHONPYCACHEPREFIX=/tmp/lunacache python3 /tmp/luna_main_prod.py
PYTHONPYCACHEPREFIX=/tmp/lunacache python3 /tmp/luna_extra_prod.py
PYTHONPYCACHEPREFIX=/tmp/lunacache python3 /tmp/luna_exchange_open.py
```

The rerun passed the complete independent production matrix. Main matrix
formatted-layer maxima were `4.363540e-39` T field absolute error and
`1.163557e-8` relative energy error; including M=4 and thin films, the overall
maxima were `4.602701e-39` T and `1.923e-8` relative energy error. The
rejection matrix, including the valid GPU-MC fixture and rebuilt fp32 binary,
also passed.

## Rejection-before-allocation matrix

The production binary was run on temporary copies of the fixture. The
Fortran/input gates rejected the following cases before the production GPU
dipole allocation path:

| Case | Result |
|---|---|
| PBC `P P P` | Reject: requires `0 0 0` |
| Mixed BC `P P 0` | Reject: requires `0 0 0` |
| MC using the valid GPU-MC fixture | Reject: GPU spin dynamics only |
| `do_dip=1` | Reject: requires `do_dip=0` |
| Non-unit block `2 1 1` | Reject: first gate requires `1 1 1` |
| Partial block, 3 cells with block 2 | Reject: blocks must divide all grid axes |
| `NA=2` | Reject: first gate requires `NA=1` |
| `gpu_dipole_alpha` or `gpu_dipole_rcut` nonzero | Reject |
| Nonzero `gpu_dipole_mesh` | Reject |
| Non-`TINFOIL` surface | Reject |
| Non-default `gpu_dipole_tol` | Reject |
| fp32 production binary | Reject: `OPEN_FFT requires fp64 GPU storage; fp32 is not accepted in the first production gate` |

The valid OPEN_FFT input reached the preflight estimate and printed
`Gpu: projected device use 2.2 kB of 13.833 GB free (16.732 GB total)`, then
reached the required geometry, memory, and operator-ready diagnostics. All
rejection cases above completed before `Initiate matrices GPU`, except the
fp32 case, whose source-level first production gate runs within GPU
initialization and rejected before any allocation.

## Memory accounting

Terra’s re-address moved the comparison to the same full-run scope as the
live tracker. For active grid 2x1x1, padded grid 3x1x1, `NA=M=1`, fp64:

```text
base matrices                         576 B
GPU measurement phase               23,184 B
integrator/thermfield                  256 B
OPEN_FFT phase                       1,040 B
full-run projected/tracker         25,056 B
release inventory                       0 B
```

The independent CUDA runs matched the full-run preflight and process tracker:

| Production input | Projected/full-run B | OPEN_FFT phase peak B | Tracker peak B | Release B | Comparison |
|---|---:|---:|---:|---:|---:|
| 2x1x1, M=1 | 25,056 | 1,040 | 25,056 | 0 | +0.0% |
| 2x2x1, M=4 | 34,400 | 6,080 | 34,400 | 0 | +0.0% |
| 4x3x1 thin film, M=1 | 37,336 | 9,440 | 37,336 | 0 | +0.0% |

The phase-local diagnostics also matched Terra’s inventory: 776/216/0 B
persistent/construction/workspace for 2x1x1, 5,048/648/0 B for M=4, and
6,632/2,520/0 B for the thin film. The former 504 B construction figure was
periodic reciprocal-alias staging that OPEN_FFT does not allocate.

CUDA memcheck passed with zero errors for both the standalone layout seam and
the actual 2x1x1 production OPEN_FFT launch:

```text
/usr/local/cuda-13.3/bin/compute-sanitizer --tool memcheck --error-exitcode=99 \
  build_gpu_wp9_fp64/bin/dipole_open_fft_layout_tests
  ERROR SUMMARY: 0 errors

/usr/local/cuda-13.3/bin/compute-sanitizer --tool memcheck --error-exitcode=99 \
  build_gpu_wp9_fp64/bin/sd.f95.cuda   # valid OPEN_FFT restart fixture
  ERROR SUMMARY: 0 errors
```

## Gate and ownership

This is an explicit **GO** for WP10.5. The independent physics, CUDA
production dispatch, memory accounting, OFF path, EWALD3D_FFT regression,
rejection gates, and CUDA memcheck are accepted. HIP numerical execution was
not available because this runner has no HIP toolchain or AMD device; that is
recorded as a limitation, not a CUDA failure.

The production gate is open only for the documented first slice: block
`1x1x1`, `NA=1`, fp64, all-open BC, GPU spin dynamics, `do_dip=0`, and default
Ewald/surface inputs. WP10.6 extensions remain gated.
