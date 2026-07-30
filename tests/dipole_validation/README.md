# Dipole validation fixtures

Generate the deterministic fixtures from the user-maintained SC bases:

```bash
python3 tests/dipole_validation/generate_cases.py --output /tmp/uppasd-dipole-cases
```

Each finite case contains `analytic_reference.json`, calculated independently
from the point-dipole tensor in SI units. These references are not copied from
UppASD CPU output. `sc2d_macro_edges` is deliberately an aggregation/layout
case rather than an atomistic analytic oracle.

The generator enables `do_prn_beff Y`, which makes UppASD write
`befftot.<simid>.out` in Tesla. After running a finite generated case, compare
that file with its analytic reference:

```bash
python3 tests/dipole_validation/check_fields.py /tmp/uppasd-dipole-cases/sc3d_pair_axial_x
```

`check_fields.py` also checks the `Dip` column in `totenergy.<simid>.out`
against the analytic energy per atom in mRy. Run the complete finite suite with:

```bash
python3 tests/dipole_validation/run.py \
  --binary build_cpu/bin/sd.f95
```

The generated finite cases explicitly set open boundaries because CPU
`do_dip=1` is a finite direct sum. Periodic 3D and periodic-slab validation
will be added separately with reciprocal-space analytic references and an
explicit `k=0` convention.

## Validation status

The finite CPU `do_dip=1` single, pair, and slab field checks pass against the
independent analytic references. This establishes the direct finite dipole
field convention, but does **not** yet validate:

- macrocell `do_dip=2`;
- CPU FFT `do_dip=3`;
- 2D/3D periodic sums; or
- GPU dipoles.

The finite direct-field and direct-energy checks are in place. The remaining
work is macrocell `do_dip=2`, CPU FFT `do_dip=3`, periodic sums, and GPU
dipoles.

## CPU captures for macrocell and FFT modes

`do_dip=2` and `do_dip=3` need regression baselines, but they must not be
presented as analytic validation.  After generating and running one of those
cases, capture the deterministic CPU result with:

```bash
python3 tests/dipole_validation/capture_references.py CASE --binary build_cpu/bin/sd.f95
```

It writes `CASE/cpu_reference.json`, containing complete printed field blocks
and the Dip energy.  A GPU run can then be compared with:

```bash
python3 tests/dipole_validation/check_captured_reference.py CASE
```

The initial GPU mode and its exact open-boundary macrocell convention are
documented in `docs/GPU_FFT_DIPOLE_DESIGN.md`.

## Complete-kernel GPU FFT validation

The complete periodic Builder A tensor is validated through a CMake/CTest
target, never by compiling a temporary CUDA driver from Python.  Configure a
double-precision GPU build, then build and run:

```bash
cmake -S . -B build_gpu -DUPPASD_GPU_BACKEND=CUDA -DUPPASD_PRECISION=DOUBLE -DBUILD_TESTING=ON
cmake --build build_gpu --target dipole_gpu_fft_tests -j1
ctest --test-dir build_gpu --output-on-failure -R '^dipole-gpu-fft-convolution$'
```

The target first runs non-physical delta convolutions over the documented
grid/basis/ensemble matrix, then compares periodic Builder A convolution
against direct host convolution.  Field and energy downloads are test-only.

The generator supplies `sc2d_macro_uniform` and
`sc2d_macro_uniform_m4` as the first capture cases.  It also keeps
`sc2d_macro_edges` as a required partial-edge gate; do not replace it with a
divisible-grid test if the CPU reference fails to run.
