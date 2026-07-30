# WP10.7b Terra handoff — OPEN_FFT uniform coarse blocks

Date: 2026-07-28

## Result for Luna

The accepted finite coarse projection from WP10.7a is now connected to the
production `OPEN_FFT` path.  The projection formula and Luna's independent
block-one oracle/goldens were not changed.

The only removed restriction is the former `block_size_x/y/z = 1 1 1`
rejection.  `OPEN_FFT` now accepts a positive, exactly divisible uniform block
shape.  The existing pre-allocation gates still require the macro-grid
quotient, the expected basis-resolved macro count, and the expected full
population in *every* active macro channel.  The Fortran macro-map histogram
is checked against those populations before GPU staging.  Partial edges,
mixed populations, out-of-range maps, and map/population disagreement remain
rejected; no coordinate-derived membership or alternate atomistic runtime was
added.

The runtime remains:

```text
existing atom -> active macro reduction -> zero-padded pack -> FFT
-> existing spectral contraction -> C2R -> existing map scatter
-> existing active-only macro energy reduction
```

Thus normalization, the single field prefactor, and per-physical-atom energy
normalization are unchanged.  The finite projected diagonal supplied by
WP10.7a is retained.

## Enabled input combinations

```text
gpu_dipole_mode  OPEN_FFT
do_dip           0
BC               0 0 0
precision        fp64 GPU spin dynamics
NA               1 or 2 (basis-resolved)
M                1 or 4
block shape      any positive (sx,sy,sz) with Ni % si == 0 on every axis
geometry         regular non-cubic and skew primitive cells are supported
```

Still rejected: fp32, MC, legacy dipoles, periodic/mixed BC, Ewald/surface
overrides, zero blocks, partial edges, active-grid mismatches, irregular maps,
and unequal or incomplete macro populations.

## Tests and accounting

`dipole_open_host_tests` retains Sol's independent `P^T K_point P` checks:
NA=1/2, block one recovery, non-cubic/skew cells, finite projected diagonal,
and nonuniform fine-moment convergence.  It passed locally, together with
Luna's unchanged oracle and host goldens.

The CUDA `dipole-open-fft-layout` seam now covers NA=1/M=1 non-cubic blocks
and NA=2/M=4 skew blocks, nonuniform macro moments, every source impulse,
block-one recovery, padded zeroing, projected-diagonal retention, scatter to
all physical member atoms, and active-only energy normalization.  Its
rejection probes execute the projection gate before stream/FFT/device
allocation for partial edges, active-grid mismatch, and zero block shape.

`dipole-open-fft-coarse-e2e` is registered for a CUDA-equipped runner.  It
uses production `OPEN_FFT` input for NA=1/M=1 and NA=2/M=4 block `2x1x1`
cases, compares coarse scatter to the conjugate average of a block-one run,
checks total/dipole energy, repeats with exchange enabled, and checks partial
edge rejection.  It derives no new golden values.

CUDA fp64 acceptance was run after rebuilding `sd.f95.cuda` from this change:

```text
ctest relevant dipole OPEN_FFT/periodic subset: 8 / 8 passed
wp10-open-coarse-e2e field error:  4.999999675228202e-39 T
wp10-open-coarse-e2e energy error: 0 mRy/atom
```

The CUDA layout/coarse seam reported `1.0658141036401503e-14` maximum field
error and `5.6843418860808015e-14` maximum energy error. CUDA Compute
Sanitizer memcheck passed it with `ERROR SUMMARY: 0 errors`.

The established allocation report remains the accounting authority:

```text
persistent + construction + workspace + 3*pme_Num_macro*M*sizeof(real)
```

and production prints all four terms plus active/padded dimensions, basis,
block shape, ensembles, and prefactor.  The GPU seam checks the persistent and
construction inventories against the estimates.  The CUDA NA=1/M=1 `N=2x1x1`,
block `2x1x1` launch reported `440 B` persistent, `72 B` construction, `0 B`
workspace, and `536 B` peak including `24 B` macro moments.

## Backend status

CUDA fp64 build, numerical E2E, and memcheck are **PASS** on the host GPU.
HIP numerical execution remains unavailable here (`hipcc` is absent). No
HIP/CUDA tolerance or expected value was altered.

## Luna acceptance request

Please independently run the registered CUDA E2E and the WP10.7c matrix,
including the explicit point-operator projection comparison, projected
diagonal, NA=1/2, M=1/4, skew/non-cubic cases, nonuniform convergence,
exchange coexistence, memory report, and all rejection paths.
