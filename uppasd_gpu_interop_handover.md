# UppASD GPU interop handover

This handover is for the next session working on the Fortran/C++ boundary in
`/home/andersb/SD/UppASD_gpu_hip_cu`.

## Current state

The recent GPU audit work is committed through:

- `18c3cf5 gpu: dispatch RNG by storage precision` (P7)
- `e44cb63 gpu: coalesce tensor exchange layout` (P4)

The existing user edit in `tests/FeCo_cuda/inpsd.dat` must be preserved. Do
not touch unrelated untracked files or build directories.

Completed work relevant to the boundary:

- C4 is complete: scalar neighbour lists use host staging uploads rather than
  mutating Fortran-owned arrays in place.
- P4 is complete: `j_tensor` is transformed during upload to the device-only
  layout `[reduced site, output axis, input axis, neighbour]`.
- P7 is complete: device RNG precision dispatch and Tensor device-indexing
  hygiene are in place.

Do not reopen the deferred MC work (C5/C6/P1/P5), broader FP32/MIXED storage
conversion (I2), or unrelated tensor layout work while doing I1.

## Goal: I1 — make setup bindings explicit and portable

`source/chelper.f90` has explicit `bind(C)` interfaces only for:

- `FortranData_setCorrelations` → `fortrandata_setcorrelations_`
- `FortranData_setGpuGeometry` → `fortrandata_setgpugeometry_`

The remaining setup calls in `FortranData_Initiate` currently rely on implicit
Fortran external resolution and the C++ trailing-underscore symbol names. This
depends on compiler naming and on hidden CHARACTER-length arguments being
appended and ignored. Replace those implicit calls with explicit `bind(C)`
interfaces, matching the two existing interfaces' style.

Convert these six calls:

| Fortran call | C++ entry point | C++ implementation |
| --- | --- | --- |
| `FortranData_setFlags` | `fortrandata_setflags_` | `source/gpu_files/fortranData.cpp` |
| `FortranData_setConstants` | `fortrandata_setconstants_` | same |
| `FortranData_setHamiltonian` | `fortrandata_sethamiltonian_` | same |
| `FortranData_setLattice` | `fortrandata_setlattice_` | same |
| `FortranData_setMeasurables` | `fortrandata_setmeasurables_` | same |
| `FortranData_setInputData` | `fortrandata_setinputdata_` | same |

`source/gpu_files/fortranData.hpp` is the type/ordering reference for the C++
side. The `extern "C"` definitions near the bottom of
`source/gpu_files/fortranData.cpp` are the definitive ABI signature to match.
The call sites, in exact argument order, are in `source/chelper.f90` around
lines 421–463.

### ABI rules for I1

1. Use `bind(C, name="fortrandata_set…_")` explicitly. Keep the current C++
   symbol spellings during this pass; changing C++ names is unnecessary risk.
2. Import and use `iso_c_binding` kinds:
   - Fortran default integer values passed to C++ `int*` should become
     `integer(c_int)` at the interface.
   - C++ `unsigned int*` parameters should also use `integer(c_int)`; the
     values are non-negative counts/flags, and C interoperability has no
     unsigned integer kind.
   - `real*` must be `real(c_double)` for the current double-precision bridge.
   - one-character flags must be `character(c_char)`.
   - correlation complex data already uses `complex(c_double)` and is the
     model to follow.
3. For arrays, use assumed-size dummies (`arg(*)`) in the `bind(C)` interface,
   not allocatable or assumed-shape dummies.
4. Mark pointer-like inputs `intent(inout)` only where the existing interface
   convention requires it. These setter calls retain addresses; the intent does
   not add runtime copying. Match the existing Correlations/Geometry style for
   consistency.
5. Do not use `integer(c_size_t)` merely because a quantity is a size in C++.
   These six setter functions currently take `int*`/`unsigned int*`; match the
   actual C++ signature. `c_size_t` belongs to the separate callback audit.
6. Keep the Fortran-to-C++ bridge double precision. `real` in C++ aliases
   `double` for this supported configuration.

### Suggested implementation sequence

1. Read the full `FortranData_Initiate` call block and the six C++ entry-point
   definitions before editing, then write all six interface declarations in
   the existing `interface` block in `source/chelper.f90`.
2. Add `import` lists only for the kinds each declaration needs. Preserve the
   argument order exactly; these are long APIs, so accidental reorderings are
   the main risk.
3. Compile after the interface addition before any cleanup/refactor. A
   compiler error is useful evidence of a mismatched Fortran declaration.
4. If compilation exposes default-`integer` variables that are not compatible
   with `integer(c_int)`, do not cast addresses. Either establish that the
   project's default integer is C-compatible or stop and decide on a narrow
   adapter strategy. Do not silently reinterpret pointers.
5. Keep changes limited to I1 files, normally `source/chelper.f90` and (only
   if a real ABI mismatch is found) the matching declarations/definitions in
   `source/gpu_files/fortranData.{hpp,cpp}`.

## Out of scope, but record for the next dedicated pass

### I2 — FP32/MIXED bridge conversion

Fortran arrays use `real(dblprec)`, while C++ `real` becomes `float` under
`SINGLE_PREC`. Therefore a SINGLE_PREC bridge currently aliases double
storage as float storage and is invalid. The eventual solution is an explicit
double-host to float-device conversion/staging layer, extending the existing
upload transforms. It must not be attempted as part of I1.

### I3 — C++ callbacks into Fortran

`fortran_print_measurables` in `source/chelper.f90` has allocatable dummies
and is declared from C++ in `source/gpu_files/c_helper.h` via `FORTNAME`.
C++ cannot safely pass Fortran allocatable descriptors as raw pointers. Audit
this call chain separately; use a `bind(C)` scalar/assumed-size interface or
the Fortran descriptor mechanism (`CFI_cdesc_t`) as appropriate. Do not fold
this callback rewrite into the setup-binding task.

The `FORTNAME`-based wrappers in `c_helper.h` (`fortran_measure*`,
`fortran_flush_measurements`, RNG, and module-variable access) are also a
separate outbound-C++-to-Fortran ABI modernization task. I1 only covers the
six Fortran-to-C++ setup setters above.

## Validation and handoff

Run at minimum:

```bash
cmake --build build_ab --parallel 1
git diff --check
git status
```

If a CUDA device is available, run the deterministic GPU regression matrix,
including `kagome_tensor`, after the build. The current host used for P4/P7
had no CUDA-capable device and no `hipcc`; state that limitation rather than
claiming runtime validation.

Before committing, ensure only focused interop files are staged. A suitable
commit message is:

```text
gpu: make Fortran data bindings explicit
```

