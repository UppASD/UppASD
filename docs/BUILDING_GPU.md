# Building UppASD with GPU support

Everything is driven by CMake presets. `cmake --list-presets` shows what is available.

Two cache options control the build:

| Option | Values | Default |
| --- | --- | --- |
| `UPPASD_GPU_BACKEND` | `OFF`, `CUDA`, `HIP` | `OFF` |
| `UPPASD_PRECISION` | `DOUBLE`, `SINGLE`, `MIXED` | `DOUBLE` |

`USE_CUDA` and `USE_HIP` still work but are deprecated aliases and warn; they map onto
`UPPASD_GPU_BACKEND`. `UPPASD_PRECISION=MIXED` is rejected at configure time — the mode is
specified against typedefs that do not exist yet.

GPU builds require **CMake >= 3.24**. CPU-only builds work with the project minimum (3.13).

## Build types

The standard CMake configurations are supported: `Debug`, `Release` (the default),
`RelWithDebInfo` and `MinSizeRel`, plus a project-specific `Coverage`. Names are accepted
case-insensitively and normalised, so the historical `RELEASE`/`DEBUG` spellings keep working.
`TESTING` is the old name for `Coverage` and still works with a deprecation warning.

Use `RelWithDebInfo` — not `Debug` — for profiling: it is optimised but keeps symbols and line
information. The `cuda-relwithdebinfo` preset also passes `-lineinfo` so Nsight can attribute
device-side samples to source lines.

`Coverage` builds with gcov instrumentation (`--coverage`) and is wired up for **gfortran
only**; other compilers configure with a warning and no instrumentation. Coverage data lands
next to the object files:

```sh
cmake --preset coverage && cmake --build --preset coverage --parallel 8
ctest --test-dir build/coverage        # or run the binary directly
gcov -o build/coverage/CMakeFiles/asdlib.dir/source/... <source>.f90
```

## NVIDIA, consumer laptop or workstation

```sh
cmake --preset cuda
cmake --build --preset cuda --parallel 8
```

`CMAKE_CUDA_ARCHITECTURES` defaults to `native`, i.e. whatever card is in the machine. On a
machine with no visible GPU — compile-only CI, a login node, cross-compilation — pass an
architecture explicitly, otherwise `native` cannot resolve:

```sh
cmake --preset cuda -DCMAKE_CUDA_ARCHITECTURES=80      # or "80;90" for a fat binary
```

## NVIDIA HPC with modules

Load a CUDA toolkit and a Fortran compiler, then use the same preset. If the site provides
MKL and you want it, add `-DUSE_MKL=ON`; BLAS/LAPACK selection is otherwise left to CMake's
own search (see `BLA_VENDOR`).

```sh
module load cuda fortran-compiler        # site-specific names
cmake --preset cuda -DUSE_MKL=ON
cmake --build --preset cuda --parallel 8
```

## AMD, RDNA laptop or desktop

ROCm packages are found through `CMAKE_PREFIX_PATH`. Point it at the ROCm install; nothing is
hardcoded in the build files.

```sh
export CMAKE_PREFIX_PATH=$ROCM_PATH
cmake --preset hip-rdna
cmake --build --preset hip-rdna --parallel 8
```

The `hip-rdna` preset targets `gfx1100;gfx1101;gfx1102`. Override with
`-DCMAKE_HIP_ARCHITECTURES=<arch>` for a different card.

## AMD CDNA (MI200 / MI300)

```sh
export CMAKE_PREFIX_PATH=$ROCM_PATH
cmake --preset hip-cdna                  # gfx90a;gfx942
cmake --build --preset hip-cdna --parallel 8
```

## LUMI

Module names change with the LUMI stack — check `module spider rocm` against what is current
rather than trusting the versions below.

```sh
module load LUMI/24.03 partition/G
module load PrgEnv-amd
module load rocm
export CMAKE_PREFIX_PATH=$ROCM_PATH

cmake --preset lumi
cmake --build --preset lumi --parallel 8
```

The `lumi` preset pins `gfx90a`, sets `CMAKE_Fortran_COMPILER=ftn` and turns on `ON_LUMI`
(Thrust complex for the HIP correlation kernels). The Cray wrappers supply libsci, so no BLAS
selection is needed — do not pass `USE_MKL`.

## Development and debugging presets

| Preset | What it is for |
| --- | --- |
| `cpu`, `cpu-debug` | No GPU backend |
| `cuda-debug` | Debug CUDA build |
| `cuda-single` | `UPPASD_PRECISION=SINGLE` |
| `cuda-fastcopy` | `USE_FAST_COPY=ON`, the pinned-transfer path |
| `cuda-ptds` | `--default-stream per-thread` |

`cuda-ptds` is the regression net for stream-ordering bugs: per-thread default stream removes
the implicit synchronisation edge between the legacy default stream and the blocking
`workStream`, so a missing explicit dependency shows up as wrong results instead of passing by
luck. Keep it building and passing.

`debug`, `release`, `debugCuda` and `releaseCuda` are kept as deprecated aliases of
`cpu-debug`, `cpu`, `cuda-debug` and `cuda`, writing to their original build directories.

## Testing a GPU build

```sh
python3 tests/gpu_regression/run.py --binary build/cuda/bin/sd.f95.cuda
python3 tests/gpu_regression/sanitize.py --binary build/cuda/bin/sd.f95.cuda
```

`sc_dm_uniaxial` currently fails in every configuration — it segfaults in
`gpusim_initiatematrices` and predates the GPU rework. Expect **11 pass, 1 fail**; treat
anything else as a regression.
