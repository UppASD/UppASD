# GPU predictor profiling

`run.py` benchmarks the standard `tests/bccFe_cuda` fixture without modifying
it.  Each warm-up and measured repetition gets a fresh private copy.  The
harness forces deterministic, low-noise GPU SD input (T=0, Depondt, disabled
GPU measurements/correlations and regular output) and prints median wall time
and range across at least three measured repetitions.

Run both the ordinary and convolution paths (the default):

```bash
python3 tests/gpu_profile/run.py --binary build_ab/sd --steps 10000 --repetitions 3
```

Use `--ncell NX NY NZ` to scale the BCC system. `--convolution Y` benchmarks
only the convolution path, `--convolution N` only the ordinary path, and
`--convolution both` is the default. Use `--keep` or `--workdir DIR` to retain
private runs; `--output DIR` writes `profile-summary.json` there. Generated
profiles and run outputs are intentionally not tracked.

CUDA timeline profiling:

```bash
nsys profile --trace=cuda,nvtx,osrt -o /tmp/uppasd-p6 \
  python3 tests/gpu_profile/run.py --binary build_ab/sd --steps 10000 --convolution Y
```

For kernel-level CUDA inspection, place the profiler before the executable
inside a retained fixture copy, or profile the harness while filtering its
child process where supported:

```bash
ncu --target-processes all --set full -o /tmp/uppasd-p6-ncu \
  python3 tests/gpu_profile/run.py --binary build_ab/sd --steps 10000 --convolution Y
```

On ROCm systems, use `rocprof` (or `rocprofv3`) similarly:

```bash
rocprof --hip-trace python3 tests/gpu_profile/run.py --binary build_hip/sd --steps 10000 --convolution Y
```

Compare before and after with exactly the same binary configuration, steps,
cell count, repetitions, and mode. The fused predictor should replace the
former local-field, damping/STT, rotation, and `emomM` elementwise launches
with one `EvolveFirst` launch; timeline/kernel metrics should therefore show
fewer predictor launches and lower global-memory traffic.
