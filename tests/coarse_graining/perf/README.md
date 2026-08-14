# RCG-09 performance fixtures (`PERF-ATOMISTIC-PROD`, `PERF-CG-SWEEP`)

Fixture provenance per remediation blueprint §6.5.

## What these cases are

Four ordinary UppASD inputs forming two matched pairs. Within a pair the two
`inpsd.dat` files are identical except for the adaptive coarse-graining block —
same `alat`, `ncell`, `cell`, boundary conditions, `posfile`, `momfile`,
`exchange`, `SDEalgh`, `Mensemble`, `temp`, `damping`, timestep, and
measurement settings. The shared `posfile`, `momfile` and `jfile` sit one level
up so that both arms provably read one tracked source rather than two copies
that could drift.

| Case | Texture | Adaptive | Role |
| --- | --- | --- | --- |
| `spiral_atomistic` | spin spiral, `initpropvec 0.125 0 0` | no | `PERF-ATOMISTIC-PROD` baseline |
| `spiral_adaptive` | same spiral | yes, `ADAPTIVE` mask | `PERF-CG-SWEEP` candidate |
| `uniform_atomistic` | uniform along the `momfile` direction | no | maximum-coarsening baseline |
| `uniform_adaptive` | same uniform state | yes, `ADAPTIVE` mask | maximum-coarsening candidate |

Physical model: bcc, two-atom basis, four isotropic exchange shells from
`jfile`; no DMI, no anisotropy, no dipole, no external field; `temp 0.0`;
single ensemble. This is the supported single-channel ferromagnetic capability
and nothing outside it.

## Geometry

`ncell 32 16 16` with a two-atom basis is 16 384 atoms. `block_size 2 2 2`
gives a 16×8×8 block grid: 1024 blocks of 16 atoms. The cell counts are exact
multiples of the block sizes, which the adaptive setup requires.

## Why the initial states are what they are

Both initial states are deterministic and reproducible from the input alone —
no opaque restart file. `Initmag 8` with `initpropvec` constructs the spiral
analytically; `Initmag 3` aligns every moment with the `momfile` direction.
Running either arm twice reproduces the same state bit for bit, which is what
makes the two arms comparable at all.

The spiral is the representative textured case. The uniform state is the
endpoint where the selector can coarsen everything: it is the most favourable
possible input for adaptive coarse graining, so it bounds what the feature can
achieve at this size rather than being a typical workload.

Neither state is a stationary point of the Hamiltonian, so both arms integrate
real dynamics. These fixtures are *performance* fixtures, however: the
correctness oracles live in the `cg13-*` suites, and a run of these cases is not
evidence that the physics is right, only that it costs what it costs.

## Units and conventions

`alat 1.0`, `temp 0.0`, `damping 0.05`, `SDEalgh 1` (Depondt) in both arms.
Note that the production Depondt integrator calls `thermfield.randomize()`
every step even at zero temperature, while the adaptive Heun path is
deterministic and does no RNG at all. That is a real difference between the two
production paths and it is disclosed in the harness output rather than
normalized away.

## Measurement scope

`do_avrg N`, `do_cumu N`, `do_tottraj N`, `plotenergy 0`,
`do_gpu_measurements N`, `do_gpu_correlations N` in every case: RCG-09 times
the Hamiltonian and the integrator, not UppASD's output machinery.

`Nstep` in the committed files is a placeholder. `run_rcg09_perf_e2e.py`
rewrites it to a short and a long value and takes the per-step steady-state
cost as the difference divided by the step difference, so process start, CUDA
context creation, input parsing, neighbour-list construction and adaptive setup
cancel exactly.

## Expected failure mode before the fix

Against `docs/RCG-09_PRODUCTION_PERFORMANCE_EVIDENCE.md`: the adaptive arm is
expected to be *slower* than the atomistic arm by orders of magnitude at this
size, dominated by the atomistic bond phase. A run in which the adaptive arm is
merely somewhat slower, or faster, contradicts the recorded evidence and should
be investigated before it is believed.
