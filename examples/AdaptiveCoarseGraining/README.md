# Adaptive coarse graining examples

These are ordinary UppASD text-input runs for the accepted deterministic,
single-ferromagnetic-channel production model. They use a 6 x 2 x 2
body-centred-cubic cell split into six 1 x 2 x 2 coarse blocks.

Run either case from its directory:

```sh
cd static_mixed
../../../build_cpu/bin/sd.f95
```

```sh
cd adaptive
../../../build_cpu/bin/sd.f95
```

The `initial_phase_texture` case starts from an `Initmag 8` spin spiral, runs
ten ordinary atomistic SD preparation steps, and only then constructs the
adaptive topology and selector:

```sh
cd initial_phase_texture
../../../build_cpu/bin/sd.f95
```

Adjust the executable path for your build. `static_mixed/mask.dat` keeps block
1 atomistic; omitted
blocks start coarse and the required interaction-width interface is added
automatically. `adaptive` starts from the selector state and permits
MAX_ANGLE decisions only after complete Heun steps.

The important physics/diagnostic keywords are:

- `block_size_x/y/z`: the physical coarse-cell partition in primitive cells.
- `cg_operator PROJECTED`: variational short-range atom/coarse coupling.
- `cg_mask_mode STATIC|ADAPTIVE`: fixed ownership or synchronized refinement.
- `cg_refine_threshold`, `cg_coarsen_threshold`: MAX_ANGLE misalignment
  thresholds in radians.
- `cg_buffer_blocks`: additional safety dilation; the interaction stencil is
  always included.
- `cg_reconstruction ALIGNED`: reconstructs each coarse atom along its block
  magnetization while retaining its atomic moment magnitude.
- `cg_energy_jump_limit`: maximum accepted resolution-change energy jump in
  joules.
- `cg_diagnostics 3`: resolved topology, ownership states, transitions,
  per-term energies, field/trajectory checksums, timings, and memory.

Supported starting-state modes are `Initmag 1` (random), `2` (random cone),
`3` (momfile direction), `5` (random signs along the momfile direction), and
`8` (spin spiral using `initpropvec`, `initrotvec`, and `initrotang`).
`Initmag 4` remains a restart and is rejected until adaptive state
serialization exists.

`ip_mode S` is also supported as an atomistic preparation stage. Its
`ip_nphase` rows use the ordinary `steps, temperature, timestep, damping`
format. Adaptive ownership is not active during this stage: the final atomic
moments are handed off at the start of the measurement phase, when CG derives
its coarse moments and initial selector state. Other `ip_mode` runners remain
outside the accepted boundary.

The current production boundary requires periodic boundaries, scalar
Heisenberg exchange, zero temperature, fixed-length deterministic Heun
dynamics, and one ferromagnetic dynamical channel. GPU runs may additionally
select `gpu_dipole_mode EWALD3D_FFT`; its basis-resolved uniform-grid field is
mapped into atomistic, interface, and coarse equations independently of the
short-range resolution mask.
