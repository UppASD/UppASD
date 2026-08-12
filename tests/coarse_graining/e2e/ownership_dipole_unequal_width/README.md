# RCG-05F: OWNERSHIP-DIPOLE-UNEQUAL-WIDTH, CPU-vs-GPU dipole ownership cross-check

**Purpose:** a tracked, git-committed, genuinely anisotropic-block-shape
(`block_size_x/y/z = 1/2/4`, `buffer_width_blocks = (2, 1, 1)`) `cg_mask_mode
STATIC` fixture, used by `run_dipole_ownership_check.py` to cross-check that
enabling GPU's `gpu_dipole_mode EWALD3D_FFT` (a genuine all-grid FFT dipole
run only GPU can perform in production; CPU rejects `do_dip`/embeds no
dipole path in `evaluate_static_hybrid_operator` at all -- see
`source/CoarseGraining/statichybridoperator.f90:115-119`) does not perturb
the fine/buffer/coarse ownership map, at RCG-05C's full per-block identity
resolution, relative to the same geometry's CPU (no-dipole) run.

## Why STATIC mode specifically

Under `cg_mask_mode STATIC` the GPU dilation kernel is never invoked at all;
GPU only ever copies CPU's setup-time `block_state` once
(`ownership_aniso_buffer/README.md` documents this same fact for the
ADAPTIVE case, by contrast). This makes STATIC mode the only mode where a
dipole-on/dipole-off ownership comparison is not confounded by dipole's
legitimate physical effect on later ADAPTIVE-mode selector decisions (a
nonzero dipole field genuinely changes the LLG trajectory, which could
legitimately change which blocks an ADAPTIVE selector later refines/
coarsens -- that would not be a bug, but it would make a literal ownership-
map equality assertion meaningless). STATIC mode's ownership is fixed at
setup, independent of the field driving the subsequent dynamics, so an exact
match is the correct, unconfounded expectation here.

## Construction

Geometry/`jfile`/`mask.dat`/`GENERATOR_MANIFEST.json` are produced by
`static_topology_oracle.unequal_width_orthogonal_fixture` (RCG-05B), reused
unchanged:

```python
unequal_width_orthogonal_fixture(
    block_size=(1, 2, 4), ncell=(6, 8, 8), bond=(2.0, 0.0, 0.0, 1.0),
    fine_block_ids={1}, cg_buffer_blocks=0, alat=1.0, simid="cg105fdipuw",
)
```

`block_grid = 6 4 2` (768 atoms, shared `../posfile`/`../momfile` 2-atom
basis). `fine=1 interface=29 coarse=18`, `buffer_width_blocks=(2,1,1)` --
genuinely anisotropic, distinct from `ownership_aniso_buffer`'s `(1,2,3)`
block shape and from `gpu_fft_static_mixed`'s `(1,2,2)`.

This `inpsd.dat` itself carries no `do_gpu`/`gpu_dipole_mode` setting (so it
runs unmodified, dipole-free, on the CPU binary as the reference). The GPU
comparand is produced at run time by `run_dipole_ownership_check.py`, which
copies this fixture and appends `run_moving_backend_parity.GPU_ENABLE_BLOCK`
plus a dipole-enabling block (`gpu_dipole_mode EWALD3D_FFT` and its
surface/tolerance/mesh settings, matching `gpu_fft_static_mixed`'s already-
validated values) -- never baking dipole into the tracked file itself, since
a CPU binary given `gpu_dipole_mode EWALD3D_FFT` without `do_gpu Y`/
`do_gpu_llg Y` is correctly rejected by
`adaptivecgproduction.f90:766-773`'s validation.
