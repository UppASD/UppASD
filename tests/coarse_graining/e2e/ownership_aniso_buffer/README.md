# RCG-05C: OWNERSHIP-ANISO-BUFFER, CPU/GPU buffer-shape defect fixture

**Purpose:** a tracked, git-committed, genuinely anisotropic-block-shape
`cg_mask_mode ADAPTIVE` fixture that engages the real GPU adaptive-transition
runtime (`GpuAdaptiveRuntime::proposeSelectorState` ->
`dilateAdaptiveState`), so `ownership_map_comparator.py`/
`run_ownership_map_comparator.py` can demonstrate, through the real
executable rather than a disposable reproduction, that the GPU/CUDA/HIP
dilation kernel collapses the CPU-correct three-component
`buffer_width_blocks(3)` to a single scalar
(`adaptivecgproduction.f90:615-616`, `int(maxval(...))`) and dilates
isotropically (`gpuAdaptiveRuntime.cpp:322-349`) instead of per axis
(`statichybridoperator.f90:196-256`). This is the permanent, reusable
successor to RCG-05A's disposable `inpsd.dat`/mask reproduction, and RCG-05D's
own negative control.

## Construction

Geometry/`jfile`/`mask.dat`/`GENERATOR_MANIFEST.json` are produced by
`static_topology_oracle.unequal_width_orthogonal_fixture` (RCG-05B), reused
unchanged rather than hand-authored:

```python
unequal_width_orthogonal_fixture(
    block_size=(1, 2, 3), ncell=(6, 10, 9), bond=(2.0, 0.0, 0.0, 1.0),
    fine_block_ids={1}, cg_buffer_blocks=0, alat=1.0, simid="cg105cani",
)
```

`block_size_x/y/z = 1/2/3` on `ncell 6 10 9` gives `block_grid = 6 5 3`
(1080 atoms, shared `../posfile`/`../momfile` 2-atom basis, single exchange
shell `dx=2.0` -- an exact same-sublattice lattice-vector neighbour, so
`neighbourmap.f90` maps it without any tolerance ambiguity). The correct,
per-axis CPU width for this shell/block-shape (`buffer_width_blocks`,
`statichybridoperator.f90:180-187`) is **`(2, 1, 1)`** -- genuinely
anisotropic. Dilating the single FINE seed (block 1) by this width gives
`fine=1 interface=44 coarse=45` (`compute_expected_topology`, matching this
fixture's own `GENERATOR_MANIFEST.json`). The isotropic collapse GPU actually
stages is `max(2,1,1)=2`, i.e. `(2,2,2)`
(`compute_isotropic_dilation_topology`), which for the *same* single seed
gives `fine=1 interface=74 coarse=15` -- already a large, analytically
predicted divergence before any dynamics run at all.

`cg_mask_mode ADAPTIVE` (not `STATIC`) is required: GPU's dilation kernel
(`dilateAdaptiveState`) is only ever launched from
`GpuSimulation::advanceAdaptiveStep`'s periodic
`proposeSelectorState`/`publishProposedState` cycle
(`gpuSimulation.cpp:1101-1112`), gated on `adaptiveMaskEnabled`
(`cg_mask_mode ADAPTIVE`) and `step % cg_update_interval == 0`. Under
`cg_mask_mode STATIC` the GPU never re-dilates at all -- it only ever copies
the CPU-computed (correct) initial `block_state` once at setup -- so a
`STATIC` fixture, however anisotropic, cannot expose this defect
operationally; this is why RCG-05C could not simply reuse RCG-05B's
generator output as a `STATIC` fixture and had to add the `ADAPTIVE`-mode
selector configuration below.

`Initmag 8` with a small `initpropvec` gives a uniform-pitch spin-spiral
modulation along `x` (the same built-in mechanism
`moving_static_mixed_bs1` etc. use, `magnetizationinit.f90:402-462`), so the
`MAX_ANGLE` selector criterion (`q = 1 - cos(angle)` between exchange-bonded
neighbours, `blockselector.f90:503`) sees a genuine, small, nonzero
misalignment near the fine/coarse interface rather than an exactly-uniform
state -- `initpropvec_x = 0.45/(4*pi) ~ 0.0358` was chosen so the per-bond
misalignment angle (`propvec_x * bond_dx = propvec_x * 2`) is `~0.45` rad,
i.e. `q = 1-cos(0.45) ~ 0.0985`, comfortably inside the
`(cg_coarsen_threshold, cg_refine_threshold) = (0.05, 0.20)` dead band for
most blocks while still crossing `cg_refine_threshold` at specific
already-atomistic boundary blocks (see `run_ownership_map_comparator.py`'s
own recorded raw output for exactly which blocks and why).

A perfectly uniform initial state (no spiral) was tried first and rejected:
with `q` exactly `0.0` everywhere, `score <= cg_coarsen_threshold` (`0 <=
0.05`) is true for *every* currently-non-coarse block, so both backends
coarsen the entire initial mask on the very first update -- vacuous (nothing
left to dilate) rather than defect-sensitive, and confirmed to also trip a
*second*, separate, already-documented capability gap (GPU's
`hardAtomisticBlockMask` is sourced only from the polarization gate, not
from `cg_static_mask_file`; see `gpuSimulation.cpp`'s own RCG-03 comment
above the `Gpu: AdaptiveCG resolved diagnostics=` print) that would have
conflated two different defects in one fixture.

## What the comparator actually finds (see `run_ownership_map_comparator.py`
## output / `docs/RCG-05_GEOMETRY_OWNERSHIP_EVIDENCE.md` RCG-05C section)

- **CPU** ends this one-step run with `fine=3` (blocks `{1, 31, 61}` --
  the seed plus its two z-periodic images, which independently crossed
  `cg_refine_threshold`), `interface=42`, `coarse=45`. This is *exactly*
  reproduced by feeding CPU's own reported FINE set back into
  `compute_expected_topology` (the correct, per-axis oracle) -- CPU's
  runtime ownership rebuild (`rebuild_static_hybrid_ownership`) is
  self-consistent with the directional formula, block for block.
- **GPU/CUDA** ends the identical input with a different FINE set (`fine=25`,
  a strict superset driven by the separate hard-mask-porting gap noted
  above, not by the buffer-width defect) and `interface=65 coarse=0`.
  Feeding GPU's own reported FINE set into `compute_expected_topology` does
  **not** reproduce GPU's actual map. Feeding the same GPU FINE set into
  `compute_isotropic_dilation_topology` (the `(2,2,2)` collapse) reproduces
  GPU's actual map **exactly**, block for block -- isolating the buffer-width
  scalarization as the specific, source-line-attributable cause of GPU's
  reported shape, independent of the (separately noted, not conflated) seed
  -selection disagreement.

No production fix is present in this fixture or its README; it exists to
demonstrate the defect for RCG-05D to repair and re-verify against.
