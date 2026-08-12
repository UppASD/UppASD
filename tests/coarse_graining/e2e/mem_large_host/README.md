# MEM-LARGE-HOST (RCG-06A)

**Purpose:** production-executable evidence that adaptive coarse graining no
longer places `O(n_atoms)` automatic (stack) arrays in the hot per-step field
path (F-13), reached through the ordinary `sd.f95` executable under a
deliberately constrained process stack limit, not through direct internal
staging.

**Geometry:** a self-contained `NA=1` simple-cubic host, `ncell 20 20 20`
(`Natom = 8000`), nearest-neighbour ferromagnetic exchange only
(`./jfile`, six directions, `J=1.0` mRy), `block_size_x/y/z 5 5 5`
(`n_spatial_blocks = 64`). `Natom=8000` is chosen as the largest size that
keeps `build_unique_bonds`'s pre-existing (out-of-RCG-06-scope, not a target
of this task) near-quadratic directed-pair deduplication tractable in CI
(`source/CoarseGraining/adaptivecgproduction.f90`) while still giving the
eliminated automatic arrays a stack footprint (see below) that reliably
exceeds a small, explicitly documented `ulimit -s`. `cg_mask_mode STATIC`
with no `cg_static_mask_file` (all-fine default) exercises
`evaluate_static_hybrid_operator`'s full per-atom field path
(`StaticHybridOperator`, including the `SmoothProjectedOperator` prolongation/
restriction it always calls regardless of mask) at the full `Natom=8000`
size, twice per step (Heun predictor and corrector).

**Stack-limit arithmetic:** before RCG-06A,
`evaluate_static_hybrid_operator` declared four `(3,n_atoms)` `real(dblprec)`
automatic arrays (`effective_direction`, `ghost_direction`, `atomistic_field`,
`ghost_reaction_field`) plus two smaller `(3,n_blocks)`-shaped ones never
reused across the call. At `Natom=8000`, the four atom-sized arrays alone
demand `4*3*8000*8 = 768000` bytes (~750 KiB). This fixture's driver
(`tests/coarse_graining/run_mem_large_host.py`) runs the real production
binary under `ulimit -s 512` (512 KiB) -- comfortably below that 750 KiB
demand -- and asserts a clean exit. The RCG-06A fix moves those arrays into
persistent, heap-allocated, operator-owned scratch
(`static_hybrid_operator_type%scratch_*`, `source/CoarseGraining/
statichybridoperator.f90`), so nothing in the per-step hot path needs stack
headroom proportional to `Natom` any more.

**Negative control (fault injection):** `tests/coarse_graining/
test_stack_overflow_fault_injection.f90` is a minimal, standalone Fortran
program reproducing exactly the eliminated automatic-array shape and size
(`Natom=8000` passed as a runtime dummy argument, not a compile-time
constant, so the compiler cannot avoid genuine dynamic stack allocation) and
writes to every element to force the pages to be touched. The same driver
runs this program under the identical `ulimit -s 512` and asserts it is
killed by `SIGSEGV` (or otherwise exits abnormally) -- confirming the chosen
`Natom`/stack-limit pair is genuinely discriminating: the eliminated pattern
fails, the accepted (fixed) implementation passes, at the same size and the
same limit, per the blueprint's `§2.3` fault-injection evidence form.

**Units and convention:** ordinary UppASD SI/mRy/lattice-unit conventions
throughout (no coarse-graining-specific unit choice is exercised by this
fixture; it deliberately uses the same physics as `static_all_fine` at a
much larger `Natom`, not a new physical scenario).

**Wired into:** a dedicated CTest target
`adaptive-cg-mem-large-host` (`ctest` label `cg13-cpu`), driven by
`tests/coarse_graining/run_mem_large_host.py`. Not part of
`run_production_e2e.py`'s fixture-dependency inventory since it is invoked
directly by its own CMake `add_test`, matching `run_setup_rejection_matrix.py`
and the other single-purpose RCG-04/05 drivers.
