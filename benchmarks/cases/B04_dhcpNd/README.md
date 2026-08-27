# B04_dhcpNd — dhcp Nd, very long-range exchange interaction cloud

Admission record for WP-08D (blueprint §3 "B04 — dhcp Nd: very long-range
interactions"). This case is the framework's deliberately interaction-heavy
neighbour-list workload: a real production dhcp Nd deck whose exchange
cutoff reaches far enough that per-atom coordination (~1338
directed-interactions/atom) is roughly 14x `B01_bccFe`'s (96) and ~223x
`B03_skyrmion3D`'s (6), while `natom` alone stays modest. Per the blueprint,
the physically intended interaction cloud is retained exactly as supplied
by the maintainer -- no shell truncation, no cutoff reduction, no discarded
weak interactions, no simplified basis, "merely to make benchmarking
easier."

## A. Template audit

### Template selection: maintainer-supplied production deck

`examples/SimpleSystems/dhcpNd` is a real maintainer production input
(present in the repository with its own historical `out.cpu`/`out.gpu`/
`meminfo`/`restart.*` run artifacts from a prior manual run at this
project's own `ncell 32 32 32` baseline) satisfying the dependency this WP
was given ("maintainer-approved dhcp Nd production input available"). It
was copied into `benchmarks/cases/B04_dhcpNd/template/`:

| File | Role |
| --- | --- |
| `inpsd.dat` | Input deck (one fix applied, see below) |
| `jfile.cut5` | Exchange list: 5,352 directed `(iatom, jatom, n1, n2, n3, J)` bonds, `maptype 3` |
| `momfile` | Four-atom basis moments (mmom=2.0 uniformly) |
| `posfile.cartesian` | Four-atom basis positions (the file `inpsd.dat` actually references; the sibling `posfile.direct` in the source example directory is unreferenced and was not copied) |
| `qpoints.bz` | q-point file kept present but inert, same finding as `B03_skyrmion3D` (see below) |

Checksums of every copied file were diffed against
`examples/SimpleSystems/dhcpNd`'s originals to confirm byte-identical
content except the one `inpsd.dat` fix below; the source example directory
itself is unmodified.

### The one required fix: a duplicate `gpu_mode` keyword bug in the maintainer's own deck

The source `inpsd.dat` declares `gpu_mode` **twice** (`gpu_mode 0` mid-file,
then `gpu_mode 1` again near the end, right before `do_gpu_measurements Y`).
`source/Input/inputhandler.f90`'s keyword loop reads sequentially and
simply overwrites the Fortran variable on each occurrence
(`case('gpu_mode'): read(ifile,*,...) gpu_mode`), so the *second*
occurrence is what the production binary actually uses -- the maintainer's
own effective default was GPU-on, which is why the source directory's own
historical run artifacts include `out.gpu`/`out.cuda*` files.

This breaks this framework's own override mechanism.
`harness/cases.py::_apply_overrides` replaces only the *first* line whose
keyword matches; a second, later, unmatched `gpu_mode` line is left
untouched. Requesting a CPU sample (`gpu_mode` not overridden, or overridden
to `0`) would silently still execute as `gpu_mode 1` at runtime -- exactly
the "CPU-only build silently no-ops the whole simulation" failure mode
`B01_bccFe/README.md`'s `gpu_mode` per-run-override rationale already warns
about, except here the danger is the *opposite* direction (a requested CPU
run silently dispatching to the GPU driver, or on a CPU-only build,
crashing/no-op'ing per that same finding).

**Fix applied**: collapsed the two `gpu_mode` lines to one, valued `0`
(matching `B01_bccFe`/`B02_skyrmion2D`/`B03_skyrmion3D`'s established
baseline-off-override-to-1-for-GPU convention). A duplicate, harmless
`do_reduced Y` line (both occurrences `Y`, no behavioural difference either
way) was also collapsed to one for the same override-mechanism-correctness
reason, though it carries no separate finding. No other line was changed;
`Nstep 10000`, `Temp 2.00000001 K`, `Initmag 1`, `ip_mode N`, `tseed 1`, and
every Hamiltonian/lattice keyword are the maintainer's own baseline values,
verbatim.

### Physical summary

| Property | Value | Source |
| --- | --- | --- |
| Lattice | dhcp (double hexagonal close-packed), four-atom basis | `posfile.cartesian` |
| Basis positions (reduced Cartesian, a=1) | `(0,0,0)`, `(0,0,1.61284)`, `(0,-0.57735,-0.80642)`, `(0,0.57735,0.80642)` | `posfile.cartesian` |
| Cell | hexagonal, `a=b`, `gamma=120deg`, `c/a=3.225675934279231` | `inpsd.dat` |
| Symmetry | `Sym 0` -- `jfile.cut5`'s bonds are the full directed list, no expansion | `inpsd.dat` |
| Hamiltonian | Exchange only (no `dm`/`anisotropy`/`dipole` keyword present) | `inpsd.dat` |
| Exchange cutoff | Long-range RKKY-type; lattice-translation offsets reach `n1,n2=+-6`, `n3=+-2` | `jfile.cut5` (measured, see §B) |
| Directed interactions/atom | **1338 mean/median exactly** (site A: 1340, site B: 1336, 50/50 split), measured directly (§D) | `struct.<simid>.out` |
| Moment | `mmom=2.0` on all four basis atoms | `momfile` |
| Boundary conditions | Fully periodic (`BC P P P`) | `inpsd.dat` |
| Initial state | `Initmag 1` (random, seeded by `tseed`) -- unlike `B01_bccFe`/`B02_skyrmion2D`/`B03_skyrmion3D`'s uniform `Initmag 3` | `inpsd.dat` |
| Temperature (maintainer default) | **2.00000001 K**, used verbatim by `dhcpNd_finite_t`, not invented | `inpsd.dat` |
| Integrator | `SDEalgh` unset -> Heun default, damping 0.5, timestep 1.0e-15 s | `inpsd.dat` |
| Initial phase | `ip_mode N`: none executes despite the `ip_mcanneal` block being present (`source/System/damping.f90:155`, same dead-input pattern `B01_bccFe/README.md` documents) | `inpsd.dat` |
| Diagnostics | `do_prnstruct 2`, `do_meminfo 1` (both maintainer baseline, not framework-added) | `inpsd.dat` |

`do_meminfo 1` being the maintainer's own baseline (not something this
admission added) is what makes real, production-measured memory-footprint
data available for this case "where possible" per blueprint B04 section B,
without any harness change (§D/§F below).

### `qfile` stays present but inert (same finding as `B03_skyrmion3D`)

`qpoints F` here too, so `source/Correlation/qvectors.f90::read_q` is
called unconditionally regardless of `do_sc N` -- removing `qfile` crashes
with `Cannot open file ''`, the same real source finding
`B03_skyrmion3D/README.md` documents. `qpoints.bz` (1.7 MB, the file the
maintainer's own deck already used) is kept byte-for-byte; nothing else in
this template reads it.

## B. The interaction cloud sets a genuine minimum legal supercell size

Unlike `B01_bccFe`/`B02_skyrmion2D`/`B03_skyrmion3D`, this case's ladder
cannot start near this project's usual ~4k-atom tier. `jfile.cut5`'s
lattice-translation offsets reach `n1,n2=+-6`, `n3=+-2` (measured directly
from the file). At a small enough replication, the same physical neighbour
is reached through two or more different periodic-boundary wraps of these
offsets -- not a crash, but a real, silent double/triple/...-counting of
that neighbour's exchange coupling (production output still "succeeds," it
is just physically wrong).

This was verified directly, not assumed: `struct.<simid>.out`'s own
`r_{ij}^x, r_{ij}^y, r_{ij}^z` columns give the real bond displacement
vector per listed neighbour, so duplicate physical neighbours inside one
atom's own listed neighbour set are directly countable (grouping atom 1's
`struct.<simid>.out` rows by rounded displacement vector, `do_prnstruct 1`,
`Nstep 2`, at each candidate `ncell n n n`):

| `n` | `natom` | neighbour-list entries (atom 1) | unique real neighbours | duplicate entries | max times one neighbour repeats |
| --- | --- | --- | --- | --- | --- |
| 4 | 256 | 1340 | 256 | 1084 | 7x |
| 6 | 864 | 1340 | 601 | 739 | 3x |
| 8 | 2048 | 1340 | 1013 | 327 | 3x |
| 10 | 5000 | 1340 | 1331 | 9 | 2x |
| 12 | 6912 | 1340 | 1340 | **0** | 1x |
| 13 | 8788 | 1340 | 1340 | 0 | 1x |
| 14 | 10976 | 1340 | 1340 | 0 | 1x |
| 16 | 16384 | 1340 | 1340 | 0 | 1x |

`n=12` is confirmed clean for every atom (spot-checked across all four
basis-atom types, not just atom 1: entries 6912/6912 atoms all show zero
duplicates at `n=12`). This is a real physical floor: below it, this case's
output is not merely "smaller," it is quietly wrong (an artificially
inflated effective exchange field on some neighbour pairs). Note that the
`directed_interactions`/`mean_neighbors` count itself does **not** change
with aliasing (each duplicate is still one line in `struct.<simid>.out`,
just pointing at a physical neighbour already listed elsewhere for the same
atom) -- naive atom/interaction-count constancy is not sufficient evidence
of correctness for this case, unlike `B01_bccFe`/`B02_skyrmion2D`/
`B03_skyrmion3D`, where it was.

**`n=16` is used as this ladder's floor** -- a comfortable margin above the
confirmed-clean `n=12`, and this project's already-established
smallest-tier convention (matches `B03_skyrmion3D`'s own floor).

## C. Size ladder

`Nx=Ny=Nz=n` cubes (four-atom dhcp basis, `natom=4*n^3`), same isotropic
methodology `B01_bccFe`/`B03_skyrmion3D` use, subject to the `n>=16`
physical floor above. Unlike those cases, the ladder's *top* is not the
usual ~4M-atom tier: per blueprint B04 section C ("stop when genuine
memory/resource limits are reached... do not require 4M atoms merely for
symmetry"), it is capped by real memory classification on this development
host (§F).

| `size_id` | `n` | `natom` | vs. nearest ~4x-ladder target | Note |
| --- | --- | --- | --- | --- |
| `16x16x16` | 16 | 16,384 | target 16,384, **exact** | floor (§B) |
| `20x20x20` | 20 | 32,000 | target 32,768, -2.3% | |
| `25x25x25` | 25 | 62,500 | target 65,536, -4.6% | |
| `32x32x32` | 32 | 131,072 | target 131,072, **exact** | maintainer's own production size (source example's baseline `ncell`) |
| `40x40x40` | 40 | 256,000 | target 262,144, -2.3% | |
| `50x50x50` | 50 | 500,000 | target 524,288, -4.6% | |
| `64x64x64` | 64 | 1,048,576 | target 1,048,576, **exact** | GPU DOUBLE becomes `unavailable_memory` here (§F) |
| `80x80x80` | 80 | 2,048,000 | target 2,097,152, -2.3% | ladder ceiling: both GPU precisions `unavailable_memory`, CPU still fits (§F) |

The next natural tier (`100x100x100`, ~4,000,000 atoms, this project's
usual top tier) is **not admitted**: §F's memory classification predicts
CPU itself becomes `unavailable_memory` there too on this development
host -- a genuine resource limit, not an arbitrary stopping point.

## D. Interaction characterization (blueprint B04 section B / master blueprint section 5)

Measured directly via `harness.cases.generate_run_directory` (not manual
file edits) with `extra_overrides={"do_prnstruct": 1, "Nstep": 20}` on
`build_cpu/bin/sd.f95`, reading `struct.<simid>.out` via this case's own
`neighbor_list_from_struct_output`:

| `size_id` | `natom` | directed interactions | mean neighbours | median neighbours | max neighbours | interaction-list memory footprint (real, `do_meminfo`) |
| --- | --- | --- | --- | --- | --- | --- |
| `16x16x16` | 16,384 | 21,921,792 | 1338.0 | 1338.0 | 1340 | 265,807,837 B (253.5 MiB) |
| `20x20x20` | 32,000 | 42,816,000 | 1338.0 | 1338.0 | 1340 | 517,975,005 B (494.0 MiB) |

Both measurements land exactly on `natom * 1338` -- the confirmed,
exact scaling law, matching the theoretical count from the raw
`jfile.cut5` template (5,352 bonds / 4 basis atoms = 1338.0) at every
tested size, with **zero** aliasing loss (§B). Median was cheap to obtain
here (not added to the shared `workload_metadata.py` parser, which still
returns `median_neighbors: null` for every case exactly as WP-02
established): `struct.<simid>.out`'s per-atom neighbour counts split
exactly 50/50 between the dhcp lattice's two crystallographic sites (8,192
atoms at 1340, 8,192 at 1336, at `16x16x16`), so median and mean coincide
exactly.

Larger sizes were not directly `do_prnstruct`-measured: a single
`do_prnstruct 1` dump at `natom=16,384` alone is a **2.4 GB**
`struct.<simid>.out` file (this case's very interaction density makes the
diagnostic dump itself expensive -- a real, disclosed cost, not
hand-waved), so tested sizes were kept to the two above, spanning a 1.95x
atom-count range with an exact match both times. Predicted
`directed_interactions` at the larger admitted sizes, from the confirmed
law (`natom * 1338`):

| `size_id` | `natom` | predicted directed interactions |
| --- | --- | --- |
| `25x25x25` | 62,500 | 83,625,000 |
| `32x32x32` | 131,072 | 175,374,336 |
| `40x40x40` | 256,000 | 342,528,000 |
| `50x50x50` | 500,000 | 669,000,000 |
| `64x64x64` | 1,048,576 | 1,402,994,688 |
| `80x80x80` | 2,048,000 | 2,740,224,000 |

## E. Sanity runs

All runs used `16x16x16` (16,384 atoms), `do_prnstruct 1`, `Nstep 50`;
CPU on `build_cpu/bin/sd.f95` (`UPPASD_GPU_BACKEND=OFF`,
`UPPASD_PRECISION=DOUBLE`), GPU on `build_gpu/bin/sd.f95.cuda`
(`UPPASD_GPU_BACKEND=CUDA`, `UPPASD_PRECISION=DOUBLE`,
`extra_overrides={"gpu_mode": 1}`).

| Run | Result |
| --- | --- |
| CPU, `dhcpNd_t0` | **PASS** -- completes; `natom=16384`, max neighbours 1340, matches §D |
| CUDA, `dhcpNd_t0` | **PASS** -- completes; identical `natom`/neighbour metadata |
| CPU, `dhcpNd_finite_t` | **PASS** -- completes, `\|m\|` conserved |
| CUDA, `dhcpNd_finite_t` | **PASS** -- completes, `\|m\|` conserved |

**`dhcpNd_t0` is a genuinely non-degenerate CPU/GPU comparison** -- a
materially different, stronger result than `B01_bccFe`/`B03_skyrmion3D`'s
own T=0 sanity checks, which are exact fixed points (uniform `Initmag 3`
states with zero net torque, so no real per-step dynamics happen at all).
Here `Initmag 1` gives a genuinely random, non-uniform initial state, so
the LLG update is doing real work every step:

* At `t=0` (`Iter 0`), CPU and GPU are **bit-identical to all 9 printed
  significant figures** in both `totenergy.<simid>.out`
  (`Exc=Tot=3.51431779E-05`, all other terms exactly `0`) and
  `averages.<simid>.out` (`<M>_x=1.58362729E-03`,
  `<M>_y=8.41882731E-05`, `<M>_z=-2.09660716E-03`,
  `<M>=2.62882572E-03`) -- confirming `Initmag 1`'s RNG-seeded random
  initial state is itself bit-reproducible across backends at the same
  `tseed`, not merely the same statistically.
* After the full 50-step relaxation, `restart.<simid>.out` shows the two
  backends have accumulated real, expected floating-point divergence from
  independently-executed per-step nonlinear dynamics: CPU
  `(M_x, M_y, M_z) = (4.09521494E-01, 4.40521385E-01, 7.98894896E-01)`
  vs. GPU `(4.09537684E-01, 4.40544966E-01, 7.98873593E-01)` -- agreeing to
  ~4 significant figures, diverging in the 4th/5th, exactly the FP
  non-associativity accumulation expected from real per-step dynamics on
  two independent code paths, not a bug. `\|Mom\|=2.66020000E-01` is
  identical on both (normalization is exact on both integrators). This is
  a **more informative** CPU/GPU LLG-update discriminator than the trivial
  fixed-point checks in `B01_bccFe`/`B03_skyrmion3D`, whose own READMEs
  note that limitation explicitly.
* **`restart.<simid>.out`'s iteration index differs by backend**
  (CPU reports `mstep=50`, the real final step; GPU reports `-1` on every
  row) and **`averages.<simid>.out`'s print cadence differs**
  (CPU prints only `Iter 0`, since `avrg_step 100 > Nstep 50`; GPU prints
  both `Iter 0` and `Iter 50`) -- both are the same pre-existing,
  general GPU backend cosmetics `B01_bccFe/README.md`/`B03_skyrmion3D/README.md`
  already document, re-observed here, not new findings specific to this
  case.
* Finite-T CPU and GPU necessarily diverge in their specific per-atom
  moment realizations (independent RNG streams per backend, same
  restriction `B01_bccFe/README.md` documents for its own
  `bcc_fe_finite_t`) -- only completion and `\|m\|` conservation were
  checked, consistent with `harness.gpu_sanity`'s restriction of strict
  moment comparison to `temperature==0`.

## F. Memory classification (blueprint B04 section D)

### GPU (`gpu_memory.py`, WP-06, reused unchanged)

`gpu_memory.classify_gpu_memory_availability` already implements exactly
the contract blueprint B04 section D asks for -- no harness change was
needed on the GPU side. Queried real device (this development host, 2x
RTX A4000): `device_memory_total_mib=16376`,
`device_memory_used_mib~=100` (idle).

### CPU (`harness/host_memory.py`, new -- this case is the first to need it)

No CPU-side equivalent existed before this WP: `B01_bccFe`/
`B02_skyrmion2D`/`B03_skyrmion3D`'s interaction lists (6-96
directed-interactions/atom) never come close to exhausting host memory at
any size on their own ladders, so nothing needed to classify CPU memory
availability. `B04_dhcpNd`'s ~1338 directed-interactions/atom changes
that -- its interaction-list footprint dominates process memory on the CPU
exactly as it does on the GPU (§D's own real `do_meminfo` measurements: 253
MiB at 16,384 atoms, growing to a predicted tens of GB by the top of this
ladder).

`harness/host_memory.py` mirrors `gpu_memory.py`'s
`classify_*_memory_availability`/`UNAVAILABLE_MEMORY` shape, but its
`BYTES_PER_DIRECTED_INTERACTION` constant (~12.07 bytes) is an **empirical
fit to two real measurements** (§D's `16x16x16`/`20x20x20` `do_meminfo`
peaks), not an architectural estimate the way `gpu_memory.py`'s GPU-side
constants necessarily are (no CUDA build prints its own real allocation
total) -- UppASD's CPU build already does, via `do_meminfo 1`
(`source/Tools/profiling.f90::memocc`), which happens to already be this
case's own maintainer-supplied baseline. `query_host_memory_mib` reads
`/proc/meminfo`'s `MemAvailable` (reclaimable-cache-aware, not naive
`MemFree`) for the real headroom side of the comparison. `detect_host_oom`
recognizes UppASD's own allocation-failure message
(`"problem of allocation of array"`, `memocc`'s real `stop` path) alongside
SIGKILL/`std::bad_alloc`/OS out-of-memory signatures, matching
`gpu_memory.detect_gpu_oom`'s evidence-based discipline.

### Classification at every ladder tier and beyond

Real queried host memory at classification time:
`MemAvailable~=49,015 MiB` (~47.9 GiB; `MemTotal~=64,096 MiB`, ~62.6 GiB
-- this development host is otherwise idle at measurement time, but see
the shared-host-contention memory on GPU availability fluctuating with
other users). `safety_margin=1.5` throughout (both modules' shared
default).

| `n` | `natom` | directed interactions (law) | CPU | GPU SINGLE | GPU DOUBLE |
| --- | --- | --- | --- | --- | --- |
| 16 | 16,384 | 21,921,792 | available | available | available |
| 20 | 32,000 | 42,816,000 | available | available | available |
| 25 | 62,500 | 83,625,000 | available | available | available |
| 32 | 131,072 | 175,374,336 | available | available | available |
| 40 | 256,000 | 342,528,000 | available | available | available |
| 50 | 500,000 | 669,000,000 | available | available | available |
| **64** | **1,048,576** | **1,402,994,688** | available | available | **`unavailable_memory`** |
| **80** | **2,048,000** | **2,740,224,000** | available | **`unavailable_memory`** | **`unavailable_memory`** |
| 100 (not admitted) | 4,000,000 | 5,352,000,000 | **`unavailable_memory`** | `unavailable_memory` | `unavailable_memory` |

This is exactly the "atom-count versus interaction-count crossover" and
"backend-dependent memory wall" characterization blueprint B04 section E
asks for: `natom` alone predicts nothing about where either backend's wall
is -- `directed_interactions` (via the confirmed `natom*1338` law) does.
GPU DOUBLE is the first to fail (`64x64x64`), then GPU SINGLE
(`80x80x80`), then, at the deliberately-unadmitted `100x100x100` tier,
CPU itself -- confirming blueprint B04 section C's "stop when genuine
memory/resource limits are reached... do not require 4M atoms merely for
symmetry" is a real constraint here, not a stylistic choice: the 4M-atom
tier every other admitted case's ladder reaches is genuinely unavailable
on every backend on this development host.

A future campaign attempting `64x64x64`/`80x80x80` GPU DOUBLE/SINGLE
samples should call `gpu_memory.classify_gpu_memory_availability` (as any
GPU campaign already must, per WP-06) and record a
`build_unavailable_memory_record` rather than attempt the run; the same
now applies to a hypothetical CPU sample at `100x100x100` or beyond via
`host_memory.classify_host_memory_availability`, never as a `FAILED`/code
error.

## Checklist

- [x] Template audited (§A: maintainer-supplied production deck
  `examples/SimpleSystems/dhcpNd`, copied verbatim except one required
  override-mechanism-correctness fix -- a real duplicate-`gpu_mode`-keyword
  bug found in the maintainer's own deck, documented and fixed).
- [x] Full intended interaction cloud retained (§A/§D: exchange cutoff,
  shell count and basis all unchanged from the maintainer deck; 5,352
  directed template bonds, confirmed present with zero loss at every
  aliasing-free tested size).
- [x] No benchmark-only truncation (§A: no cutoff reduction, no shell
  truncation, no basis simplification anywhere in the fix applied).
- [x] Ndirected measured (§D: real `struct.<simid>.out` counts at two
  sizes, exact `natom*1338` law, predicted at larger sizes from that law).
- [x] Mean/max neighbour count measured (§D: mean 1338.0, max 1340, both
  real).
- [x] Median neighbour count recorded where inexpensive (§D: 1338.0,
  exact 50/50 two-site split, computed directly from real `struct.<simid>.out`).
- [x] Memory footprint recorded where available (§D/§F: two real
  `do_meminfo` peak measurements; `harness/host_memory.py`'s empirical fit
  and the GPU/CPU classification table for every ladder tier).
- [x] Size ladder established (§C: `16x16x16` -> `80x80x80`, eight sizes,
  bounded below by a real measured periodic-aliasing floor (§B) and above
  by real memory classification (§F) -- not the usual ~4x/4M-atom
  convention, per blueprint B04 section C).
- [x] CPU sanity run passes (§E).
- [x] CUDA sanity run passes where memory allows (§E; §F documents exactly
  where it stops allowing it for larger sizes).
- [x] Memory-limit classification works (§F: `gpu_memory.py` reused
  unchanged for GPU; new `harness/host_memory.py`, empirically grounded in
  this case's own real `do_meminfo` data, for CPU -- the first case whose
  ladder needs a CPU-side classification at all).
- [x] Documentation emphasizes interaction-normalized analysis (§D/§F
  throughout: every crossover and every memory-classification result is
  presented against `directed_interactions`, not `natom` alone).
