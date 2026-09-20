# RCG-09D: Remove quadratic adaptive coarse-graining setup algorithms — evidence

**Status:** Implemented, equivalence-proven, benchmarked. Both negative controls
executed and both fail clearly on the previous algorithms.
**Date:** 2026-08-15
**Base commit:** `ad6a7a0a` (RCG-09C)
**Backend/precision:** CPU only, DOUBLE precision, GNU Fortran 13.3.0,
Release (`-O3 -Ofast -funroll-loops -finline-functions -march=native`), host
`alcazar` (11th Gen Intel Core i9-11900, 2-core cgroup-limited environment).
No GPU device is involved: both constructions are host setup code.

---

## 1. Session evidence template

```text
Task: RCG-09D -- Remove quadratic adaptive coarse-graining setup algorithms
Commit: this session's patch on ad6a7a0a
Findings addressed:
  (A) unique-bond construction scanned the already-discovered canonical-pair
      list once per directed contribution -> O(E^2);
  (B) Hamiltonian macro-block connectivity cleared or scanned an
      nblocks-length array once per source block, in five separate passes
      -> O(B^2) even for sparse physical block connectivity.
Files intentionally changed:
  source/CoarseGraining/adaptivecgproduction.f90     (build_unique_bonds ->
                                                      fold_directed_bonds)
  source/CoarseGraining/hamiltonianmacroblocks.f90   (build_macroblock_layout)
  tests/coarse_graining/test_adaptive_setup_scaling.f90  (new fixture)
  CMakeLists.txt  (register the new fixture in the cg13/cg13-cpu suite)
Equivalence: SS3 -- 54 structure dumps from 55 fixture runs, bitwise identical
  between a binary built from the pre-change sources and one built from the
  post-change sources.
Negative controls: SS5 -- both scaling assertions fail on the transcribed
  previous algorithms (172.5x and 274.6x growth against a 48x threshold).
Focused test command/result: SS4 -- ctest -L cg13-cpu 29/29 on a fresh
  out-of-tree build; ctest -R '^asd-tests$' 1/1 (feature-off regression).
Performance raw artifact: SS6.
Checklist: SS7.
```

---

## 2. What changed

### 2.1 Unique bond construction (A)

`build_unique_bonds` previously discovered canonical bonds by scanning every
canonical pair found so far for each accepted directed contribution:

```fortran
do found = 1, bond
   if (pair_i(found) == min(atom,target) .and. pair_j(found) == max(atom,target)) exit
end do
```

with the identical scan repeated in the DMI pass. For `E` accepted directed
entries this is `O(E^2)`.

The replacement, `fold_directed_bonds`, is record generation + grouping:

1. emit one record per accepted directed contribution, carrying the canonical
   key `(min(i,j), max(i,j))` and the sparse-list slot it came from — exchange
   records first, then DMI records, both in the production `atom`-major /
   `neighbour`-minor traversal order, so a record's emission index *is* its
   position in the original discovery sequence;
2. two stable counting sorts, by the high atom index then by the low one, group
   records by canonical key in `O(E + natom)`;
3. number the resulting groups by ascending first-appearance emission index;
4. fold each group in emission order, orienting the DMI contribution by
   `merge(+1,-1,atom < target)` exactly as before.

`build_unique_bonds` is now a thin production wrapper that supplies
`ham%nlist`/`ham%dmlist`/`mmom` and fills `bond_displacement_m` from
`wrapped_displacement`, which was already a pure function of the canonical
pair.

Because the sorts are stable, a group's records stay in emission order, so
step 4 replays the *same* sequence of floating-point additions per bond that
the discovered-list scan performed. Bond numbering, canonical atom pairs,
folded matrices and displacements are therefore reproduced bitwise (§3), not
merely to tolerance.

**Periodic-image key semantics are unchanged.** RCG-09A.1 established that the
production sparse lists carry target atom indices but no periodic-image
translation, so the atom index pair is the entire available key and a canonical
pair carrying more than one directed entry per direction is an alias that must
be rejected rather than merged. That rejection is still
`validate_adaptive_folded_pair_contract`, called with the same per-bond
forward/reverse counts and accumulators; it is now called from inside
`fold_directed_bonds`. The reciprocal convention (`J_ji=J_ij`, `D_ji=-D_ij`,
`K^T=-K`) and the DMI folding orientation are untouched.

### 2.2 Macro-block connectivity (B)

`build_macroblock_layout` ran five passes over the blocks, each of which
performed a full `nblocks`-length clear (`seen_target = .false.`,
`target_group = 0`) or a full `do bj = 1, nblocks` scan *per source block*.

The replacement uses one generation-stamped marker array plus, where ordering
matters, one counting sort:

- **discovery** (two passes: count, then record): a pair `(bi,bj)` is new when
  `target_stamp(bj) /= generation`; the stamp is the monotone per-source-block
  generation counter, so nothing is ever cleared;
- **CSR ordering**: the previous code obtained ascending destination blocks per
  source from the `do bj = 1, nblocks` scan. The pairs are now emitted in
  neighbour-entry order and placed by one counting sort on the destination
  block replayed into the per-source CSR cursors — `O(n_pairs + nblocks)` in
  total, not per source block;
- **the three subsequent lookup passes** stamp only their own source block's
  CSR pair range and use `target_stamp(bj) == generation` as the validity
  condition, so each is `O(entries + touched targets)`.

Every output array keeps its previous definition and ordering; the layout type
is unchanged.

### 2.3 Old implementations deleted

Neither quadratic implementation is retained in `source/`. There is no flag,
fallback or dormant copy:

```console
$ grep -n "seen_target\|do found = 1, bond\|target_group = 0" source/CoarseGraining/*.f90
(no matches)
```

The previous algorithms survive only as brute-force *references inside the test
fixture* (`tests/coarse_graining/test_adaptive_setup_scaling.f90`), which is the
pattern RCG-07 already established for `rebuild_static_hybrid_ownership`.

---

## 3. Old/new equivalence on real systems (C)

Equivalence was demonstrated **before** deletion, binary against binary.

Method: temporary instrumentation dumped every field of the built structures —
`bond_atom`, `bond_displacement_m` and `bond_matrix_j` as raw IEEE-754 bit
patterns (`transfer(x,0_int64)`, `z16.16`), and all sixteen macro-block layout
arrays as integers. Two executables were built with identical flags, one from
the pre-change sources and one from the post-change sources, and every fixture
was run under both. The instrumentation was removed afterwards and is not part
of the commit.

Corpus: all 50 `tests/coarse_graining/e2e/` fixture directories plus five
purpose-built macrocell probes (`do_macro_cells Y`), since no checked-in
fixture enables the macrocell decomposition that `build_macroblock_layout`
requires. The probes cover `block_size 2` on `ncell 12 8 4`, `block_size 1` on
`10 6 4`, anisotropic `block_size_x/y/z = 3/2/1` on `12 6 4`, an indivisible
`ncell 7 5 3` with `block_size 2` (incomplete edge blocks), and a quasi-2D
`16 16 1` with `block_size 4`.

Result:

```console
$ diff -r dumps/old dumps/new && echo IDENTICAL
IDENTICAL
```

54 dumps compared, byte for byte:

| Structure | Dumps | Range covered |
| --- | ---: | --- |
| Unique bonds (`.bonds`) | 49 | 600 to 24 000 canonical bonds; exchange-only and exchange+DMI fixtures |
| Macro-block layout (`.macro`) | 5 | 16 to 240 blocks, 80 to 12 000 block pairs, 9 660 to 38 400 Hamiltonian entries, 2 to 32 atoms per block |

This covers the contract items required by the task: same canonical bond count,
same key and bond order, same folded matrices, same block adjacency. Coarse
couplings are downstream consumers of exactly these structures and are
unchanged as a consequence; the `cg13-cpu` moving-parity and production e2e
fixtures (§4) exercise them end to end.

The permanent fixture additionally checks the constructions against the
transcribed previous algorithms on synthetic systems, including canonical bonds
that only the DMI pass discovers (the exchange and DMI neighbour shells differ
deliberately) and both rejection paths. Its matrix comparison is to roundoff
rather than bitwise — the reference is a separate translation unit and this
project builds Release with `-Ofast`, so the two compilations of the same
multiply chain legitimately differ by an ulp (observed 0.7–1.0 ulp). The
bitwise claim above rests on the binary-against-binary comparison, where both
sides are the same compilation of the same expressions.

---

## 4. Tests executed

Fresh out-of-tree build `build_rcg09d_clean`, Release, DOUBLE, GNU 13.3.0,
`cmake --build . -j8` with no new warnings from the changed files.

```console
$ ctest -L cg13-cpu
100% tests passed, 0 tests failed out of 29     (was 28; +coarse-graining-adaptive-setup-scaling)
$ ctest -R '^asd-tests$'
100% tests passed, 0 tests failed out of 1
```

Runtime physics tests are unchanged: the six `moving-parity` fixtures, the
production e2e case, the setup rejection matrix, the ownership comparator and
the transition-ownership invariants all pass without modification. No tolerance
was altered.

The new fixture `coarse-graining-adaptive-setup-scaling` reports:

```text
exchange+DMI ring(11):max relative matrix difference =  1.62836E-16
exchange-only ring(11):max relative matrix difference =  2.18899E-16
exchange+DMI ring(64):max relative matrix difference =  9.09004E-17
exchange-only ring(64):max relative matrix difference =  1.65940E-16
exchange+DMI ring(257):max relative matrix difference =  7.19688E-17
exchange-only ring(257):max relative matrix difference =  1.80784E-16
RCG-09D unique-bond fold scaling (atoms, directed entries, min-of-trials seconds):
  atoms=4096 directed=32768 seconds= 1.35577E-03
  atoms=16384 directed=131072 seconds= 5.56449E-03
  atoms=65536 directed=524288 seconds= 2.37805E-02
  growth over a 16x directed-entry increase:  1.75401E+01
RCG-09D macroblock layout scaling (blocks, atoms, min-of-trials seconds):
  blocks=1024 atoms=2048 seconds= 1.27473E-04
  blocks=4096 atoms=8192 seconds= 5.11626E-04
  blocks=16384 atoms=32768 seconds= 2.14364E-03
  growth over a 16x block-count increase:  1.68164E+01
```

---

## 5. Negative controls

Both are genuine: the fixture's own transcribed previous algorithm was
substituted for the production call under test, and both assertions failed.

| Control | Sizes | Growth over 16x | Threshold | Result |
| --- | --- | ---: | ---: | --- |
| Unique-bond fold uses the discovered-list scan | 512 / 2048 / 8192 atoms | 172.5x | 48x | **FAIL** (control works) |
| Macro-block layout uses the dense clear/scan | 1024 / 4096 / 16384 blocks | 274.6x | 48x | **FAIL** (control works) |

Raw output of the injected runs:

```text
RCG-09D unique-bond fold scaling (atoms, directed entries, min-of-trials seconds):
  atoms=512 directed=4096 seconds= 1.04795E-03
  atoms=2048 directed=16384 seconds= 1.01205E-02
  atoms=8192 directed=65536 seconds= 1.80808E-01
  growth over a 16x directed-entry increase:  1.72535E+02
FAIL: unique-bond fold growth over 16x directed entries stays far below the quadratic prediction
RCG-09D macroblock layout scaling (blocks, atoms, min-of-trials seconds):
  blocks=1024 atoms=2048 seconds= 3.86523E-04
  blocks=4096 atoms=8192 seconds= 4.62319E-03
  blocks=16384 atoms=32768 seconds= 1.06145E-01
  growth over a 16x block-count increase:  2.74615E+02
adaptive setup scaling tests failed: 2
```

The 48x threshold sits between the linear prediction (~16x) and the quadratic
one (~256x), so it discriminates the two algorithms without being
timing-flaky. Note the absolute cost as well: the previous unique-bond
construction needed 181 ms at 8 192 atoms, while the new one needs 24 ms at
65 536 atoms.

---

## 6. Setup scaling benchmark (D)

Single-threaded (`OMP_NUM_THREADS=1`), minimum of three trials per point,
fixture construction excluded from the timed region. "old" is the transcribed
previous algorithm, run only over the range where it remains affordable.

### 6.1 Unique-bond construction vs. atom count and directed bond count

Periodic ring, four exchange and four DMI neighbours per atom, so directed
entries `E = 8 N` and canonical bonds `= 6 N`.

| Atoms | Directed entries | Canonical bonds | new (s) | new ns/entry | old (s) | old ns/entry |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 024 | 8 192 | 6 144 | 4.327e-4 | 52.8 | 1.986e-3 | 242.5 |
| 2 048 | 16 384 | 12 288 | 6.367e-4 | 38.9 | 7.722e-3 | 471.3 |
| 4 096 | 32 768 | 24 576 | 1.415e-3 | 43.2 | 3.476e-2 | 1 060.7 |
| 8 192 | 65 536 | 49 152 | 2.652e-3 | 40.5 | 1.469e-1 | 2 242.2 |
| 16 384 | 131 072 | 98 304 | 5.306e-3 | 40.5 | 5.938e-1 | 4 530.2 |
| 32 768 | 262 144 | 196 608 | 1.074e-2 | 41.0 | — | — |
| 65 536 | 524 288 | 393 216 | 2.370e-2 | 45.2 | — | — |
| 131 072 | 1 048 576 | 786 432 | 4.931e-2 | 47.0 | — | — |

Cost per directed entry is flat at 39–53 ns across a 128x range. The previous
implementation's cost per entry doubles with every doubling of the system —
242 → 471 → 1 061 → 2 242 → 4 530 ns — which is the quadratic trend, and it is
removed.

### 6.2 Macro-block connectivity vs. block count

Periodic chain of blocks, 2 atoms per block, four-neighbour atomistic ring, so
block connectivity is sparse (3 destination blocks per source) and the pair
count grows linearly.

| Blocks | Atoms | new (s) | new ns/block | old (s) | old ns/block |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 256 | 512 | 3.862e-5 | 150.9 | 5.003e-5 | 195.4 |
| 512 | 1 024 | 6.680e-5 | 130.5 | 1.319e-4 | 257.6 |
| 1 024 | 2 048 | 1.324e-4 | 129.3 | 3.747e-4 | 365.9 |
| 2 048 | 4 096 | 2.562e-4 | 125.1 | 1.229e-3 | 600.3 |
| 4 096 | 8 192 | 5.152e-4 | 125.8 | 4.445e-3 | 1 085.1 |
| 8 192 | 16 384 | 1.029e-3 | 125.7 | 1.672e-2 | 2 040.8 |
| 16 384 | 32 768 | 2.147e-3 | 131.1 | — | — |
| 32 768 | 65 536 | 4.525e-3 | 138.1 | — | — |

Cost per block is flat at 125–151 ns across a 128x range; the previous
implementation's per-block cost doubles with every doubling of the block count.
At 8 192 blocks the new construction is already 16x faster, and the gap grows
without bound.

Both are setup costs, paid once, not steady-state physics costs.

---

## 7. Checklist

- [x] Unique-bond construction no longer performs pairwise discovered-list
      scans (§2.1, §2.3).
- [x] Periodic-image key semantics match RCG-09A.1 — the atom index pair is
      the whole key and multi-entry canonical pairs are still rejected with
      `periodic-image alias` (§2.1; fixture case
      `test_unique_bonds_reject_periodic_image_alias`).
- [x] DMI folding convention unchanged — same `cross_product_matrix`, same
      `merge(+1,-1,atom<target)` orientation, same reciprocal validator; proven
      by bitwise identity of `bond_matrix_j` on the DMI fixtures (§3).
- [x] Macro-block construction no longer performs avoidable `O(B^2)` scans
      (§2.2, §2.3).
- [x] Output is deterministic — no hashing and no unordered container was
      introduced; both constructions use stable counting sorts and stamped
      markers, and repeated runs of the whole corpus reproduce identical dumps
      (§3).
- [x] Old/new structure equivalence demonstrated before deletion (§3).
- [x] Old quadratic code deleted afterward (§2.3).
- [x] Setup scaling benchmark recorded (§6).
- [x] Runtime physics tests unchanged (§4).

---

## 8. Open limitations

1. **CPU only.** Both constructions are host setup code with no GPU
   counterpart, so there is no CUDA/HIP semantic divergence to align. No claim
   is made about GPU behaviour and no GPU suite was run in this session.
2. **Macro-block coverage is probe-driven.** No checked-in e2e fixture enables
   `do_macro_cells`, so the five macrocell probes used for the equivalence
   comparison (§3) are scratch inputs and are not added to the repository. The
   permanent fixture drives `build_macroblock_layout` directly instead, at six
   block-count/atoms-per-block combinations, against the dense reference.
3. **No independent reviewer.** As with the preceding RCG tasks, this evidence
   has not been reviewed by a second party.
4. The `mem_large_host` fixture contributes the largest bond dump (24 000
   canonical bonds); no larger *real* input was available on this host, so the
   million-entry end of §6.1 is synthetic.
