# ANI-NONUNIFORM-REJECT (RCG-03)

**Purpose:** production-executable evidence that
`build_production_anisotropy`'s cell-periodicity check
(`source/CoarseGraining/adaptivecgproduction.f90`) rejects spatially varying
anisotropy before integration, reached through the ordinary `sd.f95`
executable and an ordinary capability (`do_cluster`), not through direct
internal staging or source-level fault injection.

**Mechanism:** the host is the same 48-atom, `NA=2`, `ncell 6 2 2` system
used by the other `anisotropy_uniform_*` fixtures, with a uniform kfile
(`../kfile_cg_x`: basis 1 `k1=-0.002`, basis 2 `k1=-0.003`). A single-atom
`do_cluster` embedding (`ncell_clus 1 1 1`, `posfile_clus` at fractional
`(0,0,0)`) is placed exactly on top of host atom 1's own lattice site, so
`clus_expand` (`source/Clusters/clus_geometry.f90`) stays `0` and `Natom`
is unchanged (`48 = NA*N1*N2*N3`) -- the pre-existing
`Natom /= NA*N1*N2*N3` geometry rejection in `validate_configuration`
therefore never fires, unlike the `clus_expand>0` case investigated and
rejected as a candidate trigger during the original RCG-03 session (see
`docs/RCG-03_POLARIZATION_ANISOTROPY_EVIDENCE.md`). `anisotropy_clus`
(`kfile_clus`) gives that one embedded atom `k1=-0.005`: a basis-1 atom
whose anisotropy now genuinely diverges from every other basis-1 atom in
the lattice.

**Why `set_landeg 1` and `momfile_clus`'s gyromagnetic column matter:**
`do_cluster` embedding is gated behind two other pre-existing,
unrelated-to-anisotropy setup checks that a naive cluster fixture trips
first:
- `Landeg_glob` defaults to `2.0`, but the ordinary host moment path halves
  it (`Landeg(i)=Landeg_ch(...)*0.5`, `source/System/magnetizationinit.f90`)
  while `cluster_creation` (`source/Clusters/clusters.f90`) copies
  `Landeg_ch_clus` into `Landeg(iatom)` directly, with no such factor. With
  `set_landeg 0` (the default), the cluster atom would silently pick up the
  raw `2.0` default while every host atom is `1.0`, tripping the unrelated
  `Landeg/do_site_damping: ... requires uniform gamma and damping` guard
  before the anisotropy check is ever reached.
- `set_landeg 1` here makes both host and cluster momfiles' gyromagnetic
  column authoritative; `momfile_clus` therefore spells out `1.0` (not
  `2.0`) so the embedded atom's post-embedding `Landeg` matches the host's
  `1.0` exactly, given `cluster_creation` does not apply the host's `*0.5`
  convention. This is a genuine, narrow pre-existing inconsistency in the
  `do_cluster` embedding path (unrelated to RCG-03's capability-safety
  scope) that this fixture works around entirely from its own input files;
  no source change was needed or made to reach it.

**Exchange/moment plumbing:** `exchange_clus` (`jfile_clus`) is required by
`cluster_creation`'s unconditional exchange-overwrite loop even though this
fixture does not care about cluster-internal exchange; it supplies one
self-referencing bond with zero coupling so the host's real exchange is
left untouched. `momfile_clus` otherwise matches the host's basis-1
momfile line exactly (moment `2.23`, direction `+x`).

**Units and sign convention:** anisotropy values are in the same mRy
uniaxial-`k1`/`k2` convention as the ordinary (non-cluster) `kfile`
(see `source/Hamiltonian/hamiltonianinit.f90:setup_anisotropies`, which
`hamiltonianinit.f90` also calls for the cluster path with
`anisotropy_clus`/`anisotropytype_clus` as input, correctly producing
`ham_clus%kaniso_clus` etc.).

**Negative control:** replacing `kfile_clus`'s `k1=-0.005` with the host's
own uniform `k1=-0.002` (i.e. an anisotropy-matched cluster embedding)
makes the run complete normally (`AdaptiveCG: capability accepted`,
`returncode=0`) instead of rejecting -- confirming the rejection is
produced by the genuine per-atom divergence, not merely by `do_cluster`
being present. This was verified manually during construction and is
exercised as the accept path by every other `do_cluster`-free adaptive-CG
fixture in this suite.

**Expected failure mode before the RCG-03 fix:** without the
cell-periodicity check, `build_production_anisotropy` would silently
sample anisotropy from one central translated copy per basis index
(`countstart+basis`) and broadcast it to every block, homogenizing this
atom's divergent `k1=-0.005` away rather than rejecting it.

**Wired into:** `tests/coarse_graining/run_setup_rejection_matrix.py`
(`ctest` label `cg13-cpu`/`rejection`, target
`adaptive-cg-setup-rejection-matrix`), copied verbatim (not text-mutated
like the other rejection cases, since it needs its own
`posfile_clus`/`momfile_clus`/`exchange_clus`/`anisotropy_clus` files).
