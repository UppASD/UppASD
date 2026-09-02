# CPU-HAM-03A — Sparse Hamiltonian audit

Date: 2026-09-02
Repository: `UppASD_gpu_hip_cu`
Branch: `gpu_hip_cu_ab_cg`

## Outcome

Historical sparse support is **PARTIALLY_CONNECTED**. The `do_sparse` input,
setup calls, CPU spin-dynamics dispatches, CSR-like arrays, and a portable
OpenMP field apply are still present. They are not, however, the canonical
`HamiltonianActions` implementation, do not provide the canonical field/energy
contract, and are not a validated current sparse backend.

The recommendation for CPU-HAM-03B is **GO**, with a constrained redesign:

- implement a new persistent scalar-isotropic-`J` backend behind a provider
  boundary;
- use one `Natom × Natom` CSR operator and a dense three-column spin RHS;
- retain a non-MKL provider so MKL is optional;
- validate the new field against canonical `HamiltonianActions` before wiring
  it into production dispatch;
- defer DMI, tensor, anisotropy, and other terms until the scalar-J backend is
  proven.

This is an audit only. No production sparse implementation is added here.

## A. Current-source trace

### Input and setup

`do_sparse` is declared in `source/Input/inputdata.f90:147`, defaults to `N` at
`:469`, and is parsed from legacy input files by
`source/Input/inputhandler.f90:1173-1175`. Generated f90wrap accessors also
expose it in `source/f90wrap/f90wrap_inputdata.f90:1659-1673`.

During `setup_simulation`, `source/uppasd.f90:1212-1227` constructs one of
three historical representations:

| Input condition | Setup routine | Representation |
|---|---|---|
| `do_jtensor == 1` | `setupSparseTensor` | 3×3 block per neighbor |
| otherwise with `do_dm == 1` | `setupSparseBlock` | scalar `J` plus an antisymmetric DMI block |
| scalar exchange | `setupSparseScalar` | scalar CSR-like rows |

The module is compiled unconditionally by `source/System/CMakeLists.txt:12`.
Its module-level arrays are allocated during setup and are intended to remain
available for later field applications. There is no production call to
`allocateSparse(...,flag=-1)` and no sparse cleanup routine. Repeated setup
would also allocate the arrays again without an `allocated` guard.

### Runtime dispatch

The CPU SD driver imports `Sparse` and dispatches to the sparse routines in the
initial phase and field-print path (`source/sd_driver.f90:137-152` and
`:357-370`), in the ordinary SD measurement phase (`:528-541` and
`:581-599`), and in `sd_minimal` (`:1192-1207` and `:1222-1235`). The sparse
branch is outside an enclosing OpenMP parallel region; the portable sparse
routine owns its `parallel do` over target atoms for each ensemble.

`ms_driver.f90` imports `Sparse`, but no `do_sparse` conditional or sparse apply
call was found there. The other current drivers use `HamiltonianActions` and
do not dispatch to `effective_field_Sparse*`. Adaptive coarse graining rejects
`do_sparse` during setup (`source/CoarseGraining/adaptivecgproduction.f90:949-950`).

The normal current CPU path is `HamiltonianActions::effective_field`, whose
pair field is assembled by `heisenberg_field` and whose optional energy is
derived in the same per-atom implementation (`source/Hamiltonian/hamiltonianactions.f90:61-224`
and `:226-304`). Sparse is therefore a live legacy escape hatch in CPU SD, but
not the current Hamiltonian owner.

### Sparse data and apply behavior

`source/System/sparse.f90` stores:

- `Hscalar` for scalar values;
- `Hblock(3,3,*)` for block values;
- `Hcolumns` plus `HpointerB/HpointerE` for zero-based CSR-like rows.

`setupSparseScalar` fills values from `ncoup(j,i_ham,1)` and neighbor columns
from `nlist(j,i_atom)` (`:37-81`). It allocates
`max_no_neigh*Natom` slots rather than the exact nonzero count and allocates
both scalar and block storage even when only one is used (`:318-360`).
`conf_num` is accepted but the setup always selects configuration 1.

`effective_field_SparseScalar` applies the same scalar matrix to the three
Cartesian components using the component dimension as three dense RHS columns
(`:270-315`). This is conceptually the desired scalar-J operation, although
the implementation repeats the CSR traversal for every ensemble and the
portable path uses an array slice update in its inner loop.

`effective_field_SparseBlock` performs a 3×3 block multiply
(`:220-268`). `setupSparseBlock` combines diagonal `J` and an antisymmetric DMI
block, but does not use `dmlist` or `dmlistsize`; it assumes the DMI and exchange
lists agree and only warns when their maximum widths differ (`:122-218`). This
is not sufficient evidence for current DMI correctness. `setupSparseTensor`
also exists, but has no dedicated current-path oracle.

## B. Historical implementation and version history

### Located implementation

The historical implementation is `source/System/sparse.f90`, first introduced
in repository history by `67f2772` (2019-07-13), modified during the source
update at `a7432e9` (2022-10-08), and touched by `3b94ed8` (2024-02-07,
“Fix for MKL sparse routines”). The corresponding setup and dispatch wiring is
in `source/uppasd.f90` and `source/sd_driver.f90`.

The repository documentation still describes this as a historical feature:
`docs/FEATURES.md:129` says sparse effective-field support is for Jij/Dij and
hardcoded to MKL. The changelog records the same feature and its later
preprocessing change (`CHANGELOG.md:75` and `:164`). No current sparse-specific
benchmark or correctness test was found.

### Historical representations

- Scalar exchange: one scalar value per directed neighbor and a row/column
  index list.
- DMI plus exchange: a 3×3 block whose diagonal is `J` and off-diagonal part
  encodes the DMI cross-product matrix.
- Tensor exchange: one 3×3 block per directed neighbor, populated from
  `j_tens`.

The representation is a directed neighbor gather, not a proven symmetric
matrix. The scalar setup contains an unused commented alternative assignment
and the MKL descriptor is marked symmetric/lower even though the populated
rows are directed. That descriptor must not be copied into a new backend
without proving reciprocity and storage conventions.

### MKL API and obsolescence

The source conditionally includes `mkl_spblas.f90` and declares
`SPARSE_MATRIX_T`/`MATRIX_DESCR` (`sparse.f90:1-28`). Scalar setup calls
`MKL_SPARSE_D_CREATE_CSR` (`:72-78`), but the actual scalar apply calls the
legacy `mkl_dcsrmm` entry point (`:294-297`). Block apply calls legacy
`mkl_dbsrmv` (`:243-245`). The created `H_sparse` handle is not used by either
apply routine, and no destroy/release call is present.

The current CMake configuration explicitly comments that these legacy entry
points are not exported by current oneMKL releases (`CMakeLists.txt:584-592`).
The old automatic `-D__INTEL_MKL__` definition was removed by commit `6b01f8b`
(2026-07-18); current CMake instead defines
`UPPASD_MKL_LEGACY_SPARSE=0` (`CMakeLists.txt:280-285`), but the sparse source
does not consume that macro. Consequently, the `__INTEL_MKL__` branches are
not enabled by the current CMake build, even when `USE_MKL=ON`, unless a user
manually supplies the old preprocessor define. That manual route would select
the obsolete calls and is not a supported validation path.

### Threading and lifecycle

The historical non-MKL path loops over ensembles sequentially and executes an
OpenMP `parallel do` over target rows. The MKL path, if manually enabled,
would call MKL once per ensemble from the caller and would rely on MKL's own
threading. There is no explicit thread policy, no inspector state, and no
persistent provider handle used in the hot path. A new backend must establish
one owner for parallelism and avoid nested OpenMP/MKL threading.

## C. Current architecture fit

### Minimal scalar-J design

The CPU-HAM-03B design should be:

```text
CSR J:              Natom × Natom, directed rows
spin RHS X:         Natom × 3, [Mx My Mz]
field result Y:     Natom × 3, Y = J X
```

The production storage should use the exact directed nonzero count, preserve
the canonical physical atom ordering, and retain the `aham`/neighbor mapping
used to build rows. It must be constructed once after Hamiltonian setup and
before the timestep loop. CSR structure, provider handles/inspector state,
component packing, and work buffers must be persistent; no allocation or
structure analysis may occur during a timestep.

The initial backend must support only scalar isotropic `J` with the same
eligibility conditions as the canonical direct path. It must reject or fall
back for DMI, tensor exchange, rescaled exchange, LSF/configuration-dependent
couplings, and unsupported onsite or higher-order terms until each has an
explicit oracle.

### Three-RHS assessment

The three-RHS form is appropriate and should be preferred over a redundant
`3N × 3N` matrix. It removes three copies of sparse metadata and allows a
provider's SpMM operation, when efficient, to reuse the CSR traversal. A
portable fallback may process the three components in one row loop if that is
faster on the target compiler; that is an implementation measurement, not a
physics distinction.

The output must be mapped back to UppASD's existing
`beff(3,Natom,Mensemble)` layout without changing physical atom IDs, restart
ordering, or ensemble semantics.

### Energy semantics

The sparse apply routine must remain field-only and must not invent a separate
Hamiltonian energy definition. For the eligible reciprocal bilinear pair term,
the validation and optional energy path should use:

```text
E_pair = -0.5 * sum_i dot(m_i, B_pair_i)
```

with the same unit conversion used by `HamiltonianActions`. Zeeman, anisotropy,
DMI, tensor, and nonlinear terms must keep their term-specific canonical
semantics. The current sparse routines have no energy argument or reduction;
`calc_energy` separately uses the canonical Hamiltonian routines. This is a
compatibility boundary, not evidence that sparse field and sparse energy are
already equivalent.

## D. Build dependency policy

`USE_MKL` is optional in CMake. With `USE_MKL=OFF`, the project selects required
BLAS/LAPACK providers and the current sparse module uses its portable OpenMP
fallback. With `USE_MKL=ON`, CMake finds and links the optional `MKL::MKL`
target, but current configuration still does not enable the legacy sparse
branches. No sparse backend work should make MKL mandatory.

No alternative sparse library is currently selected or wrapped by the build
system. CPU-HAM-03B should therefore provide either:

1. a provider-neutral CSR implementation as the always-available baseline,
   with an optional oneMKL inspector-executor provider; or
2. an explicitly documented provider interface that initially has only the
   portable implementation.

The second option is acceptable for the first correctness milestone. A
oneMKL-specific implementation must use current APIs and must be validated on
a machine with oneMKL; it must not restore the old `mkl_dcsrmm`/`mkl_dbsrmv`
calls or define `__INTEL_MKL__` manually.

## E. Thread ownership contract for CPU-HAM-03B

The sparse backend owns the parallelism of its apply call. The caller must not
place the call inside an OpenMP loop over ensembles or target atoms. The first
portable implementation should use one OpenMP team over independent target
rows, with ensembles either included in the work partition or processed in a
controlled outer serial loop. If a provider uses internal threads, it must be
configured single-threaded whenever the backend is called from an OpenMP-owned
region; nested parallelism is prohibited. The chosen policy and environment
variables must be recorded with benchmark results.

## F. Validation required before production dispatch

Before enabling the backend for production, add a focused oracle test that:

1. builds the sparse structure from an identical current Hamiltonian;
2. applies it to random non-collinear states and multiple ensembles;
3. compares `max_abs_field_difference` and relative field error against
   canonical `HamiltonianActions` direct evaluation;
4. compares eligible pair energy using the field-derived `-1/2` relation;
5. covers zero-neighbor rows, asymmetric directed rows, multi-basis systems,
   long-range scalar `J`, and non-32-aligned atom counts;
6. confirms that ineligible terms select canonical direct evaluation; and
7. checks repeated setup/cleanup and thread-count changes.

The existing test suite has no dedicated sparse oracle. The only current
automated sparse reference is the adaptive-CG setup rejection for `do_sparse`
(`tests/coarse_graining/run_setup_rejection_matrix.py:72`), which verifies a
capability boundary rather than sparse numerical correctness.

## G. Checklist

- [x] `do_sparse` traced.
- [x] Historical implementation located if possible.
- [x] Current connectivity status classified.
- [x] Historical MKL API assessed.
- [x] Scalar-J sparse representation designed.
- [x] Three-RHS approach assessed.
- [x] Build dependency policy defined.
- [x] Thread ownership defined.
- [x] Energy semantics defined.
- [x] GO/NO_GO recommendation produced.
- [x] No production sparse implementation committed.

## Stop/Go gate

**GO for CPU-HAM-03B, constrained to a new persistent scalar-J backend and
oracle-first integration.** MKL-specific implementation and benchmark
validation remain optional and must happen on a host where oneMKL is installed;
the portable provider must remain functional without it.

## Commit

`CPU-HAM-03A: audit sparse Hamiltonian backend`
