# CPU-HAM-01A — Establish one canonical field/energy Hamiltonian contract

**Model:** Luna

## Dependency

`CPU-HAM-00` complete.

## Purpose

Preserve `effective_field()` as an energy-capable canonical Hamiltonian API while allowing callers that do not need energy to avoid unnecessary energy work.

The objective is **not** to separate field and energy physics.

The objective is:

> same term implementation → field → optional energy

## Goal

Make energy execution explicit without creating duplicate Hamiltonian kernels.

## Implementation result

The production entry points are the `effective_field` generic interface and its
`effective_field_bare`/`effective_field_full` procedures in
`source/Hamiltonian/hamiltonianactions.f90`. Before this task, the full routine
always accumulated and OpenMP-reduced `energy`, including when SD/MS predictor
and corrector callers only consumed `beff`, `beff1`, and `beff2`.

The full and bare procedures now accept an optional logical
`measure_energy`. It defaults to `.true.` for source and behavior compatibility,
so existing energy consumers remain unchanged. With `.false.`, the canonical
per-atom field assembly is still executed, but the field loop has no energy
expression or energy reduction. The dipole manager has the same optional gate
for its direct and macrocell energy reductions. The field assembly itself lives
in one shared `effective_field_atom` routine; no second Hamiltonian physics
implementation was introduced.

Ordinary CPU SD and MS initial/predictor/corrector calls in `sd_driver.f90` and
`ms_driver.f90` explicitly pass `measure_energy=.false.`. The SD field-print
and fallback CPU calls are also field-only. SLD calls retain the default energy
request because their `totenergy` value is consumed by the spin-lattice energy
path between explicit `calc_energy` calls. PT, WL, GNEB, minimization, and MC
callers that use the returned energy retain the default behavior. Global
measurement energy continues to enter through `Energy::calc_energy`, called by
the existing measurement cadence.

## Supported-term energy audit

The table records the current production term ownership and the extraction
conventions used by `Energy::calc_energy`. The factor is applied to
`-m·B_term` through `update_ene`; it is not a universal rule for every term.

| Term | Field implementation | Energy relation / existing factor | Directly from field? | Independent routine |
|---|---|---|---|---|
| Isotropic exchange `J` | `heisenberg_field` or `heisenberg_rescaling_field` | reciprocal pair: `-1/2 Σ m·B` | yes | `Energy::calc_energy` |
| Exchange tensor | `tensor_field` | reciprocal pair: `-1/2 Σ m·B` | yes | `Energy::calc_energy` |
| DMI | `dzyaloshinskii_moriya_field` | reciprocal directed-pair fold: `-1/2 Σ m·B` | yes when reciprocal contract holds | `Energy::calc_energy` |
| Symmetric anisotropic / pseudo-dipolar | `symmetric_anisotropic_field` / `pseudo_dipolar_field` | existing `update_ene` factor `1/2` | yes for the supported bilinear contract | `Energy::calc_energy` |
| BIQDM | `dzyaloshinskii_moriya_bq_field` | existing `update_ene` factor `1/2` | term-specific; not universal | `Energy::calc_energy` |
| Biquadratic exchange | `biquadratic_field` | existing `update_ene` factor `1/4` | term-specific; not universal | `Energy::calc_energy` |
| Ring exchange | `ring_field` | existing `update_ene` factor `1/4` | term-specific; not universal | `Energy::calc_energy` |
| Scalar chirality | `chirality_field` | existing `update_ene` factor `1/2` | term-specific; not universal | `Energy::calc_energy` |
| Dipole | `DipoleManager::dipolar_field`, macrocell, or FFT provider | reciprocal pair: `-1/2 Σ m·B` (macrocell uses its existing cell normalization) | yes for the supported pair operator | `Energy::calc_energy` / dipole manager |
| Uniaxial/cubic anisotropy | `uniaxial_anisotropy_field` / `cubic_anisotropy_field` | existing term-specific `update_ene` factor `1/2` | not by the pair rule | `Energy::calc_energy` |
| Zeeman / time-dependent external field | field assembly from `external_field` + `time_external_field` | `-Σ m·B_ext`, factor `1` | yes, but not halved | `Energy::calc_energy` |
| LSF | LSF field/evaluator path | `totalenergy_LSF` and its existing interpolation conventions | no | `totalenergy_LSF` |

No broad Energy-module rewrite was made. Future work can redirect global
bilinear energy reporting to cached canonical field output, but the current
term-resolved measurement behavior remains the authority for output and was
left intact in this task.

## MC exception evidence

Single-spin-flip and local-update MC kernels remain specialized by design:
`calculate_efield`/`calculate_efield_tensor` provide local fields and
`montecarlo_common::calculate_energy` evaluates local trial energy changes.
The existing `tests/coarse_graining/test_dmi_dimer_energy.f90` covers the
specialized DMI field/energy convention against a hand-derived dimer oracle.
The new CPU-HAM-01A test covers the canonical global field-energy relation for
reciprocal `J+D` and dipole pairs. A full all-term MC-before/after parity matrix
is deferred because it would broaden the MC refactor beyond this contract task.

## Field parity evidence

`tests/hamiltonian/test_field_energy_contract.f90` calls the same production
state with `measure_energy=.true.`, `.false.`, and the omitted legacy default.
It checks exact equality of total/internal/external fields for the J+D case and
the dipole case, zero energy for field-only calls, equality of explicit/default
energy, and the `-1/2 Σ m·B_pair` relation. It passes with one and four OpenMP
threads. Deterministic end-to-end T=0 restart files for bcc Fe, 2D skyrmion,
and dhcp Nd were also byte-identical between field-only and energy-enabled SD
production binaries.

## Performance evidence

The same CPU production executable and 40-step dynamics-only input were run
three times at one and eight OpenMP threads. The field-only binary used the
new SD caller gates. A temporary comparison binary used identical source and
inputs with those seven SD call sites set to `measure_energy=.true.`; it was
not retained. Values below are medians of complete-process wall time in
seconds. B02 includes its fixed triple-Q initial phase, so it is a negative
control for interpreting whole-process energy overhead.

| Workload | 1 thread field-only | 1 thread field+energy | 8 threads field-only | 8 threads field+energy | Observation |
|---|---:|---:|---:|---:|---|
| bcc Fe, 20x20x20, 16,000 atoms | 3.12 | 3.29 | 1.54 | 1.52 | low-thread reduction overhead is visible; 8-thread noise dominates |
| 2D skyrmion `J+D`, 128x128, 16,384 atoms | 0.62 | 0.62 | 0.60 | 0.57 | fixed initial phase masks the steady-state difference |
| dhcp Nd, 16x16x16, 16,384 atoms | 2.96 | 2.94 | 1.02 | 1.02 | long-range gather dominates; no material whole-process delta at this size |

The measurement establishes that field-only execution is available and removes
the reduction from ordinary LLG/MS calls; it does not claim a universal speedup
from short noisy process runs. Full metric collection (counter-level bandwidth,
IPC, LLC misses, and vectorization) belongs to CPU-HAM-00/05 and was not
repeated here.

## A. Audit current API

Identify:
- `effective_field` entry points;
- energy arguments;
- ordinary LLG callers;
- measurement callers;
- current Energy module callers;
- whether energy is currently unconditional.

## B. Introduce explicit energy request semantics

Following existing UppASD style, provide an explicit mechanism conceptually equivalent to:

```fortran
measure_energy = .false.
```

or:

```fortran
calculate_energy = .true./.false.
```

Do not invent a new global state if an existing argument mechanism is cleaner.

## C. Ordinary LLG

If CPU-HAM-00 confirms predictor/corrector LLG does not consume energy, route those calls through field-only execution without energy reduction.

Do **not** remove the ability of `effective_field()` to calculate energy.

## D. Bilinear pair energy

For eligible reciprocal bilinear pair contributions establish/test:

\[
E_{\rm pair}=-\frac12\sum_i\mathbf m_i\cdot\mathbf B_i^{\rm pair}.
\]

Use the field produced by the canonical `HamiltonianActions` pair routine.

Do not perform another neighbour traversal merely to calculate the same pair energy unless evidence demonstrates that doing so is computationally superior.

## E. Non-bilinear terms

Audit each supported term.

Create a table with:
- term;
- field routine;
- energy relation;
- whether directly derivable from field;
- required factor/convention;
- existing independent energy routine.

Do not force all terms into `-1/2 m·B`.

## F. Energy module

Do not broadly rewrite the Energy module in this task.

Document which global energy calculations could eventually be redirected to canonical `HamiltonianActions` output.

## G. MC

Document that single-flip MC/local-update `ΔE` kernels remain allowed specialized implementations.

Add or identify a parity test comparing one specialized MC `ΔE` against canonical total-energy difference if feasible without broad MC refactoring.

## H. Field parity

Field results with energy disabled/enabled must be identical within the strictest meaningful numerical criterion.

Prefer bitwise identity if only post-field energy extraction differs.

## I. Performance

Measure field-only and field+energy cost for:
- Nd;
- bcc Fe;
- short-range `J+D`.

## Checklist

- [x] Energy request semantics explicit.
- [x] `effective_field` remains energy capable.
- [x] No duplicate field implementation created.
- [x] LLG can skip energy when unneeded.
- [x] Bilinear pair energy derives from canonical field.
- [x] Term-specific energy table created.
- [x] No universal incorrect factor-1/2 assumption introduced.
- [x] Field parity with/without energy passes.
- [x] Existing measurement-energy behavior preserved.
- [x] MC exception documented.
- [x] Field-only performance measured.
- [x] Field+energy performance measured.

## Commit

`CPU-HAM-01A: unify field and energy execution semantics`
