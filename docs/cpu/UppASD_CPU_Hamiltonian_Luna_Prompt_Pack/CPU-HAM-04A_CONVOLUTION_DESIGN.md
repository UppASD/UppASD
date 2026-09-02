# CPU-HAM-04A — Establish scalar-J CPU lattice-convolution mathematics and oracle

**Model:** Luna

## Dependency

`CPU-HAM-03B` complete, or `CPU-HAM-03A` NO_GO recorded.

## Purpose

Prepare the CPU FFT-convolution backend without immediately introducing a large optimized production implementation.

## Initial scope

Only:
- scalar isotropic `J`;
- fully periodic;
- regular replicated lattice;
- reduced Hamiltonian;
- translation-invariant pair coupling.

## A. Map GPU convolution semantics

Inspect validated GPU lattice-convolution implementation.

Document lattice/basis indexing, periodic transform convention, normalization, R2C/C2R dimensions, interaction-kernel Fourier transform, basis-pair structure, and sign convention.

Use it as guidance, not something to copy blindly.

## B. Derive scalar-J CPU form

For basis `a,b`:

\[
B_{a\alpha}(q)=\sum_b J_{ab}(q)M_{b\alpha}(q),\qquad \alpha=x,y,z.
\]

Demonstrate why scalar `J` requires only `N_A^2` spectral kernels reused over x/y/z.

## C. Convolution vs dipole FFT

Document:
- lattice pair convolution → periodic circular convolution, no free-space `2N-1` padding requirement;
- dipole FFT → may use padded/free-space semantics.

Do not reuse dipole padding rules blindly.

## D. Reference implementation

Implement small correctness-first CPU reference path using an available FFT provider.

Performance is secondary.

Inefficient reference lifecycle must not be wired into production.

## E. FFT normalization

Prove normalization using:
- delta spin;
- uniform ferromagnetic state;
- single-mode spin spiral where useful.

## F. Multi-basis

Include at least one `N_A>1` fixture.

## G. Long-range oracle

Compare reference convolution against DIRECT on an Nd-derived periodic case.

## H. Ineligibility

Reject:
- nonperiodic geometry;
- non-reduced Hamiltonian;
- translationally non-invariant coupling.

## I. Design production object

Specify persistent:
- FFT plans/descriptors;
- atom↔cell/basis mapping;
- real/complex buffers;
- `J_ab(q)`;
- normalization constants.

## J. FFT provider boundary

Design minimal provider boundary so raw MKL DFTI does not spread throughout `HamiltonianActions`.

Do not build a broad generic FFT framework.

## K. Deliverable

Create `docs/CPU_HAM_04A_CONVOLUTION_DESIGN.md`.

## Checklist

- [x] GPU convolution semantics audited.
- [x] Scalar-J spectral equations documented.
- [x] `N_A^2` kernel reduction established.
- [x] Periodic normalization proven.
- [x] Dipole padding distinction documented.
- [x] Reference FFT implementation works.
- [x] Delta/uniform fixture passes.
- [x] Multi-basis fixture passes.
- [x] Nd DIRECT parity passes.
- [x] Ineligible cases rejected.
- [x] Persistent production state designed.
- [x] Minimal FFT provider boundary designed.

## Implementation result

The CPU reference is implemented in `source/Hamiltonian/cpuconvolution.f90`
with the provider-specific plan boundary in
`source/Hamiltonian/cpu_fft_provider.f90`. It is compiled only when
`USE_FFTW=ON`, and is intentionally not included in `HamiltonianActions` or
the production timestep dispatch. A future 04B implementation can reuse the
object layout and replace the provider without moving raw FFT calls into the
Hamiltonian physics module.

### GPU semantics audit

The validated GPU lattice convolution uses basis-fastest atom order:

```text
atom = basis + N_A * (cell_x + N_1 * (cell_y + N_2 * cell_z))
```

Its spin packing groups real values as `(axis, basis, ensemble, cell)`. The
R2C grid is `(N_1/2+1) x N_2 x N_3`, represented by reversed FFT dimensions
`(N_3,N_2,N_1)` for a contiguous Fortran x/cell dimension. The GPU scalar
kernel puts the same `J_ab` into the x, y, and z diagonal components and the
three components reuse those values. Tensor and DMI components are separate
GPU extensions and remain outside this scalar-J task.

The CPU reference consumes the already validated `ReducedStencil` records,
whose `delta_cell` means “source cell relative to target cell”. Since a normal
FFT product evaluates `sum_d K(d) M(R-d)`, the reference stores a record at
`-delta_cell` before transforming the kernel. This preserves the production
target-to-source convention, including negative displacements and periodic
wraps. It is the deliberate CPU provider/layout translation; GPU kernel-array
placement is not copied blindly.

### Scalar-J equations and kernel reduction

For each basis pair and Cartesian component:

```text
B_a,alpha(q) = sum_b J_ab(q) M_b,alpha(q),  alpha in {x,y,z}
```

The Cartesian components do not mix for scalar isotropic exchange, so the same
`J_ab(q)` is reused for x, y, and z. The spectral state therefore stores
exactly `N_A*N_A` complex kernels, rather than a generic `9*N_A*N_A` tensor
set. Repeated stencil records at the same wrapped displacement are summed
during setup.

### Periodic normalization and dipole distinction

The reference uses circular convolution on the original `N_1 x N_2 x N_3`
lattice. No `2N-1` free-space padding is used. FFTW's backward transform is
unnormalized, so the field is multiplied by `1/(N_1*N_2*N_3)` exactly once.
The delta fixture checks displacement/sign, and the uniform fixture checks
that a constant state returns the coupling sum without an extra grid factor.

This is different from dipole FFTs for finite/open geometry. The existing
dipole implementation uses padded/free-space dimensions (`2N-1`) where its
boundary semantics require them. Those padding rules must not be reused for a
periodic lattice pair convolution.

### Reference object and provider boundary

`cpu_convolution_t` persists:

- forward, backward, and kernel FFT plan handles;
- atom-to-cell/basis and cell/basis-to-atom mappings;
- real spin/field work storage;
- complex spin and field spectra;
- real kernel storage and the `N_A^2` complex spectral kernels;
- grid dimensions, ensemble count, and the inverse-transform scale.

`cpu_convolution_init` performs allocation and plan creation, while
`cpu_convolution_build_kernel` performs the one-time stencil-to-spectrum
construction. `cpu_convolution_apply` only packs, transforms, multiplies,
inverse-transforms, scales, and unpacks. No plan creation or allocation occurs
in apply. `cpu_fft_provider.f90` is the only file containing FFTW calls.

### Correctness evidence

`tests/hamiltonian/test_cpu_convolution.f90` covers:

- a one-basis delta displacement and explicit sign check;
- a uniform state and field-derived `-1/2 sum(m dot B)` pair energy;
- a two-basis, three-ensemble periodic fixture with boundary-crossing
  displacements;
- an Nd-shaped long-range fixture (4 basis atoms, `7 x 5 x 3` cells, seven
  records per basis) against the exact periodic DIRECT stencil oracle;
- nonperiodic, non-reduced, and deliberately non-translation-invariant input
  rejection;
- the exact `N_A^2` spectral-kernel count for the multi-basis case.

The focused test was built and run locally with FFTW:

```text
cmake -S . -B build/cpu-04a -DUSE_FFTW=ON -DBUILD_TESTING=ON \
  -DRUN_REG_TESTS=OFF -DRUN_ASD_TESTS=OFF
cmake --build build/cpu-04a --target cpu_convolution_tests -j2
OMP_NUM_THREADS=1 ./build/cpu-04a/bin/cpu_convolution_tests
CPU-HAM-04A CPU convolution tests passed
```

MKL-specific testing was not required: this task uses the local FFTW provider
and does not add an MKL code path. A separate MKL host is only needed when a
provider-specific MKL implementation or MKL threading behavior is introduced.

## Commit

`CPU-HAM-04A: validate scalar exchange CPU convolution`
