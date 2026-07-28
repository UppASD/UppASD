#pragma once

#include <array>
#include <cstddef>
#include <vector>

// Backend-neutral fp64 geometry for the finite/open displacement tensor.
// Primitive vectors are column-major [C1 C2 C3]:
// primitive_vectors[3*axis + component].
struct DipoleOpenGeometry {
   // Complete primitive-cell grid before coarse projection.  It must be
   // exactly divisible by block_shape and its quotient must equal
   // active_grid; partial edge blocks are not part of this model.
   std::array<std::size_t, 3> atomistic_grid{};
   std::array<std::size_t, 3> active_grid{};
   std::array<std::size_t, 3> fft_grid{};
   // Full primitive-cell population of every basis-resolved macro channel.
   std::array<std::size_t, 3> block_shape{{1, 1, 1}};
   std::array<double, 9> primitive_vectors{};
   unsigned int basis = 0;
   std::vector<std::array<double, 3>> basis_offsets;
};

struct DipoleOpenKernelDiagnostics {
   std::array<std::size_t, 3> atomistic_grid{};
   std::array<std::size_t, 3> active_grid{};
   std::array<std::size_t, 3> fft_grid{};
   std::array<std::size_t, 3> block_shape{{1, 1, 1}};
   std::size_t block_population = 0;
   std::size_t atomistic_cells = 0;
   std::size_t active_cells = 0;
   std::size_t fft_cells = 0;
   std::size_t kernel_batches = 0;
   std::size_t nonfinite_values = 0;
   bool all_finite = false;
   double minimum_nonself_r2 = 0.0;
   double max_reciprocity_error = 0.0;
   // True point-self contributions are excluded before projection.  This
   // remains exactly zero even though a coarse block's projected diagonal is
   // generally finite due to distinct intra-block pairs.
   double max_point_self_abs = 0.0;
   double max_projected_diagonal_abs = 0.0;
   double max_padding_gap_abs = 0.0;
};

struct DipoleOpenKernelResult {
   // [fft_cell + fft_cells*kernel_batch], with the first FFT axis fastest and
   // kernel_batch = row + 3*(column + 3*(target_basis + basis*source_basis)).
   std::vector<double> kernel;
   DipoleOpenKernelDiagnostics diagnostics;
};

// Construct the complete dimensionless finite point-dipole tensor, projected
// onto uniform block moments M_A=sum(i in A)m_i with m_i=M_A/n_A, in host
// fp64.  Thus K_AB=sum(i in A,j in B)K_ij/(n_A*n_B), including distinct
// intra-block pairs on the coarse diagonal.  This API has no Ewald settings,
// point-macrospin self term, Newell cell tensor, physical-unit prefactor, GPU
// backend, FFT dependency, or production dispatch side effect.
DipoleOpenKernelResult buildOpenDipoleDisplacementKernel(
   const DipoleOpenGeometry& geometry);
