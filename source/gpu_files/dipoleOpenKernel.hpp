#pragma once

#include <array>
#include <cstddef>
#include <vector>

// Backend-neutral fp64 geometry for the finite/open displacement tensor.
// Primitive vectors are column-major [C1 C2 C3]:
// primitive_vectors[3*axis + component].
struct DipoleOpenGeometry {
   std::array<std::size_t, 3> active_grid{};
   std::array<std::size_t, 3> fft_grid{};
   // WP10.3 accepts block one only.  The field is retained so the later
   // uniform-block projection can extend this API without changing its model.
   std::array<std::size_t, 3> block_shape{{1, 1, 1}};
   std::array<double, 9> primitive_vectors{};
   unsigned int basis = 0;
   std::vector<std::array<double, 3>> basis_offsets;
};

struct DipoleOpenKernelDiagnostics {
   std::array<std::size_t, 3> active_grid{};
   std::array<std::size_t, 3> fft_grid{};
   std::size_t active_cells = 0;
   std::size_t fft_cells = 0;
   std::size_t kernel_batches = 0;
   std::size_t nonfinite_values = 0;
   bool all_finite = false;
   double minimum_nonself_r2 = 0.0;
   double max_reciprocity_error = 0.0;
   double max_point_self_abs = 0.0;
   double max_padding_gap_abs = 0.0;
};

struct DipoleOpenKernelResult {
   // [fft_cell + fft_cells*kernel_batch], with the first FFT axis fastest and
   // kernel_batch = row + 3*(column + 3*(target_basis + basis*source_basis)).
   std::vector<double> kernel;
   DipoleOpenKernelDiagnostics diagnostics;
};

// Construct the complete dimensionless finite point-dipole tensor in host
// fp64.  This API has no Ewald settings, physical-unit prefactor, GPU backend,
// FFT dependency, or production dispatch side effect.
DipoleOpenKernelResult buildOpenDipoleDisplacementKernel(
   const DipoleOpenGeometry& geometry);
