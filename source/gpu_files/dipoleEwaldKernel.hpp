#pragma once

#include <array>
#include <cstddef>
#include <vector>

// This header deliberately has no GPU, FFT, or UppASD precision dependency.
// Builder A is an fp64 host reference used before the tensor is uploaded to a
// future device construction path.

struct DipoleKernelSettings {
   // The gpu_dipole_tol contract.  alpha and cutoffs are selected internally.
   double tolerance = 1.0e-10;
   // Test-only override used to independently exercise the Ewald split.  A
   // non-positive value retains the production automatic-alpha policy.
   double explicit_alpha_for_testing = 0.0;
};

struct DipolePeriodicGeometry {
   // H and Brec are column-major 3x3 matrices: entry 3*column + row.
   std::array<double, 9> H{};
   std::array<double, 9> Brec{};
   double volume = 0.0;
   std::array<std::size_t, 3> grid{};
   unsigned int basis = 0;
   // Cartesian offsets in the full periodic cell's length unit.
   std::vector<std::array<double, 3>> basis_offsets;
};

struct DipoleEwaldParameters {
   double alpha = 0.0;
   std::array<int, 3> real_images{};
   std::array<int, 3> reciprocal_images{};
   double tolerance = 0.0;
};

struct DipoleKernelDiagnostics {
   DipoleEwaldParameters selected{};
   // Largest final independent real/reciprocal shell changes.  They are not
   // differences of the total Ewald sum, so real/reciprocal cancellation
   // cannot make a candidate pass.
   double real_tail_residual = 0.0;
   double reciprocal_tail_residual = 0.0;
   double residual = 0.0;
   double max_alpha_difference = 0.0;
   double max_reciprocity_error = 0.0;
   double max_hermitian_error = 0.0;
   double reciprocal_identity_error = 0.0;
   std::size_t real_tensor_evaluations = 0;
   std::size_t reciprocal_tensor_evaluations = 0;
   std::size_t setup_work = 0;
   double setup_seconds = 0.0;
};

struct DipoleKernelBuildResult {
   // [cell + grid_cells * kernel_batch], with n1 fastest and
   // kernel_batch = row + 3*(column + 3*(target_basis + basis*source_basis)).
   std::vector<double> kernel;
   DipoleKernelDiagnostics diagnostics;
};

// Builds the complete, dimensionless, 3D-periodic tin-foil Ewald tensor in
// target-minus-source convention.  All arithmetic is double regardless of
// the device storage precision selected elsewhere.
DipoleKernelBuildResult buildPeriodicEwaldDisplacementKernel(
   const DipolePeriodicGeometry& geometry,
   const DipoleKernelSettings& settings = {});

// Validate target/source-basis displacement reciprocity independently of the
// construction loop.  This is also useful after a future upload/repacking.
double periodicKernelReciprocityError(const std::vector<double>& kernel,
                                      const DipolePeriodicGeometry& geometry);
