#pragma once

#include <array>
#include <cstddef>
#include <complex>
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

// Builder B keeps the short-range part in the real displacement
// representation and constructs only the reciprocal contribution directly in
// the normalized R2C spectrum.  ``real_kernel`` has the same layout as
// DipoleKernelBuildResult::kernel.  ``reciprocal_alias_spectrum`` is laid out
// [spectral_cell + spectral_cells * kernel_batch], with n1 as the retained
// R2C axis.  It is already normalized for the raw C2R contract: no 1/Ngrid
// factor is to be applied to this vector.
struct DipoleAliasSpectrumBuildResult {
   std::vector<double> real_kernel;
   std::vector<std::complex<double>> reciprocal_alias_spectrum;
   std::array<std::size_t, 3> spectral_grid{};
   DipoleKernelDiagnostics diagnostics;
};

// Builds the complete, dimensionless, 3D-periodic tin-foil Ewald tensor in
// target-minus-source convention.  All arithmetic is double regardless of
// the device storage precision selected elsewhere.
DipoleKernelBuildResult buildPeriodicEwaldDisplacementKernel(
   const DipolePeriodicGeometry& geometry,
   const DipoleKernelSettings& settings = {});

// Production Builder B.  It is algebraically equivalent to Builder A but
// avoids evaluating every reciprocal vector for every real-space
// displacement.  The real-space tensor and point self term are returned for
// one batched R2C; reciprocal aliases are added directly to that spectrum.
DipoleAliasSpectrumBuildResult buildPeriodicEwaldAliasSpectrum(
   const DipolePeriodicGeometry& geometry,
   const DipoleKernelSettings& settings = {});

// Slow, CPU-only reference transform used solely by cross-builder tests and
// benchmarks.  It returns FFT(kernel)/Ngrid in the same half-spectrum layout
// as Builder B.  Do not use it in production construction.
std::vector<std::complex<double>> referenceNormalizedKernelSpectrum(
   const std::vector<double>& kernel, const DipolePeriodicGeometry& geometry);

// Validate target/source-basis displacement reciprocity independently of the
// construction loop.  This is also useful after a future upload/repacking.
double periodicKernelReciprocityError(const std::vector<double>& kernel,
                                      const DipolePeriodicGeometry& geometry);
