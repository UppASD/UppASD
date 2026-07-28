#include "open_fft_test_seam.hpp"

#include "gpuDipoleConvolution.hpp"
#include "tensor.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace luna_open_fft_test {
namespace {

std::size_t cells(const std::array<std::size_t, 3>& grid) {
   if(grid[0] == 0 || grid[1] == 0 || grid[2] == 0 ||
      grid[0] > std::numeric_limits<std::size_t>::max() / grid[1] ||
      grid[0] * grid[1] > std::numeric_limits<std::size_t>::max() / grid[2]) {
      throw std::invalid_argument("OPEN_FFT test seam has an invalid grid extent");
   }
   return grid[0] * grid[1] * grid[2];
}

std::vector<double> doubles(const std::vector<real>& values) {
   return {values.begin(), values.end()};
}

} // namespace

Result run(const Input& input) {
   const std::size_t active_cells = cells(input.active_grid);
   const std::size_t fft_cells = cells(input.fft_grid);
   if(input.basis == 0 || input.ensembles == 0 || input.basis_offsets.size() != input.basis ||
      input.active_moments.size() != 3 * active_cells * input.basis * input.ensembles) {
      throw std::invalid_argument("OPEN_FFT test seam input has an inconsistent active macro tensor");
   }
   for(unsigned int axis = 0; axis < 3; ++axis)
      if(input.active_grid[axis] > std::numeric_limits<std::size_t>::max() / 2 + 1 ||
         input.fft_grid[axis] < 2 * input.active_grid[axis] - 1) {
         throw std::invalid_argument("OPEN_FFT test seam fft_grid is not a legal finite-convolution padding");
      }

   std::array<real, 3> c1{}, c2{}, c3{};
   for(unsigned int component = 0; component < 3; ++component) {
      c1[component] = static_cast<real>(input.primitive_vectors[component]);
      c2[component] = static_cast<real>(input.primitive_vectors[3 + component]);
      c3[component] = static_cast<real>(input.primitive_vectors[6 + component]);
   }
   GpuDipoleConvolutionDescriptor descriptor{};
   descriptor.atomistic_grid = {input.active_grid[0], input.active_grid[1], input.active_grid[2]};
   descriptor.macro_grid = descriptor.atomistic_grid;
   descriptor.fft_grid = {input.fft_grid[0], input.fft_grid[1], input.fft_grid[2]};
   descriptor.basis = input.basis;
   descriptor.ensembles = input.ensembles;
   descriptor.boundary = GpuDipoleBoundaryMode::Open;
   descriptor.discretization = GpuDipoleDiscretization::MacrospinGrid;
   descriptor.c1 = c1.data(); descriptor.c2 = c2.data(); descriptor.c3 = c3.data();
   descriptor.basis_offsets = input.basis_offsets;
   descriptor.alat = 1.0;
   descriptor.tolerance = 1.0e-10;
   descriptor.field_prefactor = 1.0;
   if(!descriptor.valid()) throw std::invalid_argument("OPEN_FFT test seam descriptor is invalid");

   const auto layout = descriptor.fftLayout();
   const std::size_t kernel_size = layout.fft_cells * layout.kernel_batches;
   if(input.padded_real_kernel.size() != kernel_size)
      throw std::invalid_argument("OPEN_FFT test seam padded real kernel has the wrong shape");

   GpuTensor<real, 3> device_moments;
   GpuDipoleConvolution solver;
   GPU_STREAM_T stream{};
   if(GPU_STREAM_CREATE(&stream) != GPU_SUCCESS) throw std::runtime_error("OPEN_FFT test seam stream creation failed");
   try {
      device_moments.Allocate(3L, static_cast<long int>(layout.active_macros), static_cast<long int>(input.ensembles));
      std::vector<real> moments(input.active_moments.begin(), input.active_moments.end());
      if(GPU_MEMCPY(device_moments.data(), moments.data(), moments.size() * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS)
         throw std::runtime_error("OPEN_FFT test seam moment upload failed");
      if(!solver.initiate(descriptor, stream)) throw std::runtime_error("OPEN_FFT test seam plan initiation failed");
      // This is deliberately a testing upload seam.  It cannot select an
      // UppASD mode and does not invoke the periodic Ewald builder.
      solver.uploadCompleteKernelForTesting(input.padded_real_kernel, false);
      solver.evaluate(device_moments);
      Result result{};
      result.active_cells = layout.active_cells;
      result.active_macros = layout.active_macros;
      result.fft_cells = layout.fft_cells;
      result.field_batches = layout.field_batches;
      result.packed_moments = doubles(solver.diagnosticPackedMomentsForTesting());
      result.padded_fields = doubles(solver.diagnosticPaddedFieldsForTesting());
      result.active_fields = doubles(solver.diagnosticFieldsForTesting());
      result.dimensionless_energy = solver.diagnosticEnergiesForTesting();
      result.persistent_bytes = GpuDipoleConvolution::estimatePersistentBytes(descriptor);
      result.construction_bytes = layout.constructionBytes();
      result.workspace_bytes = GpuDipoleConvolution::estimateWorkspaceBytes(descriptor);
      result.total_bytes = GpuDipoleConvolution::estimateBytes(descriptor);
      result.persistent_inventory_bytes =
         2 * layout.fft_cells * layout.field_batches * sizeof(real) +
         2 * layout.spectral_cells * layout.field_batches * sizeof(GpuFftComplex) +
         layout.spectral_cells * layout.kernel_batches * sizeof(GpuFftComplex) + 19 * sizeof(real);
      result.construction_inventory_bytes =
         layout.fft_cells * layout.kernel_batches * sizeof(real) +
         layout.spectral_cells * layout.kernel_batches * sizeof(GpuFftComplex);
      solver.release();
      device_moments.Free();
      if(GPU_STREAM_DESTROY(stream) != GPU_SUCCESS) throw std::runtime_error("OPEN_FFT test seam stream destruction failed");
      return result;
   } catch(...) {
      if(!device_moments.empty()) device_moments.Free();
      solver.release();
      GPU_STREAM_DESTROY(stream);
      throw;
   }
}

} // namespace luna_open_fft_test
