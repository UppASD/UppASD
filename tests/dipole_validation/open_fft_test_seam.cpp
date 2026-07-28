#include "open_fft_test_seam.hpp"

#include "dipoleOpenKernel.hpp"
#include "gpuDipoleConvolution.hpp"
#include "tensor.hpp"

#include <cmath>
#include <limits>
#include <numeric>
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

Result runImpl(const Input& input, bool production_kernel) {
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
   if(!production_kernel && input.padded_real_kernel.size() != kernel_size)
      throw std::invalid_argument("OPEN_FFT test seam padded real kernel has the wrong shape");

   GpuTensor<real, 3> device_moments;
   GpuTensor<unsigned int, 1> device_map;
   GpuTensor<real, 3> device_beff, device_eneff;
   GpuTensor<real, 2> device_energy;
   GpuDipoleConvolution solver;
   GPU_STREAM_T stream{};
   if(GPU_STREAM_CREATE(&stream) != GPU_SUCCESS) throw std::runtime_error("OPEN_FFT test seam stream creation failed");
   try {
      device_moments.Allocate(3L, static_cast<long int>(layout.active_macros), static_cast<long int>(input.ensembles));
      device_map.Allocate(static_cast<long int>(layout.active_macros));
      device_beff.Allocate(3L, static_cast<long int>(layout.active_macros), static_cast<long int>(input.ensembles));
      device_eneff.Allocate(3L, static_cast<long int>(layout.active_macros), static_cast<long int>(input.ensembles));
      device_energy.Allocate(static_cast<long int>(input.ensembles), 7L);
      std::vector<real> moments(input.active_moments.begin(), input.active_moments.end());
      std::vector<unsigned int> identity_map(layout.active_macros);
      std::iota(identity_map.begin(), identity_map.end(), 1U);
      if(GPU_MEMCPY(device_moments.data(), moments.data(), moments.size() * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS ||
         GPU_MEMCPY(device_map.data(), identity_map.data(), identity_map.size() * sizeof(unsigned int),
                    GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS)
         throw std::runtime_error("OPEN_FFT test seam moment upload failed");
      device_beff.zeros_async(stream);
      device_eneff.zeros_async(stream);
      device_energy.zeros_async(stream);
      if(!solver.initiate(descriptor, stream)) throw std::runtime_error("OPEN_FFT test seam plan initiation failed");
      if(production_kernel) {
         DipoleOpenGeometry geometry{};
         geometry.active_grid = input.active_grid;
         geometry.fft_grid = input.fft_grid;
         geometry.primitive_vectors = input.primitive_vectors;
         geometry.basis = input.basis;
         geometry.basis_offsets = input.basis_offsets;
         solver.uploadOpenKernel(buildOpenDipoleDisplacementKernel(geometry).kernel);
      } else {
         // This is deliberately a testing upload seam.  It cannot select an
         // UppASD mode and does not invoke the periodic Ewald builder.
         solver.uploadCompleteKernelForTesting(input.padded_real_kernel, false);
      }
      solver.evaluate(device_moments);
      solver.addFieldsToAtoms(device_beff, device_eneff, device_map.data(), layout.active_macros);
      solver.accumulateEnergy(device_energy, layout.active_macros);
      Result result{};
      result.active_cells = layout.active_cells;
      result.active_macros = layout.active_macros;
      result.fft_cells = layout.fft_cells;
      result.field_batches = layout.field_batches;
      result.packed_moments = doubles(solver.diagnosticPackedMomentsForTesting());
      result.padded_fields = doubles(solver.diagnosticPaddedFieldsForTesting());
      result.active_fields = doubles(solver.diagnosticFieldsForTesting());
      result.dimensionless_energy = solver.diagnosticEnergiesForTesting();
      std::vector<real> scattered(3 * layout.active_macros * input.ensembles);
      std::vector<real> energy(7 * input.ensembles);
      if(GPU_MEMCPY(scattered.data(), device_beff.data(), scattered.size() * sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS ||
         GPU_MEMCPY(energy.data(), device_energy.data(), energy.size() * sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS) {
         throw std::runtime_error("OPEN_FFT test seam production field/energy download failed");
      }
      result.scattered_fields = doubles(scattered);
      result.energy_per_atom.resize(input.ensembles);
      for(unsigned int ensemble = 0; ensemble < input.ensembles; ++ensemble)
         result.energy_per_atom[ensemble] = energy[ensemble + input.ensembles * 5];
      result.persistent_bytes = GpuDipoleConvolution::estimatePersistentBytes(descriptor);
      result.construction_bytes = GpuDipoleConvolution::estimateConstructionBytes(descriptor);
      result.workspace_bytes = GpuDipoleConvolution::estimateWorkspaceBytes(descriptor);
      result.total_bytes = GpuDipoleConvolution::estimateBytes(descriptor);
      result.persistent_inventory_bytes =
         2 * layout.fft_cells * layout.field_batches * sizeof(real) +
         2 * layout.spectral_cells * layout.field_batches * sizeof(GpuFftComplex) +
         layout.spectral_cells * layout.kernel_batches * sizeof(GpuFftComplex) + 19 * sizeof(real);
      // OPEN_FFT transforms the finite real tensor directly into the
      // persistent spectrum; unlike the periodic Builder B it has no alias
      // spectrum allocation during construction.
      result.construction_inventory_bytes =
         layout.fft_cells * layout.kernel_batches * sizeof(real);
      solver.release();
      device_energy.Free();
      device_eneff.Free();
      device_beff.Free();
      device_map.Free();
      device_moments.Free();
      if(GPU_STREAM_DESTROY(stream) != GPU_SUCCESS) throw std::runtime_error("OPEN_FFT test seam stream destruction failed");
      return result;
   } catch(...) {
      if(!device_energy.empty()) device_energy.Free();
      if(!device_eneff.empty()) device_eneff.Free();
      if(!device_beff.empty()) device_beff.Free();
      if(!device_map.empty()) device_map.Free();
      if(!device_moments.empty()) device_moments.Free();
      solver.release();
      GPU_STREAM_DESTROY(stream);
      throw;
   }
}

Result run(const Input& input) {
   return runImpl(input, false);
}

Result runProduction(const Input& input) {
   return runImpl(input, true);
}

} // namespace luna_open_fft_test
