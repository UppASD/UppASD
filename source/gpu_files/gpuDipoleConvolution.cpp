#include "gpuDipoleConvolution.hpp"

#include "gpuFftWrapper.hpp"
#include "real_type.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

bool multiply(std::size_t lhs, std::size_t rhs, std::size_t& result) {
   if(lhs == 0 || rhs == 0) {
      result = 0;
      return true;
   }
   if(lhs > std::numeric_limits<std::size_t>::max() / rhs) return false;
   result = lhs * rhs;
   return true;
}

bool add(std::size_t lhs, std::size_t rhs, std::size_t& result) {
   if(lhs > std::numeric_limits<std::size_t>::max() - rhs) return false;
   result = lhs + rhs;
   return true;
}

bool bytesFor(std::size_t elements, std::size_t element_bytes, std::size_t& bytes) {
   return multiply(elements, element_bytes, bytes);
}

bool hasNonSingularCell(const GpuDipoleConvolutionDescriptor& descriptor) {
   return descriptor.cellVolume() != static_cast<real>(0);
}

} // namespace

namespace {

__global__ void pack_macro_moments_kernel(const real* macro_moments, real* packed,
                                           std::size_t cells, unsigned int basis,
                                           unsigned int ensembles) {
   const std::size_t total = 3 * cells * basis * ensembles;
   for(std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < total; index += static_cast<std::size_t>(blockDim.x) * gridDim.x) {
      const unsigned int component = static_cast<unsigned int>(index % 3);
      const std::size_t macro_ensemble = index / 3;
      const unsigned int ensemble = static_cast<unsigned int>(macro_ensemble / (cells * basis));
      const std::size_t macro = macro_ensemble % (cells * basis);
      const unsigned int channel = static_cast<unsigned int>(macro % basis);
      const std::size_t cell = macro / basis;
      const std::size_t batch = component + 3 * (channel + basis * ensemble);
      packed[cell + cells * batch] = macro_moments[index];
   }
}

__global__ void apply_spectral_kernel(const GpuFftComplex* moments, const GpuFftComplex* kernel,
                                      GpuFftComplex* fields, std::size_t spectral_cells,
                                      unsigned int basis, unsigned int ensembles) {
   const std::size_t total = spectral_cells * 3 * basis * ensembles;
   for(std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < total; index += static_cast<std::size_t>(blockDim.x) * gridDim.x) {
      const std::size_t spectral = index % spectral_cells;
      const std::size_t batch = index / spectral_cells;
      const unsigned int row = static_cast<unsigned int>(batch % 3);
      const std::size_t channel_ensemble = batch / 3;
      const unsigned int target_basis = static_cast<unsigned int>(channel_ensemble % basis);
      const unsigned int ensemble = static_cast<unsigned int>(channel_ensemble / basis);
      real re = 0;
      real im = 0;
      for(unsigned int source_basis = 0; source_basis < basis; ++source_basis) {
         for(unsigned int column = 0; column < 3; ++column) {
            const std::size_t kernel_batch = row + 3 * (column + 3 * (target_basis + basis * source_basis));
            const std::size_t moment_batch = column + 3 * (source_basis + basis * ensemble);
            const GpuFftComplex k = kernel[spectral + spectral_cells * kernel_batch];
            const GpuFftComplex m = moments[spectral + spectral_cells * moment_batch];
            re += k.x * m.x - k.y * m.y;
            im += k.x * m.y + k.y * m.x;
         }
      }
      fields[index].x = re;
      fields[index].y = im;
   }
}

} // namespace

std::size_t GpuDipoleGridShape::cells() const {
   if(!valid() || n1 > std::numeric_limits<std::size_t>::max() / n2 ||
      n1 * n2 > std::numeric_limits<std::size_t>::max() / n3) {
      return 0;
   }
   return n1 * n2 * n3;
}

bool GpuDipoleGridShape::valid() const {
   return n1 != 0 && n2 != 0 && n3 != 0;
}

bool GpuDipoleFftLayout::valid() const {
   return real_grid.valid() && spectral_grid.valid() && real_cells != 0 && spectral_cells != 0 &&
          field_batches != 0 && kernel_batches != 0;
}

bool GpuDipoleFftLayout::fitsFftLibrary() const {
   return valid() && real_grid.n1 <= INT_MAX && real_grid.n2 <= INT_MAX && real_grid.n3 <= INT_MAX &&
          field_batches <= INT_MAX && kernel_batches <= INT_MAX && real_cells <= INT_MAX &&
          spectral_cells <= INT_MAX;
}

std::size_t GpuDipoleFftLayout::persistentBytes() const {
   if(!valid()) return 0;
   std::size_t real_fields = 0, spectral_fields = 0, kernel = 0;
   std::size_t bytes = 0, part = 0;
   if(!multiply(real_cells, field_batches, real_fields) ||
      !multiply(spectral_cells, field_batches, spectral_fields) ||
      !multiply(spectral_cells, kernel_batches, kernel) ||
      !bytesFor(real_fields, sizeof(real), bytes) ||
      !bytesFor(spectral_fields, sizeof(GpuFftComplex), part) || !add(bytes, part, bytes) ||
      !bytesFor(spectral_fields, sizeof(GpuFftComplex), part) || !add(bytes, part, bytes) ||
      !bytesFor(kernel, sizeof(GpuFftComplex), part) || !add(bytes, part, bytes) ||
      !bytesFor(19, sizeof(real), part) || !add(bytes, part, bytes)) return 0;
   return bytes;
}

std::size_t GpuDipoleFftLayout::constructionBytes() const {
   if(!valid()) return 0;
   std::size_t elements = 0, bytes = 0;
   return multiply(real_cells, kernel_batches, elements) && bytesFor(elements, sizeof(real), bytes) ? bytes : 0;
}

GpuDipoleGridShape GpuDipoleConvolutionDescriptor::activeGrid() const {
   return discretization == GpuDipoleDiscretization::MacrospinGrid ? macro_grid : atomistic_grid;
}

GpuDipoleGridShape GpuDipoleConvolutionDescriptor::paddedGrid() const {
   GpuDipoleGridShape padded = activeGrid();
   if(!padded.valid()) return {};

   const auto pad_axis = [](std::size_t extent) {
      constexpr std::size_t max = std::numeric_limits<std::size_t>::max();
      return extent > max / 2 + 1 ? 0 : 2 * extent - 1;
   };
   if(boundary == GpuDipoleBoundaryMode::Open) {
      padded.n1 = pad_axis(padded.n1);
      padded.n2 = pad_axis(padded.n2);
      padded.n3 = pad_axis(padded.n3);
   } else if(boundary == GpuDipoleBoundaryMode::Periodic2D) {
      std::size_t* extents[] = {&padded.n1, &padded.n2, &padded.n3};
      *extents[open_axis] = pad_axis(*extents[open_axis]);
   }
   return padded;
}

std::array<real, 9> GpuDipoleConvolutionDescriptor::fullCellMatrix() const {
   if(!c1 || !c2 || !c3) return {};
   const real scale[] = {static_cast<real>(atomistic_grid.n1),
                         static_cast<real>(atomistic_grid.n2),
                         static_cast<real>(atomistic_grid.n3)};
   return {scale[0] * c1[0], scale[0] * c1[1], scale[0] * c1[2],
           scale[1] * c2[0], scale[1] * c2[1], scale[1] * c2[2],
           scale[2] * c3[0], scale[2] * c3[1], scale[2] * c3[2]};
}

real GpuDipoleConvolutionDescriptor::cellVolume() const {
   const auto h = fullCellMatrix();
   return h[0] * (h[4] * h[8] - h[5] * h[7]) -
          h[1] * (h[3] * h[8] - h[5] * h[6]) +
          h[2] * (h[3] * h[7] - h[4] * h[6]);
}

std::array<real, 9> GpuDipoleConvolutionDescriptor::reciprocalCellMatrix() const {
   const auto h = fullCellMatrix();
   const real volume = cellVolume();
   if(volume == static_cast<real>(0)) return {};
   const real scale = static_cast<real>(2.0 * std::acos(-1.0)) / volume;
   // Columns are b1=2*pi*(a2 x a3)/V, b2=2*pi*(a3 x a1)/V,
   // b3=2*pi*(a1 x a2)/V, where H=[a1 a2 a3].
   return {scale * (h[4] * h[8] - h[5] * h[7]), scale * (h[5] * h[6] - h[3] * h[8]),
           scale * (h[3] * h[7] - h[4] * h[6]),
           scale * (h[7] * h[2] - h[8] * h[1]), scale * (h[8] * h[0] - h[6] * h[2]),
           scale * (h[6] * h[1] - h[7] * h[0]),
           scale * (h[1] * h[5] - h[2] * h[4]), scale * (h[2] * h[3] - h[0] * h[5]),
           scale * (h[0] * h[4] - h[1] * h[3])};
}

GpuDipoleFftLayout GpuDipoleConvolutionDescriptor::fftLayout() const {
   GpuDipoleFftLayout result{};
   result.real_grid = paddedGrid();
   if(!result.real_grid.valid()) return {};
   result.spectral_grid = result.real_grid;
   result.spectral_grid.n1 = result.real_grid.n1 / 2 + 1;
   result.real_cells = result.real_grid.cells();
   result.spectral_cells = result.spectral_grid.cells();
   if(!multiply(3, basis, result.field_batches) ||
      !multiply(result.field_batches, ensembles, result.field_batches) ||
      !multiply(9, basis, result.kernel_batches) ||
      !multiply(result.kernel_batches, basis, result.kernel_batches) || !result.valid()) return {};
   return result;
}

bool GpuDipoleConvolutionDescriptor::valid() const {
   if(basis == 0 || ensembles == 0 || !atomistic_grid.valid()) return false;
   if(discretization == GpuDipoleDiscretization::MacrospinGrid && !macro_grid.valid()) return false;
   if(boundary == GpuDipoleBoundaryMode::Periodic2D && open_axis > 2) return false;
   return hasNonSingularCell(*this) && fftLayout().fitsFftLibrary();
}

bool GpuDipoleConvolution::initiate(const GpuDipoleConvolutionDescriptor& descriptor,
                                    GPU_STREAM_T work_stream) {
   release();
   if(!descriptor.valid()) return false;
   desc = descriptor;
   layout = desc.fftLayout();
   stream = work_stream;

   // PlanMany receives the transform axes in slow-to-fast order, whereas the
   // source grid is stored with n1 as the contiguous axis.
   const int rank = 3;
   int n[] = {static_cast<int>(layout.real_grid.n3),
              static_cast<int>(layout.real_grid.n2),
              static_cast<int>(layout.real_grid.n1)};
   int inembed[] = {n[0], n[1], n[2]};
   int onembed[] = {n[0], n[1], static_cast<int>(layout.spectral_grid.n1)};
   const int idist = static_cast<int>(layout.real_cells);
   const int odist = static_cast<int>(layout.spectral_cells);
   std::size_t forward_workspace = 0, backward_workspace = 0, kernel_workspace = 0;

   try {
      // Mark the group before the first allocation so every partial failure
      // follows the same idempotent release path.
      allocated = true;
      moments_real.Allocate(static_cast<long int>(layout.real_cells),
                            static_cast<long int>(layout.field_batches));
      fields_real.Allocate(static_cast<long int>(layout.real_cells),
                           static_cast<long int>(layout.field_batches));
      moments_fft.Allocate(static_cast<long int>(layout.spectral_cells),
                          static_cast<long int>(layout.field_batches));
      fields_fft.Allocate(static_cast<long int>(layout.spectral_cells),
                         static_cast<long int>(layout.field_batches));
      kernel_fft.Allocate(static_cast<long int>(layout.spectral_cells),
                          static_cast<long int>(layout.kernel_batches));
      cell_vectors.Allocate(static_cast<long int>(3), static_cast<long int>(3));
      reciprocal_vectors.Allocate(static_cast<long int>(3), static_cast<long int>(3));
      cell_volume.Allocate(static_cast<long int>(1));
      const auto full_cell = desc.fullCellMatrix();
      const auto reciprocal_cell = desc.reciprocalCellMatrix();
      const real volume = desc.cellVolume();
      if(GPU_MEMCPY(cell_vectors.data(), full_cell.data(), 9 * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS ||
         GPU_MEMCPY(reciprocal_vectors.data(), reciprocal_cell.data(), 9 * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS ||
         GPU_MEMCPY(cell_volume.data(), &volume, sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS) {
         throw std::runtime_error("GPU dipole cell-matrix upload failed");
      }
      assertGpuFft(GPUFFT_CREATE(&forward_plan));
      forward_plan_created = true;
      assertGpuFft(GPUFFT_CREATE(&backward_plan));
      backward_plan_created = true;
      assertGpuFft(GPUFFT_CREATE(&kernel_plan));
      kernel_plan_created = true;
      assertGpuFft(GPUFFT_SET_AUTO_ALLOCATION(forward_plan, 0));
      assertGpuFft(GPUFFT_SET_AUTO_ALLOCATION(backward_plan, 0));
      assertGpuFft(GPUFFT_SET_AUTO_ALLOCATION(kernel_plan, 0));
      assertGpuFft(GPUFFT_MAKE_PLAN_MANY(forward_plan, rank, n, inembed, 1, idist,
                                          onembed, 1, odist, GPUFFT_R2C,
                                          static_cast<int>(layout.field_batches), &forward_workspace));
      assertGpuFft(GPUFFT_MAKE_PLAN_MANY(backward_plan, rank, n, onembed, 1, odist,
                                          inembed, 1, idist, GPUFFT_C2R,
                                          static_cast<int>(layout.field_batches), &backward_workspace));
      assertGpuFft(GPUFFT_MAKE_PLAN_MANY(kernel_plan, rank, n, inembed, 1, idist,
                                          onembed, 1, odist, GPUFFT_R2C,
                                          static_cast<int>(layout.kernel_batches), &kernel_workspace));
      const std::size_t workspace = std::max(forward_workspace, std::max(backward_workspace, kernel_workspace));
      if(workspace != 0) fft_workspace.Allocate(static_cast<long int>(workspace));
      assertGpuFft(GPUFFT_SET_WORK_AREA(forward_plan, fft_workspace.data()));
      assertGpuFft(GPUFFT_SET_WORK_AREA(backward_plan, fft_workspace.data()));
      assertGpuFft(GPUFFT_SET_WORK_AREA(kernel_plan, fft_workspace.data()));
      assertGpuFft(GPUFFT_SET_STREAM(forward_plan, stream));
      assertGpuFft(GPUFFT_SET_STREAM(backward_plan, stream));
      assertGpuFft(GPUFFT_SET_STREAM(kernel_plan, stream));
   } catch(...) {
      release();
      throw;
   }
   initiated = true;
   return true;
}

void GpuDipoleConvolution::release() {
   if(forward_plan_created) {
      GPUFFT_DESTROY(forward_plan);
      forward_plan_created = false;
   }
   if(backward_plan_created) {
      GPUFFT_DESTROY(backward_plan);
      backward_plan_created = false;
   }
   if(kernel_plan_created) {
      GPUFFT_DESTROY(kernel_plan);
      kernel_plan_created = false;
   }
   forward_plan = {};
   backward_plan = {};
   kernel_plan = {};
   if(allocated) {
      moments_real.Free();
      fields_real.Free();
      moments_fft.Free();
      fields_fft.Free();
      kernel_fft.Free();
      cell_vectors.Free();
      reciprocal_vectors.Free();
      cell_volume.Free();
      fft_workspace.Free();
      allocated = false;
   }
   initiated = false;
   desc = {};
   layout = {};
   stream = {};
}

GpuDipoleConvolution::~GpuDipoleConvolution() {
   release();
}

bool GpuDipoleConvolution::isInitiated() const {
   return initiated;
}

const GpuDipoleConvolutionDescriptor& GpuDipoleConvolution::descriptor() const {
   return desc;
}

const GpuDipoleFftLayout& GpuDipoleConvolution::fftLayout() const {
   return layout;
}

void GpuDipoleConvolution::packMacroMoments(const GpuTensor<real, 3>& macro_moments) {
   if(!initiated) throw std::runtime_error("GPU dipole FFT pack requested before initialization");
   const std::size_t expected = 3 * layout.real_cells * desc.basis * desc.ensembles;
   if(macro_moments.size() != expected) {
      throw std::runtime_error("GPU dipole macro-moment shape does not match the FFT grid");
   }
   constexpr unsigned int threads = 256;
   const std::size_t blocks_needed = (expected + threads - 1) / threads;
   const unsigned int blocks = static_cast<unsigned int>(std::min<std::size_t>(blocks_needed, 65535));
   pack_macro_moments_kernel<<<blocks, threads, 0, stream>>>(macro_moments.data(), moments_real.data(),
                                                               layout.real_cells, desc.basis, desc.ensembles);
   if(GPU_GET_LAST_ERROR() != GPU_SUCCESS) {
      throw std::runtime_error("GPU dipole macro-moment packing launch failed");
   }
}

void GpuDipoleConvolution::forwardTransformMoments() {
   if(!initiated) throw std::runtime_error("GPU dipole forward FFT requested before initialization");
   assertGpuFft(GPUFFT_EXEC_R2C(forward_plan, moments_real.data(), moments_fft.data()));
}

void GpuDipoleConvolution::applySpectralKernel() {
   if(!initiated) throw std::runtime_error("GPU dipole spectral contraction requested before initialization");
   const std::size_t total = layout.spectral_cells * layout.field_batches;
   constexpr unsigned int threads = 256;
   const std::size_t blocks_needed = (total + threads - 1) / threads;
   const unsigned int blocks = static_cast<unsigned int>(std::min<std::size_t>(blocks_needed, 65535));
   apply_spectral_kernel<<<blocks, threads, 0, stream>>>(moments_fft.data(), kernel_fft.data(), fields_fft.data(),
                                                           layout.spectral_cells, desc.basis, desc.ensembles);
   if(GPU_GET_LAST_ERROR() != GPU_SUCCESS) {
      throw std::runtime_error("GPU dipole spectral contraction launch failed");
   }
}

void GpuDipoleConvolution::inverseTransformFields() {
   if(!initiated) throw std::runtime_error("GPU dipole inverse FFT requested before initialization");
   assertGpuFft(GPUFFT_EXEC_C2R(backward_plan, fields_fft.data(), fields_real.data()));
}

std::size_t GpuDipoleConvolution::estimatePersistentBytes(
      const GpuDipoleConvolutionDescriptor& descriptor) {
   return descriptor.valid() ? descriptor.fftLayout().persistentBytes() : 0;
}

std::size_t GpuDipoleConvolution::estimateWorkspaceBytes(
      const GpuDipoleConvolutionDescriptor& descriptor) {
   if(!descriptor.valid()) return 0;
   const auto layout = descriptor.fftLayout();
   int n[] = {static_cast<int>(layout.real_grid.n3), static_cast<int>(layout.real_grid.n2),
              static_cast<int>(layout.real_grid.n1)};
   int inembed[] = {n[0], n[1], n[2]};
   int onembed[] = {n[0], n[1], static_cast<int>(layout.spectral_grid.n1)};
   const int idist = static_cast<int>(layout.real_cells);
   const int odist = static_cast<int>(layout.spectral_cells);
   std::size_t forward = 0, backward = 0, kernel = 0;
   if(GPUFFT_ESTIMATE_MANY(3, n, inembed, 1, idist, onembed, 1, odist, GPUFFT_R2C,
                           static_cast<int>(layout.field_batches), &forward) != GPUFFT_SUCCESS ||
      GPUFFT_ESTIMATE_MANY(3, n, onembed, 1, odist, inembed, 1, idist, GPUFFT_C2R,
                           static_cast<int>(layout.field_batches), &backward) != GPUFFT_SUCCESS ||
      GPUFFT_ESTIMATE_MANY(3, n, inembed, 1, idist, onembed, 1, odist, GPUFFT_R2C,
                           static_cast<int>(layout.kernel_batches), &kernel) != GPUFFT_SUCCESS) return 0;
   return std::max(forward, std::max(backward, kernel));
}

std::size_t GpuDipoleConvolution::estimateBytes(const GpuDipoleConvolutionDescriptor& descriptor) {
   const std::size_t buffers = estimatePersistentBytes(descriptor);
   const std::size_t workspace = estimateWorkspaceBytes(descriptor);
   std::size_t total = 0;
   return buffers != 0 && add(buffers, workspace, total) ? total : 0;
}
