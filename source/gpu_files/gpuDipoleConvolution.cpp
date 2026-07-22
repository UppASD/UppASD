#include "gpuDipoleConvolution.hpp"

#include "gpuFftWrapper.hpp"
#include "real_type.h"

#include <algorithm>
#include <climits>
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
      !bytesFor(kernel, sizeof(GpuFftComplex), part) || !add(bytes, part, bytes)) return 0;
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
   return fftLayout().fitsFftLibrary();
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

std::size_t GpuDipoleConvolution::estimatePersistentBytes(
      const GpuDipoleConvolutionDescriptor& descriptor) {
   return descriptor.valid() ? descriptor.fftLayout().persistentBytes() : 0;
}
