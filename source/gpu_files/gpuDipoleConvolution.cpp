#include "gpuDipoleConvolution.hpp"

#include "gpuFftWrapper.hpp"
#include "real_type.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

constexpr double bohr_magneton_si = 9.274009994e-24;

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
   return std::isfinite(descriptor.cellVolume()) && descriptor.cellVolume() > static_cast<real>(0);
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

__global__ void add_point_self_field_kernel(const real* moments, real* fields, std::size_t count,
                                             real prefactor) {
   for(std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < count; index += static_cast<std::size_t>(blockDim.x) * gridDim.x) {
      fields[index] += prefactor * moments[index];
   }
}

__global__ void add_real_ewald_kernel(const real* moments, real* fields, const real* centres,
   const real* h, std::size_t cells, unsigned int ensembles, real alpha, real cutoff, unsigned int extent) {
   const std::size_t total = cells * ensembles;
   for(std::size_t item = static_cast<std::size_t>(blockIdx.x)*blockDim.x+threadIdx.x; item<total;
       item += static_cast<std::size_t>(blockDim.x)*gridDim.x) {
      const std::size_t target=item%cells, ens=item/cells; real bx=0,by=0,bz=0;
      const real tx=centres[3*target],ty=centres[3*target+1],tz=centres[3*target+2];
      for(std::size_t source=0;source<cells;++source) for(int iz=-int(extent);iz<=int(extent);++iz)
      for(int iy=-int(extent);iy<=int(extent);++iy) for(int ix=-int(extent);ix<=int(extent);++ix) {
         if(target==source && ix==0 && iy==0 && iz==0) continue;
         const real x=tx-(centres[3*source]+ix*h[0]+iy*h[3]+iz*h[6]);
         const real y=ty-(centres[3*source+1]+ix*h[1]+iy*h[4]+iz*h[7]);
         const real z=tz-(centres[3*source+2]+ix*h[2]+iy*h[5]+iz*h[8]);
         const real r2=x*x+y*y+z*z; if(r2>cutoff*cutoff || r2==0) continue;
         const real r=sqrt(r2), g=exp(-alpha*alpha*r2), e=erfc(alpha*r), pi=sqrt(acos((real)-1));
         const real a=3*e/(r2*r2*r)+6*alpha*g/(pi*r2*r2)+4*alpha*alpha*alpha*g/(pi*r2);
         const real b=-e/(r2*r)-2*alpha*g/(pi*r2);
         const std::size_t base=3*(source+cells*ens); const real mx=moments[base],my=moments[base+1],mz=moments[base+2];
         const real dot=x*mx+y*my+z*mz; bx+=a*x*dot+b*mx; by+=a*y*dot+b*my; bz+=a*z*dot+b*mz;
      }
      const std::size_t base=3*(target+cells*ens); fields[base]+=bx; fields[base+1]+=by; fields[base+2]+=bz;
   }
}

__global__ void scatter_point_fields_kernel(const real* fields, real* beff, const unsigned int* cell_index,
                                             std::size_t cells, std::size_t atoms, unsigned int ensembles) {
   const std::size_t total=atoms*ensembles;
   for(std::size_t item=static_cast<std::size_t>(blockIdx.x)*blockDim.x+threadIdx.x;item<total;
       item+=static_cast<std::size_t>(blockDim.x)*gridDim.x) {
      const std::size_t atom=item%atoms, ens=item/atoms; const unsigned int one=cell_index[atom];
      if(one==0 || one>cells) continue; const std::size_t cell=one-1;
      const std::size_t out=3*(atom+atoms*ens);
      beff[out]+=fields[cell+cells*(3*ens)]; beff[out+1]+=fields[cell+cells*(1+3*ens)];
      beff[out+2]+=fields[cell+cells*(2+3*ens)];
   }
}

__global__ void reduce_point_energy_kernel(const real* moments, const real* fields, real* energy,
                                           std::size_t cells, unsigned int ensembles) {
   const std::size_t total=cells*ensembles;
   for(std::size_t item=static_cast<std::size_t>(blockIdx.x)*blockDim.x+threadIdx.x;item<total;
       item+=static_cast<std::size_t>(blockDim.x)*gridDim.x) {
      const std::size_t cell=item%cells, ens=item/cells, base=3*(cell+cells*ens);
      atomicAdd(&energy[ens], static_cast<real>(-0.5)*(moments[base]*fields[base]+moments[base+1]*fields[base+1]+moments[base+2]*fields[base+2]));
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
      !bytesFor(real_fields, sizeof(real), part) || !add(bytes, part, bytes) ||
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
   if(discretization == GpuDipoleDiscretization::MacrospinGrid && macro_count != 0) {
      std::size_t expected_count = 0;
      if(!multiply(macro_grid.cells(), basis, expected_count) || macro_count != expected_count) return false;
   }
   if(!std::isfinite(alat) || alat <= 0.0 || !std::isfinite(tolerance) || tolerance <= 0.0 ||
      !std::isfinite(field_prefactor) || field_prefactor <= 0.0) return false;
   return sizeof(real) == sizeof(double) && hasNonSingularCell(*this) && fftLayout().fitsFftLibrary();
}

bool makeEwald3dFftDipoleDescriptor(const GpuDipoleDescriptorInput& input,
                                    GpuDipoleConvolutionDescriptor& descriptor) {
   descriptor = {};
   if(input.boundaries[0] != 'P' || input.boundaries[1] != 'P' || input.boundaries[2] != 'P') return false;
   descriptor.atomistic_grid = input.atomistic_grid;
   descriptor.macro_grid = input.macro_grid;
   descriptor.basis = input.basis;
   descriptor.ensembles = input.ensembles;
   descriptor.boundary = GpuDipoleBoundaryMode::Periodic3D;
   descriptor.discretization = GpuDipoleDiscretization::MacrospinGrid;
   descriptor.c1 = input.c1;
   descriptor.c2 = input.c2;
   descriptor.c3 = input.c3;
   descriptor.macro_centers = input.macro_centers;
   descriptor.macro_count = input.macro_count;
   descriptor.alat = input.alat;
   descriptor.tolerance = input.tolerance;
   descriptor.field_prefactor = 1.0e-7 * bohr_magneton_si /
      (input.alat * input.alat * input.alat);
   return descriptor.valid();
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
      point_energy.Allocate(static_cast<long int>(desc.ensembles));
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
      point_energy.Free();
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

void GpuDipoleConvolution::buildReciprocalEwaldKernel(real alpha) {
   if(!initiated) throw std::runtime_error("GPU dipole Ewald kernel requested before initialization");
   if(desc.boundary != GpuDipoleBoundaryMode::Periodic3D || desc.basis != 1 || alpha <= static_cast<real>(0)) {
      throw std::runtime_error("reciprocal Ewald kernel currently requires EWALD3D_FFT, NA=1, and positive alpha");
   }
   const real volume = desc.cellVolume();
   if(volume <= static_cast<real>(0)) throw std::runtime_error("PME requires a right-handed positive-volume cell");
   const auto reciprocal = desc.reciprocalCellMatrix();
   const GpuFftComplex zero{};
   std::vector<GpuFftComplex> kernel(layout.spectral_cells * layout.kernel_batches, zero);
   const auto signed_index = [](std::size_t index, std::size_t extent) -> long long {
      return index <= extent / 2 ? static_cast<long long>(index) :
             static_cast<long long>(index) - static_cast<long long>(extent);
   };
   const real pi = static_cast<real>(std::acos(-1.0));
   for(std::size_t k3 = 0; k3 < layout.spectral_grid.n3; ++k3) {
      for(std::size_t k2 = 0; k2 < layout.spectral_grid.n2; ++k2) {
         for(std::size_t k1 = 0; k1 < layout.spectral_grid.n1; ++k1) {
            const long long n1 = static_cast<long long>(k1);
            const long long n2 = signed_index(k2, layout.real_grid.n2);
            const long long n3 = signed_index(k3, layout.real_grid.n3);
            if(n1 == 0 && n2 == 0 && n3 == 0) continue; // tin-foil k=0
            const real wave[] = {
               reciprocal[0] * n1 + reciprocal[3] * n2 + reciprocal[6] * n3,
               reciprocal[1] * n1 + reciprocal[4] * n2 + reciprocal[7] * n3,
               reciprocal[2] * n1 + reciprocal[5] * n2 + reciprocal[8] * n3};
            const real wave2 = wave[0] * wave[0] + wave[1] * wave[1] + wave[2] * wave[2];
            const real prefactor = -static_cast<real>(4) * pi * std::exp(-wave2 / (static_cast<real>(4) * alpha * alpha)) /
                                   (volume * wave2);
            const std::size_t spectral = k1 + layout.spectral_grid.n1 * (k2 + layout.spectral_grid.n2 * k3);
            for(unsigned int row = 0; row < 3; ++row) {
               for(unsigned int column = 0; column < 3; ++column) {
                  const std::size_t batch = row + 3 * column;
                  kernel[spectral + layout.spectral_cells * batch].x = prefactor * wave[row] * wave[column];
               }
            }
         }
      }
   }
   if(GPU_MEMCPY(kernel_fft.data(), kernel.data(), kernel.size() * sizeof(GpuFftComplex),
                 GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS) {
      throw std::runtime_error("GPU reciprocal Ewald kernel upload failed");
   }
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

void GpuDipoleConvolution::addPointSelfField(real alpha) {
   if(!initiated) throw std::runtime_error("GPU dipole self field requested before initialization");
   if(alpha <= static_cast<real>(0)) throw std::runtime_error("GPU dipole self field requires positive alpha");
   const std::size_t count = layout.real_cells * layout.field_batches;
   constexpr unsigned int threads = 256;
   const std::size_t blocks_needed = (count + threads - 1) / threads;
   const unsigned int blocks = static_cast<unsigned int>(std::min<std::size_t>(blocks_needed, 65535));
   const real prefactor = static_cast<real>(4) * alpha * alpha * alpha /
                          (static_cast<real>(3) * static_cast<real>(std::sqrt(std::acos(-1.0))));
   add_point_self_field_kernel<<<blocks, threads, 0, stream>>>(moments_real.data(), fields_real.data(), count, prefactor);
   if(GPU_GET_LAST_ERROR() != GPU_SUCCESS) {
      throw std::runtime_error("GPU dipole self-field launch failed");
   }
}

void GpuDipoleConvolution::addRealSpaceField(real alpha, real cutoff, unsigned int image_extent) {
   if(!initiated || desc.basis != 1 || !desc.macro_centers || desc.macro_count != layout.real_cells || alpha<=0 || cutoff<=0)
      throw std::runtime_error("real Ewald primitive requires NA=1 centres, positive alpha and cutoff");
   constexpr unsigned int threads=128; const std::size_t total=layout.real_cells*desc.ensembles;
   const unsigned int blocks=static_cast<unsigned int>(std::min<std::size_t>((total+threads-1)/threads,65535));
   add_real_ewald_kernel<<<blocks,threads,0,stream>>>(moments_real.data(),fields_real.data(),desc.macro_centers,
      cell_vectors.data(),layout.real_cells,desc.ensembles,alpha,cutoff,image_extent);
   if(GPU_GET_LAST_ERROR()!=GPU_SUCCESS) throw std::runtime_error("GPU real Ewald launch failed");
}

void GpuDipoleConvolution::evaluatePointEwald(const GpuTensor<real, 3>& macro_moments, real alpha,
                                              real cutoff, unsigned int image_extent) {
   if(desc.basis != 1 || desc.boundary != GpuDipoleBoundaryMode::Periodic3D) {
      throw std::runtime_error("point Ewald evaluation currently requires EWALD3D_FFT with NA=1");
   }
   buildReciprocalEwaldKernel(alpha);
   packMacroMoments(macro_moments);
   forwardTransformMoments();
   applySpectralKernel();
   inverseTransformFields();
   addRealSpaceField(alpha, cutoff, image_extent);
   addPointSelfField(alpha);
   reducePointEwaldEnergy();
}

void GpuDipoleConvolution::scatterPointFields(GpuTensor<real, 3>& beff, const unsigned int* one_based_cell_index,
                                              std::size_t atom_count) {
   if(!initiated || desc.basis!=1 || !one_based_cell_index || beff.size()!=3*atom_count*desc.ensembles)
      throw std::runtime_error("invalid NA=1 dipole field scatter");
   constexpr unsigned int threads=256; const std::size_t total=atom_count*desc.ensembles;
   const unsigned int blocks=static_cast<unsigned int>(std::min<std::size_t>((total+threads-1)/threads,65535));
   scatter_point_fields_kernel<<<blocks,threads,0,stream>>>(fields_real.data(),beff.data(),one_based_cell_index,
      layout.real_cells,atom_count,desc.ensembles);
   if(GPU_GET_LAST_ERROR()!=GPU_SUCCESS) throw std::runtime_error("GPU dipole field scatter launch failed");
}

void GpuDipoleConvolution::reducePointEwaldEnergy() {
   if(!initiated || desc.basis!=1) throw std::runtime_error("point Ewald energy requires initialized NA=1 grid");
   GPU_MEMSET_ASYNC(point_energy.data(),0,desc.ensembles*sizeof(real),stream);
   constexpr unsigned int threads=256; const std::size_t total=layout.real_cells*desc.ensembles;
   const unsigned int blocks=static_cast<unsigned int>(std::min<std::size_t>((total+threads-1)/threads,65535));
   reduce_point_energy_kernel<<<blocks,threads,0,stream>>>(moments_real.data(),fields_real.data(),point_energy.data(),layout.real_cells,desc.ensembles);
   if(GPU_GET_LAST_ERROR()!=GPU_SUCCESS) throw std::runtime_error("GPU point Ewald energy reduction launch failed");
}

std::vector<real> GpuDipoleConvolution::pointEwaldEnergies() const {
   if(!initiated) throw std::runtime_error("point Ewald energy requested before initialization");
   if(GPU_STREAM_SYNC(stream) != GPU_SUCCESS) {
      throw std::runtime_error("GPU point Ewald energy stream synchronization failed");
   }
   std::vector<real> result(desc.ensembles);
   if(GPU_MEMCPY(result.data(),point_energy.data(),result.size()*sizeof(real),GPU_MEMCPY_DEVICE_TO_HOST)!=GPU_SUCCESS)
      throw std::runtime_error("GPU point Ewald energy download failed");
   return result;
}

std::vector<real> GpuDipoleConvolution::pointEwaldFields() const {
   if(!initiated || desc.basis != 1) throw std::runtime_error("point Ewald field requested from an invalid grid");
   if(GPU_STREAM_SYNC(stream) != GPU_SUCCESS) throw std::runtime_error("GPU point Ewald field stream synchronization failed");
   std::vector<real> packed(layout.real_cells * layout.field_batches);
   if(GPU_MEMCPY(packed.data(), fields_real.data(), packed.size()*sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS)
      throw std::runtime_error("GPU point Ewald field download failed");
   std::vector<real> result(3 * layout.real_cells * desc.ensembles);
   for(unsigned int ensemble=0; ensemble<desc.ensembles; ++ensemble)
      for(std::size_t cell=0; cell<layout.real_cells; ++cell)
         for(unsigned int component=0; component<3; ++component)
            result[component + 3*(cell + layout.real_cells*ensemble)] =
               packed[cell + layout.real_cells*(component + 3*ensemble)];
   return result;
}

std::size_t GpuDipoleConvolution::estimatePersistentBytes(
      const GpuDipoleConvolutionDescriptor& descriptor) {
   if(!descriptor.valid()) return 0;
   std::size_t energy_bytes=0,total=0;
   return bytesFor(descriptor.ensembles,sizeof(real),energy_bytes) &&
          add(descriptor.fftLayout().persistentBytes(),energy_bytes,total) ? total : 0;
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
   const std::size_t construction = descriptor.fftLayout().constructionBytes();
   std::size_t total = 0;
   return buffers != 0 && construction != 0 && add(buffers, workspace, total) &&
          add(total, construction, total) ? total : 0;
}
