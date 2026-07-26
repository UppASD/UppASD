#include "gpuDipoleConvolution.hpp"

#include "gpuFftWrapper.hpp"
#include "real_type.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
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

struct HostComplex {
   double real = 0.0;
   double imag = 0.0;
};

std::size_t spectralIndex(const GpuDipoleFftLayout& layout, std::size_t q1, std::size_t q2, std::size_t q3) {
   return q1 + layout.spectral_grid.n1 * (q2 + layout.real_grid.n2 * q3);
}

HostComplex fullSpectrumValue(const std::vector<GpuFftComplex>& spectrum, const GpuDipoleFftLayout& layout,
                              std::size_t q1, std::size_t q2, std::size_t q3, std::size_t batch) {
   bool conjugate = false;
   if(q1 >= layout.spectral_grid.n1) {
      q1 = layout.real_grid.n1 - q1;
      q2 = (layout.real_grid.n2 - q2) % layout.real_grid.n2;
      q3 = (layout.real_grid.n3 - q3) % layout.real_grid.n3;
      conjugate = true;
   }
   const auto value = spectrum[spectralIndex(layout, q1, q2, q3) + layout.spectral_cells * batch];
   return {value.x, conjugate ? -value.y : value.y};
}

double completeKernelReciprocityError(const std::vector<double>& kernel, const GpuDipoleFftLayout& layout,
                                      unsigned int basis) {
   double maximum = 0.0;
   for(std::size_t d3 = 0; d3 < layout.real_grid.n3; ++d3) for(std::size_t d2 = 0; d2 < layout.real_grid.n2; ++d2)
      for(std::size_t d1 = 0; d1 < layout.real_grid.n1; ++d1) {
         const std::size_t cell = d1 + layout.real_grid.n1 * (d2 + layout.real_grid.n2 * d3);
         const std::size_t reverse = (layout.real_grid.n1 - d1) % layout.real_grid.n1 +
            layout.real_grid.n1 * ((layout.real_grid.n2 - d2) % layout.real_grid.n2 +
            layout.real_grid.n2 * ((layout.real_grid.n3 - d3) % layout.real_grid.n3));
         for(unsigned int target = 0; target < basis; ++target) for(unsigned int source = 0; source < basis; ++source)
            for(unsigned int row = 0; row < 3; ++row) for(unsigned int column = 0; column < 3; ++column) {
               const std::size_t batch = row + 3 * (column + 3 * (target + basis * source));
               const std::size_t paired = column + 3 * (row + 3 * (source + basis * target));
               maximum = std::max(maximum, std::abs(kernel[cell + layout.real_cells * batch] -
                                                     kernel[reverse + layout.real_cells * paired]));
            }
      }
   return maximum;
}

GpuDipoleSpectrumDiagnostics validateSpectrum(const std::vector<double>& kernel,
                                              const std::vector<GpuFftComplex>& spectrum,
                                              const GpuDipoleFftLayout& layout, unsigned int basis) {
   GpuDipoleSpectrumDiagnostics result{};
   result.max_reciprocity_error = completeKernelReciprocityError(kernel, layout, basis);
   for(const auto& value : spectrum) {
      if(!std::isfinite(value.x) || !std::isfinite(value.y))
         throw std::runtime_error("complete periodic dipole spectrum contains a non-finite value");
   }
   for(std::size_t q3 = 0; q3 < layout.real_grid.n3; ++q3) for(std::size_t q2 = 0; q2 < layout.real_grid.n2; ++q2)
      for(std::size_t q1 = 0; q1 < layout.real_grid.n1; ++q1) {
         const std::size_t r1 = (layout.real_grid.n1 - q1) % layout.real_grid.n1;
         const std::size_t r2 = (layout.real_grid.n2 - q2) % layout.real_grid.n2;
         const std::size_t r3 = (layout.real_grid.n3 - q3) % layout.real_grid.n3;
         for(std::size_t batch = 0; batch < layout.kernel_batches; ++batch) {
            const auto value = fullSpectrumValue(spectrum, layout, q1, q2, q3, batch);
            const auto reverse = fullSpectrumValue(spectrum, layout, r1, r2, r3, batch);
            result.max_conjugacy_error = std::max(result.max_conjugacy_error,
               std::max(std::abs(reverse.real - value.real), std::abs(reverse.imag + value.imag)));
         }
         for(unsigned int target = 0; target < basis; ++target) for(unsigned int source = 0; source < basis; ++source)
            for(unsigned int row = 0; row < 3; ++row) for(unsigned int column = 0; column < 3; ++column) {
               const std::size_t batch = row + 3 * (column + 3 * (target + basis * source));
               const std::size_t paired = column + 3 * (row + 3 * (source + basis * target));
               const auto value = fullSpectrumValue(spectrum, layout, q1, q2, q3, batch);
               const auto transpose = fullSpectrumValue(spectrum, layout, q1, q2, q3, paired);
               result.max_hermitian_error = std::max(result.max_hermitian_error,
                  std::max(std::abs(value.real - transpose.real), std::abs(value.imag + transpose.imag)));
            }
      }
   return result;
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

// cuFFT/hipFFT leave C2R unnormalised.  The normalization contract is kept
// explicit here, at kernel construction, so raw C2R applies K exactly once.
__global__ void normalize_kernel_spectrum(GpuFftComplex* kernel, std::size_t count, real inverse_grid) {
   for(std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < count; index += static_cast<std::size_t>(blockDim.x) * gridDim.x) {
      kernel[index].x *= inverse_grid;
      kernel[index].y *= inverse_grid;
   }
}

// macro_cell_index retains Fortran's one-based map.  WP5 admits precisely one
// atom per macro cell and one basis channel, but retaining the map here makes
// the scatter independent of atom ordering and keeps the invariant explicit.
__global__ void add_dipole_fields(const real* fields, real* beff, real* eneff,
                                  const unsigned int* macro_cell_index,
                                  std::size_t cells, std::size_t atoms,
                                  unsigned int ensembles, real prefactor,
                                  bool shared_field_storage) {
   const std::size_t total = 3 * atoms * ensembles;
   for(std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < total; index += static_cast<std::size_t>(blockDim.x) * gridDim.x) {
      const unsigned int component = static_cast<unsigned int>(index % 3);
      const std::size_t atom_ensemble = index / 3;
      const std::size_t ensemble = atom_ensemble / atoms;
      const std::size_t atom = atom_ensemble % atoms;
      const unsigned int one_based_cell = macro_cell_index[atom];
      // Initialization validates every map entry, so this is defensive only.
      if(one_based_cell == 0 || one_based_cell > cells) continue;
      const std::size_t cell = one_based_cell - 1;
      const std::size_t field_index = cell + cells * (component + 3 * ensemble);
      const real value = prefactor * fields[field_index];
      beff[index] += value;
      if(!shared_field_storage) eneff[index] += value;
   }
}

// energyM is [ensemble, energy_column]; columns 5 and 6 are total and Dip.
// The physical mu_B/mRy conversion remains in the existing measurement path.
__global__ void accumulate_dipole_energy(const real* moments, const real* fields,
                                         real* energy, std::size_t cells,
                                         unsigned int ensembles, real prefactor,
                                         real atoms_inverse) {
   const std::size_t total = cells * 3 * ensembles;
   for(std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < total; index += static_cast<std::size_t>(blockDim.x) * gridDim.x) {
      const std::size_t batch = index / cells;
      const unsigned int ensemble = static_cast<unsigned int>(batch / 3);
      const real value = static_cast<real>(-0.5) * moments[index] * (prefactor * fields[index]) * atoms_inverse;
      atomicAdd(&energy[ensemble + ensembles * 5], value);
      atomicAdd(&energy[ensemble + ensembles * 6], value);
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
   kernel_ready = false;
   spectrum_diagnostics = {};

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
      if(kernel_real_allocated) {
         kernel_real.Free();
         kernel_real_allocated = false;
      }
      cell_vectors.Free();
      reciprocal_vectors.Free();
      cell_volume.Free();
      fft_workspace.Free();
      allocated = false;
   }
   initiated = false;
   kernel_ready = false;
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

bool GpuDipoleConvolution::kernelReady() const {
   return kernel_ready;
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

void GpuDipoleConvolution::uploadCompleteKernelForTesting(const std::vector<double>& complete_kernel,
                                                           bool validate_physics) {
   if(!initiated) throw std::runtime_error("GPU dipole kernel upload requested before plan initialization");
   if(sizeof(real) != sizeof(double)) throw std::runtime_error("complete-kernel GPU validation requires fp64 storage");
   const std::size_t expected = layout.real_cells * layout.kernel_batches;
   if(complete_kernel.size() != expected) throw std::runtime_error("complete host tensor does not match GPU kernel layout");

   kernel_ready = false;
   spectrum_diagnostics = {};
   try {
      for(const double value : complete_kernel) {
         if(!std::isfinite(value)) throw std::runtime_error("complete periodic dipole tensor contains a non-finite value");
      }
      kernel_real.Allocate(static_cast<long int>(layout.real_cells), static_cast<long int>(layout.kernel_batches));
      kernel_real_allocated = true;
      if(GPU_MEMCPY(kernel_real.data(), complete_kernel.data(), expected * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS)
         throw std::runtime_error("complete periodic dipole tensor upload failed");
      assertGpuFft(GPUFFT_EXEC_R2C(kernel_plan, kernel_real.data(), kernel_fft.data()));
      constexpr unsigned int threads = 256;
      const std::size_t count = layout.spectral_cells * layout.kernel_batches;
      const unsigned int blocks = static_cast<unsigned int>(std::min<std::size_t>((count + threads - 1) / threads, 65535));
      normalize_kernel_spectrum<<<blocks, threads, 0, stream>>>(kernel_fft.data(), count,
                                                                  static_cast<real>(1) / static_cast<real>(layout.real_cells));
      if(GPU_GET_LAST_ERROR() != GPU_SUCCESS) throw std::runtime_error("GPU dipole spectrum-normalization launch failed");
      // A kernel cannot become observable until its R2C and normalization are complete.
      if(GPU_STREAM_SYNC(stream) != GPU_SUCCESS) throw std::runtime_error("GPU dipole kernel construction synchronization failed");
      std::vector<GpuFftComplex> spectrum(layout.spectral_cells * layout.kernel_batches);
      if(GPU_MEMCPY(spectrum.data(), kernel_fft.data(), spectrum.size() * sizeof(GpuFftComplex),
                    GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS) {
         throw std::runtime_error("GPU dipole spectrum diagnostic download failed");
      }
      spectrum_diagnostics = validateSpectrum(complete_kernel, spectrum, layout, desc.basis);
      constexpr double validation_tolerance = 1.0e-12;
      if(validate_physics && (spectrum_diagnostics.max_reciprocity_error > validation_tolerance ||
                               spectrum_diagnostics.max_conjugacy_error > validation_tolerance ||
                               spectrum_diagnostics.max_hermitian_error > validation_tolerance)) {
         throw std::runtime_error("complete periodic dipole kernel violates fp64 validation: reciprocity=" +
                                  std::to_string(spectrum_diagnostics.max_reciprocity_error) +
                                  " conjugacy=" + std::to_string(spectrum_diagnostics.max_conjugacy_error) +
                                  " hermitian=" + std::to_string(spectrum_diagnostics.max_hermitian_error));
      }
      kernel_real.Free();
      kernel_real_allocated = false;
      kernel_ready = true;
   } catch(...) {
      if(kernel_real_allocated) {
         kernel_real.Free();
         kernel_real_allocated = false;
      }
      throw;
   }
}

void GpuDipoleConvolution::uploadCompleteKernelForTesting(const DipoleKernelBuildResult& complete_kernel) {
   uploadCompleteKernelForTesting(complete_kernel.kernel);
}

void GpuDipoleConvolution::buildPeriodicKernel() {
   if(!initiated) throw std::runtime_error("GPU periodic dipole kernel requested before plan initialization");
   if(desc.boundary != GpuDipoleBoundaryMode::Periodic3D || desc.discretization != GpuDipoleDiscretization::MacrospinGrid ||
      desc.basis != 1 || sizeof(real) != sizeof(double)) {
      throw std::runtime_error("GPU production dipole construction is limited to fp64 NA=1 periodic macrospin grids");
   }
   DipolePeriodicGeometry geometry{};
   const auto h = desc.fullCellMatrix();
   const auto brec = desc.reciprocalCellMatrix();
   for(unsigned int index = 0; index < 9; ++index) {
      geometry.H[index] = static_cast<double>(h[index]);
      geometry.Brec[index] = static_cast<double>(brec[index]);
   }
   geometry.volume = static_cast<double>(desc.cellVolume());
   const auto grid = desc.activeGrid();
   geometry.grid = {grid.n1, grid.n2, grid.n3};
   geometry.basis = 1;
   geometry.basis_offsets = {{0.0, 0.0, 0.0}};
   DipoleKernelSettings settings{};
   settings.tolerance = desc.tolerance;
   const auto complete = buildPeriodicEwaldDisplacementKernel(geometry, settings);
   uploadCompleteKernelForTesting(complete);
}

void GpuDipoleConvolution::evaluate(const GpuTensor<real, 3>& macro_moments) {
   if(!initiated) throw std::runtime_error("GPU dipole evaluation requested before plan initialization");
   if(!kernel_ready) throw std::runtime_error("GPU dipole evaluation requested before complete kernel is ready");
   packMacroMoments(macro_moments);
   forwardTransformMoments();
   applySpectralKernel();
   inverseTransformFields();
}

void GpuDipoleConvolution::addFieldsToAtoms(GpuTensor<real, 3>& beff, GpuTensor<real, 3>& eneff,
                                             const unsigned int* macro_cell_index,
                                             std::size_t atom_count) {
   if(!initiated || !kernel_ready) throw std::runtime_error("GPU dipole field scatter requested before a ready kernel");
   if(desc.basis != 1 || layout.real_cells != atom_count || !macro_cell_index ||
      beff.size() != 3 * atom_count * desc.ensembles || eneff.size() != beff.size()) {
      throw std::runtime_error("GPU dipole production scatter does not match the accepted NA=1 block-one layout");
   }
   constexpr unsigned int threads = 256;
   const std::size_t count = 3 * atom_count * desc.ensembles;
   const unsigned int blocks = static_cast<unsigned int>(std::min<std::size_t>((count + threads - 1) / threads, 65535));
   add_dipole_fields<<<blocks, threads, 0, stream>>>(fields_real.data(), beff.data(), eneff.data(),
                                                       macro_cell_index, layout.real_cells, atom_count,
                                                       desc.ensembles, static_cast<real>(desc.field_prefactor),
                                                       beff.data() == eneff.data());
   if(GPU_GET_LAST_ERROR() != GPU_SUCCESS) throw std::runtime_error("GPU dipole field-scatter launch failed");
}

void GpuDipoleConvolution::accumulateEnergy(GpuTensor<real, 2>& energyM, std::size_t atom_count) {
   if(!initiated || !kernel_ready) throw std::runtime_error("GPU dipole energy requested before a ready kernel");
   if(desc.basis != 1 || layout.real_cells != atom_count || energyM.extent(0) != desc.ensembles || energyM.extent(1) < 7) {
      throw std::runtime_error("GPU dipole energy buffer does not match the accepted NA=1 block-one layout");
   }
   constexpr unsigned int threads = 256;
   const std::size_t count = layout.real_cells * 3 * desc.ensembles;
   const unsigned int blocks = static_cast<unsigned int>(std::min<std::size_t>((count + threads - 1) / threads, 65535));
   accumulate_dipole_energy<<<blocks, threads, 0, stream>>>(moments_real.data(), fields_real.data(), energyM.data(),
                                                              layout.real_cells, desc.ensembles,
                                                              static_cast<real>(desc.field_prefactor),
                                                              static_cast<real>(1) / static_cast<real>(atom_count));
   if(GPU_GET_LAST_ERROR() != GPU_SUCCESS) throw std::runtime_error("GPU dipole energy-reduction launch failed");
   // Acceptance-only trace seam.  Normal field evaluation never enters this
   // branch and therefore remains free of host synchronization.
   if(std::getenv("UPPASD_GPU_TEST_TRACE_DIP_ENERGY")) {
      if(GPU_STREAM_SYNC(stream) != GPU_SUCCESS)
         throw std::runtime_error("GPU dipole energy trace synchronization failed");
      std::vector<real> trace(desc.ensembles * 7);
      if(GPU_MEMCPY(trace.data(), energyM.data(), trace.size() * sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS)
         throw std::runtime_error("GPU dipole energy trace download failed");
      for(unsigned int ensemble = 0; ensemble < desc.ensembles; ++ensemble)
         std::printf("GPU_TEST_DIP_ENERGY ensemble=%u total=%.17e dip=%.17e\n", ensemble,
                     static_cast<double>(trace[ensemble + desc.ensembles * 5]),
                     static_cast<double>(trace[ensemble + desc.ensembles * 6]));
   }
}

std::vector<real> GpuDipoleConvolution::diagnosticFieldsForTesting() const {
   if(!initiated || !kernel_ready) throw std::runtime_error("GPU dipole diagnostic field requested before a ready evaluation");
   if(GPU_STREAM_SYNC(stream) != GPU_SUCCESS) throw std::runtime_error("GPU dipole diagnostic field synchronization failed");
   std::vector<real> packed(layout.real_cells * layout.field_batches);
   if(GPU_MEMCPY(packed.data(), fields_real.data(), packed.size() * sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS)
      throw std::runtime_error("GPU dipole diagnostic field download failed");
   std::vector<real> result(3 * layout.real_cells * desc.basis * desc.ensembles);
   for(unsigned int ensemble = 0; ensemble < desc.ensembles; ++ensemble)
      for(std::size_t cell = 0; cell < layout.real_cells; ++cell)
         for(unsigned int basis = 0; basis < desc.basis; ++basis)
            for(unsigned int component = 0; component < 3; ++component)
               result[component + 3 * (basis + desc.basis * (cell + layout.real_cells * ensemble))] =
                  packed[cell + layout.real_cells * (component + 3 * (basis + desc.basis * ensemble))];
   return result;
}

std::vector<real> GpuDipoleConvolution::diagnosticEnergiesForTesting() const {
   if(!initiated || !kernel_ready) throw std::runtime_error("GPU dipole diagnostic energy requested before a ready evaluation");
   if(GPU_STREAM_SYNC(stream) != GPU_SUCCESS) throw std::runtime_error("GPU dipole diagnostic energy synchronization failed");
   std::vector<real> moments(layout.real_cells * layout.field_batches);
   std::vector<real> fields(layout.real_cells * layout.field_batches);
   if(GPU_MEMCPY(moments.data(), moments_real.data(), moments.size() * sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS ||
      GPU_MEMCPY(fields.data(), fields_real.data(), fields.size() * sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS)
      throw std::runtime_error("GPU dipole diagnostic energy download failed");
   std::vector<real> result(desc.ensembles, static_cast<real>(0));
   for(unsigned int ensemble = 0; ensemble < desc.ensembles; ++ensemble)
      for(std::size_t cell = 0; cell < layout.real_cells; ++cell)
         for(unsigned int basis = 0; basis < desc.basis; ++basis)
            for(unsigned int component = 0; component < 3; ++component) {
               const std::size_t index = cell + layout.real_cells * (component + 3 * (basis + desc.basis * ensemble));
               result[ensemble] -= static_cast<real>(0.5) * moments[index] * fields[index];
            }
   return result;
}

bool GpuDipoleConvolution::diagnosticConstructionStorageAllocatedForTesting() const {
   return kernel_real_allocated;
}

GpuDipoleSpectrumDiagnostics GpuDipoleConvolution::diagnosticSpectrumForTesting() const {
   if(!initiated || !kernel_ready) throw std::runtime_error("GPU dipole spectrum diagnostic requested before a ready kernel");
   return spectrum_diagnostics;
}

std::size_t GpuDipoleConvolution::estimatePersistentBytes(
      const GpuDipoleConvolutionDescriptor& descriptor) {
   if(!descriptor.valid()) return 0;
   return descriptor.fftLayout().persistentBytes();
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
