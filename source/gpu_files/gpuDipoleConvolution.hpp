#pragma once

#include "gpu_wrappers.h"
#include "gpuFftWrapper.hpp"
#include "real_type.h"
#include "tensor.hpp"

#include <cstddef>

// CV6 deliberately separates the grid used by the dipole solver from the
// exchange-convolution grid.  Dipole boundary conditions define the physical
// long-range sum and must not be inferred from exchange's current PBC-only
// convolution eligibility.
enum class GpuDipoleBoundaryMode {
   Open,
   Periodic3D,
   Periodic2D
};

enum class GpuDipoleDiscretization {
   Atomistic,
   MacrospinGrid
};

struct GpuDipoleGridShape {
   std::size_t n1 = 0;
   std::size_t n2 = 0;
   std::size_t n3 = 0;

   std::size_t cells() const;
   bool valid() const;
};

// Buffer dimensions for the regular-grid R2C/C2R path.  Keeping this as a
// value type makes the FFT-library integer limits and the memory contract
// testable before a plan or a device allocation exists.
struct GpuDipoleFftLayout {
   GpuDipoleGridShape real_grid{};
   GpuDipoleGridShape spectral_grid{};
   std::size_t real_cells = 0;
   std::size_t spectral_cells = 0;
   std::size_t field_batches = 0;   // 3 * basis * ensembles
   std::size_t kernel_batches = 0;  // 9 * basis * basis

   bool valid() const;
   bool fitsFftLibrary() const;
   std::size_t persistentBytes() const;
   std::size_t constructionBytes() const;
};

struct GpuDipoleConvolutionDescriptor {
   // Atomistic grid dimensions.  For MacrospinGrid these describe the source
   // lattice; macro_grid describes the actual FFT grid.
   GpuDipoleGridShape atomistic_grid{};
   GpuDipoleGridShape macro_grid{};
   unsigned int basis = 0;
   unsigned int ensembles = 0;
   GpuDipoleBoundaryMode boundary = GpuDipoleBoundaryMode::Open;
   GpuDipoleDiscretization discretization = GpuDipoleDiscretization::Atomistic;

   // Columns of the physical cell matrix H=[C1 C2 C3].  The periodic
   // reciprocal tensor must be formed from this full (possibly skew) matrix;
   // Cartesian grid extents alone are never a valid substitute.
   const real* c1 = nullptr;
   const real* c2 = nullptr;
   const real* c3 = nullptr;

   // Only used for Periodic2D.  The selected axis is the non-periodic slab
   // normal and is linearly padded; the remaining axes are periodic.
   unsigned int open_axis = 2;

   GpuDipoleGridShape activeGrid() const;
   GpuDipoleGridShape paddedGrid() const;
   GpuDipoleFftLayout fftLayout() const;
   bool valid() const;
};

// Lifecycle and geometry shell for CV6.2+.  This intentionally owns no FFT
// plans or tensors yet: CV6.0 must first select and document the atomistic or
// macrospin physics contract.  Once initiated, all future plans must use the
// supplied work stream.
class GpuDipoleConvolution {
public:
   GpuDipoleConvolution() = default;
   ~GpuDipoleConvolution();
   GpuDipoleConvolution(const GpuDipoleConvolution&) = delete;
   GpuDipoleConvolution& operator=(const GpuDipoleConvolution&) = delete;

   bool initiate(const GpuDipoleConvolutionDescriptor& descriptor, GPU_STREAM_T work_stream);
   void release();

   bool isInitiated() const;
   const GpuDipoleConvolutionDescriptor& descriptor() const;
   const GpuDipoleFftLayout& fftLayout() const;

   // Persistent field and spectral buffers required by the eventual regular
   // grid solver.  Tensor construction staging is deliberately separate
   // because it is released after the kernel FFT.
   static std::size_t estimatePersistentBytes(const GpuDipoleConvolutionDescriptor& descriptor);
   static std::size_t estimateWorkspaceBytes(const GpuDipoleConvolutionDescriptor& descriptor);
   static std::size_t estimateBytes(const GpuDipoleConvolutionDescriptor& descriptor);

private:
   GpuDipoleConvolutionDescriptor desc{};
   GpuDipoleFftLayout layout{};
   GPU_STREAM_T stream{};
   bool initiated = false;
   bool allocated = false;
   bool forward_plan_created = false;
   bool backward_plan_created = false;
   bool kernel_plan_created = false;
   GpuFftHandle forward_plan{};
   GpuFftHandle backward_plan{};
   GpuFftHandle kernel_plan{};
   GpuTensor<real, 2> moments_real;
   GpuTensor<real, 2> fields_real;
   GpuTensor<GpuFftComplex, 2> moments_fft;
   GpuTensor<GpuFftComplex, 2> fields_fft;
   GpuTensor<GpuFftComplex, 2> kernel_fft;
   GpuTensor<real, 2> cell_vectors;
   GpuTensor<unsigned char, 1> fft_workspace;
};
