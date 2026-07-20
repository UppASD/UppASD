#pragma once

#include "gpuFftWrapper.hpp"
#include "gpuStructures.hpp"
#include "gpu_wrappers.h"
#include "tensor.hpp"

struct GpuLatticeConvolutionDescriptor {
   unsigned int n1 = 0;
   unsigned int n2 = 0;
   unsigned int n3 = 0;
   unsigned int basis = 0;
   unsigned int ensembles = 0;
   char bc1 = '0';
   char bc2 = '0';
   char bc3 = '0';
   const real* c1 = nullptr;
   const real* c2 = nullptr;
   const real* c3 = nullptr;
   const real* basis_positions = nullptr;

   unsigned int cells() const {
      return n1 * n2 * n3;
   }

   unsigned int atoms() const {
      return cells() * basis;
   }

   bool valid() const {
      return n1 > 0 && n2 > 0 && n3 > 0 && basis > 0 && ensembles > 0;
   }

   bool matches(unsigned int atom_count) const {
      return valid() && atoms() == atom_count;
   }
};

class GpuLatticeConvolutionHamiltonian {
public:
   GpuLatticeConvolutionHamiltonian();
   ~GpuLatticeConvolutionHamiltonian();

   bool initiate(const GpuLatticeConvolutionDescriptor& descriptor, GPU_STREAM_T stream);
   // Projected device bytes for the FFT buffers, derived from the cell grid
   // (N1,N2,N3), basis (NA) and ensembles (M). Mirror the allocations in initiate().
   static std::size_t estimateBytes(const SimulationParameters& SimParam);
   void release();
   bool isInitiated() const;
   bool supports(const GpuLatticeConvolutionDescriptor& descriptor, unsigned int atom_count) const;

   GpuTensor<GpuFftComplex, 2>& spectralKernel();
   const GpuTensor<GpuFftComplex, 2>& spectralKernel() const;
   const GpuTensor<real, 2>& interactionRij() const;
   bool buildIsotropicDmKernel(const GpuTensor<real, 2>& exchange,
                               const GpuTensor<unsigned int, 2>& exchange_pos,
                               unsigned int exchange_mnn,
                               const GpuTensor<real, 3>& dm,
                               const GpuTensor<unsigned int, 2>& dm_pos,
                               unsigned int dm_mnn,
                               bool include_dm,
                               const GpuTensor<real, 2>& anisotropy_k,
                               const GpuTensor<real, 2>& anisotropy_axis,
                               const GpuTensor<unsigned int, 1>& anisotropy_type,
                               bool include_uniaxial_anisotropy,
                               GPU_STREAM_T stream);

   bool buildTensorKernel(const GpuTensor<real, 4>& exchange_tensor,
                          const GpuTensor<unsigned int, 2>& exchange_pos,
                          unsigned int exchange_mnn,
                          const GpuTensor<real, 2>& anisotropy_k,
                          const GpuTensor<real, 2>& anisotropy_axis,
                          const GpuTensor<unsigned int, 1>& anisotropy_type,
                          bool include_uniaxial_anisotropy,
                          GPU_STREAM_T stream);

   void apply(deviceLattice& gpuLattice, const GpuTensor<real, 3>& external_field, GPU_STREAM_T stream);

private:
   GpuLatticeConvolutionDescriptor desc;
   bool initiated;
   bool allocated;
   bool geometry_allocated;
   bool rij_allocated;
   GpuFftHandle forward_plan;
   GpuFftHandle backward_plan;
   GpuFftHandle kernel_plan;
   GpuTensor<real, 2> field_real;
   GpuTensor<GpuFftComplex, 2> spin_fft;
   GpuTensor<GpuFftComplex, 2> field_fft;
   GpuTensor<GpuFftComplex, 2> kernel_fft;
   GpuTensor<real, 2> kernel_real;
   GpuTensor<real, 2> cell_vectors;
   GpuTensor<real, 2> basis_positions;
   GpuTensor<real, 2> interaction_rij;
};
