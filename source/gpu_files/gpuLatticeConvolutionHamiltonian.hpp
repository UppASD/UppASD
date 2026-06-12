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

   unsigned int cells() const {
      return n1 * n2 * n3;
   }

   bool oneSpinPerCell(unsigned int atom_count) const {
      return basis == 1 && cells() == atom_count;
   }
};

class GpuLatticeConvolutionHamiltonian {
public:
   GpuLatticeConvolutionHamiltonian();
   ~GpuLatticeConvolutionHamiltonian();

   bool initiate(const GpuLatticeConvolutionDescriptor& descriptor, GPU_STREAM_T stream);
   void release();
   bool isInitiated() const;
   bool supports(const GpuLatticeConvolutionDescriptor& descriptor, unsigned int atom_count) const;

   GpuTensor<GpuFftComplex, 2>& spectralKernel();
   const GpuTensor<GpuFftComplex, 2>& spectralKernel() const;
   void buildIsotropicDmKernel(const GpuTensor<real, 2>& exchange,
                               const GpuTensor<unsigned int, 2>& exchange_pos,
                               unsigned int exchange_mnn,
                               const GpuTensor<real, 3>& dm,
                               const GpuTensor<unsigned int, 2>& dm_pos,
                               unsigned int dm_mnn,
                               bool include_dm,
                               GPU_STREAM_T stream);

   void apply(deviceLattice& gpuLattice, const GpuTensor<real, 3>& external_field, GPU_STREAM_T stream);

private:
   GpuLatticeConvolutionDescriptor desc;
   bool initiated;
   bool allocated;
   GpuFftHandle plan;
   GpuFftHandle kernel_plan;
   GpuTensor<GpuFftComplex, 2> spin_fft;
   GpuTensor<GpuFftComplex, 2> field_fft;
   GpuTensor<GpuFftComplex, 2> kernel_fft;
};
