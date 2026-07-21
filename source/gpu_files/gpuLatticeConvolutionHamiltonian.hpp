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

// On-site anisotropy, evaluated directly in the unpack/energy kernel rather than
// folded into the spectral exchange kernel. Kept out of the tensor because the
// on-site anisotropy is biquadratic: its field prefactor (2*k1) differs from its
// energy prefactor (k1), so folding it into the (bilinear) convolution kernel
// would corrupt the exchange energy. Raw device pointers mirror the sparse
// GpuHamiltonianCalculations::Heisge indexing (per site: kaniso[2], eaniso[3]).
struct GpuLatticeConvolutionAnisotropy {
   const real* kaniso = nullptr;
   const real* eaniso = nullptr;
   const unsigned int* taniso = nullptr;
   const real* sb = nullptr;
   bool present = false;
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
                               GPU_STREAM_T stream);

   bool buildTensorKernel(const GpuTensor<real, 4>& exchange_tensor,
                          const GpuTensor<unsigned int, 2>& exchange_pos,
                          unsigned int exchange_mnn,
                          GPU_STREAM_T stream);

   // Evaluate the convolution field into gpuLattice.beff/eneff (adding the
   // external field and on-site anisotropy). When energies != nullptr the same
   // kernel also reduces the per-ensemble energy columns of energyM: the bilinear
   // (exchange+DM) energy into column energy_col (0 for isotropic exchange, 3 for
   // tensor exchange, matching the CPU do_jtensor convention), plus the on-site
   // anisotropy (col 1), external (col 4) and total (col 5). The caller must have
   // zeroed energyM beforehand.
   void apply(deviceLattice& gpuLattice, const GpuTensor<real, 3>& external_field,
              const GpuLatticeConvolutionAnisotropy& anisotropy, bool includeAnisotropy,
              deviceEnergies* energies, unsigned int energy_col, GPU_STREAM_T stream);

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
