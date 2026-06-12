#include "gpuLatticeConvolutionHamiltonian.hpp"

#include "base.hpp"
#include "real_type.h"

namespace {

__device__ inline GpuFftComplex make_fft_complex(real re, real im) {
   GpuFftComplex z{};
   z.x = re;
   z.y = im;
   return z;
}

__device__ inline GpuFftComplex complex_add(GpuFftComplex a, GpuFftComplex b) {
   return make_fft_complex(a.x + b.x, a.y + b.y);
}

__device__ inline GpuFftComplex complex_mul(GpuFftComplex a, GpuFftComplex b) {
   return make_fft_complex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

__global__ void pack_spins_to_fft(const real* __restrict__ emomM,
                                  GpuFftComplex* __restrict__ spin_fft,
                                  unsigned int grid_size,
                                  unsigned int ensembles) {
   const unsigned int total = grid_size * ensembles * 3;
   const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
   if(id >= total) return;

   const unsigned int k = id % grid_size;
   const unsigned int axis_ensemble = id / grid_size;
   const unsigned int axis = axis_ensemble % 3;
   const unsigned int ensemble = axis_ensemble / 3;
   const unsigned int atom = ensemble * grid_size + k;

   spin_fft[id] = make_fft_complex(emomM[atom * 3 + axis], (real)0.0);
}

__global__ void multiply_spectral_kernel(const GpuFftComplex* __restrict__ spin_fft,
                                         GpuFftComplex* __restrict__ field_fft,
                                         const GpuFftComplex* __restrict__ kernel_fft,
                                         unsigned int grid_size,
                                         unsigned int ensembles) {
   const unsigned int total = grid_size * ensembles;
   const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
   if(id >= total) return;

   const unsigned int k = id % grid_size;
   const unsigned int ensemble = id / grid_size;
   const unsigned int spin_base = k + grid_size * (3 * ensemble);

   const GpuFftComplex sx = spin_fft[spin_base + grid_size * 0];
   const GpuFftComplex sy = spin_fft[spin_base + grid_size * 1];
   const GpuFftComplex sz = spin_fft[spin_base + grid_size * 2];

   const GpuFftComplex kxx = kernel_fft[k + grid_size * 0];
   const GpuFftComplex kxy = kernel_fft[k + grid_size * 1];
   const GpuFftComplex kxz = kernel_fft[k + grid_size * 2];
   const GpuFftComplex kyx = kernel_fft[k + grid_size * 3];
   const GpuFftComplex kyy = kernel_fft[k + grid_size * 4];
   const GpuFftComplex kyz = kernel_fft[k + grid_size * 5];
   const GpuFftComplex kzx = kernel_fft[k + grid_size * 6];
   const GpuFftComplex kzy = kernel_fft[k + grid_size * 7];
   const GpuFftComplex kzz = kernel_fft[k + grid_size * 8];

   field_fft[spin_base + grid_size * 0] = complex_add(complex_add(complex_mul(kxx, sx), complex_mul(kxy, sy)), complex_mul(kxz, sz));
   field_fft[spin_base + grid_size * 1] = complex_add(complex_add(complex_mul(kyx, sx), complex_mul(kyy, sy)), complex_mul(kyz, sz));
   field_fft[spin_base + grid_size * 2] = complex_add(complex_add(complex_mul(kzx, sx), complex_mul(kzy, sy)), complex_mul(kzz, sz));
}

__global__ void fill_isotropic_dm_kernel(GpuFftComplex* __restrict__ kernel_fft,
                                         const real* __restrict__ exchange,
                                         const unsigned int* __restrict__ exchange_pos,
                                         unsigned int exchange_mnn,
                                         const real* __restrict__ dm,
                                         const unsigned int* __restrict__ dm_pos,
                                         unsigned int dm_mnn,
                                         bool include_dm,
                                         unsigned int grid_size) {
   const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
   if(i >= exchange_mnn + (include_dm ? dm_mnn : 0)) return;

   if(i < exchange_mnn) {
      const unsigned int site = exchange_pos[i * grid_size];
      const real j = exchange[i];
      atomicAdd(&kernel_fft[site + grid_size * 0].x, j);
      atomicAdd(&kernel_fft[site + grid_size * 4].x, j);
      atomicAdd(&kernel_fft[site + grid_size * 8].x, j);
      return;
   }

   const unsigned int dm_i = i - exchange_mnn;
   const unsigned int site = dm_pos[dm_i];
   const real dx = dm[0 + 3 * dm_i];
   const real dy = dm[1 + 3 * dm_i];
   const real dz = dm[2 + 3 * dm_i];

   atomicAdd(&kernel_fft[site + grid_size * 1].x, -dz);
   atomicAdd(&kernel_fft[site + grid_size * 2].x, dy);
   atomicAdd(&kernel_fft[site + grid_size * 3].x, dz);
   atomicAdd(&kernel_fft[site + grid_size * 5].x, -dx);
   atomicAdd(&kernel_fft[site + grid_size * 6].x, -dy);
   atomicAdd(&kernel_fft[site + grid_size * 7].x, dx);
}

__global__ void unpack_fft_field(const GpuFftComplex* __restrict__ field_fft,
                                 const real* __restrict__ ext_f,
                                 real* __restrict__ beff,
                                 real* __restrict__ eneff,
                                 unsigned int grid_size,
                                 unsigned int ensembles,
                                 real scale) {
   const unsigned int total = grid_size * ensembles * 3;
   const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
   if(id >= total) return;

   const unsigned int k = id % grid_size;
   const unsigned int axis_ensemble = id / grid_size;
   const unsigned int axis = axis_ensemble % 3;
   const unsigned int ensemble = axis_ensemble / 3;
   const unsigned int atom = ensemble * grid_size + k;
   const unsigned int element = atom * 3 + axis;
   const unsigned int spectral_id = k + grid_size * axis_ensemble;
   const real value = field_fft[spectral_id].x * scale + ext_f[element];

   beff[element] = value;
   eneff[element] = value;
}

} // namespace

GpuLatticeConvolutionHamiltonian::GpuLatticeConvolutionHamiltonian()
   : initiated(false), allocated(false), plan{}, kernel_plan{} {
}

GpuLatticeConvolutionHamiltonian::~GpuLatticeConvolutionHamiltonian() {
   release();
}

bool GpuLatticeConvolutionHamiltonian::supports(const GpuLatticeConvolutionDescriptor& descriptor,
                                                unsigned int atom_count) const {
   return descriptor.n1 > 0 && descriptor.n2 > 0 && descriptor.n3 > 0 &&
          descriptor.ensembles > 0 && descriptor.oneSpinPerCell(atom_count);
}

bool GpuLatticeConvolutionHamiltonian::initiate(const GpuLatticeConvolutionDescriptor& descriptor,
                                                GPU_STREAM_T stream) {
   release();
   if(!supports(descriptor, descriptor.cells())) {
      return false;
   }

   desc = descriptor;
   const long int grid_size = static_cast<long int>(desc.cells());
   const long int batches = static_cast<long int>(3 * desc.ensembles);

   spin_fft.Allocate(grid_size, batches);
   field_fft.Allocate(grid_size, batches);
   kernel_fft.Allocate(grid_size, static_cast<long int>(9));
   allocated = true;
   kernel_fft.zeros();

   int n[] = {static_cast<int>(desc.n3), static_cast<int>(desc.n2), static_cast<int>(desc.n1)};
   const int rank = 3;
   const int istride = 1;
   const int ostride = 1;
   const int idist = static_cast<int>(grid_size);
   const int odist = static_cast<int>(grid_size);
   int* inembed = n;
   int* onembed = n;

   assertGpuFft(GPUFFT_PLAN_MANY(&plan, rank, n, inembed, istride, idist,
                                 onembed, ostride, odist, GPUFFT_C2C,
                                 static_cast<int>(batches)));
   assertGpuFft(GPUFFT_SET_STREAM(plan, stream));
   assertGpuFft(GPUFFT_PLAN_MANY(&kernel_plan, rank, n, inembed, istride, idist,
                                 onembed, ostride, odist, GPUFFT_C2C, 9));
   assertGpuFft(GPUFFT_SET_STREAM(kernel_plan, stream));

   initiated = true;
   return true;
}

void GpuLatticeConvolutionHamiltonian::release() {
   if(initiated) {
      GPUFFT_DESTROY(plan);
      GPUFFT_DESTROY(kernel_plan);
      initiated = false;
      plan = {};
      kernel_plan = {};
   }
   if(allocated) {
      spin_fft.Free();
      field_fft.Free();
      kernel_fft.Free();
      allocated = false;
   }
}

bool GpuLatticeConvolutionHamiltonian::isInitiated() const {
   return initiated;
}

GpuTensor<GpuFftComplex, 2>& GpuLatticeConvolutionHamiltonian::spectralKernel() {
   return kernel_fft;
}

const GpuTensor<GpuFftComplex, 2>& GpuLatticeConvolutionHamiltonian::spectralKernel() const {
   return kernel_fft;
}

void GpuLatticeConvolutionHamiltonian::buildIsotropicDmKernel(const GpuTensor<real, 2>& exchange,
                                                              const GpuTensor<unsigned int, 2>& exchange_pos,
                                                              unsigned int exchange_mnn,
                                                              const GpuTensor<real, 3>& dm,
                                                              const GpuTensor<unsigned int, 2>& dm_pos,
                                                              unsigned int dm_mnn,
                                                              bool include_dm,
                                                              GPU_STREAM_T stream) {
   if(!initiated) return;

   kernel_fft.zeros_async(stream);
   const unsigned int total = exchange_mnn + (include_dm ? dm_mnn : 0);
   if(total == 0) return;

   const unsigned int threads = 256;
   const dim3 grid((total + threads - 1) / threads);
   fill_isotropic_dm_kernel<<<grid, threads, 0, stream>>>(kernel_fft.data(),
                                                          exchange.data(),
                                                          exchange_pos.data(),
                                                          exchange_mnn,
                                                          dm.data(),
                                                          dm_pos.data(),
                                                          dm_mnn,
                                                          include_dm,
                                                          desc.cells());
   assertGpuFft(GPUFFT_EXEC_C2C(kernel_plan, kernel_fft.data(), kernel_fft.data(), GPUFFT_FORWARD));
}

void GpuLatticeConvolutionHamiltonian::apply(deviceLattice& gpuLattice,
                                             const GpuTensor<real, 3>& external_field,
                                             GPU_STREAM_T stream) {
   if(!initiated) return;

   const unsigned int grid_size = desc.cells();
   const unsigned int threads = 256;
   const unsigned int packed_total = grid_size * desc.ensembles * 3;
   const unsigned int vector_total = grid_size * desc.ensembles;
   const dim3 packed_grid((packed_total + threads - 1) / threads);
   const dim3 vector_grid((vector_total + threads - 1) / threads);

   pack_spins_to_fft<<<packed_grid, threads, 0, stream>>>(gpuLattice.emomM.data(), spin_fft.data(),
                                                          grid_size, desc.ensembles);
   assertGpuFft(GPUFFT_EXEC_C2C(plan, spin_fft.data(), spin_fft.data(), GPUFFT_FORWARD));

   multiply_spectral_kernel<<<vector_grid, threads, 0, stream>>>(spin_fft.data(), field_fft.data(),
                                                                 kernel_fft.data(), grid_size,
                                                                 desc.ensembles);
   assertGpuFft(GPUFFT_EXEC_C2C(plan, field_fft.data(), field_fft.data(), GPUFFT_BACKWARD));

   const real scale = (real)1.0 / static_cast<real>(grid_size);
   unpack_fft_field<<<packed_grid, threads, 0, stream>>>(field_fft.data(), external_field.data(),
                                                         gpuLattice.beff.data(), gpuLattice.eneff.data(),
                                                         grid_size, desc.ensembles, scale);
}
