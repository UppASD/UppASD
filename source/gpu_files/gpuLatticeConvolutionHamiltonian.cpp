#include "gpuLatticeConvolutionHamiltonian.hpp"

#include "base.hpp"
#include "real_type.h"

namespace {

struct DeviceCellGeometry {
   const real* cell_vectors = nullptr;
   const real* basis_positions = nullptr;
   unsigned int n1 = 0;
   unsigned int n2 = 0;
   unsigned int n3 = 0;
   unsigned int basis_count = 0;
   char bc1 = '0';
   char bc2 = '0';
   char bc3 = '0';
};

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

__device__ inline real sq_norm3(real x, real y, real z) {
   return x * x + y * y + z * z;
}

__device__ inline void cell_index3(unsigned int cell,
                                   const DeviceCellGeometry& geom,
                                   unsigned int& i1,
                                   unsigned int& i2,
                                   unsigned int& i3) {
   i1 = cell % geom.n1;
   const unsigned int t = cell / geom.n1;
   i2 = t % geom.n2;
   i3 = t / geom.n2;
}

__device__ inline void wrapped_coord_diff(const DeviceCellGeometry& geom,
                                          unsigned int basis_out,
                                          unsigned int neighbor_atom,
                                          real& rx,
                                          real& ry,
                                          real& rz) {
   if(!geom.cell_vectors || !geom.basis_positions || geom.basis_count == 0) {
      rx = (real)0.0;
      ry = (real)0.0;
      rz = (real)0.0;
      return;
   }

   const unsigned int cell = neighbor_atom / geom.basis_count;
   const unsigned int basis_in = neighbor_atom - cell * geom.basis_count;
   unsigned int i1 = 0;
   unsigned int i2 = 0;
   unsigned int i3 = 0;
   cell_index3(cell, geom, i1, i2, i3);

   const real* c1 = geom.cell_vectors;
   const real* c2 = geom.cell_vectors + 3;
   const real* c3 = geom.cell_vectors + 6;
   const real* bin = geom.basis_positions + 3 * basis_in;
   const real* bout = geom.basis_positions + 3 * basis_out;

   const real odiff_x = static_cast<real>(i1) * c1[0] + static_cast<real>(i2) * c2[0] +
                        static_cast<real>(i3) * c3[0] + bin[0] - bout[0];
   const real odiff_y = static_cast<real>(i1) * c1[1] + static_cast<real>(i2) * c2[1] +
                        static_cast<real>(i3) * c3[1] + bin[1] - bout[1];
   const real odiff_z = static_cast<real>(i1) * c1[2] + static_cast<real>(i2) * c2[2] +
                        static_cast<real>(i3) * c3[2] + bin[2] - bout[2];

   rx = odiff_x;
   ry = odiff_y;
   rz = odiff_z;
   real best_norm = sq_norm3(rx, ry, rz);

   const int xmin = geom.bc1 == 'P' ? -1 : 0;
   const int xmax = geom.bc1 == 'P' ? 1 : 0;
   const int ymin = geom.bc2 == 'P' ? -1 : 0;
   const int ymax = geom.bc2 == 'P' ? 1 : 0;
   const int zmin = geom.bc3 == 'P' ? -1 : 0;
   const int zmax = geom.bc3 == 'P' ? 1 : 0;

   for(int z = zmin; z <= zmax; ++z) {
      for(int y = ymin; y <= ymax; ++y) {
         for(int x = xmin; x <= xmax; ++x) {
            const real sx = odiff_x + static_cast<real>(x * static_cast<int>(geom.n1)) * c1[0] +
                            static_cast<real>(y * static_cast<int>(geom.n2)) * c2[0] +
                            static_cast<real>(z * static_cast<int>(geom.n3)) * c3[0];
            const real sy = odiff_y + static_cast<real>(x * static_cast<int>(geom.n1)) * c1[1] +
                            static_cast<real>(y * static_cast<int>(geom.n2)) * c2[1] +
                            static_cast<real>(z * static_cast<int>(geom.n3)) * c3[1];
            const real sz = odiff_z + static_cast<real>(x * static_cast<int>(geom.n1)) * c1[2] +
                            static_cast<real>(y * static_cast<int>(geom.n2)) * c2[2] +
                            static_cast<real>(z * static_cast<int>(geom.n3)) * c3[2];
            const real trial_norm = sq_norm3(sx, sy, sz);
            if(trial_norm < best_norm) {
               best_norm = trial_norm;
               rx = sx;
               ry = sy;
               rz = sz;
            }
         }
      }
   }
}

__global__ void pack_spins_to_real(const real* __restrict__ emomM,
                                  real* __restrict__ spin_real,
                                  unsigned int grid_size,
                                  unsigned int basis_count,
                                  unsigned int ensembles) {
   const unsigned int total = grid_size * basis_count * ensembles * 3;
   const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
   if(id >= total) return;

   const unsigned int k = id % grid_size;
   const unsigned int component = id / grid_size;
   const unsigned int axis = component % 3;
   const unsigned int basis = (component / 3) % basis_count;
   const unsigned int ensemble = component / (3 * basis_count);
   const unsigned int atoms = grid_size * basis_count;
   const unsigned int atom = ensemble * atoms + k * basis_count + basis;

   spin_real[id] = emomM[atom * 3 + axis];
}

__global__ void multiply_spectral_kernel(const GpuFftComplex* __restrict__ spin_fft,
                                         GpuFftComplex* __restrict__ field_fft,
                                         const GpuFftComplex* __restrict__ kernel_fft,
                                         unsigned int spectral_size,
                                         unsigned int basis_count,
                                         unsigned int ensembles) {
   const unsigned int total = spectral_size * basis_count * ensembles * 3;
   const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
   if(id >= total) return;

   const unsigned int k = id % spectral_size;
   const unsigned int component = id / spectral_size;
   const unsigned int axis_out = component % 3;
   const unsigned int basis_out = (component / 3) % basis_count;
   const unsigned int ensemble = component / (3 * basis_count);

   GpuFftComplex value = make_fft_complex((real)0.0, (real)0.0);
   for(unsigned int basis_in = 0; basis_in < basis_count; ++basis_in) {
      for(unsigned int axis_in = 0; axis_in < 3; ++axis_in) {
         const unsigned int spin_component = axis_in + 3 * (basis_in + basis_count * ensemble);
         const unsigned int kernel_component = axis_out + 3 * (axis_in + 3 * (basis_out + basis_count * basis_in));
         value = complex_add(value,
                             complex_mul(kernel_fft[k + spectral_size * kernel_component],
                                         spin_fft[k + spectral_size * spin_component]));
      }
   }

   field_fft[id] = value;
}

__global__ void fill_isotropic_dm_kernel(real* __restrict__ kernel_real,
                                         const real* __restrict__ exchange,
                                         const unsigned int* __restrict__ exchange_pos,
                                         unsigned int exchange_mnn,
                                         const real* __restrict__ dm,
                                         const unsigned int* __restrict__ dm_pos,
                                         unsigned int dm_mnn,
                                         bool include_dm,
                                         unsigned int grid_size,
                                         unsigned int basis_count,
                                         unsigned int atom_count,
                                         DeviceCellGeometry geom,
                                         real* __restrict__ rij_cache,
                                         unsigned int rij_stride) {
   const unsigned int basis_out = blockIdx.y;
   const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
   if(basis_out >= basis_count || i >= exchange_mnn + (include_dm ? dm_mnn : 0)) return;

   if(i < exchange_mnn) {
      const unsigned int neighbor_atom = exchange_pos[basis_out + i * atom_count];
      const unsigned int cell = neighbor_atom / basis_count;
      const unsigned int basis_in = neighbor_atom - cell * basis_count;
      if(cell >= grid_size || basis_in >= basis_count) return;

      real rx = (real)0.0;
      real ry = (real)0.0;
      real rz = (real)0.0;
      wrapped_coord_diff(geom, basis_out, neighbor_atom, rx, ry, rz);
      if(rij_cache) {
         const unsigned int rij_entry = basis_out * rij_stride + i;
         rij_cache[0 + 3 * rij_entry] = rx;
         rij_cache[1 + 3 * rij_entry] = ry;
         rij_cache[2 + 3 * rij_entry] = rz;
      }

      const real j = exchange[basis_out + i * basis_count];
      const unsigned int xx = 0 + 3 * (0 + 3 * (basis_out + basis_count * basis_in));
      const unsigned int yy = 1 + 3 * (1 + 3 * (basis_out + basis_count * basis_in));
      const unsigned int zz = 2 + 3 * (2 + 3 * (basis_out + basis_count * basis_in));
      atomicAdd(&kernel_real[cell + grid_size * xx], j);
      atomicAdd(&kernel_real[cell + grid_size * yy], j);
      atomicAdd(&kernel_real[cell + grid_size * zz], j);
      return;
   }

   const unsigned int dm_i = i - exchange_mnn;
   const unsigned int neighbor_atom = dm_pos[basis_out * dm_mnn + dm_i];
   const unsigned int cell = neighbor_atom / basis_count;
   const unsigned int basis_in = neighbor_atom - cell * basis_count;
   if(cell >= grid_size || basis_in >= basis_count) return;

   real rx = (real)0.0;
   real ry = (real)0.0;
   real rz = (real)0.0;
   wrapped_coord_diff(geom, basis_out, neighbor_atom, rx, ry, rz);
   if(rij_cache) {
      const unsigned int rij_entry = basis_out * rij_stride + i;
      rij_cache[0 + 3 * rij_entry] = rx;
      rij_cache[1 + 3 * rij_entry] = ry;
      rij_cache[2 + 3 * rij_entry] = rz;
   }

   const real dx = dm[0 + 3 * (dm_i + dm_mnn * basis_out)];
   const real dy = dm[1 + 3 * (dm_i + dm_mnn * basis_out)];
   const real dz = dm[2 + 3 * (dm_i + dm_mnn * basis_out)];

   const unsigned int xy = 0 + 3 * (1 + 3 * (basis_out + basis_count * basis_in));
   const unsigned int xz = 0 + 3 * (2 + 3 * (basis_out + basis_count * basis_in));
   const unsigned int yx = 1 + 3 * (0 + 3 * (basis_out + basis_count * basis_in));
   const unsigned int yz = 1 + 3 * (2 + 3 * (basis_out + basis_count * basis_in));
   const unsigned int zx = 2 + 3 * (0 + 3 * (basis_out + basis_count * basis_in));
   const unsigned int zy = 2 + 3 * (1 + 3 * (basis_out + basis_count * basis_in));

   atomicAdd(&kernel_real[cell + grid_size * xy], -dz);
   atomicAdd(&kernel_real[cell + grid_size * xz], dy);
   atomicAdd(&kernel_real[cell + grid_size * yx], dz);
   atomicAdd(&kernel_real[cell + grid_size * yz], -dx);
   atomicAdd(&kernel_real[cell + grid_size * zx], -dy);
   atomicAdd(&kernel_real[cell + grid_size * zy], dx);
}

__global__ void fill_tensor_kernel(real* __restrict__ kernel_real,
                                   const real* __restrict__ exchange_tensor,
                                   const unsigned int* __restrict__ exchange_pos,
                                   unsigned int exchange_mnn,
                                   unsigned int grid_size,
                                   unsigned int basis_count,
                                   DeviceCellGeometry geom,
                                   real* __restrict__ rij_cache) {
   const unsigned int total = basis_count * exchange_mnn * 9;
   const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
   if(id >= total) return;

   const unsigned int axis_out = id % 3;
   const unsigned int axis_in = (id / 3) % 3;
   const unsigned int neighbor = (id / 9) % exchange_mnn;
   const unsigned int basis_out = id / (9 * exchange_mnn);
   const unsigned int neighbor_atom = exchange_pos[basis_out * exchange_mnn + neighbor];
   const unsigned int cell = neighbor_atom / basis_count;
   const unsigned int basis_in = neighbor_atom - cell * basis_count;
   if(cell >= grid_size || basis_in >= basis_count) return;

   real rx = (real)0.0;
   real ry = (real)0.0;
   real rz = (real)0.0;
   wrapped_coord_diff(geom, basis_out, neighbor_atom, rx, ry, rz);
   if(rij_cache && axis_out == 0 && axis_in == 0) {
      const unsigned int rij_entry = basis_out * exchange_mnn + neighbor;
      rij_cache[0 + 3 * rij_entry] = rx;
      rij_cache[1 + 3 * rij_entry] = ry;
      rij_cache[2 + 3 * rij_entry] = rz;
   }

   const real j = exchange_tensor[axis_out + 3 * (axis_in + 3 * (neighbor + exchange_mnn * basis_out))];
   const unsigned int kernel_component = axis_out + 3 * (axis_in + 3 * (basis_out + basis_count * basis_in));
   atomicAdd(&kernel_real[cell + grid_size * kernel_component], j);
}

__global__ void add_uniaxial_anisotropy_kernel(real* __restrict__ kernel_real,
                                               const real* __restrict__ anisotropy_k,
                                               const real* __restrict__ anisotropy_axis,
                                               const unsigned int* __restrict__ anisotropy_type,
                                               unsigned int grid_size,
                                               unsigned int basis_count,
                                               unsigned int* __restrict__ unsupported) {
   const unsigned int basis = blockIdx.x * blockDim.x + threadIdx.x;
   if(basis >= basis_count) return;

   const unsigned int type = anisotropy_type[basis];
   const real k1 = anisotropy_k[0 + 2 * basis];
   const real k2 = anisotropy_k[1 + 2 * basis];
   if(type == 0) return;
   if(type != 1 || k2 != (real)0.0) {
      *unsupported = 1;
      return;
   }

   const real ex = anisotropy_axis[0 + 3 * basis];
   const real ey = anisotropy_axis[1 + 3 * basis];
   const real ez = anisotropy_axis[2 + 3 * basis];
   const real pref = (real)-2.0 * k1;
   const real e[3] = {ex, ey, ez};

   for(unsigned int axis_out = 0; axis_out < 3; ++axis_out) {
      for(unsigned int axis_in = 0; axis_in < 3; ++axis_in) {
         const unsigned int kernel_component = axis_out + 3 * (axis_in + 3 * (basis + basis_count * basis));
         atomicAdd(&kernel_real[grid_size * kernel_component], pref * e[axis_out] * e[axis_in]);
      }
   }
}

__global__ void unpack_real_field(const real* __restrict__ field_real,
                                 const real* __restrict__ ext_f,
                                 real* __restrict__ beff,
                                 real* __restrict__ eneff,
                                 unsigned int grid_size,
                                 unsigned int basis_count,
                                 unsigned int ensembles,
                                 real scale) {
   const unsigned int total = grid_size * basis_count * ensembles * 3;
   const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
   if(id >= total) return;

   const unsigned int k = id % grid_size;
   const unsigned int component = id / grid_size;
   const unsigned int axis = component % 3;
   const unsigned int basis = (component / 3) % basis_count;
   const unsigned int ensemble = component / (3 * basis_count);
   const unsigned int atoms = grid_size * basis_count;
   const unsigned int atom = ensemble * atoms + k * basis_count + basis;
   const unsigned int element = atom * 3 + axis;
   const unsigned int field_id = k + grid_size * component;
   const real value = field_real[field_id] * scale + ext_f[element];

   beff[element] = value;
   eneff[element] = value;
}

} // namespace

GpuLatticeConvolutionHamiltonian::GpuLatticeConvolutionHamiltonian()
   : initiated(false), allocated(false), geometry_allocated(false), rij_allocated(false),
     forward_plan{}, backward_plan{}, kernel_plan{} {
}

GpuLatticeConvolutionHamiltonian::~GpuLatticeConvolutionHamiltonian() {
   release();
}

bool GpuLatticeConvolutionHamiltonian::supports(const GpuLatticeConvolutionDescriptor& descriptor,
                                                unsigned int atom_count) const {
   return descriptor.matches(atom_count);
}

bool GpuLatticeConvolutionHamiltonian::initiate(const GpuLatticeConvolutionDescriptor& descriptor,
                                                GPU_STREAM_T stream) {
   release();
   if(!descriptor.valid()) {
      return false;
   }

   desc = descriptor;
   const long int grid_size = static_cast<long int>(desc.cells());
   const long int field_batches = static_cast<long int>(3 * desc.basis * desc.ensembles);
   const long int kernel_batches = static_cast<long int>(9 * desc.basis * desc.basis);
   const long int spectral_size = static_cast<long int>((desc.n1 / 2 + 1) * desc.n2 * desc.n3);

   field_real.Allocate(grid_size, field_batches);
   spin_fft.Allocate(spectral_size, field_batches);
   field_fft.Allocate(spectral_size, field_batches);
   kernel_fft.Allocate(spectral_size, kernel_batches);
   kernel_real.Allocate(grid_size, kernel_batches);
   allocated = true;
   kernel_real.zeros();

   if(desc.c1 && desc.c2 && desc.c3 && desc.basis_positions) {
      cell_vectors.Allocate(static_cast<long int>(3), static_cast<long int>(3));
      basis_positions.Allocate(static_cast<long int>(3), static_cast<long int>(desc.basis));
      geometry_allocated = true;
      ASSERT_GPU(GPU_MEMCPY(cell_vectors.data(), desc.c1, 3 * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE));
      ASSERT_GPU(GPU_MEMCPY(cell_vectors.data() + 3, desc.c2, 3 * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE));
      ASSERT_GPU(GPU_MEMCPY(cell_vectors.data() + 6, desc.c3, 3 * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE));
      ASSERT_GPU(GPU_MEMCPY(basis_positions.data(), desc.basis_positions,
                            3 * desc.basis * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE));
   }

   int n[] = {static_cast<int>(desc.n3), static_cast<int>(desc.n2), static_cast<int>(desc.n1)};
   const int rank = 3;
   const int istride = 1;
   const int ostride = 1;
   const int idist = static_cast<int>(grid_size);
   const int odist = static_cast<int>(spectral_size);
   int* inembed = n;
   int onembed[] = {static_cast<int>(desc.n3), static_cast<int>(desc.n2),
                    static_cast<int>(desc.n1 / 2 + 1)};

   assertGpuFft(GPUFFT_PLAN_MANY(&forward_plan, rank, n, inembed, istride, idist,
                                 onembed, ostride, odist, GPUFFT_R2C,
                                 static_cast<int>(field_batches)));
   assertGpuFft(GPUFFT_SET_STREAM(forward_plan, stream));
   assertGpuFft(GPUFFT_PLAN_MANY(&backward_plan, rank, n, onembed, ostride, odist,
                                 inembed, istride, idist, GPUFFT_C2R,
                                 static_cast<int>(field_batches)));
   assertGpuFft(GPUFFT_SET_STREAM(backward_plan, stream));
   assertGpuFft(GPUFFT_PLAN_MANY(&kernel_plan, rank, n, inembed, istride, idist,
                                 onembed, ostride, odist, GPUFFT_R2C,
                                 static_cast<int>(kernel_batches)));
   assertGpuFft(GPUFFT_SET_STREAM(kernel_plan, stream));

   initiated = true;
   return true;
}

void GpuLatticeConvolutionHamiltonian::release() {
   if(initiated) {
      GPUFFT_DESTROY(forward_plan);
      GPUFFT_DESTROY(backward_plan);
      GPUFFT_DESTROY(kernel_plan);
      initiated = false;
      forward_plan = {};
      backward_plan = {};
      kernel_plan = {};
   }
   if(allocated) {
      field_real.Free();
      spin_fft.Free();
      field_fft.Free();
      kernel_fft.Free();
      kernel_real.Free();
      allocated = false;
   }
   if(geometry_allocated) {
      cell_vectors.Free();
      basis_positions.Free();
      geometry_allocated = false;
   }
   if(rij_allocated) {
      interaction_rij.Free();
      rij_allocated = false;
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

const GpuTensor<real, 2>& GpuLatticeConvolutionHamiltonian::interactionRij() const {
   return interaction_rij;
}

bool GpuLatticeConvolutionHamiltonian::buildIsotropicDmKernel(const GpuTensor<real, 2>& exchange,
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
                                                              GPU_STREAM_T stream) {
   if(!initiated) return false;

   kernel_real.zeros_async(stream);
   if(rij_allocated) {
      interaction_rij.Free();
      rij_allocated = false;
   }
   const unsigned int total = exchange_mnn + (include_dm ? dm_mnn : 0);
   const unsigned int threads = 256;
   DeviceCellGeometry geom{};
   geom.cell_vectors = cell_vectors.data();
   geom.basis_positions = basis_positions.data();
   geom.n1 = desc.n1;
   geom.n2 = desc.n2;
   geom.n3 = desc.n3;
   geom.basis_count = desc.basis;
   geom.bc1 = desc.bc1;
   geom.bc2 = desc.bc2;
   geom.bc3 = desc.bc3;

   if(total > 0) {
      interaction_rij.Allocate(static_cast<long int>(3), static_cast<long int>(desc.basis * total));
      rij_allocated = true;
      const dim3 blocks((total + threads - 1) / threads, desc.basis);
      fill_isotropic_dm_kernel<<<blocks, threads, 0, stream>>>(kernel_real.data(),
                                                               exchange.data(),
                                                               exchange_pos.data(),
                                                               exchange_mnn,
                                                               dm.data(),
                                                               dm_pos.data(),
                                                               dm_mnn,
                                                               include_dm,
                                                               desc.cells(),
                                                               desc.basis,
                                                               desc.atoms(),
                                                               geom,
                                                               interaction_rij.data(),
                                                               total);
   }

   if(include_uniaxial_anisotropy) {
      GpuTensor<unsigned int, 1> unsupported;
      unsupported.Allocate(static_cast<long int>(1));
      unsupported.zeros_async(stream);
      const dim3 blocks((desc.basis + threads - 1) / threads);
      add_uniaxial_anisotropy_kernel<<<blocks, threads, 0, stream>>>(kernel_real.data(),
                                                                     anisotropy_k.data(),
                                                                     anisotropy_axis.data(),
                                                                     anisotropy_type.data(),
                                                                     desc.cells(),
                                                                     desc.basis,
                                                                     unsupported.data());
      ASSERT_GPU(GPU_STREAM_SYNC(stream));
      unsigned int host_unsupported = 0;
      ASSERT_GPU(GPU_MEMCPY(&host_unsupported, unsupported.data(), sizeof(unsigned int),
                            GPU_MEMCPY_DEVICE_TO_HOST));
      unsupported.Free();
      if(host_unsupported != 0) {
         return false;
      }
   }

   assertGpuFft(GPUFFT_EXEC_R2C(kernel_plan, kernel_real.data(), kernel_fft.data()));
   return true;
}

bool GpuLatticeConvolutionHamiltonian::buildTensorKernel(const GpuTensor<real, 4>& exchange_tensor,
                                                         const GpuTensor<unsigned int, 2>& exchange_pos,
                                                         unsigned int exchange_mnn,
                                                         const GpuTensor<real, 2>& anisotropy_k,
                                                         const GpuTensor<real, 2>& anisotropy_axis,
                                                         const GpuTensor<unsigned int, 1>& anisotropy_type,
                                                         bool include_uniaxial_anisotropy,
                                                         GPU_STREAM_T stream) {
   if(!initiated || exchange_mnn == 0) return false;

   kernel_real.zeros_async(stream);
   if(rij_allocated) {
      interaction_rij.Free();
      rij_allocated = false;
   }
   interaction_rij.Allocate(static_cast<long int>(3), static_cast<long int>(desc.basis * exchange_mnn));
   rij_allocated = true;
   const unsigned int total = desc.basis * exchange_mnn * 9;
   const unsigned int threads = 256;
   const dim3 grid((total + threads - 1) / threads);
   DeviceCellGeometry geom{};
   geom.cell_vectors = cell_vectors.data();
   geom.basis_positions = basis_positions.data();
   geom.n1 = desc.n1;
   geom.n2 = desc.n2;
   geom.n3 = desc.n3;
   geom.basis_count = desc.basis;
   geom.bc1 = desc.bc1;
   geom.bc2 = desc.bc2;
   geom.bc3 = desc.bc3;
   fill_tensor_kernel<<<grid, threads, 0, stream>>>(kernel_real.data(),
                                                    exchange_tensor.data(),
                                                    exchange_pos.data(),
                                                    exchange_mnn,
                                                    desc.cells(),
                                                    desc.basis,
                                                    geom,
                                                    interaction_rij.data());
   if(include_uniaxial_anisotropy) {
      GpuTensor<unsigned int, 1> unsupported;
      unsupported.Allocate(static_cast<long int>(1));
      unsupported.zeros_async(stream);
      const dim3 ani_grid((desc.basis + threads - 1) / threads);
      add_uniaxial_anisotropy_kernel<<<ani_grid, threads, 0, stream>>>(kernel_real.data(),
                                                                       anisotropy_k.data(),
                                                                       anisotropy_axis.data(),
                                                                       anisotropy_type.data(),
                                                                       desc.cells(),
                                                                       desc.basis,
                                                                       unsupported.data());
      ASSERT_GPU(GPU_STREAM_SYNC(stream));
      unsigned int host_unsupported = 0;
      ASSERT_GPU(GPU_MEMCPY(&host_unsupported, unsupported.data(), sizeof(unsigned int),
                            GPU_MEMCPY_DEVICE_TO_HOST));
      unsupported.Free();
      if(host_unsupported != 0) {
         return false;
      }
   }
   assertGpuFft(GPUFFT_EXEC_R2C(kernel_plan, kernel_real.data(), kernel_fft.data()));
   return true;
}

void GpuLatticeConvolutionHamiltonian::apply(deviceLattice& gpuLattice,
                                             const GpuTensor<real, 3>& external_field,
                                             GPU_STREAM_T stream) {
   if(!initiated) return;

   const unsigned int grid_size = desc.cells();
   const unsigned int threads = 256;
   const unsigned int packed_total = grid_size * desc.basis * desc.ensembles * 3;
   const dim3 packed_grid((packed_total + threads - 1) / threads);

   // This real buffer is reusable: R2C has consumed the packed spins before
   // the inverse transform writes the real field back into it.
   pack_spins_to_real<<<packed_grid, threads, 0, stream>>>(gpuLattice.emomM.data(), field_real.data(),
                                                          grid_size, desc.basis, desc.ensembles);
   assertGpuFft(GPUFFT_EXEC_R2C(forward_plan, field_real.data(), spin_fft.data()));

   const unsigned int spectral_size = (desc.n1 / 2 + 1) * desc.n2 * desc.n3;
   const dim3 spectral_grid((spectral_size * desc.basis * desc.ensembles * 3 + threads - 1) / threads);
   multiply_spectral_kernel<<<spectral_grid, threads, 0, stream>>>(spin_fft.data(), field_fft.data(),
                                                                 kernel_fft.data(), spectral_size,
                                                                 desc.basis, desc.ensembles);
   assertGpuFft(GPUFFT_EXEC_C2R(backward_plan, field_fft.data(), field_real.data()));

   const real scale = (real)1.0 / static_cast<real>(grid_size);
   unpack_real_field<<<packed_grid, threads, 0, stream>>>(field_real.data(), external_field.data(),
                                                         gpuLattice.beff.data(), gpuLattice.eneff.data(),
                                                         grid_size, desc.basis, desc.ensembles, scale);
}
