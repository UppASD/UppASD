#pragma once


#include <numeric>

#include "tensor.cuh"


template <typename T>
__device__ __forceinline__ T tnorm(const CudaMatrix<T>& M, index_t i) {
   return sqrt(M(i, 0) * M(i, 0) + M(i, 1) * M(i, 1) + M(i, 2) * M(i, 2));
}


// function parameters have shallow copy semantics
// cannot be passed by reference on GPU(in our case), might crash
template <typename T, unsigned threads_pb>
__global__ void GPUMeasurements(const CudaMatrix<T> M, CudaVector<T> Ms, CudaVector<T> Mss) {
   __shared__ T loc_ms[threads_pb];
   __shared__ T loc_mss[threads_pb];

   index_t ind = blockIdx.x * blockDim.x + threadIdx.x;
   T mt{}, mtt{};
   for(; ind < M.extent(0); ind += blockDim.x * gridDim.x) {
      T x = tnorm(M, ind);
      mt += x;
      mtt += x * x;
   }
   loc_ms[threadIdx.x] = mt;
   loc_mss[threadIdx.x] = mtt;

   __syncthreads();

   auto it = blockDim.x / 2;
   for(; it != 0; it /= 2) {
      if(threadIdx.x < it) {
         loc_ms[threadIdx.x] += loc_ms[threadIdx.x + it];
         loc_mss[threadIdx.x] += loc_mss[threadIdx.x + it];
      }
      __syncthreads();
   }

   if(threadIdx.x == 0) {
      Ms[blockIdx.x] = loc_ms[0];
      Mss[blockIdx.x] = loc_mss[0];
   }
}


template <typename T>
void AllocateMeasurements(Vector<T>& Ms, Vector<T>& Mss, CudaVector<T>& Ms_gpu, CudaVector<T>& Mss_gpu,
                          index_t size) {
   Ms.AllocateHost(size);
   Mss.AllocateHost(size);

   Ms_gpu.Allocate(size);
   Mss_gpu.Allocate(size);
}


template <typename T>
void FreeMeasurements(Vector<T>& Ms, Vector<T>& Mss, CudaVector<T>& Ms_gpu, CudaVector<T>& Mss_gpu) {
   Ms.FreeHost();
   Mss.FreeHost();

   Ms_gpu.Free();
   Mss_gpu.Free();
}


// sum is constrained to Vector<T>, but can be made more flexible in the future if needed
template <typename T>
T sum(const Vector<T>& x) {
   return std::accumulate(x.begin(), x.end(), T{});
}


// function parameters have shallow copy semantics; you can pass it by constant copy if you prefer
template <typename T>
T Measurements(const Matrix<T>& M) {
   // assume that  M is of size Nt x 3 rather than 3 x Nt
   assert(M.extent(1) == 3);

   // M_gpu is not allocated in AllocateMeasurements because how it is allocated and computed might change
   CudaMatrix<T> M_gpu;
   M_gpu.Allocate(M.extents());
   M_gpu.copy_sync(M);

   constexpr dim3 threads_pb(256);
   const dim3 blocks((static_cast<unsigned>(M.extent(0)) + threads_pb.x - 1) / threads_pb.x);

   // Ms, Mss can be allocated static as a fixed-size array(depending on number of blocks) or preallocated
   Vector<T> Ms, Mss;
   CudaVector<T> Ms_gpu, Mss_gpu;
   AllocateMeasurements(Ms, Mss, Ms_gpu, Mss_gpu, blocks.x);

   GPUMeasurements<T, threads_pb.x><<<blocks, threads_pb>>>(M_gpu, Ms_gpu, Mss_gpu);

   Ms.copy_sync(Ms_gpu);
   Mss.copy_sync(Mss_gpu);

   T m_sum = sum(Ms) / static_cast<T>(M.extent(0));
   T m_sum_squared = sum(Mss) / static_cast<T>(M.extent(0));
   T delta_m = m_sum_squared - m_sum * m_sum;

   FreeMeasurements(Ms, Mss, Ms_gpu, Mss_gpu);

   // Similarly, not freed in FreeMeasurements
   M_gpu.Free();

   return delta_m;
}
