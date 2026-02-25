#pragma once

#include "tensor.hpp"
#include "real_type.h"
#include <cmath>
#include "gpu_wrappers.h"
#if defined(HIP_V)
#include <hip/hip_runtime.h>
#include <hip/hip_cooperative_groups.h>
#elif defined(CUDA_V)
#include <cuda_runtime.h>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <cuda.h>
#endif
namespace cg = cooperative_groups;



__inline__ __device__
void warpReduceScalar(real &sum) {
    const unsigned int FULL_MASK = 0xffffffff;
    
    #pragma unroll
    for (unsigned int offset = warpSize / 2; offset > 0; offset >>= 1) {
#if defined (CUDA_V)        
#if CUDA_VERSION < 9000
        real shfl = __shfl_down(sum, offset);
#else
        real shfl = __shfl_down_sync(FULL_MASK, sum, offset);
#endif
        sum += shfl;
#elif defined(HIP_V)
        sum += __shfl_down(sum, offset);

#endif
    }
}


__global__ void fill_spinwait(GpuTensor<real, 4> spinwait, const GpuTensor<real, 3> emom, const int taskNum, const int swIdx);

__global__ void calc_autocorr_block(GpuTensor<real, 2> ac_block, const GpuTensor<real, 4> spinwait, const GpuTensor<real, 3> emom);


__global__ void calc_autocorr_final(GpuTensor<real, 2> ac_block, GpuTensor<real, 2> ac, const real norm, const int ac_count, const int numBlocks);
