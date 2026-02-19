#pragma once

#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include <numeric>
#include <hip/hip_runtime.h>
#include <hip/hip_cooperative_groups.h>
#include <hip/hip_complex.h>

namespace cg = cooperative_groups;
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// Optimized register-based warp reduction for raw real/imaginary components
// Modifies sum_re and sum_im in-place, avoiding complex object overhead

__inline__ __device__
void warpReduceSum(real &sum_re, real &sum_im)
{
    #pragma unroll
    for (int offset = warpSize / 2; offset > 0; offset >>= 1)
    {
        sum_re += __shfl_down(sum_re, offset);
        sum_im += __shfl_down(sum_im, offset);
    }
}

// Complex multiply-accumulate: accumulate (a + bi) * (c + di) into (sum_re + sum_im*i)
// Result: sum_re' = a*c - b*d + sum_re, sum_im' = a*d + b*c + sum_im
__inline__ __device__
void complexMulAdd(real &sum_re, real &sum_im, real a_re, real a_im, real b_re, real b_im) {
    sum_re += a_re * b_re - a_im * b_im;
    sum_im += a_re * b_im + a_im * b_re;
}

// Complex scalar multiply-accumulate (for when imaginary part is zero)
// Accumulate (a + 0*i) * (c + d*i) = a*c + a*d*i
__inline__ __device__
void complexScalarMulAdd(real &sum_re, real &sum_im, real a, real c, real d) {
    sum_re += a * c;
    sum_im += a * d;
}

__device__ real sc_window_fac(int sc_window_fun, unsigned int step, unsigned int nstep);

template <size_t dim>
__global__ void setZero(GpuTensor<gpu_complex, dim> sc, unsigned int size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        sc[idx] = gpu_complex(0, 0);;
    }
}

__global__ void GPUSqSum(const GpuTensor<real, 3> spin, const GpuTensor<real, 2> coord, const GpuTensor<real, 2> q, const GpuTensor<real, 1> r_mid, GpuTensor<gpu_complex, 2> scblock, int tasks, unsigned int N); 

__global__ void GPUSqFinalSum_stat(GpuTensor<gpu_complex, 2> scblock, GpuTensor<gpu_complex, 2> scsum, int numBlocks);

__global__ void GPUSqFinalSum_dyn(GpuTensor<gpu_complex, 2> scblock, GpuTensor<gpu_complex, 3> scsum, int numBlocks, unsigned int t_cur);

__global__ void GPUSqFinalSum_both(GpuTensor<gpu_complex, 2> scblock, GpuTensor<gpu_complex, 2> scsum_q, GpuTensor<gpu_complex, 3> scsum_qt, int numBlocks, unsigned int t_cur, unsigned int both_flag);

__global__ void GPUSwSum(const GpuTensor<gpu_complex, 3> sq, const GpuTensor<real, 1> dt, const GpuTensor<real, 1> w, GpuTensor<gpu_complex, 3> scblock, int tasks, unsigned int tSize, unsigned int nq, int sc_max_nstep, int sc_window_fun);

__global__ void GPUSwFinalSum(GpuTensor<gpu_complex, 3> scblock, GpuTensor<gpu_complex, 3> scsum, int numBlocks, int nq);
