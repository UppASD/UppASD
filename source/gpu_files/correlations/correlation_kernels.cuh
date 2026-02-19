#pragma once

#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include <numeric>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#if defined(USE_CUCOMPLEX_CORRELATIONS)
#include <cuComplex.h>
#else
#include <thrust/complex.h>
#endif
#include <curand.h>
#include <cuda.h>

// Type alias for correlation complex numbers (must match gpuStructures.hpp)
#if defined(USE_CUCOMPLEX_CORRELATIONS)
using complex_type = cuDoubleComplex;
#else
using complex_type = thrust::complex<real>;
#endif

namespace cg = cooperative_groups;
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// Optimized register-based warp reduction for raw real/imaginary components
// Modifies sum_re and sum_im in-place, avoiding complex object overhead
__inline__ __device__
void warpReduceSum(real &sum_re, real &sum_im) {
    const unsigned int FULL_MASK = 0xffffffff;
    
    #pragma unroll
    for (int offset = warpSize / 2; offset > 0; offset /= 2) {
#if CUDA_VERSION < 9000
        real shfl_re = __shfl_down(sum_re, offset);
        real shfl_im = __shfl_down(sum_im, offset);
#else
        real shfl_re = __shfl_down_sync(FULL_MASK, sum_re, offset);
        real shfl_im = __shfl_down_sync(FULL_MASK, sum_im, offset);
#endif
        sum_re += shfl_re;
        sum_im += shfl_im;
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
__global__ void setZero(GpuTensor<complex_type, dim> sc, unsigned int size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
#if defined(USE_CUCOMPLEX_CORRELATIONS)
        sc[idx] = make_cuDoubleComplex(0.0, 0.0);
#else
        sc[idx] = thrust::complex<real>(0, 0);
#endif
    }
}

__global__ void GPUSqSum(const GpuTensor<real, 3> spin, const GpuTensor<real, 2> coord, const GpuTensor<real, 2> q, const GpuTensor<real, 1> r_mid, GpuTensor<complex_type, 2> scblock, int tasks, unsigned int N); 

__global__ void GPUSqFinalSum_stat(GpuTensor<complex_type, 2> scblock, GpuTensor<complex_type, 2> scsum, int numBlocks);

__global__ void GPUSqFinalSum_dyn(GpuTensor<complex_type, 2> scblock, GpuTensor<complex_type, 3> scsum, int numBlocks, unsigned int t_cur);

__global__ void GPUSqFinalSum_both(GpuTensor<complex_type, 2> scblock, GpuTensor<complex_type, 2> scsum_q, GpuTensor<complex_type, 3> scsum_qt, int numBlocks, unsigned int t_cur, unsigned int both_flag);

__global__ void GPUSwSum(const GpuTensor<complex_type, 3> sq, const GpuTensor<real, 1> dt, const GpuTensor<real, 1> w, GpuTensor<complex_type, 3> scblock, int tasks, unsigned int tSize, unsigned int nq, int sc_max_nstep, int sc_window_fun);

__global__ void GPUSwFinalSum(GpuTensor<complex_type, 3> scblock, GpuTensor<complex_type, 3> scsum, int numBlocks, int nq);
