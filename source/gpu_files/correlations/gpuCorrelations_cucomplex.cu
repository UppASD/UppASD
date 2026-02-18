#pragma once

#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "gpuStructures.hpp"
#include "gpuCorrelations.cuh"
#include "fortranData.hpp"
#include <numeric>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <thrust/complex.h>
#include <curand.h>
#include <cuda.h>

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


__device__ real sc_window_fac(int sc_window_fun, unsigned int step, unsigned int nstep) {
    real dum = 1.0;
    switch (sc_window_fun) {
        //Hann
    case 2:
        dum = (0.5 - 0.5 * cos(2.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)));
            //Hamming
    case 3:
        dum = (0.54 - 0.46 * cos(2.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)));
            //Hamming v2
    case 32:
        dum = (0.53836 - 0.46164 * cos(2.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)));
            //Blackman - Harris
    case 4:
        dum =
            (0.35785 - 0.48829 * cos(2.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)) +
                0.14128 * cos(4.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)) -
                0.01168 * cos(6.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)));
            //Nuttal
    case 5:
        dum =
            (0.355768 - 0.478396 * cos(2.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)) +
                0.144232 * cos(4.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)) -
                0.012604 * cos(6.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)));
    }
 
        
    return dum;
}

template <size_t dim>
__global__ void setZero(GpuTensor<thrust::complex<real>, dim> sc, unsigned int size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        sc[idx] = thrust::complex<real>(0, 0);;
    }
}

__global__ void GPUSqSum(const GpuTensor<real, 3> spin, const GpuTensor<real, 2> coord, const GpuTensor<real, 2> q, const GpuTensor<real, 1> r_mid, GpuTensor<thrust::complex<real>, 2> scblock, int tasks, unsigned int N) {
    auto grid = cg::this_grid();
    auto block = cg::this_thread_block();
    auto warp = cg::tiled_partition<32>(block);

    int lane = warp.thread_rank();
    int wid = warp.meta_group_rank();
    int wSize = warp.size();
    int wNum = warp.meta_group_size();
    int tid = grid.thread_rank();
    int tid_in_block = block.thread_rank();

    int qInd = grid.block_index().y;
    int tid_in_X = grid.block_index().x * block.num_threads() + tid_in_block;
    int stride = grid.dim_blocks().x * block.num_threads();

    // Register-based accumulation: store real and imaginary parts separately (reduced register pressure)
    real sum_re[3] = {0.0, 0.0, 0.0};
    real sum_im[3] = {0.0, 0.0, 0.0};
    
    unsigned int rInd, mInd, cInd, ii;
    real inv_N = 1.0 / N;
    __shared__ real shared_re[3][32];
    __shared__ real shared_im[3][32];

    // Main computation loop over tasks (unroll for better performance)
    #pragma unroll 4
    for (int id = tid_in_X; id < tasks; id += stride) {
        ii = id / 3;
        cInd = id % 3;
        rInd = ii % N;
        mInd = ii / N; 

        // Calculate phase: 2 * PI * q · (r - r_mid)
        real phase = 2.0 * M_PI * (q(0, qInd) * (coord(0, rInd) - r_mid(0)) + 
                                   q(1, qInd) * (coord(1, rInd) - r_mid(1)) + 
                                   q(2, qInd) * (coord(2, rInd) - r_mid(2)));

        real s, c;
        sincos(phase, &s, &c);

        real spin_val = spin(cInd, rInd, mInd) * inv_N;

        // Accumulate using optimized scalar multiply: (spin_val + 0i) * (c + si)
        complexScalarMulAdd(sum_re[cInd], sum_im[cInd], spin_val, c, s);
    }
    // #pragma unroll 4
    // for (int id = tid_in_X; id < tasks; id += stride) {
    //     ii = id / 3;
    //     cInd = id % 3;
    //     rInd = ii % N;
    //     mInd = ii / N; 
    //     qdr = q(0, qInd) * (coord(0, rInd) - r_mid(0)) + q(1, qInd) * (coord(1, rInd) - r_mid(1)) + q(2, qInd) * (coord(2, rInd) - r_mid(2));



    //     mySum[cInd] += thrust::complex<real>(spin(cInd, rInd, mInd), 0) * thrust::exp(iqfac * thrust::complex<real>(qdr, 0))/N;
    //     //printf("qInd = %i, Re = %.3lf, Im = %.3lf,qdr = %.3lf\n", qInd, thrust::complex<real>(qdr, 0).real(), thrust::complex<real>(qdr, 0).imag(), qdr);
    //     //printf("0:rmid = %i, q = %.3lf, Im = %.3lf,qdr = %.3lf\n", qInd, thrust::complex<real>(qdr, 0).real(), thrust::complex<real>(qdr, 0).imag(), qdr);

    //     //  printf("tid = %i, mInd = %i, stride = %i, data_id = %i, mySum = %.3f\n", tid, mInd , stride, id + offsetM, mySum[id % 3]);
    // }
    warp.sync();

    // Warp-level reduction using register-based functions
    warpReduceSum(sum_re[0], sum_im[0]);
    warpReduceSum(sum_re[1], sum_im[1]);
    warpReduceSum(sum_re[2], sum_im[2]);

    // Store warp results to shared memory
    if (lane == 0) {
        shared_re[0][wid] = sum_re[0];
        shared_im[0][wid] = sum_im[0];
        shared_re[1][wid] = sum_re[1];
        shared_im[1][wid] = sum_im[1];
        shared_re[2][wid] = sum_re[2];
        shared_im[2][wid] = sum_im[2];
    }

    __syncthreads();              // Wait for all partial reductions
    
    // Load results from shared memory for final warp reduction
    sum_re[0] = (tid_in_block < wNum) ? shared_re[0][lane] : 0;
    sum_im[0] = (tid_in_block < wNum) ? shared_im[0][lane] : 0;
    sum_re[1] = (tid_in_block < wNum) ? shared_re[1][lane] : 0;
    sum_im[1] = (tid_in_block < wNum) ? shared_im[1][lane] : 0;
    sum_re[2] = (tid_in_block < wNum) ? shared_re[2][lane] : 0;
    sum_im[2] = (tid_in_block < wNum) ? shared_im[2][lane] : 0;

    // Final reduction in first warp
    if (wid == 0) {
        warpReduceSum(sum_re[0], sum_im[0]);
        warpReduceSum(sum_re[1], sum_im[1]);
        warpReduceSum(sum_re[2], sum_im[2]);
    }

    // Reconstruct complex objects and write only at final step
    if (tid_in_block == 0) {
        scblock(3 * block.group_index().x + 0, qInd) = thrust::complex<real>(sum_re[0], sum_im[0]);
        scblock(3 * block.group_index().x + 1, qInd) = thrust::complex<real>(sum_re[1], sum_im[1]);
        scblock(3 * block.group_index().x + 2, qInd) = thrust::complex<real>(sum_re[2], sum_im[2]);
    }
}
__global__ void GPUSqFinalSum_stat(GpuTensor<thrust::complex<real>, 2> scblock, GpuTensor<thrust::complex<real>, 2> scsum, int numBlocks)
{
    auto grid = cg::this_grid();
    auto block = cg::this_thread_block();
    auto warp = cg::tiled_partition<32>(block);

    int lane = warp.thread_rank();
    int wid = warp.meta_group_rank();
    int wSize = warp.size();
    int wNum = warp.meta_group_size();
    int tid = grid.thread_rank();
    int tNum = block.size();
    int tid_in_block = block.thread_rank();

    int qInd = grid.block_index().x;
    int tid_in_Q = tid_in_block;

    // Register-based accumulators
    real sum_re[3] = {0.0, 0.0, 0.0};
    real sum_im[3] = {0.0, 0.0, 0.0};
    
    __shared__ real shared_re[3][32];
    __shared__ real shared_im[3][32];

    if (tid_in_Q < numBlocks) {
        for (int k = 0; k < 3; k++) {
            thrust::complex<real> val = scblock(3 * tid_in_Q + k, qInd);
            sum_re[k] += val.real();
            sum_im[k] += val.imag();
        }
    }

    warp.sync();

    // Warp-level reduction
    warpReduceSum(sum_re[0], sum_im[0]);
    warpReduceSum(sum_re[1], sum_im[1]);
    warpReduceSum(sum_re[2], sum_im[2]);

    if (lane == 0) {
        shared_re[0][wid] = sum_re[0];
        shared_im[0][wid] = sum_im[0];
        shared_re[1][wid] = sum_re[1];
        shared_im[1][wid] = sum_im[1];
        shared_re[2][wid] = sum_re[2];
        shared_im[2][wid] = sum_im[2];
    }

    __syncthreads();
    sum_re[0] = (tid_in_block < wNum) ? shared_re[0][lane] : 0;
    sum_im[0] = (tid_in_block < wNum) ? shared_im[0][lane] : 0;
    sum_re[1] = (tid_in_block < wNum) ? shared_re[1][lane] : 0;
    sum_im[1] = (tid_in_block < wNum) ? shared_im[1][lane] : 0;
    sum_re[2] = (tid_in_block < wNum) ? shared_re[2][lane] : 0;
    sum_im[2] = (tid_in_block < wNum) ? shared_im[2][lane] : 0;
    
    if (wid == 0) {
        warpReduceSum(sum_re[0], sum_im[0]);
        warpReduceSum(sum_re[1], sum_im[1]);
        warpReduceSum(sum_re[2], sum_im[2]);
    }

    if (tid_in_block == 0) {
        scsum(0, qInd) += thrust::complex<real>(sum_re[0], sum_im[0]);
        scsum(1, qInd) += thrust::complex<real>(sum_re[1], sum_im[1]);
        scsum(2, qInd) += thrust::complex<real>(sum_re[2], sum_im[2]);

        /*mblock_gpu[block.group_index().x] += mySum[0];
        mblock_gpu[block.group_index().x + grid.group_dim().x] += mySum[1];
        mblock_gpu[block.group_index().x + 2 * grid.group_dim().x] += mySum[2];*/
       // printf("qInd = %i, mblock0 = %lf, mblock1 = %lf, mblock2 = %lf\n", mInd, msum(0, curstep, mInd), msum(1, curstep, mInd), msum(2, curstep, mInd));
    }
}

__global__ void GPUSqFinalSum_dyn(GpuTensor<thrust::complex<real>, 2> scblock, GpuTensor<thrust::complex<real>, 3> scsum, int numBlocks, unsigned int t_cur)
{
    auto grid = cg::this_grid();
    auto block = cg::this_thread_block();
    auto warp = cg::tiled_partition<32>(block);

    int lane = warp.thread_rank();
    int wid = warp.meta_group_rank();
    int wSize = warp.size();
    int wNum = warp.meta_group_size();
    int tid = grid.thread_rank();
    int tNum = block.size();
    int tid_in_block = block.thread_rank();

    int qInd = grid.block_index().x;
    int tid_in_Q = tid_in_block;

    // Register-based accumulators
    real sum_re[3] = {0.0, 0.0, 0.0};
    real sum_im[3] = {0.0, 0.0, 0.0};
    
    __shared__ real shared_re[3][32];
    __shared__ real shared_im[3][32];

    if (tid_in_Q < numBlocks) {
        for (int k = 0; k < 3; k++) {
            thrust::complex<real> val = scblock(3 * tid_in_Q + k, qInd);
            sum_re[k] += val.real();
            sum_im[k] += val.imag();
        }
    }

    warp.sync();

    // Warp-level reduction
    warpReduceSum(sum_re[0], sum_im[0]);
    warpReduceSum(sum_re[1], sum_im[1]);
    warpReduceSum(sum_re[2], sum_im[2]);

    if (lane == 0) {
        shared_re[0][wid] = sum_re[0];
        shared_im[0][wid] = sum_im[0];
        shared_re[1][wid] = sum_re[1];
        shared_im[1][wid] = sum_im[1];
        shared_re[2][wid] = sum_re[2];
        shared_im[2][wid] = sum_im[2];
    }

    __syncthreads();
    sum_re[0] = (tid_in_block < wNum) ? shared_re[0][lane] : 0;
    sum_im[0] = (tid_in_block < wNum) ? shared_im[0][lane] : 0;
    sum_re[1] = (tid_in_block < wNum) ? shared_re[1][lane] : 0;
    sum_im[1] = (tid_in_block < wNum) ? shared_im[1][lane] : 0;
    sum_re[2] = (tid_in_block < wNum) ? shared_re[2][lane] : 0;
    sum_im[2] = (tid_in_block < wNum) ? shared_im[2][lane] : 0;
    
    if (wid == 0) {
        warpReduceSum(sum_re[0], sum_im[0]);
        warpReduceSum(sum_re[1], sum_im[1]);
        warpReduceSum(sum_re[2], sum_im[2]);
    }

    if (tid_in_block == 0) {
            // sc_qt_gpu is now (3, nq, sc_max_nstep) so index is (component, qInd, t_cur)
            scsum(0, qInd, t_cur) += thrust::complex<real>(sum_re[0], sum_im[0]);
            scsum(1, qInd, t_cur) += thrust::complex<real>(sum_re[1], sum_im[1]);
            scsum(2, qInd, t_cur) += thrust::complex<real>(sum_re[2], sum_im[2]);
        
        //printf("qInd = %i, t_cur = %i, Re = %.3lf, Im = %.3lf\n", qInd, t_cur, mySum[2].real(), mySum[2].imag());

        /*mblock_gpu[block.group_index().x] += mySum[0];
        mblock_gpu[block.group_index().x + grid.group_dim().x] += mySum[1];
        mblock_gpu[block.group_index().x + 2 * grid.group_dim().x] += mySum[2];*/
        // printf("qInd = %i, mblock0 = %lf, mblock1 = %lf, mblock2 = %lf\n", mInd, msum(0, curstep, mInd), msum(1, curstep, mInd), msum(2, curstep, mInd));
    }
}

__global__ void GPUSqFinalSum_both(GpuTensor<thrust::complex<real>, 2> scblock, GpuTensor<thrust::complex<real>, 2> scsum_q, GpuTensor<thrust::complex<real>, 3> scsum_qt, int numBlocks, unsigned int t_cur, unsigned int both_flag)
{
    auto grid = cg::this_grid();
    auto block = cg::this_thread_block();
    auto warp = cg::tiled_partition<32>(block);

    int lane = warp.thread_rank();
    int wid = warp.meta_group_rank();
    int wSize = warp.size();
    int wNum = warp.meta_group_size();
    int tid = grid.thread_rank();
    int tNum = block.size();
    int tid_in_block = block.thread_rank();

    int qInd = grid.block_index().x;
    int tid_in_Q = tid_in_block;

    // Register-based accumulators
    real sum_re[3] = {0.0, 0.0, 0.0};
    real sum_im[3] = {0.0, 0.0, 0.0};
    
    __shared__ real shared_re[3][32];
    __shared__ real shared_im[3][32];

    if (tid_in_Q < numBlocks) {
        for (int k = 0; k < 3; k++) {
            thrust::complex<real> val = scblock(3 * tid_in_Q + k, qInd);
            sum_re[k] += val.real();
            sum_im[k] += val.imag();
        }
    }

    warp.sync();

    // Warp-level reduction
    warpReduceSum(sum_re[0], sum_im[0]);
    warpReduceSum(sum_re[1], sum_im[1]);
    warpReduceSum(sum_re[2], sum_im[2]);

    if (lane == 0) {
        shared_re[0][wid] = sum_re[0];
        shared_im[0][wid] = sum_im[0];
        shared_re[1][wid] = sum_re[1];
        shared_im[1][wid] = sum_im[1];
        shared_re[2][wid] = sum_re[2];
        shared_im[2][wid] = sum_im[2];
    }

    __syncthreads();
    sum_re[0] = (tid_in_block < wNum) ? shared_re[0][lane] : 0;
    sum_im[0] = (tid_in_block < wNum) ? shared_im[0][lane] : 0;
    sum_re[1] = (tid_in_block < wNum) ? shared_re[1][lane] : 0;
    sum_im[1] = (tid_in_block < wNum) ? shared_im[1][lane] : 0;
    sum_re[2] = (tid_in_block < wNum) ? shared_re[2][lane] : 0;
    sum_im[2] = (tid_in_block < wNum) ? shared_im[2][lane] : 0;
    
    if (wid == 0) {
        warpReduceSum(sum_re[0], sum_im[0]);
        warpReduceSum(sum_re[1], sum_im[1]);
        warpReduceSum(sum_re[2], sum_im[2]);
    }

    if (tid_in_block == 0) {
        if ((both_flag == 1) || (both_flag == 2)) {
            scsum_qt(0, qInd, t_cur) += thrust::complex<real>(sum_re[0], sum_im[0]);
            scsum_qt(1, qInd, t_cur) += thrust::complex<real>(sum_re[1], sum_im[1]);
            scsum_qt(2, qInd, t_cur) += thrust::complex<real>(sum_re[2], sum_im[2]);
        }

        if ((both_flag == 0) || (both_flag == 2)) {
            scsum_q(0, qInd) += thrust::complex<real>(sum_re[0], sum_im[0]);
            scsum_q(1, qInd) += thrust::complex<real>(sum_re[1], sum_im[1]);
            scsum_q(2, qInd) += thrust::complex<real>(sum_re[2], sum_im[2]);
        }

        //printf("Re = %.3lf, Im = %.3lf\n", mySum[1].real(), mySum[1].imag());

        /*mblock_gpu[block.group_index().x] += mySum[0];
        mblock_gpu[block.group_index().x + grid.group_dim().x] += mySum[1];
        mblock_gpu[block.group_index().x + 2 * grid.group_dim().x] += mySum[2];*/
        // printf("qInd = %i, mblock0 = %lf, mblock1 = %lf, mblock2 = %lf\n", mInd, msum(0, curstep, mInd), msum(1, curstep, mInd), msum(2, curstep, mInd));
    }
}

__global__ void GPUSwSum(const GpuTensor<thrust::complex<real>, 3> sq, const GpuTensor<real, 1> dt, const GpuTensor<real, 1> w, GpuTensor<thrust::complex<real>, 3> scblock, int tasks, unsigned int tSize, unsigned int nq, int sc_max_nstep, int sc_window_fun) {
    auto grid = cg::this_grid();
    auto block = cg::this_thread_block();
    auto warp = cg::tiled_partition<32>(block);

    int lane = warp.thread_rank();
    int wid = warp.meta_group_rank();
    int wSize = warp.size();
    int wNum = warp.meta_group_size();
    int tid = grid.thread_rank();
    int tid_in_block = block.thread_rank();

    int qInd = grid.block_index().y%nq;
    int wInd = grid.block_index().y/nq;

    int tid_in_X = grid.block_index().x * block.num_threads() + tid_in_block;
    int stride = grid.dim_blocks().x * block.num_threads();

    // Register-based accumulation: store real and imaginary parts separately
    real sum_re[3] = {0.0, 0.0, 0.0};
    real sum_im[3] = {0.0, 0.0, 0.0};
    
    unsigned int tInd, cInd, ii;
    __shared__ real shared_re[3][32];
    __shared__ real shared_im[3][32];

    // Fourier transform loop: unroll for performance
    #pragma unroll 2
    for (int id = tid_in_X; id < tasks; id += stride) {
        tInd = id / 3;
        cInd = id % 3;

        // 1. Calculate the real-valued phase: Phase = t * dt(t) * w(w)
        real phase = (real)tInd * dt(tInd) * w(wInd);

        // 2. Hardware-accelerated sine and cosine
        real s, c;
        sincos(phase, &s, &c);

        // 3. Extract S(q,t) values
        thrust::complex<real> sq_val = sq(cInd, qInd, tInd);
        real sq_re = sq_val.real();
        real sq_im = sq_val.imag();
        
        // 4. Windowing function
        real win = sc_window_fac(sc_window_fun, (tInd - 1), sc_max_nstep);

        // 5. Complex multiply-accumulate: (sq_re + sq_im*i) * (c + s*i) * win
        // Result = ((sq_re*c - sq_im*s) + i*(sq_re*s + sq_im*c)) * win
        complexMulAdd(sum_re[cInd], sum_im[cInd], sq_re, sq_im, c*win, s*win);
    }
    // #pragma unroll 2
    // for (int id = tid_in_X; id < tasks; id += stride) {
    //     tInd = id / 3;
    //     cInd = id % 3;
    //     //qdr = q(0, qInd) * (coord(0, rInd) - r_mid(0)) + q(1, qInd) * (coord(1, rInd) - r_mid(1)) + q(2, qInd) * (coord(2, rInd) - r_mid(2));
    //     //wfac = 1.0_dblprec
    //     //corr_kw = 0.0_dblprec
    //     //i = (0.0_dblprec, 1.0_dblprec)
    //      //tidx = cc % sc_max_nstep
    //     //tt = i * wfac * (step - 1) * dt(step)            
    //     //epowwt = exp(cc % w(iw) * tt) * sc_window_fac(sc_window_fun, step, tidx)
    //     tw = thrust::complex<real>(0, 1) * tInd * dt(tInd)*w(wInd);
    //     mySum[cInd] += thrust::exp(tw) * sc_window_fac(sc_window_fun, (tInd-1), sc_max_nstep)*sq(cInd, qInd, tInd);
    //     //mySum[cInd] += thrust::complex<real>(spin(cInd, rInd, mInd), 0) * thrust::exp(iqfac * thrust::complex<real>(qdr, 0)) / N;
    //     //printf("qInd = %i, Re = %.3lf, Im = %.3lf,qdr = %.3lf\n", qInd, thrust::complex<real>(qdr, 0).real(), thrust::complex<real>(qdr, 0).imag(), qdr);
    //     //printf("0:rmid = %i, q = %.3lf, Im = %.3lf,qdr = %.3lf\n", qInd, thrust::complex<real>(qdr, 0).real(), thrust::complex<real>(qdr, 0).imag(), qdr);

    //     //  printf("tid = %i, mInd = %i, stride = %i, data_id = %i, mySum = %.3f\n", tid, mInd , stride, id + offsetM, mySum[id % 3]);
    // }
    warp.sync();

    // Warp-level reduction using register-based functions
    warpReduceSum(sum_re[0], sum_im[0]);
    warpReduceSum(sum_re[1], sum_im[1]);
    warpReduceSum(sum_re[2], sum_im[2]);

    // Store warp results to shared memory
    if (lane == 0) {
        shared_re[0][wid] = sum_re[0];
        shared_im[0][wid] = sum_im[0];
        shared_re[1][wid] = sum_re[1];
        shared_im[1][wid] = sum_im[1];
        shared_re[2][wid] = sum_re[2];
        shared_im[2][wid] = sum_im[2];
    }

    __syncthreads();              // Wait for all partial reductions
    
    // Load results from shared memory for final warp reduction
    sum_re[0] = (tid_in_block < wNum) ? shared_re[0][lane] : 0;
    sum_im[0] = (tid_in_block < wNum) ? shared_im[0][lane] : 0;
    sum_re[1] = (tid_in_block < wNum) ? shared_re[1][lane] : 0;
    sum_im[1] = (tid_in_block < wNum) ? shared_im[1][lane] : 0;
    sum_re[2] = (tid_in_block < wNum) ? shared_re[2][lane] : 0;
    sum_im[2] = (tid_in_block < wNum) ? shared_im[2][lane] : 0;

    // Final reduction in first warp
    if (wid == 0) {
        warpReduceSum(sum_re[0], sum_im[0]);
        warpReduceSum(sum_re[1], sum_im[1]);
        warpReduceSum(sum_re[2], sum_im[2]);
    }

    // Reconstruct complex objects and write only at final step
    if (tid_in_block == 0) {
        scblock(3 * block.group_index().x + 0, qInd, wInd) = thrust::complex<real>(sum_re[0], sum_im[0]);
        scblock(3 * block.group_index().x + 1, qInd, wInd) = thrust::complex<real>(sum_re[1], sum_im[1]);
        scblock(3 * block.group_index().x + 2, qInd, wInd) = thrust::complex<real>(sum_re[2], sum_im[2]);
    }
}

__global__ void GPUSwFinalSum(GpuTensor<thrust::complex<real>, 3> scblock, GpuTensor<thrust::complex<real>, 3> scsum, int numBlocks, int nq)
{
    auto grid = cg::this_grid();
    auto block = cg::this_thread_block();
    auto warp = cg::tiled_partition<32>(block);

    int lane = warp.thread_rank();
    int wid = warp.meta_group_rank();
    int wSize = warp.size();
    int wNum = warp.meta_group_size();
    int tid = grid.thread_rank();
    int tNum = block.size();
    int tid_in_block = block.thread_rank();

    int qInd = (grid.block_index().x) % nq;
    int wInd = (grid.block_index().x) / nq;
    int tid_in_Q = tid_in_block;

    // Register-based accumulators for intermediate block results
    real sum_re[3] = {0.0, 0.0, 0.0};
    real sum_im[3] = {0.0, 0.0, 0.0};
    
    __shared__ real shared_re[3][32];
    __shared__ real shared_im[3][32];

    if (tid_in_Q < numBlocks) {
        // Load and sum block results
        for (int k = 0; k < 3; k++) {
            thrust::complex<real> val = scblock(3 * tid_in_Q + k, qInd, wInd);
            sum_re[k] += val.real();
            sum_im[k] += val.imag();
        }
    }

    warp.sync();

    // Warp-level reduction using register-based functions
    warpReduceSum(sum_re[0], sum_im[0]);
    warpReduceSum(sum_re[1], sum_im[1]);
    warpReduceSum(sum_re[2], sum_im[2]);

    // Store warp results to shared memory
    if (lane == 0) {
        shared_re[0][wid] = sum_re[0];
        shared_im[0][wid] = sum_im[0];
        shared_re[1][wid] = sum_re[1];
        shared_im[1][wid] = sum_im[1];
        shared_re[2][wid] = sum_re[2];
        shared_im[2][wid] = sum_im[2];
    }

    __syncthreads();              // Wait for all partial reductions
    
    // Load results from shared memory for final warp reduction
    sum_re[0] = (tid_in_block < wNum) ? shared_re[0][lane] : 0;
    sum_im[0] = (tid_in_block < wNum) ? shared_im[0][lane] : 0;
    sum_re[1] = (tid_in_block < wNum) ? shared_re[1][lane] : 0;
    sum_im[1] = (tid_in_block < wNum) ? shared_im[1][lane] : 0;
    sum_re[2] = (tid_in_block < wNum) ? shared_re[2][lane] : 0;
    sum_im[2] = (tid_in_block < wNum) ? shared_im[2][lane] : 0;

    // Final reduction in first warp
    if (wid == 0) {
        warpReduceSum(sum_re[0], sum_im[0]);
        warpReduceSum(sum_re[1], sum_im[1]);
        warpReduceSum(sum_re[2], sum_im[2]);
    }

    // Accumulate results (reconstruct complex only at final write)
    if (tid_in_block == 0) {
        scsum(0, qInd, wInd) += thrust::complex<real>(sum_re[0], sum_im[0]);
        scsum(1, qInd, wInd) += thrust::complex<real>(sum_re[1], sum_im[1]);
        scsum(2, qInd, wInd) += thrust::complex<real>(sum_re[2], sum_im[2]);
    }
}
__global__ void GPUSqAvrg(GpuTensor<thrust::complex<real>, 2> sc, int n_steps, int tasks, int M) {
    auto grid = cg::this_grid();
    int tid = grid.thread_rank();
    if (tid < tasks) {
        // Do NOT normalize here - let Fortran handle normalization
        // sc[tid] = sc[tid] / (n_steps*M);
    }
}
// Constructor
GpuCorrelations::GpuCorrelations(const Flag Flags, const SimulationParameters SimParam, const deviceLattice& gpuLattice, const hostCorrelations& cpuCorrelations)
: emomM(gpuLattice.emomM)
, emom(gpuLattice.emom)
, mmom(gpuLattice.mmom) {

    isallocated = 0; 
    if(!initiate(Flags, SimParam, cpuCorrelations)) {  
      std::fprintf(stderr, "GpuCorrelations: correlations failed to initiate!\n");
      return;
   }
}
// Destructor
GpuCorrelations::~GpuCorrelations() {
    release();
}
// Initiator
bool GpuCorrelations::initiate(const Flag Flags, const SimulationParameters SimParam, const hostCorrelations& cpuCorrelations) {
    // Assert that we're not already initialized
    //release();

    // Parameters
    if(Flags.do_gpu_correlations){

        N = SimParam.N;
        M = SimParam.M;
        nq = SimParam.nq;
        sc_max_nstep = SimParam.sc_max_nstep;
        sc_window_fun = SimParam.sc_window_fun;
        nw = SimParam.nw;
        delta_t = SimParam.delta_t;
        t_cur = 0;
        n_samples = 0;
        do_sc = Flags.do_sc;
        sc_sep = SimParam.sc_sep;
        sc_step = SimParam.sc_step;
        // nainv = 1 / N;
        // Blocks and threads
        maxThreads = 512;
        maxBlocks = 1024; 
        tasksTot_q = 3 * N * M;
        tasksTot_w = 3 * sc_max_nstep;
        // maxBlocks = 1023; //must be devidable by 3, less than 1024
        numThreads = maxThreads;
        //numBlocks = std::min(((3 * ((spinTot + 2) / 3) + numThreads - 1) / numThreads), maxBlocks);
        numBlocksX_q = std::min(((tasksTot_q + numThreads - 1) / numThreads), maxBlocks);
        numBlocksX_w = std::min(((tasksTot_w + numThreads - 1) / numThreads), maxBlocks);
        numBlocksY_q = nq;
        numBlocksY_w = nq*nw;//TODO
        blocks_q = { numBlocksX_q, numBlocksY_q, 1 };
        blocks_w = { numBlocksX_w, numBlocksY_w, 1 };//TODO
        threads = { numThreads, 1, 1 };
        //printf("numBlocks = %i\n", numBlocksX_q);

        //iqfac = thrust::complex<real>(0, 2 * M_PI);
        r_mid.Allocate(static_cast <long int>(3));
        q.Allocate(static_cast <long int>(3), static_cast <long int>(nq));
        coord.Allocate(static_cast <long int>(3), static_cast <long int>(N));

        r_mid.copy_sync(cpuCorrelations.r_mid);
        q.copy_sync(cpuCorrelations.q);
        coord.copy_sync(cpuCorrelations.coord);
        int bl;

        sc_block_gpu.Allocate(static_cast <long int>(3 * numBlocksX_q), static_cast <long int>(nq));
        if ((do_sc == 'C') || (do_sc == 'Y')) {
            sc_q_gpu.Allocate(static_cast <long int>(3), static_cast <long int>(nq));
            sc_q_cpu.AllocateHost(static_cast <long int>(3), static_cast <long int>(nq));
            bl = (3 * nq + numThreads - 1) / numThreads;
            setZero<2> << <bl, numThreads >> > (sc_q_gpu, 3 * nq);

        }
        if ((do_sc == 'Q') || (do_sc == 'Y')) {
            // CRITICAL: Match Fortran memory layout: (component, q, time) NOT (component, time, q)
            sc_qt_gpu.Allocate(static_cast <long int>(3), static_cast <long int>(nq), static_cast <long int>(sc_max_nstep));
            sc_qw_gpu.Allocate(static_cast <long int>(3), static_cast <long int>(nq), static_cast <long int>(nw));
            sc_block_w_gpu.Allocate(static_cast <long int>(3 * numBlocksX_w), static_cast <long int>(nq), static_cast <long int>(nw));
            dt.Allocate(static_cast <long int>(sc_max_nstep));
            dt_cpu.AllocateHost(static_cast <long int>(sc_max_nstep));
            sc_step_arr_cpu.AllocateHost(static_cast <long int>(sc_max_nstep));
            w.Allocate(static_cast <long int>(nw));
            //dt.copy_sync(cpuCorrelations.dt);const deviceLattice& gpuLattice, const int curstep
            w.copy_sync(cpuCorrelations.w);


            bl = (3 * nq* sc_max_nstep + numThreads - 1) / numThreads;
            setZero<3> << <bl, numThreads >> > (sc_qt_gpu, 3 * nq* sc_max_nstep);
            bl = (3 * nq * nw + numThreads - 1) / numThreads;
            setZero<3> << <bl, numThreads >> > (sc_qw_gpu, 3 * nq*nw);
            //bl = (3 * nq * sc_max_nstep + numThreads - 1) / numThreads;
        }

        //mbuff_gpu.Allocate(static_cast <long int>(3), static_cast <long int>(avrg_buff), static_cast <long int>(M));
        isallocated = 1;
        bl = (3 * numBlocksX_q * nq + numThreads - 1) / numThreads;
        setZero<2> << <bl, numThreads >> > (sc_block_gpu, 3 * numBlocksX_q * nq);
        //sc_block_gpu.zeros();
        //sc_gpu.zeros();
    }

    // All initialized?
    if (cudaDeviceSynchronize() != cudaSuccess) {
        release();
        return false;
    }

    return true;
}
void GpuCorrelations::release() {
    if (isallocated) {
        r_mid.Free();
        coord.Free();
        q.Free();
        sc_block_gpu.Free();
        if ((do_sc == 'C') || (do_sc == 'Y')) {
            sc_q_gpu.Free();
            sc_q_cpu.FreeHost();
        }
        if ((do_sc == 'Q') || (do_sc == 'Y')) {
            sc_qt_gpu.Free();
            sc_qw_gpu.Free();
            sc_block_w_gpu.Free();
            w.Free();
            dt.Free();
            dt_cpu.FreeHost();
            sc_step_arr_cpu.FreeHost();
        }
        isallocated = 0;
    }

}

void GpuCorrelations::measure(std::size_t mstep) {
    
    std::size_t curstep = mstep;
    switch (do_sc) {
    case 'C':
        if ((curstep % sc_sep) == 0) {
            GPUSqSum << <blocks_q, threads >> > (emomM, coord, q, r_mid, sc_block_gpu, tasksTot_q, N);
            GPUSqFinalSum_stat << <nq, 1024 >> > (sc_block_gpu, sc_q_gpu, numBlocksX_q);
            cudaDeviceSynchronize();
            n_samples++;
        }
        break;

    case 'Q':
        if ((curstep % sc_step) == 0 && t_cur < static_cast<unsigned int>(sc_max_nstep)) {
            // printf("[GPU-SAMPLE Q] mstep=%zu, curstep=%zu, condition: (curstep %% sc_step)==0, t_cur=%u (max=%u)\n", 
            //        mstep, curstep, t_cur, sc_max_nstep);
            
            // Record metadata BEFORE kernel call (mirrors Fortran: increment then record)
            if (t_cur < static_cast<unsigned int>(dt_cpu.extent(0))) {
                dt_cpu(int(t_cur)) = delta_t * sc_step;  // Record time step at current index

            }
            if (t_cur < static_cast<unsigned int>(sc_step_arr_cpu.extent(0))) {
                sc_step_arr_cpu(int(t_cur)) = static_cast<real>(sc_step);  // Record step width

            }
            
            // Kernel writes to m_kt(:,:,t_cur)
            GPUSqSum << <blocks_q, threads >> > (emomM, coord, q, r_mid, sc_block_gpu, tasksTot_q, N);
            GPUSqFinalSum_dyn << <nq, 1024 >> > (sc_block_gpu, sc_qt_gpu, numBlocksX_q, t_cur);
            cudaDeviceSynchronize();
            t_cur++;  // Increment AFTER writing to that time slice
        } else {
            if ((curstep % sc_step) == 0) {
            }
        }
        break;

    case 'Y':
        if (((curstep % sc_step) == 0) && ((curstep % sc_sep) == 0) && t_cur < static_cast<unsigned int>(sc_max_nstep)) {
            both_flag = 2;

            // Record metadata for time-domain sample
            if (t_cur < static_cast<unsigned int>(dt_cpu.extent(0))) {
                dt_cpu(int(t_cur)) = delta_t * sc_step;
            }
            if (t_cur < static_cast<unsigned int>(sc_step_arr_cpu.extent(0))) {
                sc_step_arr_cpu(int(t_cur)) = static_cast<real>(sc_step);
            }
            
            GPUSqSum << <blocks_q, threads >> > (emomM, coord, q, r_mid, sc_block_gpu, tasksTot_q, N);
            GPUSqFinalSum_both << <nq, 1024 >> > (sc_block_gpu, sc_q_gpu, sc_qt_gpu, numBlocksX_q, t_cur, both_flag);
            cudaDeviceSynchronize();
            t_cur++;
            n_samples++;

        }
        else if ((curstep % sc_step) == 0 && t_cur < static_cast<unsigned int>(sc_max_nstep)) {
            both_flag = 1;

            // Record metadata for time-domain sample
            if (t_cur < static_cast<unsigned int>(dt_cpu.extent(0))) {
                dt_cpu(int(t_cur)) = delta_t * sc_step;
            }
            if (t_cur < static_cast<unsigned int>(sc_step_arr_cpu.extent(0))) {
                sc_step_arr_cpu(int(t_cur)) = static_cast<real>(sc_step);
            }
            
            GPUSqSum << <blocks_q, threads >> > (emomM, coord, q, r_mid, sc_block_gpu, tasksTot_q, N);
            GPUSqFinalSum_both << <nq, 1024 >> > (sc_block_gpu, sc_q_gpu, sc_qt_gpu, numBlocksX_q, t_cur, both_flag);
            cudaDeviceSynchronize();
            t_cur++;
        }
        else if ((curstep % sc_sep) == 0) {
            both_flag = 0;

            GPUSqSum << <blocks_q, threads >> > (emomM, coord, q, r_mid, sc_block_gpu, tasksTot_q, N);
            GPUSqFinalSum_both << <nq, 1024 >> > (sc_block_gpu, sc_q_gpu, sc_qt_gpu, numBlocksX_q, t_cur, both_flag);
            cudaDeviceSynchronize();
            n_samples++;
        }
        break;

    }    

}

void GpuCorrelations::flushCorrelations(hostCorrelations& cpuCorrelations, std::size_t mstep) {
    //TODO cpuCorrelations.m_kt.copy_sync(sc_qt_gpu);
    thrust::complex<real> sc_cur;
    
       /* for (int k = 0; k < nq; k++) {
            for (int i = 0; i < sc_max_nstep; i++) {
                for (int j = 0; j < 3; j++) {
                sc_cur = cpuCorrelations.m_kt(j, i, k);
                printf("iq = %i, it  =%i, is = %i, Re = %.3lf, Im = %.3lf\n", k, i, j, sc_cur.real(), sc_cur.imag());
            }
        }
    }*/
    int tasks; int bl;
    switch (do_sc) {
    case 'C': {
        tasks = 3 * nq;
        bl = (tasks + maxThreads - 1) / maxThreads;
        GPUSqAvrg << <bl, maxThreads >> > (sc_q_gpu, n_samples, tasks, M);  
        cudaDeviceSynchronize();
        
        // Transfer to host CPU tensor first
        sc_q_cpu.copy_sync(sc_q_gpu);
        
        // Transfer to cpuCorrelations for Fortran
        cpuCorrelations.m_k.copy_sync(sc_q_gpu);
        break;
    }

    case 'Q': {
        // Copy time step data to GPU
        dt.copy_sync(dt_cpu);
        
        // Zero intermediate and output buffers
        int bl_w = (3 * numBlocksX_w * nq * nw + numThreads - 1) / numThreads;
        setZero<3> << <bl_w, numThreads >> > (sc_block_w_gpu, 3 * numBlocksX_w * nq * nw);
        bl_w = (3 * nq * nw + numThreads - 1) / numThreads;
        setZero<3> << <bl_w, numThreads >> > (sc_qw_gpu, 3 * nq * nw);
        cudaDeviceSynchronize();
        
        // Compute partial S(q,ω) from S(q,t) using Fourier transform
        GPUSwSum << <blocks_w, threads >> > (sc_qt_gpu, dt, w, sc_block_w_gpu, tasksTot_w, sc_max_nstep, nq, sc_max_nstep, sc_window_fun);
        cudaDeviceSynchronize();
        
        // Reduce block results to final S(q,ω)
        GPUSwFinalSum << <nq * nw, 1024 >> > (sc_block_w_gpu, sc_qw_gpu, numBlocksX_w, nq);
        cudaDeviceSynchronize();
        
        // Transfer time-domain correlations for Fortran reference
        if (sc_qt_gpu.extent(0) == cpuCorrelations.m_kt.extent(0) &&
            sc_qt_gpu.extent(1) == cpuCorrelations.m_kt.extent(1) &&
            sc_qt_gpu.extent(2) == cpuCorrelations.m_kt.extent(2)) {
            cpuCorrelations.m_kt.copy_sync(sc_qt_gpu);
            cpuCorrelations.sc_tidx = static_cast<int>(sc_qt_gpu.extent(2));
        }
        
        // Transfer frequency-domain correlations (m_kw)
        if (sc_qw_gpu.extent(0) == cpuCorrelations.m_kw.extent(0) &&
            sc_qw_gpu.extent(1) == cpuCorrelations.m_kw.extent(1) &&
            sc_qw_gpu.extent(2) == cpuCorrelations.m_kw.extent(2)) {
            cpuCorrelations.m_kw.copy_sync(sc_qw_gpu);
        }
        break;
    }

    case 'Y': {
        // Copy time step data to GPU
        dt.copy_sync(dt_cpu);
        
        // Zero intermediate and output FFT buffers
        int bl_w = (3 * numBlocksX_w * nq * nw + numThreads - 1) / numThreads;
        setZero<3> << <bl_w, numThreads >> > (sc_block_w_gpu, 3 * numBlocksX_w * nq * nw);
        bl_w = (3 * nq * nw + numThreads - 1) / numThreads;
        setZero<3> << <bl_w, numThreads >> > (sc_qw_gpu, 3 * nq * nw);
        cudaDeviceSynchronize();
        
        // Average static S(q) correlations
        tasks = 3 * nq;
        bl = (tasks + maxThreads - 1) / maxThreads;
        GPUSqAvrg << <bl, maxThreads >> > (sc_q_gpu, n_samples, tasks, M);
        cudaDeviceSynchronize();
        
        // Compute partial S(q,ω) from S(q,t) using Fourier transform
        GPUSwSum << <blocks_w, threads >> > (sc_qt_gpu, dt, w, sc_block_w_gpu, tasksTot_w, sc_max_nstep, nq, sc_max_nstep, sc_window_fun);
        cudaDeviceSynchronize();
        
        // Reduce block results to final S(q,ω)
        GPUSwFinalSum << <nq * nw, 1024 >> > (sc_block_w_gpu, sc_qw_gpu, numBlocksX_w, nq);
        cudaDeviceSynchronize();
        
        // Transfer static S(q)
        cpuCorrelations.m_k.copy_sync(sc_q_gpu);
        
        // Transfer time-domain S(q,t)
        if (sc_qt_gpu.extent(0) == cpuCorrelations.m_kt.extent(0) &&
            sc_qt_gpu.extent(1) == cpuCorrelations.m_kt.extent(1) &&
            sc_qt_gpu.extent(2) == cpuCorrelations.m_kt.extent(2)) {
            cpuCorrelations.m_kt.copy_sync(sc_qt_gpu);
            cpuCorrelations.sc_tidx = static_cast<int>(sc_qt_gpu.extent(2));
        }
        
        // Transfer frequency-domain S(q,ω)
        if (sc_qw_gpu.extent(0) == cpuCorrelations.m_kw.extent(0) &&
            sc_qw_gpu.extent(1) == cpuCorrelations.m_kw.extent(1) &&
            sc_qw_gpu.extent(2) == cpuCorrelations.m_kw.extent(2)) {
            cpuCorrelations.m_kw.copy_sync(sc_qw_gpu);
        }
        break;
    }
    
    default: {
        // Case C: S(q) static correlations only
        tasks = 3 * nq;
        bl = (tasks + maxThreads - 1) / maxThreads;
        GPUSqAvrg << <bl, maxThreads >> > (sc_q_gpu, n_samples, tasks, M);
        cudaDeviceSynchronize();
        
        // Transfer to cpuCorrelations for Fortran
        cpuCorrelations.m_k.copy_sync(sc_q_gpu);
        break;
    }
    }  // end switch
    
    cudaDeviceSynchronize();
    
    // Publish sampling info to Fortran (CRITICAL: updates sc_nsamp and sc_tidx)
    publishSamplingInfo(cpuCorrelations);

}

void GpuCorrelations::recordSample() {
    // Record current time and step size when a sample is taken
    if (n_samples < static_cast<std::size_t>(dt_cpu.extent(0))) {
        dt_cpu(int(n_samples)) = delta_t;
    }
    // Track which step this sample corresponds to
    if (t_cur < static_cast<unsigned int>(sc_step_arr_cpu.extent(0))) {
        sc_step_arr_cpu(int(t_cur)) = static_cast<real>(sc_step);
    }
    t_cur++;
    n_samples++;
}

void GpuCorrelations::publishSamplingInfo(hostCorrelations& cpuCorrelations) {
    // Copy the recorded delta_t values back to Fortran's deltat_corr array
    if (FortranData::deltat_corr != nullptr && dt_cpu.extent(0) > 0) {
        std::size_t n_copy = std::min(static_cast<std::size_t>(dt_cpu.extent(0)), n_samples);
        for (std::size_t i = 0; i < n_copy; i++) {
            FortranData::deltat_corr[i] = dt_cpu(int(i));
        }
    }
    
    // Copy the recorded sc_step values back to Fortran's scstep_arr array
    if (FortranData::scstep_arr != nullptr && sc_step_arr_cpu.extent(0) > 0) {
        std::size_t n_copy = std::min(static_cast<std::size_t>(sc_step_arr_cpu.extent(0)), n_samples);
        for (std::size_t i = 0; i < n_copy; i++) {
            FortranData::scstep_arr[i] = sc_step_arr_cpu(int(i));
        }
    }
    
    // Update cpuCorrelations with the sample count and time index
    // GPU does NOT normalize - Fortran will handle normalization using sc_nsamp
    cpuCorrelations.sc_nsamp = static_cast<int>(n_samples);
    cpuCorrelations.sc_tidx = static_cast<int>(t_cur);
    
    // Write sample count and time index back to Fortran through FortranData pointers
    if (FortranData::sc_nsamp_ptr != nullptr) {
        *FortranData::sc_nsamp_ptr = static_cast<int>(n_samples);
    }
    if (FortranData::sc_tidx_ptr != nullptr) {
        *FortranData::sc_tidx_ptr = static_cast<int>(t_cur);
    }
}






