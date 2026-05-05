#pragma once

#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include <numeric>
#include "gpu_wrappers.h"
#include <thrust/complex.h>
#include <correlation_kernels.cuh>

namespace cg = cooperative_groups;
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif


__device__ real sc_window_fac(int sc_window_fun, unsigned int step, unsigned int nstep) {
    real dum = 1.0;
    switch (sc_window_fun) {
        //Hann
    case 2:
        dum = (0.5 - 0.5 * cos(2.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)));
        break;
        //Hamming
    case 3:
        dum = (0.54 - 0.46 * cos(2.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)));
        break;
            //Hamming v2
    case 32:
        dum = (0.53836 - 0.46164 * cos(2.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)));
        break;
        //Blackman - Harris
    case 4:
        dum =
            (0.35785 - 0.48829 * cos(2.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)) +
                0.14128 * cos(4.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)) -
                0.01168 * cos(6.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)));
        break;
        //Nuttal
    case 5:
        dum =
            (0.355768 - 0.478396 * cos(2.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)) +
                0.144232 * cos(4.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)) -
                0.012604 * cos(6.0 * M_PI * ((real)step - 1.0) / ((real)nstep - 1.0)));
        break;
    }
 
        
    return dum;
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

    int qInd = grid.block_index().y;
    int wInd = grid.block_index().z;

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
        real win = sc_window_fac(sc_window_fun, (tInd + 1), sc_max_nstep);

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

__global__ void GPUSqProjSum(const GpuTensor<real, 3> spin, const GpuTensor<real, 2> coord, const GpuTensor<real, 2> q, const GpuTensor<real, 1> r_mid, 
                             GpuVector<int>aproj, GpuTensor<thrust::complex<real>, 3> scblock, int tasks, unsigned int N) {
    auto grid = cg::this_grid();
    auto block = cg::this_thread_block();
    auto warp = cg::tiled_partition<32>(block);

    int lane = warp.thread_rank();
    int wid = warp.meta_group_rank();
    int wSize = warp.size();
    int wNum = warp.meta_group_size();
    int tid = grid.thread_rank();
    int tid_in_block = block.thread_rank();

    int qInd = grid.block_index().z;
    int pInd = grid.block_index().y;
    int tid_in_X = grid.block_index().x * block.num_threads() + tid_in_block;
    int stride = grid.dim_blocks().x * block.num_threads();
    int nproj = aproj.extent(0);

    // Register-based accumulation: store real and imaginary parts separately (reduced register pressure)
    real sum_re[3] = {0.0, 0.0, 0.0};
    real sum_im[3] = {0.0, 0.0, 0.0};
    
    unsigned int rInd, mInd, cInd, ii;
    real inv_N = 1.0 / N;
    __shared__ real shared_re[3][32];
    __shared__ real shared_im[3][32];

    int it = 0;

    // Main computation loop over tasks (unroll for better performance)
    #pragma unroll 4
    for (int id = tid_in_X; id < tasks; id += stride) {
        ii = id / 3;
        cInd = id % 3;
        rInd = ii % N;
        mInd = ii / N; 
        it = aproj(rInd) - 1;
        if(it == pInd){
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
        scblock(3 * block.group_index().x + 0, pInd, qInd) = thrust::complex<real>(sum_re[0], sum_im[0]);
        scblock(3 * block.group_index().x + 1, pInd, qInd) = thrust::complex<real>(sum_re[1], sum_im[1]);
        scblock(3 * block.group_index().x + 2, pInd, qInd) = thrust::complex<real>(sum_re[2], sum_im[2]);
    }
}

__global__ void GPUSqProjFinalSum_stat(GpuTensor<thrust::complex<real>, 3> scblock, GpuTensor<thrust::complex<real>, 3> scsum, int numBlocks)
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

    int qInd = grid.block_index().y;
    int pInd = grid.block_index().x;
    int tid_in_Q = tid_in_block;

    // Register-based accumulators
    real sum_re[3] = {0.0, 0.0, 0.0};
    real sum_im[3] = {0.0, 0.0, 0.0};
    
    __shared__ real shared_re[3][32];
    __shared__ real shared_im[3][32];

    if (tid_in_Q < numBlocks) {
        for (int k = 0; k < 3; k++) {
            thrust::complex<real> val = scblock(3 * tid_in_Q + k, pInd, qInd);
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
        scsum(0, pInd, qInd) += thrust::complex<real>(sum_re[0], sum_im[0]);
        scsum(1, pInd, qInd) += thrust::complex<real>(sum_re[1], sum_im[1]);
        scsum(2, pInd, qInd) += thrust::complex<real>(sum_re[2], sum_im[2]);

        /*mblock_gpu[block.group_index().x] += mySum[0];
        mblock_gpu[block.group_index().x + grid.group_dim().x] += mySum[1];
        mblock_gpu[block.group_index().x + 2 * grid.group_dim().x] += mySum[2];*/
       // printf("qInd = %i, mblock0 = %lf, mblock1 = %lf, mblock2 = %lf\n", mInd, msum(0, curstep, mInd), msum(1, curstep, mInd), msum(2, curstep, mInd));
    }
}

__global__ void GPUSqProjFinalSum_dyn(GpuTensor<thrust::complex<real>, 3> scblock, GpuTensor<thrust::complex<real>, 4> scsum, int numBlocks, unsigned int t_cur)
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

    int qInd = grid.block_index().y;
    int pInd = grid.block_index().x;
    int tid_in_Q = tid_in_block;

    // Register-based accumulators
    real sum_re[3] = {0.0, 0.0, 0.0};
    real sum_im[3] = {0.0, 0.0, 0.0};
    
    __shared__ real shared_re[3][32];
    __shared__ real shared_im[3][32];

    if (tid_in_Q < numBlocks) {
        for (int k = 0; k < 3; k++) {
            thrust::complex<real> val = scblock(3 * tid_in_Q + k, pInd, qInd);
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
            scsum(0, pInd, qInd, t_cur) += thrust::complex<real>(sum_re[0], sum_im[0]);
            scsum(1, pInd, qInd, t_cur) += thrust::complex<real>(sum_re[1], sum_im[1]);
            scsum(2, pInd, qInd, t_cur) += thrust::complex<real>(sum_re[2], sum_im[2]);
        
        //printf("qInd = %i, t_cur = %i, Re = %.3lf, Im = %.3lf\n", qInd, t_cur, mySum[2].real(), mySum[2].imag());

        /*mblock_gpu[block.group_index().x] += mySum[0];
        mblock_gpu[block.group_index().x + grid.group_dim().x] += mySum[1];
        mblock_gpu[block.group_index().x + 2 * grid.group_dim().x] += mySum[2];*/
        // printf("qInd = %i, mblock0 = %lf, mblock1 = %lf, mblock2 = %lf\n", mInd, msum(0, curstep, mInd), msum(1, curstep, mInd), msum(2, curstep, mInd));
    }
}

__global__ void GPUSqProjFinalSum_both(GpuTensor<thrust::complex<real>, 3> scblock, GpuTensor<thrust::complex<real>, 3> scsum_q, GpuTensor<thrust::complex<real>, 4> scsum_qt, int numBlocks, unsigned int t_cur, unsigned int both_flag)
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

    int qInd = grid.block_index().y;
    int pInd = grid.block_index().x;
    int tid_in_Q = tid_in_block;

    // Register-based accumulators
    real sum_re[3] = {0.0, 0.0, 0.0};
    real sum_im[3] = {0.0, 0.0, 0.0};
    
    __shared__ real shared_re[3][32];
    __shared__ real shared_im[3][32];

    if (tid_in_Q < numBlocks) {
        for (int k = 0; k < 3; k++) {
            thrust::complex<real> val = scblock(3 * tid_in_Q + k, pInd, qInd);
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
            scsum_qt(0, pInd, qInd, t_cur) += thrust::complex<real>(sum_re[0], sum_im[0]);
            scsum_qt(1, pInd, qInd, t_cur) += thrust::complex<real>(sum_re[1], sum_im[1]);
            scsum_qt(2, pInd, qInd, t_cur) += thrust::complex<real>(sum_re[2], sum_im[2]);
        }

        if ((both_flag == 0) || (both_flag == 2)) {
            scsum_q(0, pInd, qInd) += thrust::complex<real>(sum_re[0], sum_im[0]);
            scsum_q(1, pInd, qInd) += thrust::complex<real>(sum_re[1], sum_im[1]);
            scsum_q(2, pInd, qInd) += thrust::complex<real>(sum_re[2], sum_im[2]);
        }

        //printf("Re = %.3lf, Im = %.3lf\n", mySum[1].real(), mySum[1].imag());

        /*mblock_gpu[block.group_index().x] += mySum[0];
        mblock_gpu[block.group_index().x + grid.group_dim().x] += mySum[1];
        mblock_gpu[block.group_index().x + 2 * grid.group_dim().x] += mySum[2];*/
        // printf("qInd = %i, mblock0 = %lf, mblock1 = %lf, mblock2 = %lf\n", mInd, msum(0, curstep, mInd), msum(1, curstep, mInd), msum(2, curstep, mInd));
    }
}

__global__ void GPUSwProjSum(const GpuTensor<thrust::complex<real>, 4> sq, const GpuTensor<real, 1> dt, const GpuTensor<real, 1> w, GpuTensor<thrust::complex<real>, 4> scblock,
                             unsigned int blockN, int tasks, unsigned int tSize, unsigned int nq, int sc_max_nstep, int sc_window_fun) {
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
    int wInd = grid.block_index().z;

    int pInd = grid.block_index().x / blockN;
    int bInd = grid.block_index().x % blockN;
    
    int tid_in_X = bInd * blockDim.x + threadIdx.x;
    int stride   = blockN * blockDim.x;
    //int tid_in_X = grid.block_index().x * block.num_threads() + tid_in_block;
    //int stride = grid.dim_blocks().x * block.num_threads();

    // Register-based accumulation: store real and imaginary parts separately
    real sum_re[3] = {0.0, 0.0, 0.0};
    real sum_im[3] = {0.0, 0.0, 0.0};
    
    unsigned int tInd, cInd, ii;
    __shared__ real shared_re[3][32];
    __shared__ real shared_im[3][32];


    thrust::complex<real> sq_val;
    // Fourier transform loop: unroll for performance
    #pragma unroll 2
    for (int id = tid_in_X; id < tasks; id += stride) {
        cInd = id % 3;
        ii = id / 3;
        tInd = ii % sc_max_nstep;

        //printf("c = %i, p = %i, q = %i, t = %i\n", cInd, pInd, qInd, tInd);
       // printf("id = %i, tasks = %i\n",id, tasks);


        // 1. Calculate the real-valued phase: Phase = t * dt(t) * w(w)
        real phase = (real)tInd * dt(tInd) * w(wInd);

        // 2. Hardware-accelerated sine and cosine
        real s, c;
        sincos(phase, &s, &c);

        // 3. Extract S(q,t) values
        sq_val = sq(cInd, pInd, qInd, tInd);
       //sq_val = thrust::complex<real>(0,0);
        real sq_re = sq_val.real();
        real sq_im = sq_val.imag();
        
        // 4. Windowing function
        real win = sc_window_fac(sc_window_fun, (tInd + 1), sc_max_nstep);

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
       scblock(3 * bInd + 0, pInd, qInd, wInd) = thrust::complex<real>(sum_re[0], sum_im[0]);
       scblock(3 * bInd + 1, pInd, qInd, wInd) = thrust::complex<real>(sum_re[1], sum_im[1]);
       scblock(3 * bInd + 2, pInd, qInd, wInd) = thrust::complex<real>(sum_re[2], sum_im[2]);
       //printf("pInd = %i, scblock = %.4lf\n", pInd, scblock(3 * bInd + 0, pInd, qInd, wInd).real());
    }
}

__global__ void GPUSwProjFinalSum(GpuTensor<thrust::complex<real>, 4> scblock, GpuTensor<thrust::complex<real>, 4> scsum, int numBlocks, int nq)
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

    int qInd = (grid.block_index().y);
    int wInd = (grid.block_index().z);
    int pInd = (grid.block_index().x);

    int tid_in_Q = tid_in_block;

    // Register-based accumulators for intermediate block results
    real sum_re[3] = {0.0, 0.0, 0.0};
    real sum_im[3] = {0.0, 0.0, 0.0};
    
    __shared__ real shared_re[3][32];
    __shared__ real shared_im[3][32];

    if (tid_in_Q < numBlocks) {
        // Load and sum block results
        for (int k = 0; k < 3; k++) {
            thrust::complex<real> val = scblock(3 * tid_in_Q + k, pInd, qInd, wInd);
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
        scsum(0, pInd, qInd, wInd) += thrust::complex<real>(sum_re[0], sum_im[0]);
        scsum(1, pInd, qInd, wInd) += thrust::complex<real>(sum_re[1], sum_im[1]);
        scsum(2, pInd, qInd, wInd) += thrust::complex<real>(sum_re[2], sum_im[2]);
       // printf("p = %i, scsum = %.4lf, %.4lf\n", pInd, scsum(0, pInd, qInd, wInd).real(), scsum(0, pInd, qInd, wInd).imag());
    }
}
