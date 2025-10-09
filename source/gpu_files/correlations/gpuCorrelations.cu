#pragma once

#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "gpuStructures.hpp"
#include "gpuCorrelations.cuh"
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

__inline__ __device__
thrust::complex<real> warpReduceSum(thrust::complex<real> val) {
    real valr = val.real();
    real vali = val.imag();
#if CUDA_VERSION < 9000
    for (int offset = warpSize / 2; offset > 0; offset /= 2) {

        valr += __shfl_down(valr, offset);
        vali += __shfl_down(vali, offset);
    }
#else
    for (int offset = warpSize / 2; offset > 0; offset /= 2) {

        valr += __shfl_down_sync(0xffffffffffff, valr, offset);
        vali += __shfl_down_sync(0xffffffffffff, vali, offset);
    }
#endif

    val = thrust::complex<real>(valr, vali);
    return val;
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
    //int offsetM = mInd * tasks;
    int tid_in_X = grid.block_index().x * block.num_threads() + tid_in_block;
    int stride = grid.dim_blocks().x * block.num_threads();

    thrust::complex<real>iqfac = thrust::complex<real>(0, 2 * M_PI);
    unsigned int rInd, mInd, cInd, ii;
    real qdr;
    thrust::complex<real> mySum[3] = {thrust::complex < real>(0.0, 0.0), thrust::complex < real>(0.0, 0.0), thrust::complex < real>(0.0, 0.0) };
    static __shared__ thrust::complex<real> shared0[32];
    static __shared__ thrust::complex<real> shared1[32];
    static __shared__ thrust::complex<real> shared2[32];

   // qdr = q(1, iq) * (coord(1, r) - r_mid(1)) + q(2, iq) * (coord(2, r) - r_mid(2)) + q(3, iq) * (coord(3, r) - r_mid(3))
   // epowqr = exp(iqfac * qdr) * nainv!*sc_window_fac(sc_window_fun, iq, nq)
   // wA = wA + epowqr * SA(:, r, l)
    for (int id = tid_in_X; id < tasks; id += stride) {
        ii = id / 3;
        cInd = id % 3;
        rInd = ii % N;
        mInd = ii / N; 
        qdr = q(0, qInd) * (coord(0, rInd) - r_mid(0)) + q(1, qInd) * (coord(1, rInd) - r_mid(1)) + q(2, qInd) * (coord(2, rInd) - r_mid(2));



        mySum[cInd] += thrust::complex<real>(spin(cInd, rInd, mInd), 0) * thrust::exp(iqfac * thrust::complex<real>(qdr, 0))/N;
        //printf("qInd = %i, Re = %.3lf, Im = %.3lf,qdr = %.3lf\n", qInd, thrust::complex<real>(qdr, 0).real(), thrust::complex<real>(qdr, 0).imag(), qdr);
        //printf("0:rmid = %i, q = %.3lf, Im = %.3lf,qdr = %.3lf\n", qInd, thrust::complex<real>(qdr, 0).real(), thrust::complex<real>(qdr, 0).imag(), qdr);

        //  printf("tid = %i, mInd = %i, stride = %i, data_id = %i, mySum = %.3f\n", tid, mInd , stride, id + offsetM, mySum[id % 3]);
    }
    warp.sync();

    mySum[0] = warpReduceSum(mySum[0]);
    mySum[1] = warpReduceSum(mySum[1]);
    mySum[2] = warpReduceSum(mySum[2]);

    if (lane == 0) {
        shared0[wid] = mySum[0];
        shared1[wid] = mySum[1];
        shared2[wid] = mySum[2];
    }

    __syncthreads();              // Wait for all partial reductions
    mySum[0] = (tid_in_block < wNum) ? shared0[lane] : 0;
    mySum[1] = (tid_in_block < wNum) ? shared1[lane] : 0;
    mySum[2] = (tid_in_block < wNum) ? shared2[lane] : 0;

    if (wid == 0) {
        mySum[0] = warpReduceSum(mySum[0]); //Final reduce within first warp
        mySum[1] = warpReduceSum(mySum[1]);
        mySum[2] = warpReduceSum(mySum[2]);
    }
   // printf("Re = %.3lf, Im = %.3lf\n", mySum[1].real(), mySum[1].imag());

    if (tid_in_block == 0) {
        scblock(block.group_index().x, qInd) = mySum[0];
        scblock(block.group_index().x + grid.group_dim().x, qInd) = mySum[1];
        scblock(block.group_index().x + 2 * grid.group_dim().x, qInd) = mySum[2];
       // printf("qInd = %i, sInd = %i, Re = %.3lf, Im = %.3lf, q = %.3lf\n", qInd, block.group_index().x,  mySum[1].real(), mySum[1].imag(), q(1, qInd));
        // printf("mInd = %i, bid = %i, mblock0 = %.3f\n", mInd, block.group_index().x, mblock(block.group_index().x, mInd));
         //printf("tid = %i, mInd = %i, mblock0 = %lf, mblock1 = %lf, mblock2 = %lf\n", tid, mInd, mblock(block.group_index().x, mInd), mblock(block.group_index().x + grid.group_dim().x, mInd), mblock(block.group_index().x + 2 * grid.group_dim().x, mInd));
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
    //int offsetM = mInd * numBlocks;
    int tid_in_Q = tid_in_block;
    //printf("numblocks = %i\n", numBlocks);

    thrust::complex<real> mySum[3] = {thrust::complex < real>(0.0, 0.0), thrust::complex < real>(0.0, 0.0), thrust::complex < real>(0.0, 0.0) };
    static __shared__ thrust::complex<real> shared0[32];
    static __shared__ thrust::complex<real> shared1[32];
    static __shared__ thrust::complex<real> shared2[32];

    if (tid_in_Q < numBlocks) {

        mySum[0] += scblock(tid_in_Q, qInd);
        mySum[1] += scblock(tid_in_Q + numBlocks, qInd);
        mySum[2] += scblock(tid_in_Q + 2 * numBlocks, qInd);
        //printf("tid_in_m = %i, mInd = %i, mblock = %.3f\n", tid_in_M, mInd, mblock(tid_in_M, mInd));
    }

    warp.sync();

    mySum[0] = warpReduceSum(mySum[0]); //Final reduce within first warp
    mySum[1] = warpReduceSum(mySum[1]); //Final reduce within first warp
    mySum[2] = warpReduceSum(mySum[2]); //Final reduce within first warp

    if (lane == 0) {
        shared0[wid] = mySum[0];
        shared1[wid] = mySum[1];
        shared2[wid] = mySum[2];
    }

    __syncthreads();              // Wait for all partial reductions
    mySum[0] = (tid_in_block < wNum) ? shared0[lane] : 0;
    mySum[1] = (tid_in_block < wNum) ? shared1[lane] : 0;
    mySum[2] = (tid_in_block < wNum) ? shared2[lane] : 0;
    if (wid == 0) mySum[0] = warpReduceSum(mySum[0]); //Final reduce within first warp
    if (wid == 0) mySum[1] = warpReduceSum(mySum[1]); //Final reduce within first warp
    if (wid == 0) mySum[2] = warpReduceSum(mySum[2]); //Final reduce within first warp

    if (tid_in_block == 0) {
        scsum(0, qInd) += mySum[0];
        scsum(1, qInd) += mySum[1];
        scsum(2, qInd) += mySum[2];
        //printf("Re = %.3lf, Im = %.3lf\n", mySum[1].real(), mySum[1].imag());

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
    //int offsetM = mInd * numBlocks;
    int tid_in_Q = tid_in_block;
    //printf("numblocks = %i\n", numBlocks);

    thrust::complex<real> mySum[3] = { thrust::complex < real>(0.0, 0.0), thrust::complex < real>(0.0, 0.0), thrust::complex < real>(0.0, 0.0) };
    static __shared__ thrust::complex<real> shared0[32];
    static __shared__ thrust::complex<real> shared1[32];
    static __shared__ thrust::complex<real> shared2[32];

    if (tid_in_Q < numBlocks) {

        mySum[0] += scblock(tid_in_Q, qInd);
        mySum[1] += scblock(tid_in_Q + numBlocks, qInd);
        mySum[2] += scblock(tid_in_Q + 2 * numBlocks, qInd);
        //printf("tid_in_m = %i, mInd = %i, mblock = %.3f\n", tid_in_M, mInd, mblock(tid_in_M, mInd));
    }

    warp.sync();

    mySum[0] = warpReduceSum(mySum[0]); //Final reduce within first warp
    mySum[1] = warpReduceSum(mySum[1]); //Final reduce within first warp
    mySum[2] = warpReduceSum(mySum[2]); //Final reduce within first warp

    if (lane == 0) {
        shared0[wid] = mySum[0];
        shared1[wid] = mySum[1];
        shared2[wid] = mySum[2];
    }

    __syncthreads();              // Wait for all partial reductions
    mySum[0] = (tid_in_block < wNum) ? shared0[lane] : 0;
    mySum[1] = (tid_in_block < wNum) ? shared1[lane] : 0;
    mySum[2] = (tid_in_block < wNum) ? shared2[lane] : 0;
    if (wid == 0) mySum[0] = warpReduceSum(mySum[0]); //Final reduce within first warp
    if (wid == 0) mySum[1] = warpReduceSum(mySum[1]); //Final reduce within first warp
    if (wid == 0) mySum[2] = warpReduceSum(mySum[2]); //Final reduce within first warp

    if (tid_in_block == 0) {
            scsum(0, t_cur, qInd) += mySum[0];
            scsum(1, t_cur, qInd) += mySum[1];
            scsum(2, t_cur, qInd) += mySum[2];
        
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
    //int offsetM = mInd * numBlocks;
    int tid_in_Q = tid_in_block;
    //printf("numblocks = %i\n", numBlocks);

    thrust::complex<real> mySum[3] = { thrust::complex < real>(0.0, 0.0), thrust::complex < real>(0.0, 0.0), thrust::complex < real>(0.0, 0.0) };
    static __shared__ thrust::complex<real> shared0[32];
    static __shared__ thrust::complex<real> shared1[32];
    static __shared__ thrust::complex<real> shared2[32];

    if (tid_in_Q < numBlocks) {

        mySum[0] += scblock(tid_in_Q, qInd);
        mySum[1] += scblock(tid_in_Q + numBlocks, qInd);
        mySum[2] += scblock(tid_in_Q + 2 * numBlocks, qInd);
        //printf("tid_in_m = %i, mInd = %i, mblock = %.3f\n", tid_in_M, mInd, mblock(tid_in_M, mInd));
    }

    warp.sync();

    mySum[0] = warpReduceSum(mySum[0]); //Final reduce within first warp
    mySum[1] = warpReduceSum(mySum[1]); //Final reduce within first warp
    mySum[2] = warpReduceSum(mySum[2]); //Final reduce within first warp

    if (lane == 0) {
        shared0[wid] = mySum[0];
        shared1[wid] = mySum[1];
        shared2[wid] = mySum[2];
    }

    __syncthreads();              // Wait for all partial reductions
    mySum[0] = (tid_in_block < wNum) ? shared0[lane] : 0;
    mySum[1] = (tid_in_block < wNum) ? shared1[lane] : 0;
    mySum[2] = (tid_in_block < wNum) ? shared2[lane] : 0;
    if (wid == 0) mySum[0] = warpReduceSum(mySum[0]); //Final reduce within first warp
    if (wid == 0) mySum[1] = warpReduceSum(mySum[1]); //Final reduce within first warp
    if (wid == 0) mySum[2] = warpReduceSum(mySum[2]); //Final reduce within first warp

    if (tid_in_block == 0) {
        if ((both_flag == 1) || (both_flag == 2)) {
            scsum_qt(0, t_cur, qInd) += mySum[0];
            scsum_qt(1, t_cur, qInd) += mySum[1];
            scsum_qt(2, t_cur, qInd) += mySum[2];
        }

        if ((both_flag == 0) || (both_flag == 2)) {
            scsum_q(0, qInd) += mySum[0];
            scsum_q(1, qInd) += mySum[1];
            scsum_q(2, qInd) += mySum[2];
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

    int qInd = grid.block_index().y%nq;//TODO
    int wInd = grid.block_index().y/nq;//TODO

    int tid_in_X = grid.block_index().x * block.num_threads() + tid_in_block;
    int stride = grid.dim_blocks().x * block.num_threads();

   // thrust::complex<real>iqfac = thrust::complex<real>(0, 2 * M_PI);
    unsigned int tInd, cInd, ii;
    thrust::complex<real> tw;
    thrust::complex<real> mySum[3] = { thrust::complex < real>(0.0, 0.0), thrust::complex < real>(0.0, 0.0), thrust::complex < real>(0.0, 0.0) };
    static __shared__ thrust::complex<real> shared0[32];
    static __shared__ thrust::complex<real> shared1[32];
    static __shared__ thrust::complex<real> shared2[32];

    // qdr = q(1, iq) * (coord(1, r) - r_mid(1)) + q(2, iq) * (coord(2, r) - r_mid(2)) + q(3, iq) * (coord(3, r) - r_mid(3))
    // epowqr = exp(iqfac * qdr) * nainv!*sc_window_fac(sc_window_fun, iq, nq)
    // wA = wA + epowqr * SA(:, r, l)
    for (int id = tid_in_X; id < tasks; id += stride) {
        tInd = id / 3;
        cInd = id % 3;
        //qdr = q(0, qInd) * (coord(0, rInd) - r_mid(0)) + q(1, qInd) * (coord(1, rInd) - r_mid(1)) + q(2, qInd) * (coord(2, rInd) - r_mid(2));
        //wfac = 1.0_dblprec
        //corr_kw = 0.0_dblprec
        //i = (0.0_dblprec, 1.0_dblprec)
         //tidx = cc % sc_max_nstep
        //tt = i * wfac * (step - 1) * dt(step)            
        //epowwt = exp(cc % w(iw) * tt) * sc_window_fac(sc_window_fun, step, tidx)
        tw = thrust::complex<real>(0, 1) * tInd * dt(tInd)*w(wInd);
        mySum[cInd] += thrust::exp(tw) * sc_window_fac(sc_window_fun, (tInd-1), sc_max_nstep)*sq(cInd, tInd, qInd);
        //mySum[cInd] += thrust::complex<real>(spin(cInd, rInd, mInd), 0) * thrust::exp(iqfac * thrust::complex<real>(qdr, 0)) / N;
        //printf("qInd = %i, Re = %.3lf, Im = %.3lf,qdr = %.3lf\n", qInd, thrust::complex<real>(qdr, 0).real(), thrust::complex<real>(qdr, 0).imag(), qdr);
        //printf("0:rmid = %i, q = %.3lf, Im = %.3lf,qdr = %.3lf\n", qInd, thrust::complex<real>(qdr, 0).real(), thrust::complex<real>(qdr, 0).imag(), qdr);

        //  printf("tid = %i, mInd = %i, stride = %i, data_id = %i, mySum = %.3f\n", tid, mInd , stride, id + offsetM, mySum[id % 3]);
    }
    warp.sync();

    mySum[0] = warpReduceSum(mySum[0]);
    mySum[1] = warpReduceSum(mySum[1]);
    mySum[2] = warpReduceSum(mySum[2]);

    if (lane == 0) {
        shared0[wid] = mySum[0];
        shared1[wid] = mySum[1];
        shared2[wid] = mySum[2];
    }

    __syncthreads();              // Wait for all partial reductions
    mySum[0] = (tid_in_block < wNum) ? shared0[lane] : 0;
    mySum[1] = (tid_in_block < wNum) ? shared1[lane] : 0;
    mySum[2] = (tid_in_block < wNum) ? shared2[lane] : 0;

    if (wid == 0) {
        mySum[0] = warpReduceSum(mySum[0]); //Final reduce within first warp
        mySum[1] = warpReduceSum(mySum[1]);
        mySum[2] = warpReduceSum(mySum[2]);
    }
    // printf("Re = %.3lf, Im = %.3lf\n", mySum[1].real(), mySum[1].imag());

    if (tid_in_block == 0) {
        scblock(block.group_index().x, qInd, wInd) = mySum[0];
        scblock(block.group_index().x + grid.group_dim().x, qInd, wInd) = mySum[1];
        scblock(block.group_index().x + 2 * grid.group_dim().x, qInd, wInd) = mySum[2];
    //    printf("wInd = %i, qInd = %i, Re = %.3lf, Im = %.3lf\n", wInd, qInd, mySum[2].real(), mySum[2].imag());

        // printf("qInd = %i, sInd = %i, Re = %.3lf, Im = %.3lf, q = %.3lf\n", qInd, block.group_index().x,  mySum[1].real(), mySum[1].imag(), q(1, qInd));
         // printf("mInd = %i, bid = %i, mblock0 = %.3f\n", mInd, block.group_index().x, mblock(block.group_index().x, mInd));
          //printf("tid = %i, mInd = %i, mblock0 = %lf, mblock1 = %lf, mblock2 = %lf\n", tid, mInd, mblock(block.group_index().x, mInd), mblock(block.group_index().x + grid.group_dim().x, mInd), mblock(block.group_index().x + 2 * grid.group_dim().x, mInd));
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

    int qInd = (grid.block_index().x) % nq;//TODO
    int wInd = (grid.block_index().x) / nq;//TODO
    //int offsetM = mInd * numBlocks;
    int tid_in_Q = tid_in_block;
    //printf("numblocks = %i\n", numBlocks);
     //   printf("wInd = %i, qInd = %i\n", wInd, qInd);

    thrust::complex<real> mySum[3] = { thrust::complex < real>(0.0, 0.0), thrust::complex < real>(0.0, 0.0), thrust::complex < real>(0.0, 0.0) };
    static __shared__ thrust::complex<real> shared0[32];
    static __shared__ thrust::complex<real> shared1[32];
    static __shared__ thrust::complex<real> shared2[32];

    if (tid_in_Q < numBlocks) {

        mySum[0] += scblock(tid_in_Q, qInd, wInd);
        mySum[1] += scblock(tid_in_Q + numBlocks, qInd, wInd);
        mySum[2] += scblock(tid_in_Q + 2 * numBlocks, qInd, wInd);
        //printf("tid_in_m = %i, mInd = %i, mblock = %.3f\n", tid_in_M, mInd, mblock(tid_in_M, mInd));
    }

    warp.sync();

    mySum[0] = warpReduceSum(mySum[0]); //Final reduce within first warp
    mySum[1] = warpReduceSum(mySum[1]); //Final reduce within first warp
    mySum[2] = warpReduceSum(mySum[2]); //Final reduce within first warp

    if (lane == 0) {
        shared0[wid] = mySum[0];
        shared1[wid] = mySum[1];
        shared2[wid] = mySum[2];
    }

    __syncthreads();              // Wait for all partial reductions
    mySum[0] = (tid_in_block < wNum) ? shared0[lane] : 0;
    mySum[1] = (tid_in_block < wNum) ? shared1[lane] : 0;
    mySum[2] = (tid_in_block < wNum) ? shared2[lane] : 0;
    if (wid == 0) mySum[0] = warpReduceSum(mySum[0]); //Final reduce within first warp
    if (wid == 0) mySum[1] = warpReduceSum(mySum[1]); //Final reduce within first warp
    if (wid == 0) mySum[2] = warpReduceSum(mySum[2]); //Final reduce within first warp

    if (tid_in_block == 0) {
        scsum(0, qInd, wInd) += mySum[0];
        scsum(1, qInd, wInd) += mySum[1];
        scsum(2, qInd, wInd) += mySum[2];
        printf("wInd = %i, qInd = %i, Re = %.3lf, Im = %.3lf\n", wInd, qInd, mySum[2].real(), mySum[2].imag());


        //printf("Re = %.3lf, Im = %.3lf\n", mySum[1].real(), mySum[1].imag());

        /*mblock_gpu[block.group_index().x] += mySum[0];
        mblock_gpu[block.group_index().x + grid.group_dim().x] += mySum[1];
        mblock_gpu[block.group_index().x + 2 * grid.group_dim().x] += mySum[2];*/
        // printf("qInd = %i, mblock0 = %lf, mblock1 = %lf, mblock2 = %lf\n", mInd, msum(0, curstep, mInd), msum(1, curstep, mInd), msum(2, curstep, mInd));
    }
}
__global__ void GPUSqAvrg(GpuTensor<thrust::complex<real>, 2> sc, int n_steps, int tasks, int M) {
    auto grid = cg::this_grid();
    int tid = grid.thread_rank();
    if (tid < tasks) {
        sc[tid] = sc[tid] / (n_steps*M);
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
            bl = (3 * nq + numThreads - 1) / numThreads;
            setZero<2> << <bl, numThreads >> > (sc_q_gpu, 3 * nq);

        }
        if ((do_sc == 'Q') || (do_sc == 'Y')) {
            sc_qt_gpu.Allocate(static_cast <long int>(3), static_cast <long int>(sc_max_nstep), static_cast <long int>(nq));
            //sc_qt_gpu.Allocate(static_cast <long int>(3), static_cast <long int>(sc_max_nstep), static_cast <long int>(nq));
            sc_qw_gpu.Allocate(static_cast <long int>(3), static_cast <long int>(nq), static_cast <long int>(nw));
            sc_block_w_gpu.Allocate(static_cast <long int>(3 * numBlocksX_w), static_cast <long int>(nq), static_cast <long int>(nw));
            dt.Allocate(static_cast <long int>(sc_max_nstep));
            dt_cpu.AllocateHost(static_cast <long int>(sc_max_nstep));
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
        }
        if ((do_sc == 'Q') || (do_sc == 'Y')) {
            sc_qt_gpu.Free();
            sc_qw_gpu.Free();
            sc_block_w_gpu.Free();
            w.Free();
            dt.Free();
            dt_cpu.FreeHost();
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
        if ((curstep % sc_step) == 0) {
            GPUSqSum << <blocks_q, threads >> > (emomM, coord, q, r_mid, sc_block_gpu, tasksTot_q, N);
            GPUSqFinalSum_dyn << <nq, 1024 >> > (sc_block_gpu, sc_qt_gpu, numBlocksX_q, t_cur);
            cudaDeviceSynchronize();
            dt_cpu[t_cur] = delta_t * sc_step;
            t_cur++;

        }   
        break;

    case 'Y':
        if (((curstep % sc_step) == 0) && ((curstep % sc_sep) == 0)) {
            both_flag = 2;
            GPUSqSum << <blocks_q, threads >> > (emomM, coord, q, r_mid, sc_block_gpu, tasksTot_q, N);
            GPUSqFinalSum_both << <nq, 1024 >> > (sc_block_gpu, sc_q_gpu, sc_qt_gpu, numBlocksX_q, t_cur, both_flag);
            cudaDeviceSynchronize();
            dt_cpu[t_cur] = delta_t * sc_step;
            t_cur++;
            n_samples++;

        }
        else if ((curstep % sc_step) == 0) {
            both_flag = 1;
            GPUSqSum << <blocks_q, threads >> > (emomM, coord, q, r_mid, sc_block_gpu, tasksTot_q, N);
            GPUSqFinalSum_both << <nq, 1024 >> > (sc_block_gpu, sc_q_gpu, sc_qt_gpu, numBlocksX_q, t_cur, both_flag);
            cudaDeviceSynchronize();
            dt_cpu[t_cur] = delta_t * sc_step;
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
    case 'C':
        tasks = 3 * nq;
        bl = (tasks + maxThreads - 1) / maxThreads;
        GPUSqAvrg << <bl, maxThreads >> > (sc_q_gpu, n_samples, tasks, M);  
        //TODO cpuCorrelations.m_k.copy_sync(sc_q_gpu);
        break;

    case 'Q':
        dt.copy_sync(dt_cpu);
        GPUSwSum << <blocks_w, maxThreads >> > (sc_qt_gpu, dt, w, sc_block_w_gpu, tasksTot_w, sc_max_nstep, nq, sc_max_nstep, sc_window_fun);//TODO
        GPUSwFinalSum << <nw * nq, 1024 >> > (sc_block_w_gpu, sc_qw_gpu, numBlocksX_w, nq);
        //TODO cpuCorrelations.m_kw.copy_sync(sc_qw_gpu);
        break;

    case 'Y':
        dt.copy_sync(dt_cpu);
        tasks = 3 * nq;
        bl = (tasks + maxThreads - 1) / maxThreads;
        GPUSqAvrg << <bl, maxThreads >> > (sc_q_gpu, n_samples, tasks, M);
         //TODO cpuCorrelations.m_k.copy_sync(sc_q_gpu); //TODO: async
        GPUSwSum << <blocks_w, maxThreads >> > (sc_qt_gpu, dt, w, sc_block_w_gpu, tasksTot_w, sc_max_nstep, nq, sc_max_nstep, sc_window_fun);//TODO
        GPUSwFinalSum << <nw*nq, 1024 >> > (sc_block_w_gpu, sc_qw_gpu, numBlocksX_w, nq);//TODO blocks
        //TODO cpuCorrelations.m_kw.copy_sync(sc_qw_gpu);
        break;

    }
    cudaDeviceSynchronize();

}





