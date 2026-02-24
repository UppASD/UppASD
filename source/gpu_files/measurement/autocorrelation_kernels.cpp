#pragma once
#include "autocorrelation_kernels.hpp"

__global__ void fill_spinwait(GpuTensor<real, 4> spinwait, const GpuTensor<real, 3> emom, const int taskNum, const int swIdx){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if(idx < taskNum){
    const uint N = static_cast<uint>(spinwait.extent(1));
    const uint M = static_cast<uint>(spinwait.extent(2));

    int sp = 3 * N;
    int mInd = idx / sp;
    int nn = (idx % sp);
    int nInd = nn / 3;
    int cInd = nn % 3;

    spinwait(cInd, nInd, mInd, swIdx) = emom(cInd, nInd, mInd);
    } 
}

__global__ void calc_autocorr_block(GpuTensor<real, 2> ac_block, const GpuTensor<real, 4> spinwait, const GpuTensor<real, 3> emom, const int taskNum, const int swIdx){
    //int idx = blockIdx.x*blockDim.x + threadIdx.x;

}

__global__ void calc_autocorr_final(GpuTensor<real, 2> ac_block, GpuTensor<real, 3> ac, const int taskNum, const int ac_count){

    
};
