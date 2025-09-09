#pragma once

#include "gpuHamiltonianCalculations.hpp"
#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "gpuStructures.hpp"
#include "gpuMetropolis_bruteforce.hpp"
#include "thrust/device_vector.h"
#include "thrust/host_vector.h"

#include "gpu_wrappers.h"
#if defined(HIP_V)
#include <hiprand/hiprand.h>
#elif defined(CUDA_V)
#include <curand.h>
#endif

__global__ void moms_bf(int tasks, GpuTensor<real, 2> mmom, GpuTensor<real, 3> emomM, GpuTensor<real, 3> emom, GpuTensor<real, 3> emom2,GpuTensor<real, 2> mmom0, GpuTensor<real, 2> mmom2, GpuTensor<real, 2> mmomi) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (idx < tasks){
    unsigned int mInd = blockIdx.y;
    unsigned int nInd = idx;   
    //printf("mInd = %i, nInd = %i\n", mInd, nInd);
    emom(0, nInd, mInd) = emomM(0, nInd, mInd)/mmom(nInd, mInd);
    emom(1, nInd, mInd) = emomM(1, nInd, mInd)/mmom(nInd, mInd);
    emom(2, nInd, mInd) = emomM(2, nInd, mInd)/mmom(nInd, mInd);
    emom2(0, nInd, mInd) = emom(0, nInd, mInd);
    emom2(1, nInd, mInd) = emom(1, nInd, mInd);
    emom2(2, nInd, mInd) = emom(2, nInd, mInd);
    mmom0(nInd, mInd) = mmom(nInd, mInd);
    mmom2(nInd, mInd) = mmom(nInd, mInd);
    mmomi(nInd, mInd) = 1/mmom(nInd, mInd);

    }

}

__device__ void GetRandomXYZ_bf(real& x, real& y, real& z, real magnitude, GPU_RAND_STATE* d_state, int& idx) {
    real norm;
    x = GPU_NORMAL_DOUBLE(d_state + idx);
    y = GPU_NORMAL_DOUBLE(d_state + idx);
    z = GPU_NORMAL_DOUBLE(d_state + idx);

    norm = magnitude / sqrt(x * x + y * y + z * z);
    x *= norm;
    y *= norm;
    z *= norm;
}

__global__ void InitGenerator_bf(GPU_RAND_STATE* state, unsigned long long seed, unsigned int taskMax) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (idx < taskMax)
        GPU_RAND_INIT(seed, idx, 0, &state[idx]);
}

__global__ void MCSweep_bf(GpuTensor<GPU_RAND_STATE, 2> d_state, GpuTensor<real, 2> mmom, GpuTensor<real, 3> emomM, GpuTensor<real, 3> emom, GpuTensor<real, 3> eneff, real beta, unsigned int N, unsigned int tasknum, real k_bolt, real mub) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (idx < N) {
        unsigned int mInd = blockIdx.y;
        unsigned int nInd = idx;
        //unsigned int offset = 3 * (nInd + mInd * N);
        //unsigned int actId_x = offset;
        //unsigned int actId_y = 1 + offset;
        //unsigned int actId_z = 2 + offset;
        real xNew, yNew, zNew;
        //printf("subL = %i, idx = %i, nInd = %i, mInd = %i\n", subL, idx, nInd, mInd);
        real fx = eneff(0, nInd, mInd);
        real fy = eneff(1, nInd, mInd);
        real fz = eneff(2, nInd, mInd);
        real magnitude = mmom(nInd, mInd);
        int dIdx = mInd * N + idx;
       // printf("dInd = %i, beta = %lf\n",dIdx,  beta) ;
        GetRandomXYZ_bf(xNew, yNew, zNew, magnitude, d_state.data(), dIdx);

        real hamOld = fx * emomM(0, nInd, mInd) + fy * emomM(1, nInd, mInd) + fz * emomM(2, nInd, mInd);
        real hamNew = fx * xNew + fy * yNew + fz * zNew;
       // printf("exp = %lf\n", -beta * mub*(hamNew - hamOld));
        if (exp(beta * mub*(hamNew - hamOld)) < curand_uniform_double(d_state.data() + dIdx)) {
          //  printf("NOT sweeped\n");

            return;
        }
        //printf("sweeped\n");
        emomM(0, nInd, mInd) = xNew;
        emomM(1, nInd, mInd) = yNew;
        emomM(2, nInd, mInd) = zNew;

        emom(0, nInd, mInd) = xNew/magnitude;
        emom(1, nInd, mInd) = yNew/magnitude;
        emom(2, nInd, mInd) = zNew/magnitude;
        //printf("magnitude_real = %lf, magnitude_cur = %lf\n", magnitude, sqrt(xNew*xNew + yNew*yNew+ zNew*zNew));
    }
}

// Constructor
GpuMetropolis_bruteforce::GpuMetropolis_bruteforce() {
    isallocated = false;
    thread_num = 256;
}
// Destructor
GpuMetropolis_bruteforce::~GpuMetropolis_bruteforce() {
    if (isallocated && !isfreed) { release(); }
}

void GpuMetropolis_bruteforce::rnd_init() {
    srand(time(NULL));
    unsigned long long seed = (unsigned long long)rand();
    InitGenerator_bf << <taskMax, 1 >> > (d_state.data(), seed, taskMax);
}

bool GpuMetropolis_bruteforce::initiate(const SimulationParameters SimParam) {

    // Assert that we're not already initialized
    release();
    N = SimParam.N;
    M = SimParam.M;
    k_bolt = SimParam.k_bolt;
    mnn = SimParam.mnn;
    mub = SimParam.mub;

    d_state.Allocate(static_cast <long int>(N), static_cast <long int>(M));

    isallocated = false;
    isfreed = true;
    taskMax = M * N;

    rnd_init();

   if(d_state.empty()) {
      isallocated = false;
       std::printf("Unable to allocate d_state in brute force Metropolis\n");  
      return false;
   }
    return true;

}

void GpuMetropolis_bruteforce::release() {
    if (isallocated && !isfreed) {
        d_state.Free();
        isfreed = true;
        isallocated = false;
    }
    // Assert that we're not already initialized

}

void GpuMetropolis_bruteforce::MCrun(deviceLattice& gpuLattice, real beta) {
    threads = { thread_num, 1, 1 };
      
    blocks = { static_cast <unsigned int>((N + thread_num - 1)/thread_num),  static_cast <unsigned int>(M), static_cast <unsigned int>(1) };
    //printf("blocks = %i, M = %i\n", static_cast <unsigned int>(block_subL_cpu(i)),  static_cast <unsigned int>(M));
    MCSweep_bf << <blocks, threads >> > (d_state, gpuLattice.mmom,gpuLattice.emomM,gpuLattice.emom, gpuLattice.eneff, beta, N, taskMax, k_bolt, mub);
    

    
}

void GpuMetropolis_bruteforce::mom_update(deviceLattice& gpuLattice){
    threads = {thread_num, 1, 1};
    blocks = {(N + thread_num - 1)/thread_num, M, 1};
     moms_bf << <blocks, threads >> > (N, gpuLattice.mmom, gpuLattice.emomM, gpuLattice.emom, gpuLattice.emom2, gpuLattice.mmom0, gpuLattice.mmom2, gpuLattice.mmomi);
}