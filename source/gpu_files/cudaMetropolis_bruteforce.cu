#pragma once

#include <curand.h>
#include "cudaHamiltonianCalculations.hpp"

#include "c_headers.hpp"
#include "tensor.cuh"
#include "real_type.h"
#include "cudaStructures.hpp"
#include "cudaMetropolis_bruteforce.cuh"
#include "thrust/device_vector.h"
#include "thrust/host_vector.h"

__global__ void moms_bf(int tasks, CudaTensor<real, 2> mmom, CudaTensor<real, 3> emomM, CudaTensor<real, 3> emom, CudaTensor<real, 3> emom2,CudaTensor<real, 2> mmom0, CudaTensor<real, 2> mmom2, CudaTensor<real, 2> mmomi) {
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

__device__ void GetRandomXYZ_bf(real& x, real& y, real& z, real magnitude, curandState* d_state, int& idx) {
    real norm;
    x = curand_normal_double(d_state + idx);
    y = curand_normal_double(d_state + idx);
    z = curand_normal_double(d_state + idx);

    norm = magnitude / sqrt(x * x + y * y + z * z);
    x *= norm;
    y *= norm;
    z *= norm;
}

__global__ void InitGenerator_bf(curandState* state, unsigned long long seed, unsigned int taskMax) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (idx < taskMax)
        curand_init(seed, idx, 0, &state[idx]);
}

__global__ void MCSweep_bf(CudaTensor<curandState, 2> d_state, CudaTensor<real, 2> mmom, CudaTensor<real, 3> emomM, CudaTensor<real, 3> emom, CudaTensor<real, 3> eneff, real beta, unsigned int N, unsigned int tasknum, real k_bolt, real mub) {
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
CudaMetropolis_bruteforce::CudaMetropolis_bruteforce() {
    isallocated = false;
    thread_num = 256;
}
// Destructor
CudaMetropolis_bruteforce::~CudaMetropolis_bruteforce() {
    if (isallocated && !isfreed) { release(); }
}

void CudaMetropolis_bruteforce::rnd_init() {
    srand(time(NULL));
    unsigned long long seed = (unsigned long long)rand();
    InitGenerator_bf << <taskMax, 1 >> > (d_state.data(), seed, taskMax);
}

bool CudaMetropolis_bruteforce::initiate(const SimulationParameters SimParam) {

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

void CudaMetropolis_bruteforce::release() {
    if (isallocated && !isfreed) {
        d_state.Free();
        isfreed = true;
        isallocated = false;
    }
    // Assert that we're not already initialized

}

void CudaMetropolis_bruteforce::MCrun(cudaLattice& gpuLattice, real beta) {
    threads = { thread_num, 1, 1 };
      
    blocks = { static_cast <unsigned int>((N + thread_num - 1)/thread_num),  static_cast <unsigned int>(M), static_cast <unsigned int>(1) };
    //printf("blocks = %i, M = %i\n", static_cast <unsigned int>(block_subL_cpu(i)),  static_cast <unsigned int>(M));
    MCSweep_bf << <blocks, threads >> > (d_state, gpuLattice.mmom,gpuLattice.emomM,gpuLattice.emom, gpuLattice.eneff, beta, N, taskMax, k_bolt, mub);
    

    
}

void CudaMetropolis_bruteforce::mom_update(cudaLattice& gpuLattice){
    threads = {thread_num, 1, 1};
    blocks = {(N + thread_num - 1)/thread_num, M, 1};
     moms_bf << <blocks, threads >> > (N, gpuLattice.mmom, gpuLattice.emomM, gpuLattice.emom, gpuLattice.emom2, gpuLattice.mmom0, gpuLattice.mmom2, gpuLattice.mmomi);
}