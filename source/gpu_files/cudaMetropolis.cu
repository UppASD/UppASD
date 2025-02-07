#pragma once

#include <curand.h>

#include "c_headers.hpp"
#include "tensor.cuh"
#include "real_type.h"
#include "cudaStructures.hpp"
#include "cudaMetropolis.cuh"
#include "thrust/device_vector.h"
#include "thrust/host_vector.h"

__global__ void moms(int tasks, CudaTensor<real, 2> mmom, CudaTensor<real, 3> emomM, CudaTensor<real, 3> emom, CudaTensor<real, 3> emom2,CudaTensor<real, 2> mmom0, CudaTensor<real, 2> mmom2, CudaTensor<real, 2> mmomi) {
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


__device__ void GetRandomXYZ(real& x, real& y, real& z, real magnitude, curandState* d_state, int& idx) {
    real norm;
    x = curand_normal_double(d_state + idx);
    y = curand_normal_double(d_state + idx);
    z = curand_normal_double(d_state + idx);

    norm = magnitude / sqrt(x * x + y * y + z * z);
    x *= norm;
    y *= norm;
    z *= norm;
}

__global__ void InitGenerator(curandState* state, unsigned long long seed, unsigned int taskMax) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (idx < taskMax)
        curand_init(seed, idx, 0, &state[idx]);
}

__global__ void MCSweep(CudaTensor<curandState, 2> d_state, CudaTensor<real, 2> mmom, CudaTensor<real, 3> emomM, CudaTensor<real, 3> emom, CudaTensor<real, 3> eneff, CudaTensor<unsigned int, 2> subLIdx, real beta, unsigned int N, unsigned int tasknum, unsigned int subL, unsigned int max_spins, real k_bolt, real mub) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (idx < tasknum) {
        unsigned int mInd = blockIdx.y;
        unsigned int nInd = subLIdx(idx, subL);
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
        int dIdx = mInd * max_spins + idx;
       // printf("dInd = %i, beta = %lf\n",dIdx,  beta) ;
        GetRandomXYZ(xNew, yNew, zNew, magnitude, d_state.data(), dIdx);

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
        //printf("magbnitude_real = %lf, magnitude_cur = %lf\n", magnitude, sqrt(xNew*xNew + yNew*yNew+ zNew*zNew));
    }
}
// Constructor
CudaMetropolis::CudaMetropolis() {
    isallocated = false;
    thread_num = 256;
}
// Destructor
CudaMetropolis::~CudaMetropolis() {
    if (isallocated && !isfreed) { release(); }
}
void CudaMetropolis::fillneighbours_in_play(unsigned int* neigbours_in_play, Tensor<unsigned int, 2> nlist, int i) {
    for (int k = 0; k < mnn; k++) {
        neigbours_in_play[k] = nlist(k, i);
        //printf("%i ", neigbours_in_play[k]);
    }
    //printf("\n");
}
unsigned int CudaMetropolis::increaseneighbours_in_play(unsigned int* neigbours_in_play, Tensor<unsigned int, 2> nlist, unsigned int mnn_cur, int i) {
    int nli = 0;
    int npi = mnn_cur;
    for (int k = 0; k < mnn; k++) {
        if (nlist(k, i) == 0) {
            nli = k + 1;
            break;
        }
    }
    int ext_nl = (nli == 0) ? mnn : nli;


    for (int k = 0; k < mnn_cur; k++) {
        if (neigbours_in_play[k] == 0) {
            npi = k;
            break;
        }
    }

    int ext_np = mnn_cur - npi;
    int increase = (ext_np > ext_nl) ? 0 : (ext_nl - ext_np);//?????????? +1???????
    //printf("k = %i, increase = %i\n", i, increase);
    return (mnn_cur + increase);
}
void CudaMetropolis::refillneigbours_in_play(unsigned int* neigbours_in_play, Tensor<unsigned int, 2> nlist, int i) {
    int k = 0;
    int l = 0;
    while (neigbours_in_play[k] != 0) k++;
    //k++;
    while ((nlist(l, i) != 0) && (l < mnn)) {
        neigbours_in_play[k] = nlist(l, i);
        k++; l++;
    }

}
void CudaMetropolis::split_lattice(const Tensor<unsigned int, 2> nlist) {
    Tensor<unsigned int, 1> used_sites;
    used_sites.AllocateHost(static_cast <long int>(N));
    used_sites.zeros();
    unsigned int rw = 0;
    int flag = 0;
    int flag2 = 0;
    unsigned int spins_in_subL = 0;
    unsigned int used_num = 0;
    unsigned int mnn_cur = mnn;
    thrust::host_vector<unsigned int> neigbours_in_play(mnn);
    thrust::fill(neigbours_in_play.begin(), neigbours_in_play.end(), 0);
    thrust::host_vector<unsigned int> m_ind(N);
    thrust::fill(m_ind.begin(), m_ind.end(), 0);
    printf("\n");

    for (int i = 0; i < N; i++) {

        flag = 0;
        spins_in_subL = 0;
        neigbours_in_play.resize(mnn);
        thrust::fill(neigbours_in_play.begin(), neigbours_in_play.end(), 0);
        fillneighbours_in_play(neigbours_in_play.data(), nlist, i);
        //for (int lll = 0; lll < neigbours_in_play.size(); lll++) printf("%i ", neigbours_in_play[lll]);
       // printf("\n\n");
        if (used_sites[i] != 1) {
            flag2 = 1;
            used_sites[i] = 1;
            m_ind[rw * N + i] = 1;
            used_num++;
            num_subL++;
            spins_in_subL++;
        }
        else {
            //rw++;
            continue;
        }
        for (int k = i + 1; k < N; k++) {
            //printf("i = %i, k = %i\n", i, k);
            flag = 0;
            for (int ll = 0; ll < neigbours_in_play.size(); ll++) {
                if (k == (neigbours_in_play[ll] - 1)) {
                    flag = 1;
                    //printf("i = %i, k = %i\n", i, k);
                    break;
                }
            }
            if ((used_sites(k) != 1) && (!flag)) {
                m_ind[rw * N + k] = 1;
                mnn_cur = increaseneighbours_in_play(neigbours_in_play.data(), nlist, neigbours_in_play.size(), k);
                neigbours_in_play.resize(mnn_cur);
                //for (int lll = 0; lll < neigbours_in_play.size(); lll++) printf("%i ", neigbours_in_play[lll]);
               // printf("\n\n");
                refillneigbours_in_play(neigbours_in_play.data(), nlist, k);
                //for (int lll = 0; lll < neigbours_in_play.size(); lll++) printf("%i ", neigbours_in_play[lll]);
                //printf("\n\n");
                used_sites[k] = 1;
                used_num++;
                spins_in_subL++;
                //printf("i = %i, k = %i, used = %i, spins_in_subL = %i\n", i, k, used_num, spins_in_subL);

                if (used_num == N) break;

            }
            max_spins = std::max(max_spins, spins_in_subL);
            if (used_num == N) break;
        }
        if (flag2) {
            rw++;
            m_ind.resize(m_ind.size() + N);
            //printf("\nHERe\n");
        }
        //neigbours_in_play.zeros();
        if (used_num == N) break;
    }


    printf("number of sublattices = %i, max spins = %i\n", num_subL, max_spins);
    //for (int i = 0; i < num_subL; i++) {
    //    for (int k = 0; k < N; k++) {
    //        printf("%i ", m_ind[i * N + k]);
    //    }
    //    printf("\n");
    //}

    //subIdx_cpu.resize(num_subL * max_spins);
    subIdx_cpu.resize(Extents<2>{static_cast <long int>(max_spins), static_cast <long int>(num_subL)}, false);
    int si = 0; int sk = 0; flag = 0;
    subIdx_cpu.zeros();

    for (int i = 0; i < num_subL; i++) {
        sk = 0; flag = 0;
        for (int k = 0; k < N; k++) {
            if (m_ind[i * N + k] != 0) {
                subIdx_cpu[si * max_spins + sk] = k;
                //printf("i = %i, k = %i, si = %i, sk = %i, subL = %i, idx = %i \n", i, k, si, sk, subIndex(si, sk), i * N + k);
                //printf("i = %i,  k = %i, si = %i, sk = %i, subL = %i, \n", i, k, si, sk, subIndex(0, 1));
                sk++;
                flag = 1;
            }

            if (sk == max_spins) break;
        }
        si++;

        //if (flag == 1) si++;
        //printf("i = %i,  si = %i, sk = %i, subL = %i, \n", i, si, sk, subIndex(0, 1));
        //if (si == num_subL) break;
    }

     for (int i = 0; i < num_subL; i++) {
         for (int k = 0; k < max_spins; k++) {
             printf("%i ", subIdx_cpu[i * max_spins + k]);
         }
         printf("\n");
     }

    used_sites.FreeHost();
    // return subIdx_cpu;


}
void CudaMetropolis::count_spins() {
    unsigned int countL = subIdx_cpu.extent(1);
    unsigned int max_spins = subIdx_cpu.extent(0);
    unsigned int cur_spins = 1;
    for (int i = 0; i < countL; i++) {
        cur_spins = 1;
        for (int k = 1; k < max_spins; k++) {
            // printf("i = %i, k = %i, subL_idx_cpu(k, i) = %i\n", i, k, subL_idx_cpu(k, i));

            if (subIdx_cpu(k, i) == 0) break;
            cur_spins++;

        }
        subL_spnum_cpu(i) = cur_spins;
        block_subL_cpu(i) = (cur_spins + thread_num - 1) / thread_num;
    }

}

void CudaMetropolis::rnd_init() {
    srand(time(NULL));
    unsigned long long seed = (unsigned long long)rand();
    InitGenerator << <taskMax, 1 >> > (d_state.data(), seed, taskMax);
};
bool CudaMetropolis::initiate(const SimulationParameters SimParam, const hostHamiltonian& cpuHam, const hostLattice& cpuLattice) {

    // Assert that we're not already initialized
    release();
    num_subL = 0;
    max_spins = 0;
    N = SimParam.N;
    M = SimParam.M;
    k_bolt = SimParam.k_bolt;
    mnn = SimParam.mnn;
    mub = SimParam.mub;


    subIdx_cpu.AllocateHost(static_cast <long int>(1), static_cast <long int>(1));
    split_lattice(cpuHam.nlist);

    subIdx_gpu.Allocate(static_cast <long int>(max_spins), static_cast <long int>(num_subL));
    subL_spnum_cpu.AllocateHost(static_cast <long int>(num_subL));
    subL_spnum_gpu.Allocate(static_cast <long int>(num_subL));
    d_state.Allocate(static_cast <long int>(max_spins), static_cast <long int>(M));
    block_subL_cpu.AllocateHost(static_cast <long int>(num_subL));


    isallocated = false;
    isfreed = true;
    taskMax = M * max_spins;

    subIdx_gpu.copy_sync(subIdx_cpu);

    count_spins();
    rnd_init();
    subL_spnum_gpu.copy_sync(subL_spnum_cpu);

   if(subL_spnum_gpu.empty()) {
      isallocated = false;
       //std::printf("HERE - 1\n");  
      return false;
   }
    //printf("\n\n");
    //for (int i = 0; i < num_subL; i++) {
    //    for (int k = 0; k < max_spins; k++) {
    //        printf("%i ", subIdx_cpu(k, i));
    //    }
    //    printf("\n");
    //}
    return true;

}
void CudaMetropolis::release() {
    if (isallocated && !isfreed) {
        subIdx_gpu.Free();
        subL_spnum_cpu.FreeHost();
        subL_spnum_gpu.Free();
        d_state.Free();
        block_subL_cpu.FreeHost();
        subIdx_cpu.FreeHost();
        isfreed = true;
        isallocated = false;
    }
    // Assert that we're not already initialized

}

void CudaMetropolis::MCrun(cudaLattice& gpuLattice, real beta) {
    threads = { thread_num, 1, 1 };

    for (unsigned int i = 0; i < num_subL; i++) {
        
        blocks = { static_cast <unsigned int>(block_subL_cpu(i)),  static_cast <unsigned int>(M), static_cast <unsigned int>(1) };
        task_num = subL_spnum_cpu(i);
        //printf("blocks = %i, M = %i\n", static_cast <unsigned int>(block_subL_cpu(i)),  static_cast <unsigned int>(M));
        MCSweep << <blocks, threads >> > (d_state, gpuLattice.mmom,gpuLattice.emomM,gpuLattice.emom, gpuLattice.eneff, subIdx_gpu, beta, N, task_num, i, max_spins, k_bolt, mub);
    }
}

void CudaMetropolis::mom_update(cudaLattice& gpuLattice){
    threads = {thread_num, 1, 1};
    blocks = {(N + thread_num - 1)/thread_num, M, 1};
     moms << <blocks, threads >> > (N, gpuLattice.mmom, gpuLattice.emomM, gpuLattice.emom, gpuLattice.emom2, gpuLattice.mmom0, gpuLattice.mmom2, gpuLattice.mmomi);
}