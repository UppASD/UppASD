#pragma once
#include "gpuHamiltonianCalculations.hpp"

#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "gpuStructures.hpp"

#include "device_launch_parameters.h"
#include "thrust/host_vector.h"

#include "gpu_wrappers.h"
#if defined(HIP_V)
#include <hip/hip_runtime.h>
#include <hiprand/hiprand.h>
#include <hip/hiprand_kernel.h>
#include <hip/hip_cooperative_groups.h>
//reduce??
#elif defined(CUDA_V)
#include "cuda_runtime.h"
#include <curand.h>
#include <curand_kernel.h>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#endif

class GpuMetropolis {
private:
    // System parameters
    real k_bolt;

    unsigned int num_subL;
    unsigned int max_spins;
    std::size_t N;
    std::size_t M;
    std::size_t mnn;
    unsigned int taskMax;
    bool isallocated;
    bool isfreed;
    unsigned int thread_num;
    dim3 blocks;
    dim3 threads;
    unsigned int task_num;
    real beta;
    real mub;


    // Class local matrices
    GpuTensor<unsigned int, 2> subIdx_gpu;
    Tensor<unsigned int, 1> subL_spnum_cpu;
    GpuTensor<unsigned int, 1> subL_spnum_gpu;
    Tensor<unsigned int, 1> block_subL_cpu;
    GpuTensor<GPU_RAND_STATE, 2> d_state;
    Tensor<unsigned int, 2> subIdx_cpu;


    // Timer
   // StopwatchDeviceSync stopwatch;

    // RNG initialization & splitting into sublattices 
    void fillneighbours_in_play(unsigned int* neigbours_in_play, Tensor<unsigned int, 2> nlist, int i);
    void split_lattice(const Tensor<unsigned int, 2> nlist);
    unsigned int increaseneighbours_in_play(unsigned int* neigbours_in_play, Tensor<unsigned int, 2> nlist, unsigned int mnn_cur, int i);
    void refillneigbours_in_play(unsigned int* neigbours_in_play, Tensor<unsigned int, 2> nlist, int i);
    void rnd_init();
    void count_spins();

public:

    // Constructor
    GpuMetropolis();

    // Destructor
    ~GpuMetropolis();

    // Initiator
    unsigned initiate(const SimulationParameters SimParam, const hostHamiltonian& cpuHam, const hostLattice& cpuLattice);
   
    // Releaser
    void release();

    // Algorithm
    void MCrun(deviceLattice& gpuLattice, real beta, unsigned int sub);
    void mom_update(deviceLattice& gpuLattice);

};

