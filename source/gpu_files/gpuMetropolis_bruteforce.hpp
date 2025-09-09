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
class GpuMetropolis_bruteforce {
private:
    // System parameters
    real k_bolt;
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
     GpuTensor<GPU_RAND_STATE, 2> d_state;

    // Timer
   // StopwatchDeviceSync stopwatch;

    // RNG initialization & splitting into sublattices 
     void rnd_init();

public:

    // Constructor
    GpuMetropolis_bruteforce();

    // Destructor
    ~GpuMetropolis_bruteforce();

    // Initiator
    bool initiate(const SimulationParameters SimParam);
   
    // Releaser
    void release();

    // Algorithm
    void MCrun(deviceLattice& gpuLattice, real beta);
    void mom_update(deviceLattice& gpuLattice);

};

