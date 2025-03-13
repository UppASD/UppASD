#pragma once
#include "cudaHamiltonianCalculations.hpp"

#include <curand.h>
#include <curand_kernel.h>
#include "c_headers.hpp"
#include "tensor.cuh"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "cudaStructures.hpp"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include "thrust/host_vector.h"

class CudaMetropolis_bruteforce {
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
     CudaTensor<curandState, 2> d_state;

    // Timer
   // StopwatchDeviceSync stopwatch;

    // RNG initialization & splitting into sublattices 
     void rnd_init();

public:

    // Constructor
    CudaMetropolis_bruteforce();

    // Destructor
    ~CudaMetropolis_bruteforce();

    // Initiator
    bool initiate(const SimulationParameters SimParam);
   
    // Releaser
    void release();

    // Algorithm
    void MCrun(cudaLattice& gpuLattice, real beta);
    void mom_update(cudaLattice& gpuLattice);

};

