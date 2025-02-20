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

class CudaMetropolis {
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
    CudaTensor<unsigned int, 2> subIdx_gpu;
    Tensor<unsigned int, 1> subL_spnum_cpu;
    CudaTensor<unsigned int, 1> subL_spnum_gpu;
    Tensor<unsigned int, 1> block_subL_cpu;
    CudaTensor<curandState, 2> d_state;
    Tensor<unsigned int, 2> subIdx_cpu;


    // Timer
   // StopwatchDeviceSync stopwatch;

    // RNG initialization & splitting into sublattices 
    void fillneighbours_in_play(unsigned int* neigbours_in_play, Tensor<unsigned int, 2> nlist, int i);
    void split_lattice(const Tensor<unsigned int, 2> nlist);
    unsigned int increaseneighbours_in_play(unsigned int* neigbours_in_play, Tensor<unsigned int, 2> nlist, unsigned int mnn_cur, int i);
    void refillneigbours_in_play(unsigned int* neigbours_in_play, Tensor<unsigned int, 2> nlist, int i);
    void count_spins();
    void rnd_init();

public:

    // Constructor
    CudaMetropolis();

    // Destructor
    ~CudaMetropolis();

    // Initiator
    bool initiate(const SimulationParameters SimParam, const hostHamiltonian& cpuHam, const hostLattice& cpuLattice);

    // Releaser
    void release();

    // Algorithm
    void MCrun(cudaLattice& gpuLattice, real beta, CudaHamiltonianCalculations& hamCalc);
    void mom_update(cudaLattice& gpuLattice);

};

