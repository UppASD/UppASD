#pragma once


#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "gpuStructures.hpp"
#include <numeric>
//#include <cooperative_groups.h>
//#include <cooperative_groups/reduce.h>
//#include <thrust/complex.h>


#include <hiprand/hiprand.h>


//#include <complex>

class GpuCorrelations {
private:
    unsigned int isallocated;

    unsigned int numThreads;
    unsigned int numBlocksX_q;
    unsigned int numBlocksX_w;
    unsigned int numBlocksY_q;
    unsigned int numBlocksY_w;
    //std::size_t step_num;
    //std::size_t cc_step;
    std::size_t N;
    std::size_t M;
    std::size_t nq;
    //std::size_t tidx;
    std::size_t sc_step;
    std::size_t sc_sep;

    std::size_t n_samples;
    char do_sc;
    std::size_t sc_max_nstep;
    std::size_t sc_window_fun;
    std::size_t nw;
    std::size_t both_flag;
    real delta_t;

    unsigned int  t_cur;
    //unsigned int  spinTot;
    unsigned int  tasksTot_q;
    unsigned int  tasksTot_w;
    unsigned int maxThreads;
    unsigned int maxBlocks;
    dim3 blocks_q;
    dim3 blocks_w;
    dim3 threads;

    //real nainv;
    //thrust::complex<real> iqfac;

    //Block variables TODO
    //GpuTensor<thrust::complex<real>, 2> sc_block_gpu;
    //GpuTensor<thrust::complex<real>, 3> sc_block_w_gpu;
    //GpuTensor<thrust::complex<real>, 2> sc_q_gpu;
    //GpuTensor<thrust::complex<real>, 3> sc_qt_gpu;
    //GpuTensor<thrust::complex<real>, 3> sc_qw_gpu;


    // Buffer variables 
    //Tensor<real, 2> SC;     // 3 x buff x M
    GpuTensor<real, 1> r_mid;
    GpuTensor<real, 2> q;
    GpuTensor<real, 2> coord;
    GpuTensor<real, 1> dt;
    GpuTensor<real, 1> w;
    Tensor<real, 1> dt_cpu;

public:
    // Constructor
    GpuCorrelations();
    // Destructor
    ~GpuCorrelations();

    // Initiator
    bool initiate(const Flag Flags, const SimulationParameters SimParam, const hostCorrelations& cpuCorrelations);
    // Releaser
    void release();
    // Measurements
    void runCorrelation_spin(const deviceLattice& gpuLattice, const int curstep);
    void finalCorrelation(hostCorrelations& cpuCorrelations);
};

