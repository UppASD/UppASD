#pragma once


#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "gpuStructures.hpp"
#include <numeric>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <thrust/complex.h>
#include <curand.h>
#include "correlation.hpp"
#include "correlation_types.h"


//#include <complex>

class GpuCorrelations : public Correlation{
private:
    unsigned int isallocated;

    const GpuTensor<real, 3>& emomM;
    const GpuTensor<real, 3>& emom;
    const GpuTensor<real, 2>& mmom;

    const std::size_t N;
    const std::size_t M;
    const std::size_t nq;
    //std::size_t tidx;
    const std::size_t sc_step;
    const std::size_t sc_sep;

    std::size_t n_samples;
    const char do_sc;
    const char do_proj;
    const char do_projch;
    const std::size_t sc_max_nstep;
    const std::size_t sc_window_fun;
    const std::size_t nw;

    const std::size_t Nchmax;
    const std::size_t NT;

    std::size_t both_flag;
    real delta_t;

    unsigned int  t_cur;

    unsigned int maxThreads;
    unsigned int maxBlocks;
    unsigned int numThreads;


    blocksQW blQ;
    blocksQW blW;
    blocksQWproj blQproj;
    blocksQWproj blWproj;
    blocksQWproj blQprojch;
    blocksQWproj blWprojch;
    dim3 threads;

    GpuTensor<real, 1> r_mid;
    GpuTensor<real, 2> q;
    GpuTensor<real, 2> coord;
    GpuTensor<real, 1> dt;
    GpuTensor<real, 1> w;
    Tensor<real, 1> dt_cpu;
    Tensor<real, 1> sc_step_arr_cpu;  // Host buffer for sc_step array bookkeeping

    SC_proj sc_proj;
    SC_proj sc_projch;
    SC sc;

public:
    // Constructor
    GpuCorrelations(const Flag Flags, const SimulationParameters SimParam, const deviceLattice& gpuLattice, const hostCorrelations& cpuCorrelations);
    // Destructor
    ~GpuCorrelations() override;;

    // Initiator
    bool initiate(const Flag Flags, const SimulationParameters SimParam, const hostCorrelations& cpuCorrelations);
    // Releaser
    void release();
    // Measurements
    void measure(std::size_t mstep) override;
    void flushCorrelations(hostCorrelations& cpuCorrelations, std::size_t mstep) override;
    void recordSample();
    void publishSamplingInfo(hostCorrelations& cpuCorrelations);

private:
    void measure_SC(std::size_t mstep);
    void measure_SC_proj(std::size_t mstep, SC_proj& scp, blocksQWproj blQp, char sc_type);

    void flush_SC(std::size_t mstep, hostCorrelations& cpuCorrelations);
    void flush_SC_proj(std::size_t mstep, char p, hostCorrelations& cpuCorrelations, SC_proj& scp, blocksQWproj blWp, char sc_type);

};

