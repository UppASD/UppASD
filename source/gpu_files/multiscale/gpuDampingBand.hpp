#pragma once

#if defined (CUDA_V)
#include <curand.h>
#elif defined (HIP_V)
 #include <hiprand/hiprand.h>
#endif
#include "c_headers.hpp"
#include "tensor.hpp"
#include "gpuStructures.hpp"
#include "real_type.h"
#include "gpu_wrappers.h"
#include "interpolation_kernels.hpp"

class GpuDampingBand {
private:

    unsigned int maxThreads;
    unsigned int maxBlocks;
    dim3 threads;
    dim3 blocks;
    const unsigned int N;
    const unsigned int M;
    GpuTensor<real, 4>& backbuffer;
    unsigned int backbufferHead;
    deviceDampingBand& gpuDampingBand;
    deviceLattice& gpuLattice;
    


public:
    GpuDampingBand(const unsigned int p_N, const unsigned int p_M, deviceDampingBand& p_gpuDampingBand, 
                     deviceLattice& p_gpuLattice, GpuTensor<real, 4> p_backbuffer, unsigned int p_backbufferHead);
    ~GpuDampingBand();

    void preInterpolation();  
    void corrPreInterpolation();  
    void corrPreInterpolationAvg();  
};

