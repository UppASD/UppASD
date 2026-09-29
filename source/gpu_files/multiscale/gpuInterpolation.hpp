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

class GpuInterpolation {
private:

    unsigned int maxThreads;
    unsigned int maxBlocks;
    dim3 threads;
    dim3 blocks;
    const unsigned int N;
    const unsigned int M;
    deviceInterpolationInfo& gpuInterpolationInfo;
    GpuTensor<real, 4> backbuffer;
    int backbufferHead;
    deviceLattice& gpuLattice;
    


public:
    GpuInterpolation(const unsigned int p_N, const unsigned int p_M, deviceInterpolationInfo& p_gpuInterpolationInfo, 
                     deviceLattice& p_gpuLattice);
    ~GpuInterpolation();

    void interpolate();  
};

