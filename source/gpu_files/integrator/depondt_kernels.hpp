#pragma once

#include "gpu_wrappers.h"
#include "real_type.h"
#include "tensor.hpp"


#if defined(HIP_V)
#include <hiprand/hiprand_kernel.h>
#elif defined(CUDA_V)
#include <curand_kernel.h>
#endif

using GpuRngState = GPU_RAND_STATE;

__global__ void InitializeDeviceRngStates(GpuRngState* states,                                       
                                         unsigned long long seed,
                                         unsigned int count,);

__global__ void InitializeThermalSigmaFactor(GpuTensor<real, 1> sigmaFactor,
                                             const GpuTensor<real, 1> temperature,
                                             unsigned int count,
                                             real damping);

__global__ void DepondtEvolveFirst();

__global__ void DepondtEvolveSecond(); 