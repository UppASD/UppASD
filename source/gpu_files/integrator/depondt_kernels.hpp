#pragma once

#include "gpu_wrappers.h"
#include "real_type.h"

#if defined(HIP_V)
#include <hiprand/hiprand_kernel.h>
#elif defined(CUDA_V)
#include <curand_kernel.h>
#endif

using GpuRngState = GPU_RAND_STATE;

__global__ void InitializeDeviceRngStates(GpuRngState* states,
                                         unsigned int count,
                                         unsigned long long seed);

__global__ void InitializeThermalSigmaFactor(real* sigmaFactor,
                                             unsigned int count,
                                             real damping);

__global__ void DepondtEvolveFirst();

__global__ void DepondtEvolveSecond(); 