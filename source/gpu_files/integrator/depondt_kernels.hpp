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

__global__ void DepondtEvolveFirst(GpuTensor<real, 3> thermfield, const GpuTensor<real, 1> sigmaFactor, GpuTensor<GPU_RAND_STATE, 1> d_state,
                                  const GpuTensor<real, 3> beff, GpuTensor<real, 3> emom, GpuTensor<real, 3> emom2, GpuTensor<real, 2> mmom, 
                                  GpuTensor<real, 3> emomM, const GpuTensor<real, 3> btorque, int N, int M, bool applyTorque);

__global__ void DepondtEvolveSecond(); 