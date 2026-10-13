#include "depondt_kernals.hpp"


__global__ void InitializeDeviceRngStates(GpuRngState* states,
                                          unsigned int count,
                                          unsigned long long seed) {
   const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
   if (id < count) {
      GPU_RAND_INIT(seed, static_cast<unsigned long long>(id), 0ULL, &states[id]);
   }
}

__global__ void InitializeThermalSigmaFactor(real* sigmaFactor,
                                             unsigned int count,
                                             real damping) {
   const unsigned int site = blockIdx.x * blockDim.x + threadIdx.x;
   if (site < count) {
      sigmaFactor[site] = sqrt(damping * sigmaFactor[site]);
   }
}

__device__ inline void rotateDepondtVector(){}

__device__ inline void BuildEffectiveField(){}

_global__ void DepondtEvolveFirst(){};

_global__ void DepondtEvolveSecond(){}; 

