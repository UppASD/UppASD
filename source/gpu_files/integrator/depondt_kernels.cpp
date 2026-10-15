#include "depondt_kernals.hpp"


__global__ void InitializeDeviceRngStates(GpuRngState* states,
                                          unsigned long long seed,
                                          unsigned int count) {
   const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
   if (id < count) {
      GPU_RAND_INIT(seed, static_cast<unsigned long long>(id), 0ULL, &states[id]);
   }
}

__global__ void InitializeThermalSigmaFactor(GpuTensor<real, 1> sigmaFactor,
                                             const GpuTensor<real, 1> temperature,
                                             unsigned int count,
                                             real damping
                                          ) {
   const unsigned int site = blockIdx.x * blockDim.x + threadIdx.x;
   if (site < count) {
      sigmaFactor[site] = sqrt(damping * temperature[site]);
   }
}

__device__ inline void rotateDepondtVector(){}

__device__ inline void BuildEffectiveField(){}

__device__ inline void BuildThermalField(){


}

__global__ void DepondtEvolveFirst(){


   
};

__global__ void DepondtEvolveSecond(){}; 

