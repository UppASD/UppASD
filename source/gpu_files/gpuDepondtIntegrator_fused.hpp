#pragma once

#include "c_headers.hpp"

#include "tensor.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "gpuStructures.hpp"
#include "gpu_wrappers.h"
#include "depondt_kernels.hpp"
#if defined(HIP_V)
#include <hiprand/hiprand.h>
#elif defined(CUDA_V)
#include <curand.h>
#endif


class GpuDepondtIntegrator_fused {
private:
   // System parameters
   real gamma;
   real k_bolt;
   real mub;
   real damping;
   real dp_factor;

   // Integrator parameters
   char stt;  // Spin transfer torque mode
   real timestep;

   unsigned int N;
   unsigned int M;

   // Class local matrices
   deviceLattice& gpuLattice;
   GpuTensor<real, 3> thermfield;
   GpuTensor<real, 1> sigmaFactor;
   GpuTensor<GPU_RAND_STATE, 1> d_state;

   bool dataInitiated{false};
   bool constantsInitiated{false};
   bool applyTorque;

   dim3 threads;
   dim3 blocks;
   unsigned int tasks3NM;
   unsigned int tasksNM;
   

   // Timer
   StopwatchDeviceSync stopwatch;


public:

   // Constructor
   GpuDepondtIntegrator_fused(const SimulationParameters SimParam, deviceLattice& p_gpuLattice, const Tensor<real, 1> temperature);

   // Destructor
   ~GpuDepondtIntegrator_fused();

   void resetConstants(const Tensor<real, 1> temperature, real phaseTimestep, real phaseDamping);

   // Releaser
   void release();

   // Algorithm
   void evolveFirst(deviceLattice& gpuLattice);

   void evolveSecond(deviceLattice& gpuLattice);
};

