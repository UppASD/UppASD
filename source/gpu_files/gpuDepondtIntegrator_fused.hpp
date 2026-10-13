#pragma once

#include "c_headers.hpp"

#include "tensor.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "gpuStructures.hpp"
#include "gpu_wrappers.h"
#include "depondt_kernels.h"
#if defined(HIP_V)
#include <hiprand/hiprand.h>
#elif defined(CUDA_V)
#include <curand.h>
#endif


class GpuDepondtIntegrator {
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

   // Class local matrices
   GpuTensor<real, 3> mrod;
   GpuTensor<real, 3> blocal;
   GpuTensor<real, 3> bdup;
ß

   // Timer
   StopwatchDeviceSync stopwatch;


   // Algorithm
   void rotate(const GpuTensor<real, 3>& emom, real delta_t);
   void buildbeff(const GpuTensor<real, 3>& emom, const GpuTensor<real, 3>& btorque);

public:

   // Constructor
   GpuDepondtIntegrator();

   // Destructor
   ~GpuDepondtIntegrator();

   // Initiator
   bool initiate(const SimulationParameters SimParam);

   // Set up constants
   bool initiateConstants(const SimulationParameters SimParam, const Tensor<real, 1>temperature);
   void resetConstants(const Tensor<real, 1> temperature, real phaseTimestep, real phaseDamping);

   // Releaser
   void release();

   // Algorithm
   void evolveFirst(deviceLattice& gpuLattice);

   void evolveSecond(deviceLattice& gpuLattice);
};

