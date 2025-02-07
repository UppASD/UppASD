#pragma once

#include <curand.h>

#include "c_headers.hpp"
#include "cudaParallelizationHelper.hpp"
#include "cudaThermfield.hpp"
#include "tensor.cuh"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "cudaStructures.hpp"

class CudaDepondtIntegrator {
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
   CudaTensor<real, 3> mrod;
   CudaTensor<real, 3> blocal;
   CudaTensor<real, 3> bdup;

   // Thermfield
   CudaThermfield thermfield;

   // Timer
   StopwatchDeviceSync stopwatch;

   // Parallelization helper
   CudaParallelizationHelper& parallel;

   // Algorithm
   void rotate(const CudaTensor<real, 3>& emom, real delta_t);
   void buildbeff(const CudaTensor<real, 3>& emom, const CudaTensor<real, 3>& btorque);

public:
   // Parallelization helpers
   class Rotate;
   class BuildEffectiveField;

   // Constructor
   CudaDepondtIntegrator();

   // Destructor
   ~CudaDepondtIntegrator();

   // Initiator
   bool initiate(const SimulationParameters SimParam);

   // Set up constants
   bool initiateConstants(const SimulationParameters SimParam, const hostLattice& cpuLattice);

   // Releaser
   void release();

   // Algorithm
   void evolveFirst(cudaLattice& gpuLattice);

   void evolveSecond(cudaLattice& gpuLattice);
};

