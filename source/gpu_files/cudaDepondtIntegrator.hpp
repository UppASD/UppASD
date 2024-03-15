#pragma once

#include <curand.h>

#include "cudaMatrix.hpp"
#include "cudaParallelizationHelper.hpp"
#include "cudaThermfield.hpp"
#include "fortMatrix.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"

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
   cudaMatrix<real, 3, 3> mrod;
   cudaMatrix<real, 3, 3> blocal;
   cudaMatrix<real, 3, 3> bdup;

   // Thermfield
   CudaThermfield thermfield;

   // Timer
   StopwatchDeviceSync stopwatch;

   // Parallelization helper
   CudaParallelizationHelper& parallel;

   // Algorithm
   void rotate(const cudaMatrix<real, 3, 3>& emom, real delta_t);
   void buildbeff(const cudaMatrix<real, 3, 3>& emom, const cudaMatrix<real, 3, 3>& btorque);

public:
   // Parallelization helpers
   class Rotate;
   class BuildEffectiveField;

   // Constructor
   CudaDepondtIntegrator();

   // Destructor
   ~CudaDepondtIntegrator();

   // Initiator
   bool initiate(std::size_t N, std::size_t M, char _stt, real _timestep, curandRngType_t rng,
                 unsigned long long seed);

   // Set up constants
   bool initiateConstants(const fortMatrix<real, 1>& temperature, real timestep, real gamma_const,
                          real k_bolt_const, real mub_const, real damping_const);

   // Releaser
   void release();

   // Algorithm
   void evolveFirst(const cudaMatrix<real, 3, 3>& beff, cudaMatrix<real, 3, 3>& b2eff,
                    const cudaMatrix<real, 3, 3>& btorque, cudaMatrix<real, 3, 3>& emom,
                    cudaMatrix<real, 3, 3>& emom2, cudaMatrix<real, 3, 3>& emomM,
                    const cudaMatrix<real, 2>& mmom);

   void evolveSecond(const cudaMatrix<real, 3, 3>& beff, const cudaMatrix<real, 3, 3>& b2eff,
                     const cudaMatrix<real, 3, 3>& btorque, cudaMatrix<real, 3, 3>& emom,
                     cudaMatrix<real, 3, 3>& emom2);
};

