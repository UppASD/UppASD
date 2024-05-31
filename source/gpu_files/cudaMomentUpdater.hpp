#pragma once

#include "tensor.cuh"
#include "cudaParallelizationHelper.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "cudaStructures.hpp"

class CudaMomentUpdater {
private:
   // Moments to update
   CudaTensor<real, 2>& mmom;
   CudaTensor<real, 2>& mmom0;
   CudaTensor<real, 2>& mmom2;
   CudaTensor<real, 3>& emom;
   CudaTensor<real, 3>& emom2;
   CudaTensor<real, 3>& emomM;
   CudaTensor<real, 2>& mmomi;

   // Parameters
   int mompar;
   char initexc;

   // Timer
   StopwatchDeviceSync stopwatch;

   // Parallelization helper
   CudaParallelizationHelper& parallel;

public:
   // Parallelization classes
   class Mompar1;
   class Mompar2;
   class Copy1;
   class Copy2;

   // Constructor
   CudaMomentUpdater(cudaLattice& gpuLattice, int mompar, char initexc);

   // Updater
   void update();
};


