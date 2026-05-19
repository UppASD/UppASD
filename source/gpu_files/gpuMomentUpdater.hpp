#pragma once

#include "tensor.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "gpuStructures.hpp"
#include "gpuParallelizationHelper.hpp"
#include "gpu_wrappers.h"
using ParallelizationHelper = GpuParallelizationHelper;

class GpuMomentUpdater {
private:
   // Moments to update
   GpuTensor<real, 2>& mmom;
   GpuTensor<real, 2>& mmom0;
   GpuTensor<real, 2>& mmom2;
   GpuTensor<real, 3>& emom;
   GpuTensor<real, 3>& emom2;
   GpuTensor<real, 3>& emomM;
   GpuTensor<real, 2>& mmomi;

   // Parameters
   int mompar;
   char initexc;

   // Timer
   StopwatchDeviceSync stopwatch;

   // Parallelization helper
   ParallelizationHelper& parallel;

public:
   // Parallelization classes
   class Mompar1;
   class Mompar2;
   class Copy1;
   class Copy2;

   // Constructor
   GpuMomentUpdater(deviceLattice& gpuLattice, int mompar, char initexc);

   // Updater
   void update();
};


