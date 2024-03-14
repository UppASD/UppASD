#pragma once

#include "cudaMatrix.hpp"
#include "cudaParallelizationHelper.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"

class CudaMomentUpdater {
private:
   // Moments to update
   cudaMatrix<real, 2>& mmom;
   cudaMatrix<real, 2>& mmom0;
   cudaMatrix<real, 2>& mmom2;
   cudaMatrix<real, 3, 3>& emom;
   cudaMatrix<real, 3, 3>& emom2;
   cudaMatrix<real, 3, 3>& emomM;
   cudaMatrix<real, 2>& mmomi;

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
   CudaMomentUpdater(cudaMatrix<real, 2>& mmom, cudaMatrix<real, 2>& mmom0, cudaMatrix<real, 2>& mmom2,
                     cudaMatrix<real, 3, 3>& emom, cudaMatrix<real, 3, 3>& emom2,
                     cudaMatrix<real, 3, 3>& emomM, cudaMatrix<real, 2>& mmomi, int mompar, char initexc);

   // Updater
   void update();
};

