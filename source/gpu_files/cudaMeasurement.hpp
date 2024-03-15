#pragma once

#include <cuda_runtime.h>

#include "c_headers.hpp"
#include "cudaEventPool.hpp"
#include "cudaMatrix.hpp"
#include "cudaParallelizationHelper.hpp"
#include "hostMatrix.hpp"
#include "measurementQueue.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"

#if defined(USE_FAST_COPY)
#define DEFAULT_FAST_COPY true
#else
#define DEFAULT_FAST_COPY false
#endif

class CudaMeasurement {
   // Queue callback data struct
   struct queue_callback_data {
      queue_callback_data(CudaMeasurement* m, std::size_t s) : me(m), step(s) {
      }

      CudaMeasurement* me;
      std::size_t step;
   };

   // Queue callback
   static void queue_callback(cudaStream_t, cudaError_t, void* data);

   // Temporary device storage vectors
   cudaMatrix<real, 3, 3> tmp_emomM;
   cudaMatrix<real, 3, 3> tmp_emom;
   cudaMatrix<real, 2> tmp_mmom;

   // Temporary host storage (pinned memory)
   real* pinned_emomM;
   real* pinned_emom;
   real* pinned_mmom;

   // Vectors to copy
   const cudaMatrix<real, 3, 3>& emomM;
   const cudaMatrix<real, 3, 3>& emom;
   const cudaMatrix<real, 2>& mmom;

   // Event stack
   CudaEventPool eventPool;

   // Measure queue
   MeasurementQueue measurementQueue;

   // Timer
   StopwatchDeviceSync stopwatch;

   // Parallelization helper
   CudaParallelizationHelper& parallel;

   // Control flags
   bool alwaysCopy;
   bool fastCopy;

   // Helpers
   void queueMeasurement(std::size_t mstep);
   void copyQueueFast(std::size_t mstep);
   void copyQueueSlow(std::size_t mstep);

public:
   // Constructor / deconstructor ???
   // TODO add flag for fast_copy
   CudaMeasurement(const cudaMatrix<real, 3, 3>& emomM, const cudaMatrix<real, 3, 3>& emom,
                   const cudaMatrix<real, 2>& mmom, bool fastCopy = DEFAULT_FAST_COPY,
                   bool alwaysCopy = false);
   ~CudaMeasurement();

   // Access methods
   void measure(std::size_t mstep);
   void flushMeasurements(std::size_t mstep);
};

