#pragma once

#include "c_headers.hpp"
#include "gpuEventPool.hpp"
#include "tensor.hpp"
#include "measurementQueue.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "gpu_wrappers.h"
#include "gpuParallelizationHelper.hpp"
#if defined (HIP_V)
#include <hip/hip_runtime.h>
#elif defined(CUDA_V)
#include <cuda_runtime.h>
#endif
using ParallelizationHelper = GpuParallelizationHelper;


#if defined(USE_FAST_COPY)
#define DEFAULT_FAST_COPY true
#else
#define DEFAULT_FAST_COPY false
#endif

class CpuRestMeasurement {
   // Queue callback data struct
   struct queue_callback_data {
      queue_callback_data(CpuRestMeasurement* m, std::size_t s) : me(m), step(s) {
      }

      CpuRestMeasurement* me;
      std::size_t step;
   };

   // Queue callback
   static void queue_callback(void* data);

   //ext_emomM, ext_emom, ext_mmom, ext_beff, ext_mstep

   // Temporary device storage vectors
   GpuTensor<real, 3> tmp_emomM;
   GpuTensor<real, 3> tmp_emom;
   GpuTensor<real, 2> tmp_mmom;
   GpuTensor<real, 3> tmp_beff;


   // Temporary host storage (pinned memory)
   Tensor<real, 3> pinned_emomM;
   Tensor<real, 3> pinned_emom;
   Tensor<real, 2> pinned_mmom;
   Tensor<real, 3> pinned_beff;

   // Vectors to copy
   const GpuTensor<real, 3>& emomM;
   const GpuTensor<real, 3>& emom;
   const GpuTensor<real, 2>& mmom;
   const GpuTensor<real, 3>& beff;

   Tensor<real, 3>& fortran_emomM;
   Tensor<real, 3>& fortran_emom;
   Tensor<real, 2>& fortran_mmom;
   Tensor<real, 3>& fortran_beff;

   // Event stack
   GpuEventPool eventPool;

   // Measure queue
   MeasurementQueue &measurementQueue;

   // Timer
   StopwatchDeviceSync stopwatch;

   // Parallelization helper
   ParallelizationHelper& parallel;

   // Control flags
   bool alwaysCopy;
   bool fastCopy;

   // Helpers
   void queueMeasurement(std::size_t mstep);
   void copyQueueFast(std::size_t mstep);
   void copyQueueSlow(std::size_t mstep);

public:
   // TODO add flag for fast_copy
   CpuRestMeasurement(const GpuTensor<real, 3>& emomM, const GpuTensor<real, 3>& emom,
                   const GpuTensor<real, 2>& mmom, const GpuTensor<real, 3>& beff, Tensor<real, 3>& f_emomM, Tensor<real, 3>& f_emom,
                   Tensor<real, 2>& f_mmom, Tensor<real, 3>& f_beff, MeasurementQueue& mq, bool fastCopy = DEFAULT_FAST_COPY,
                   bool alwaysCopy = false);
   ~CpuRestMeasurement(); 

   // Access methods
   void measure(std::size_t mstep);
   void flushMeasurements(std::size_t mstep);
};
