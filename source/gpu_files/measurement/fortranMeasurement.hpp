#pragma once

#include "c_headers.hpp"
#include "measurable.hpp"
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

class FortranMeasurement : public Measurable {
   // Queue callback data struct
   struct queue_callback_data {
      queue_callback_data(FortranMeasurement* m, std::size_t s) : me(m), step(s) {
      }

      FortranMeasurement* me;
      std::size_t step;
   };

   // Queue callback
   static void queue_callback(GPU_STREAM_T, GPU_ERROR_T, void* data);

   // Temporary device storage vectors
   GpuTensor<real, 3> tmp_emomM;
   GpuTensor<real, 3> tmp_emom;
   GpuTensor<real, 2> tmp_mmom;

   // Temporary host storage (pinned memory)
   Tensor<real, 3> pinned_emomM;
   Tensor<real, 3> pinned_emom;
   Tensor<real, 2> pinned_mmom;

   // Vectors to copy
   const GpuTensor<real, 3>& emomM;
   const GpuTensor<real, 3>& emom;
   const GpuTensor<real, 2>& mmom;

   Tensor<real, 3>& fortran_emomM;
   Tensor<real, 3>& fortran_emom;
   Tensor<real, 2>& fortran_mmom;

   // Event stack
   GpuEventPool eventPool;

   // Measure queue
   MeasurementQueue measurementQueue;

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
   FortranMeasurement(const GpuTensor<real, 3>& emomM, const GpuTensor<real, 3>& emom,
                   const GpuTensor<real, 2>& mmom, Tensor<real, 3>& f_emomM, Tensor<real, 3>& f_emom,
                   Tensor<real, 2>& f_mmom, bool fastCopy = DEFAULT_FAST_COPY,
                   bool alwaysCopy = false);
   ~FortranMeasurement() override;

   // Access methods
   void measure(std::size_t mstep) override;
   void flushMeasurements(std::size_t mstep) override;
};

