/*
 * Cuda event pool class
 *  Niklas Fejes 2013
 */
#pragma once

#include <atomic>
#include <vector>
#include "gpu_wrappers.h"
#if defined(HIP_V)
#include <hip/hip_runtime.h>
#elif defined(CUDA_V)
#include <cuda_runtime.h>
#endif


class GpuEventPool {
public:
   // Event wrapper class
   class Event {
      friend GpuEventPool;
      // Cleared from the GPU host-callback thread, read/set from the host
      // thread that calls get(), so it must be atomic.
      std::atomic<bool> active;
      GPU_EVENT_T _event;
      static void deactivate_callback(void *e);

   public:
      Event();
      GPU_EVENT_T event();
      void deactivate();
      void addDeactivateCallback(GPU_STREAM_T stream);
   };

private:
   std::vector<Event *> stack;

public:
   // Get an activated event from the pool
   Event &get();

   // Destructor
   ~GpuEventPool();
};

