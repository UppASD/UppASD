/*
 * Cuda event pool class
 *  Niklas Fejes 2013
 */

#include "c_headers.hpp"
#include <vector>

#include "gpuEventPool.hpp"
#include "gpu_wrappers.h"
#if defined(HIP_V)
#include <hip/hip_runtime.h>
#elif defined(CUDA_V)
#include <cuda.h>
#endif

#ifdef NVPROF
#include <nvToolsExtCuda.h>
#endif

// Event class methods
GpuEventPool::Event::Event() {
   GPU_EVENT_CREATE(&_event);
   active = true;
}

GPU_EVENT_T GpuEventPool::Event::event() {
   return _event;
}

void GpuEventPool::Event::deactivate() {
   active.store(false, std::memory_order_release);
}

void GpuEventPool::Event::deactivate_callback(void *e) {
#ifdef NVPROF
   nvtxRangePush("deactivate_callback");
#endif
   ((Event *)e)->active.store(false, std::memory_order_release);
#ifdef NVPROF
   nvtxRangePop();
#endif
}

void GpuEventPool::Event::addDeactivateCallback(GPU_STREAM_T stream) {
   GPU_STREAM_ADD_CALLBACK(stream, deactivate_callback, this);
}

// Pool class methods
GpuEventPool::Event &GpuEventPool::get() {
   std::vector<Event *>::iterator it;
   for(it = stack.begin(); it != stack.end(); it++) {
      Event &e = (**it);
      if(e.active.load(std::memory_order_acquire) == false) {
         e.active.store(true, std::memory_order_relaxed);
         return e;
      }
   }

   stack.push_back(new Event());
   return *stack.back();
}

// Destroy all events in pool when done
GpuEventPool::~GpuEventPool() {
   //	std::printf("Event stack size: %ld\n", stack.size());
   std::vector<Event *>::iterator it;
   for(it = stack.begin(); it != stack.end(); it++) {
      delete *it;
   }
}
