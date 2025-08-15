/*
 * Cuda event pool class
 *  Niklas Fejes 2013
 */

#include <cuda.h>

#include "c_headers.hpp"
#include <vector>

#include "cudaEventPool.hpp"

#ifdef NVPROF
#include <nvToolsExtCuda.h>
#endif


// Event class methods
CudaEventPool::Event::Event() {
   cudaEventCreate(&_event);
   active = true;
}


cudaEvent_t CudaEventPool::Event::event() {
   return _event;
}


void CudaEventPool::Event::deactivate() {
   active = false;
}


void CudaEventPool::Event::deactivate_callback(cudaStream_t, cudaError_t, void *e) {
#ifdef NVPROF
   nvtxRangePush("deactivate_callback");
#endif
   ((Event *)e)->active = false;
#ifdef NVPROF
   nvtxRangePop();
#endif
}


void CudaEventPool::Event::addDeactivateCallback(cudaStream_t stream) {
   cudaStreamAddCallback(stream, deactivate_callback, this, 0);
}


// Pool class methods
CudaEventPool::Event &CudaEventPool::get() {
   std::vector<Event *>::iterator it;
   for(it = stack.begin(); it != stack.end(); it++) {
      Event &e = (**it);
      if(e.active == false) {
         e.active = true;
         return e;
      }
   }

   stack.push_back(new Event());
   return *stack.back();
}


// Destroy all events in pool when done
CudaEventPool::~CudaEventPool() {
   //	std::printf("Event stack size: %ld\n", stack.size());
   std::vector<Event *>::iterator it;
   for(it = stack.begin(); it != stack.end(); it++) {
      delete *it;
   }
}


