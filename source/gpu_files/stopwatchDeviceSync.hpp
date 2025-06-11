// Stopwatch Device Sync class
// Niklas Fejes 2012-2013

#pragma once

#include <cuda_runtime.h>

#include "c_headers.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "gpu_wrappers.h"
#if defined(HIP_V)
#include <hip/hip_runtime.h>
#elif defined(CUDA_V)
#include <cuda_runtime.h>
#endif

class StopwatchDeviceSync {
#if defined(DUMMY_STOPWATCH) || defined(ASYNC_STOPWATCH)
   inline void sync() {
   }
#else
   inline void sync() {
      GPU_DEVICE_SYNCHRONIZE();
   }
#endif
public:
   // Parent stopwatch
   Stopwatch &parent;

   // Constructor
   StopwatchDeviceSync(Stopwatch &p) : parent(p) {
   }

   // Wrappers
   void startPoint() {
      sync();
      parent.startPoint();
   }

   void skip() {
      sync();
      parent.skip();
   }

   void add(const char *name) {
      sync();
      parent.add(name);
   }

   void add(const char *name, std::size_t len) {
      sync();
      parent.add(name, len);
   }

   void add(const std::string &name) {
      sync();
      parent.add(name);
   }
};

