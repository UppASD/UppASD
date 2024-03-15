// Stopwatch Device Sync class
// Niklas Fejes 2012-2013

#pragma once

#include <cuda_runtime.h>

#include "c_headers.hpp"
#include "real_type.h"
#include "stopwatch.hpp"

class StopwatchDeviceSync {
#if defined(DUMMY_STOPWATCH) || defined(ASYNC_STOPWATCH)
   inline void sync() {
   }
#else
   inline void sync() {
      cudaDeviceSynchronize();
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

   void add(const char *name, usd_int len) {
      sync();
      parent.add(name, len);
   }

   void add(const std::string &name) {
      sync();
      parent.add(name);
   }
};

