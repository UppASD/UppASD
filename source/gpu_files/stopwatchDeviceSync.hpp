// Host wall-clock and asynchronous device-event stopwatch.
#pragma once

#include <string>
#include <vector>

#include "gpu_wrappers.h"
#include "stopwatch.hpp"

class StopwatchDeviceSync {
   struct PendingSample {
      std::string name;
      GPU_EVENT_T start;
      GPU_EVENT_T end;
   };

   struct StreamState {
      GPU_STREAM_T stream;
      GPU_EVENT_T start;
      bool active;
   };

   GPU_STREAM_T defaultStream;
   std::vector<PendingSample> pending;
   std::vector<StreamState> streams;
   std::vector<GPU_EVENT_T> freeEvents;

   GPU_EVENT_T acquireEvent() {
      if(!freeEvents.empty()) {
         GPU_EVENT_T event = freeEvents.back();
         freeEvents.pop_back();
         return event;
      }
      GPU_EVENT_T event;
      GPU_EVENT_CREATE(&event);
      return event;
   }

   void releaseEvent(GPU_EVENT_T event) {
      freeEvents.push_back(event);
   }

   StreamState& state(GPU_STREAM_T stream) {
      for(auto& entry : streams) {
         if(entry.stream == stream) return entry;
      }
      streams.push_back({stream, GPU_EVENT_T(), false});
      return streams.back();
   }

   void begin(GPU_STREAM_T stream) {
      if(!Stopwatch::timingEnabled()) return;
      StreamState& entry = state(stream);
      // A coarse timer spans the complete lifetime of the stopwatch on each
      // stream, so intermediate phase boundaries must not reset it.
      if(entry.active && Stopwatch::coarseTiming()) return;
      if(entry.active) releaseEvent(entry.start);
      entry.start = acquireEvent();
      GPU_EVENT_RECORD(entry.start, stream);
      entry.active = true;
   }

   void collect(bool wait) {
      for(auto it = pending.begin(); it != pending.end();) {
         const GPU_ERROR_T status = wait ? GPU_EVENT_SYNCHRONIZE(it->end) : GPU_EVENT_QUERY(it->end);
         if(status != GPU_SUCCESS) {
            ++it;
            continue;
         }
         float elapsed = 0.0f;
         if(GPU_EVENT_ELAPSED_TIME(&elapsed, it->start, it->end) == GPU_SUCCESS) {
            parent.addGpu(it->name, elapsed);
         }
         releaseEvent(it->start);
         releaseEvent(it->end);
         it = pending.erase(it);
      }
   }

   void addGpu(const std::string& name, GPU_STREAM_T stream) {
      if(!Stopwatch::timingEnabled()) return;
      if(Stopwatch::coarseTiming()) return;
      StreamState& entry = state(stream);
      if(!entry.active) begin(stream);
      GPU_EVENT_T end = acquireEvent();
      GPU_EVENT_RECORD(end, stream);
      pending.push_back({name, entry.start, end});
      entry.start = acquireEvent();
      GPU_EVENT_RECORD(entry.start, stream);
      entry.active = true;
      collect(false);
   }

 public:
   Stopwatch& parent;

   explicit StopwatchDeviceSync(Stopwatch& p, GPU_STREAM_T stream = 0)
      : defaultStream(stream), parent(p) {
      if(Stopwatch::coarseTiming()) {
         parent.startPoint();
         begin(defaultStream);
      }
   }

   ~StopwatchDeviceSync() {
      if(!Stopwatch::timingEnabled()) return;
      if(Stopwatch::coarseTiming()) {
         for(auto& entry : streams) {
            if(!entry.active) continue;
            GPU_EVENT_T end = acquireEvent();
            GPU_EVENT_RECORD(end, entry.stream);
            if(GPU_EVENT_SYNCHRONIZE(end) == GPU_SUCCESS) {
               float elapsed = 0.0f;
               if(GPU_EVENT_ELAPSED_TIME(&elapsed, entry.start, end) == GPU_SUCCESS) {
                  parent.addGpuTotal(elapsed);
               }
            }
            GPU_EVENT_DESTROY(entry.start);
            releaseEvent(end);
         }
         for(const auto& event : freeEvents) GPU_EVENT_DESTROY(event);
         return;
      }
      collect(true);
      for(auto& entry : streams) {
         if(entry.active) GPU_EVENT_DESTROY(entry.start);
      }
      for(const auto& event : freeEvents) GPU_EVENT_DESTROY(event);
   }

   void startPoint() { startPoint(defaultStream); }
   void startPoint(GPU_STREAM_T stream) {
      parent.startPoint();
      begin(stream);
   }

   void skip() { skip(defaultStream); }
   void skip(GPU_STREAM_T stream) {
      parent.skip();
      begin(stream);
   }

   void add(const char* name) { add(std::string(name), defaultStream); }
   void add(const char* name, std::size_t len) { add(std::string(name, len), defaultStream); }
   void add(const std::string& name) { add(name, defaultStream); }
   void add(const char* name, GPU_STREAM_T stream) { add(std::string(name), stream); }
   void add(const std::string& name, GPU_STREAM_T stream) {
      parent.add(name);
      addGpu(name, stream);
   }
};
