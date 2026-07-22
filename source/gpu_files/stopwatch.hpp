// Stopwatch class
// Niklas Fejes 2012-2013

#pragma once

#include <sys/time.h>

#include <ctime>
#include <list>
#include <string>

#include "c_headers.hpp"
#include "real_type.h"

class Stopwatch {
   // Allow StopwatchPool to access private methods
   friend class StopwatchPool;
   friend class StopwatchDeviceSync;

private:
   // Timepoint struct
   typedef struct {
      std::string name;
      float time;
      float gpu_time;
   } timepoint;

   // Members
   std::list<timepoint> time_list;
   struct timeval start;
   struct timeval last;
   bool initialized = false;
   bool coarse_timing_recorded = false;
   float gpu_total = 0.0f;

public:
   // N disables timing, C records totals only, and Y records all categories.
   static void setTimingMode(char mode) {
      timingModeStorage() = mode == 'Y' ? TimingMode::Detailed :
                            mode == 'C' ? TimingMode::Coarse : TimingMode::Disabled;
   }

   static bool timingEnabled() {
      return timingModeStorage() != TimingMode::Disabled;
   }

   static bool detailedTiming() {
      return timingModeStorage() == TimingMode::Detailed;
   }

   static bool coarseTiming() {
      return timingModeStorage() == TimingMode::Coarse;
   }

   Stopwatch() {
      reset();
   }

   void startPoint() {
      if(!timingEnabled()) return;
      gettimeofday(&last, 0);
      if(!initialized) {
         start = last;
         initialized = true;
      }
   }

   void skip() {
      if(!timingEnabled()) return;
      add("-");
   }

   void add(const char *name) {
      if(!timingEnabled()) return;
      add(std::string(name));
   }

   void add(const char *name, std::size_t len) {
      if(!timingEnabled()) return;
      add(std::string(name, len));
   }

   void add(const std::string& name) {
      if(!timingEnabled()) return;
      struct timeval now;
      gettimeofday(&now, 0);
      if(!initialized) {
         start = now;
         last = now;
         initialized = true;
      }
      float time = diff(&last, &now);
      last = now;

      if(coarseTiming()) {
         if(name == "-") {
            addTime(name, time);
         } else {
            coarse_timing_recorded = true;
         }
         return;
      }

      addTime(name, time);
   }

   void addGpuTotal(float time) {
      if(timingEnabled()) {
         gpu_total += time;
         if(coarseTiming()) coarse_timing_recorded = true;
      }
   }

   void reset() {
      time_list.clear();
      initialized = false;
      coarse_timing_recorded = false;
      gpu_total = 0.0f;
      if(!timingEnabled()) return;
      gettimeofday(&start, 0);
      gettimeofday(&last, 0);
      initialized = true;
   }

   void print() {
      if(!timingEnabled()) return;
      print(2, 0);
   }

   bool empty() {
      if(!timingEnabled()) return true;
      if(coarseTiming()) return !coarse_timing_recorded;
      if(time_list.empty()) {
         return true;
      }
      if(time_list.size() == 1) {
         return (time_list.front().name.compare("-") == 0);
      }
      return false;
   }

private:
   void addTime(const std::string& name, float time) {

      // Check if the name is already in the list
      std::list<timepoint>::iterator it;
      for(it = time_list.begin(); it != time_list.end(); it++) {
         if(name.compare(it->name) == 0) {
            it->time += time;
            return;
         }
      }

      // Otherwise add new
      timepoint p;
      p.name = name;
      p.time = time;
      p.gpu_time = 0.0f;
      time_list.push_back(p);
   }

   // Add elapsed device time to an existing wall-clock category.  GPU event
   // samples are collected asynchronously by StopwatchDeviceSync.
   void addGpu(const std::string& name, float time) {
      if(!timingEnabled()) return;
      for(auto it = time_list.begin(); it != time_list.end(); ++it) {
         if(name == it->name) {
            it->gpu_time += time;
            return;
         }
      }
      timepoint p;
      p.name = name;
      p.time = 0.0f;
      p.gpu_time = time;
      time_list.push_back(p);
   }

   enum class TimingMode { Disabled, Coarse, Detailed };

   static TimingMode& timingModeStorage() {
      static TimingMode mode = TimingMode::Disabled;
      return mode;
   }

   inline float diff(const struct timeval* a, const struct timeval* b) {
      return (float)((1000.0 * (b->tv_sec - a->tv_sec)) + (0.001 * (b->tv_usec - a->tv_usec)));
   }

   // Get the length of the longest name
   int minNameLen() {
      std::list<timepoint>::iterator it;
      int minlen = 5;  // = strlen("Total");

      // Find the length of the longest name
      for(it = time_list.begin(); it != time_list.end(); it++) {
         if(minlen < it->name.size()) {
            minlen = it->name.size();
         }
      }
      return minlen;
   }

   void indent(int n) {
      while(n-- > 0) {
         std::fputc(' ', stdout);
      }
   }

   // Prints the result
   void print(int ind, int minlen) {
      std::list<timepoint>::iterator it;
      float total = diff(&start, &last);
      float total_gpu = gpu_total;

      // Minimum 5 chars (= strlen("Total"))
      if(minlen < 5) {
         minlen = 5;
      }

      // Find the length of the longest name
      for(it = time_list.begin(); it != time_list.end(); it++) {
         if(minlen < it->name.size()) {
            minlen = it->name.size();
         }
      }

      // Excluded time?
      const std::string exclude_name("-");
      for(it = time_list.begin(); it != time_list.end(); it++) {
         if(exclude_name.compare(it->name) == 0) {
            total -= it->time;
            break;
         }
      }

      for(it = time_list.begin(); it != time_list.end(); it++) {
         if(exclude_name.compare(it->name) != 0) total_gpu += it->gpu_time;
      }

      // Print the total time
      indent(2);
      std::printf("%-*s  %9.3f ms wall  %9.3f ms GPU\n", minlen, "Total", total, total_gpu);

      if(coarseTiming()) return;

      //		float total_percent = 0.0f;

      for(it = time_list.begin(); it != time_list.end(); it++) {
         if(exclude_name.compare(it->name) != 0) {
            float time = it->time;
            float percent = (time / total) * 100.0f;
            //			total_percent += percent;
            indent(2);
            std::printf("%-*s  %9.3f ms wall  %9.3f ms GPU  (%5.2f%%)\n", minlen, it->name.c_str(), time,
                        it->gpu_time, percent);
         }
      }
      std::printf("\n");
   }
};
