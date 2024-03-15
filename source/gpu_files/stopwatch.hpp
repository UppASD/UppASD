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

private:
   // Timepoint struct
   typedef struct {
      std::string name;
      float time;
   } timepoint;

   // Members
   std::list<timepoint> time_list;
   struct timeval start;
   struct timeval last;

public:
#ifdef DUMMY_STOPWATCH
   // Dummy functions
   Stopwatch() {
   }

   void startPoint() {
   }

   void skip() {
   }

   void add(const char *name) {
   }

   void add(const char *name, usd_int len) {
   }

   void add(const std::string &name) {
   }

   void reset() {
   }

   void print() {
   }

   bool empty() {
      return true;
   }
#else

   // Constructor
   Stopwatch() {
      reset();
   }

   // Use to force start a new timing event (may cause the
   // percentage to sum together to less than 100%)
   void startPoint() {
      gettimeofday(&last, 0);
   }

   // Adds the time since the last timing event to a category that
   // will not be included in the results (wont affect percentage)
   void skip() {
      add("-");
   }

   // Add a new time point or add time to the last point with the
   // same name. The time added is the time elapsed since the last
   // called add / startPoint / reset.
   void add(const char* name) {
      add(std::string(name));
   }

   void add(const char* name, usd_int len) {
      add(std::string(name, len));
   }

   void add(const std::string& name) {
      struct timeval now;
      gettimeofday(&now, 0);
      float time = diff(&last, &now);
      last = now;

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
      time_list.push_back(p);
   }

   // Reset the Stopwatch object
   void reset() {
      time_list.clear();
      gettimeofday(&start, 0);
      gettimeofday(&last, 0);
   }

   // Prints the result
   // TODO: ask is this should be kept
   void print() {
      print(2, 0);
   }

   // Empty stopwatch?
   bool empty() {
      if(time_list.empty()) {
         return true;
      }
      if(time_list.size() == 1) {
         return (time_list.front().name.compare("-") == 0);
      }
      return false;
   }

   // Helpers
private:
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
         fputc(' ', stdout);
      }
   }

   // Prints the result
   void print(int ind, int minlen) {
      std::list<timepoint>::iterator it;
      float total = diff(&start, &last);

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

      // Print the total time
      indent(2);
      std::printf("%-*s  %9.3f ms\n", minlen, "Total", total);

      //		float total_percent = 0.0f;

      for(it = time_list.begin(); it != time_list.end(); it++) {
         if(exclude_name.compare(it->name) != 0) {
            float time = it->time;
            float percent = (time / total) * 100.0f;
            //			total_percent += percent;
            indent(2);
            std::printf("%-*s  %9.3f ms (%5.2f%%)\n", minlen, it->name.c_str(), time, percent);
         }
      }
      std::printf("\n");
   }
#endif  // #ifndef DUMMY_STOPWATCH
};

