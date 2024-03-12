#include <stdio.h>
#include <string.h>

#include <string>
#include <vector>

using namespace std;

#include "stopwatch.hpp"
#include "stopwatchPool.hpp"

// Reset a timer
extern "C" void stopwatch_reset_(const char* category, size_t strlen) {
   GlobalStopwatchPool::get(category, strlen).reset();
}

// Exclude the time since the last add from this category
extern "C" void stopwatch_skip_(const char* category, size_t strlen) {
   GlobalStopwatchPool::get(category, strlen).skip();
}

// extern "C" void stopwatch_add_(const char * category, const char * event, size_t strlen, size_t strlen2);
extern "C" void stopwatch_add_(const char* category, size_t strlen, const char* event, size_t strlen2) {
   // The ordering of the parameters may differ between compilers.
   // If strlen > 255 it is probably not the length but the event string pointer
   // as FORTRAN makes the function call as:
   // void stopwatch_add_(const char * category, size_t strlen, const char * event, size_t strlen2)
   if(strlen > 0xff) {
      size_t tmp = (size_t)event;
      event = (const char*)strlen;
      strlen = tmp;
   }
   // printf("add (%s, %ld, %s, %ld)\n", category, strlen, event, strlen2);
   GlobalStopwatchPool::get(category, strlen).add(event, strlen2);
}

// Print category
extern "C" void stopwatch_print_(const char* category, size_t strlen) {
   GlobalStopwatchPool::get(category, strlen).print();
}

