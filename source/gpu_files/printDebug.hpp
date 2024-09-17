#pragma once

#include "c_headers.hpp"

#define _dpr \
   { printDebug(__LINE__, __FILE__); }

inline void printDebug(const int line, const char* file) {
   std::printf("Reached: %d in %s\n", line, file);
}
