#pragma once

#define _dpr \
   { printDebug(__LINE__, __FILE__); }

#include <stdio.h>

inline void printDebug(const int line, const char* file) {
   printf("Reached: %d in %s\n", line, file);
}
