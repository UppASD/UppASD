#pragma once

#include "real_type.h"
#include <cstdio>

class RandomNumbers {
public:
   RandomNumbers();
   ~RandomNumbers();

   void fillArray(real* v, std::size_t len);
};


