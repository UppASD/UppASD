#pragma once

#include "real_type.h"
#include "c_headers.hpp"

class RandomNumbers {
public:
   RandomNumbers();
   ~RandomNumbers();

   void fillArray(real* v, usd_int len);
};


