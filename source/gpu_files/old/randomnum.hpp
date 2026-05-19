#pragma once

#include "c_headers.hpp"
#include "real_type.h"

class RandomNumbers {
public:
   RandomNumbers();
   ~RandomNumbers();

   void fillArray(real* v, std::size_t len);
};

