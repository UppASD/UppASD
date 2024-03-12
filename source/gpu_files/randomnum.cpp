#include "randomnum.hpp"

#include <stddef.h>
#include <stdio.h>
#include <time.h>

#include "c_helper.h"
#include "real_type.h"

RandomNumbers::RandomNumbers() {
}

RandomNumbers::~RandomNumbers() {
}

void RandomNumbers::fillArray(real* v, std::size_t len) {
   fill_rngarray(v, len);
}

