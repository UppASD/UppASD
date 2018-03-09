#include <time.h>
#include <stddef.h>
#include <stdio.h>

#include "real_type.h"

#include "c_helper.h"

#include "randomnum.hpp"

RandomNumbers::RandomNumbers() {}
RandomNumbers::~RandomNumbers() {}

void RandomNumbers::fillArray(real * v, size_t len) {
	fill_rngarray(v, len);
}

