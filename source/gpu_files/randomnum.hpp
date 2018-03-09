#ifndef __RANDOMNUM_HPP__
#define __RANDOMNUM_HPP__

#include "real_type.h"

class RandomNumbers {
public:
	RandomNumbers();
	~RandomNumbers();

	void fillArray(real * v, size_t len);
};


#endif

