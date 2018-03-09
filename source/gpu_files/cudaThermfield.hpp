#ifndef __CUDA_THERMFIELD_HPP__
#define __CUDA_THERMFIELD_HPP__

#include "real_type.h"

#include <cuda.h>
#include <curand.h>

#include "cudaParallelizationHelper.hpp"

#include "fortMatrix.hpp"
#include "cudaMatrix.hpp"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"

// CUDA thermfield class
// Optimized for simulations with constant temperatures and timestep.

class CudaThermfield {
private:
	// Generator
	curandGenerator_t gen;

	// Initiation flags
	bool dataInitiated;
	bool constantsInitiated;

	// Data
	cudaMatrix<real,3,3> field;
	cudaMatrix<real,1>   sigmaFactor; // = sqrt(Dp*temperature(i))

	// Timer
	StopwatchDeviceSync stopwatch;

	// Parallelization helper
	CudaParallelizationHelper &parallel;

public:
	// Parallelization helpers
	class SetupSigmaFactor;
	class SetupField;

	// Constructor / destructor
	CudaThermfield();
	~CudaThermfield();

	// Initiate
	bool initiate(size_t N, size_t M, 
		curandRngType_t rngType = CURAND_RNG_PSEUDO_DEFAULT,
		unsigned long long seed = 0);
	bool initiateConstants(const fortMatrix<real,1> &temperature, 
		real timestep, real gamma, real k_bolt, real mub,
		real damping);

	// Initiated?
	inline bool initiated() {return dataInitiated && constantsInitiated;}

	// Get field
	inline const cudaMatrix<real,3,3> &getField() {return field;}

	// Randomize
	void randomize(const cudaMatrix<real,2> &mmom);
};

#endif

