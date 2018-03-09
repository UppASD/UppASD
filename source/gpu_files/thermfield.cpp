#include <cstddef>
#include <cmath>

using namespace std;

#include "real_type.h"

#include "thermfield.hpp"

#include "randomnum.hpp"
#include "matrix.hpp"
#include "fortMatrix.hpp"
#include "stopwatch.hpp"
#include "stopwatchPool.hpp"


// Constructor
Thermfield::Thermfield() : stopwatch(GlobalStopwatchPool::get("Thermfield")) {
}

// Destructor
Thermfield::~Thermfield() {
	delete field.get_data();
	delete sigmaFactor.get_data();

	field.set(NULL, 0, 0, 0);
	sigmaFactor.set(NULL, 0);
}

// Initiate
bool Thermfield::initiate(size_t N, size_t M) {
	field.set(new real[3*N*M], 3, N, M);
	sigmaFactor.set(new real[N], N);
	dataInitiated = true;
	return true;
}

// Initiate constants
bool Thermfield::initiateConstants(const matrix<real,1> &temperature, 
	real timestep, real gamma, real k_bolt, real mub,
	real damping) {

	// Initiated?
	if (!dataInitiated) return false;

	// Damping parameter
	real dp = (2.0*damping*k_bolt) / (timestep*gamma*mub*(1+damping*damping));

	// Size
	size_t N = temperature.dimension_size(0);

	// Set up sigmaFactor
	#pragma omp parallel for
	for (size_t i = 0; i < N; i++) {
		sigmaFactor(i) = sqrt(dp * temperature(i));
	}

	// Flag that we're initiated
	constantsInitiated = true;
	return true;
}

// Randomize
void Thermfield::randomize(const matrix<real,2> &mmom) {
	// Initiated?
	if (!initiated()) return;

	// Sizes
	size_t M = field.dimension_size(2);
	size_t N = field.dimension_size(1);

	// Timing
	stopwatch.skip();

	// Randomize
	rng.fillArray(field.get_data(), field.size());
	stopwatch.add("RNG");		

	// Expand
	#pragma omp parallel for collapse(2)
	for (size_t k = 0; k < M; k++) {
		for (size_t i = 0; i < N; i++) {
			real sigma = sigmaFactor(i) / sqrt(mmom(i,k));
			field(0,i,k) *= sigma;
			field(1,i,k) *= sigma;
			field(2,i,k) *= sigma;
		}
	}
	stopwatch.add("loop");	
}

