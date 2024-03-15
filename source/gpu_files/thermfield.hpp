#pragma once

#include "c_headers.hpp"
#include "fortMatrix.hpp"
#include "randomnum.hpp"
#include "real_type.h"
#include "stopwatch.hpp"

class Thermfield {
private:
   bool dataInitiated;
   bool constantsInitiated;

   // Data
   fortMatrix<real, 3, 3> field;
   fortMatrix<real, 1> sigmaFactor;  // = sqrt(Dp*temperature(i))

   // RNG
   RandomNumbers rng;

   Stopwatch &stopwatch;

public:
   // Constructor / destructor
   Thermfield();
   ~Thermfield();

   // Initiate
   bool initiate(std::size_t N, std::size_t M);
   bool initiateConstants(const matrix<real, 1>& temperature, real timestep, real gamma, real k_bolt,
                          real mub, real damping);

   // Initiated?
   inline bool initiated() {
      return dataInitiated && constantsInitiated;
   }

   // Get field
   inline const fortMatrix<real, 3, 3> &getField() {
      return field;
   }

   // Randomize
   void randomize(const matrix<real, 2> &mmom);
};

