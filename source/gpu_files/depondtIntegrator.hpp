#pragma once

#include "hostMatrix.hpp"
#include "matrix.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "thermfield.hpp"

class DepondtIntegrator {
private:
   // System parameters
   usd_int Natom;
   usd_int Mensemble;
   real gamma;
   real damping;
   real timestep;

   // Integrator parameters
   char stt;  // Spin transfer torque

   // Class local matrices
   // <real,3,3> -- real precision, 3 dimensions, first dim always 3
   fortMatrix<real, 3, 3> mrod;    // Rotated magnetic moment
   fortMatrix<real, 3, 3> blocal;  // Local effective field
   fortMatrix<real, 3, 3> bdup;    // Resulting effective field

   // Thermfield
   Thermfield tfield;

   // Timer
   Stopwatch& stopwatch;

   // Algorithm
   bool rotate(const hostMatrix<real, 3, 3>& emom, real delta_t);
   void thermfield(hostMatrix<real, 2>& mmom, real deltat, const hostMatrix<real, 1>& temperature);
   void buildbeff(const hostMatrix<real, 3, 3>& emom, const hostMatrix<real, 3, 3>& btorque);

public:
   // Constructor
   DepondtIntegrator();

   // Destructor
   ~DepondtIntegrator();

   // Initiator
   bool initiate(usd_int Natom, usd_int Mensemble, char stt);

   // Set up constants
   bool initiateConstants(real gamma_const, real k_bolt_const, real mub_const, real damping_const,
                          const hostMatrix<real, 1>& temp, real timestep);

   void setDamping(real d);

   // Releaser
   void release();

   // Algorithm
   void evolveFirst(const hostMatrix<real, 3, 3>& beff, hostMatrix<real, 3, 3>& b2eff,
                    const hostMatrix<real, 3, 3>& btorque, hostMatrix<real, 3, 3>& emom,
                    hostMatrix<real, 3, 3>& emom2, hostMatrix<real, 3, 3>& emomM,
                    const hostMatrix<real, 2>& mmom);

   void evolveSecond(const hostMatrix<real, 3, 3>& beff, const hostMatrix<real, 3, 3>& b2eff,
                     const hostMatrix<real, 3, 3>& btorque, hostMatrix<real, 3, 3>& emom,
                     hostMatrix<real, 3, 3>& emom2);
};

