#pragma once

#include "c_headers.hpp"

#include "tensor.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "gpuStructures.hpp"
#include "gpuParallelizationHelper.hpp"
#include "gpuThermfield.hpp"
using ParallelizationHelper = GpuParallelizationHelper;

class GpuDepondtIntegrator {
private:
   // System parameters
   real gamma;
   real k_bolt;
   real mub;
   real damping;
   real dp_factor;

   // Integrator parameters
   char stt;  // Spin transfer torque mode
   real timestep;

   // Class local matrices
   GpuTensor<real, 3> mrod;
   GpuTensor<real, 3> blocal;
   GpuTensor<real, 3> bdup;

   // Thermfield
   GpuThermfield thermfield;

   // Timer
   StopwatchDeviceSync stopwatch;

   // Parallelization helper
   ParallelizationHelper& parallel;

   // Algorithm
   void rotate(const GpuTensor<real, 3>& emom, real delta_t);
   void buildbeff(const GpuTensor<real, 3>& emom, const GpuTensor<real, 3>& btorque);

public:
   // Parallelization helpers
   class Rotate;
   class BuildEffectiveField;
   class EvolveFirst;
   class CommitActivePredictor;
   class RestoreActiveInitial;
   class CommitActiveCorrector;

   // Constructor
   GpuDepondtIntegrator();

   // Destructor
   ~GpuDepondtIntegrator();

   // Initiator
   bool initiate(const SimulationParameters SimParam);

   // Projected device bytes (mrod/blocal/bdup + owned thermfield); mirror initiate().
   static std::size_t estimateBytes(const SimulationParameters& SimParam);

   // Set up constants
   bool initiateConstants(const SimulationParameters SimParam, const Tensor<real, 1>temperature);
   void resetConstants(const Tensor<real, 1> temperature, real phaseTimestep, real phaseDamping);

   // Releaser
   void release();

   // Algorithm
   void evolveFirst(deviceLattice& gpuLattice);

   // Execute the production Depondt predictor for a compact set of one-based
   // physical atom ids in device memory. Each id is applied in every ensemble;
   // its compact position never enters state or RNG indexing. The thermal generator
   // deliberately remains the production full-state generator, so a site's
   // draw is tied to its ordinary (site, ensemble) field location even when
   // this list is reordered or compacted.
   void evolveFirst(deviceLattice& gpuLattice, const int* oneBasedAtoms,
                    std::size_t activeAtomCount);

   void evolveSecond(deviceLattice& gpuLattice);

   // Correct the matching active-subset predictor.  Call with precisely the
   // same list and count supplied to evolveFirst().
   void evolveSecond(deviceLattice& gpuLattice, const int* oneBasedAtoms,
                     std::size_t activeAtomCount);
};
