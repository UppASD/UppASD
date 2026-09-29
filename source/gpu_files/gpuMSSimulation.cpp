
#include "c_headers.hpp"
#include "c_helper.h"
#include "gpuDepondtIntegrator_MS.hpp"

#include "gpuHamiltonianCalculations.hpp"
#include "gpuMeasurement.hpp"
#include "gpuMomentUpdater.hpp"
#include "gpuSimulation.hpp"
#include "gpuStructures.hpp"
#include "gpuStructures.hpp"
#include "fortranData.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "stopwatchPool.hpp"
#include "tensor.hpp"
#include "gpuParallelizationHelper.hpp"
#include "measurementFactory.hpp"
#include "correlationFactory.hpp"
#include "measurementQueue.hpp"
#include "cpuRestMeasurement.hpp"
#include "gpuInterpolation.hpp"

#include "gpu_wrappers.h"
#include "gpuCorrelations.hpp"

using ParallelizationHelper = GpuParallelizationHelper;

GpuSimulation::GpuMSSimulation::GpuMSSimulation() {
   // isInitiatedSD = false;
}

GpuSimulation::GpuMSSimulation::~GpuMSSimulation() {
}

// Printing simulation status
// Added copy to fortran line so that simulation status is printed correctly > Thomas Nystrand 14/09/09
void GpuSimulation::GpuMSSimulation::printMdStatus(std::size_t mstep, GpuSimulation& gpuSim) {

}

// Multiscale initial phase
void GpuSimulation::GpuMSSimulation::MSiphase(GpuSimulation& gpuSim) {
   
}

// Multiscale measurement phase
void GpuSimulation::GpuMSSimulation::MSmphase(GpuSimulation& gpuSim) {

      // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);
   std::printf("GpuSDSimulation: SD measurement phase starting\n");

   // Initiated?
   if(!gpuSim.isInitiated) {
      std::fprintf(stderr, "GpuSimulation: not initiated!\n");
      return;
   }

   // Reload lattice state from Fortran to guarantee phase handover continuity,
   // including all moment arrays that may have been updated in SDiphase.
   gpuSim.copyFromFortran();

   // Make phase boundary explicit: start SDmphase with emom2 synchronized to emom.
   // This removes any dependence on historical emom2 content from prior phase bookkeeping.
   gpuSim.gpuLattice.emom2.copy_sync(gpuSim.gpuLattice.emom);
   gpuSim.copyToFortran();

   // Timer
   StopwatchDeviceSync stopwatch = StopwatchDeviceSync(GlobalStopwatchPool::get("GPU measurement phase"));

   // Initiate default parallelization helper
    ParallelizationHelperInstance.initiate(gpuSim.SimParam.N, gpuSim.SimParam.M, gpuSim.SimParam.NH);
   // Depontd integrator
   GpuDepondtIntegrator_MS integrator;

   // Hamiltonian calculations
   GpuHamiltonianCalculations hamCalc;

   // Moment updater
   GpuMomentUpdater momUpdater(gpuSim.gpuLattice, gpuSim.SimParam.mompar, gpuSim.SimParam.initexc);
   //Queue
   MeasurementQueue mqueue;
   // Measurement
   const auto measurement = MeasurementFactory::create(gpuSim.gpuLattice, gpuSim.cpuLattice, gpuSim.gpuEnergies, mqueue, gpuSim.Flags.do_jtensor);
   //CPU residing measurements
   //CpuRestMeasurement cpuMeas(gpuSim.gpuLattice.emomM, gpuSim.gpuLattice.emom, gpuSim.gpuLattice.mmom, 
   //                gpuSim.gpuLattice.beff, gpuSim.cpuLattice.emomM, gpuSim.cpuLattice.emom,
    //               gpuSim.cpuLattice.mmom, gpuSim.cpuLattice.beff, mqueue);
   //Corrrelations
   const auto correlation = CorrelationFactory::create(gpuSim.gpuLattice, gpuSim.cpuLattice, 
            gpuSim.Flags, gpuSim.SimParam, gpuSim.cpuCorrelations, mqueue);

   // Initiate integrator and Hamiltonian
   if(!integrator.initiate(gpuSim.SimParam)) {  // TODO
      std::fprintf(stderr, "GpuSDSimulation: integrator failed to initiate!\n");
      return;
   }

   if(!hamCalc.initiate(gpuSim.Flags, gpuSim.SimParam, gpuSim.gpuHamiltonian)) {  // TODO
      std::fprintf(stderr, "GpuSDSimulation: Hamiltonian failed to initiate!\n");
      return;
   }

   GpuInterpolation interpolation(gpuSim.SimParam.N, gpuSim.SimParam.M, gpuSim.gpuInterpolationInfo, 
                    gpuSim.gpuLattice);

   int mnn = gpuSim.cpuHamiltonian.j_tensor.extent(2);
   int l = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   std::printf("_______________________________________________\n");

   // Initiate constants for integrator
   integrator.initiateConstants(gpuSim.SimParam, gpuSim.cpuLattice.temperature);

   // Timing
   stopwatch.add("initiate");

   size_t nstep = gpuSim.SimParam.nstep;
   size_t rstep = gpuSim.SimParam.rstep;

   bool measure_ene;

   // Time step loop
   for(std::size_t mstep = rstep + 1; mstep <= rstep + nstep; mstep++) {
      // Measure
      measurement->measure(mstep);
      correlation->measure(mstep);

      stopwatch.add("measurement");

      // Print simulation status for each 5% of the simulation length
      printMdStatus(mstep, gpuSim);

      // Apply Hamiltonian to obtain effective field
      hamCalc.heisge(gpuSim.gpuLattice, gpuSim.gpuEnergies, false);
      stopwatch.add("hamiltonian");

      // Perform first step of SDE solver
      integrator.evolveFirst(gpuSim.gpuLattice); //TODO
      stopwatch.add("evolution");

      measure_ene = ((gpuSim.Flags.do_ene > 0 ) && (gpuSim.Flags.do_gpu_measurements)&&
            (((mstep-1)%gpuSim.SimParam.ene_step == 0)||((gpuSim.Flags.do_cumu)&&((mstep-1)%gpuSim.SimParam.cumu_step == 0))));

      // Apply Hamiltonian to obtain effective field
      hamCalc.heisge(gpuSim.gpuLattice, gpuSim.gpuEnergies, measure_ene);
      stopwatch.add("hamiltonian");
  

      // Perform second (corrector) step of SDE solver
      integrator.evolveSecond(gpuSim.gpuLattice); //TODO
      stopwatch.add("evolution");
      // Update magnetic moments after time evolution step
      momUpdater.update();
      stopwatch.add("moments");

      measurement->updateAC(mstep);

      // Check for error
      GPU_ERROR_T e = GPU_GET_LAST_ERROR();
      if(e != GPU_SUCCESS) {
         std::printf("Uncaught GPU error %d: %s\n", e, GPU_GET_ERROR_STRING(e));
         GPU_DEVICE_RESET();
         std::exit(EXIT_FAILURE);
      }    real cv{};            // Specific heat


   }  // End loop over simulation steps

   // Final measure and print remaining measurements to file

   measure_ene = ((gpuSim.Flags.do_ene > 0 ) && (gpuSim.Flags.do_gpu_measurements));
   hamCalc.heisge(gpuSim.gpuLattice, gpuSim.gpuEnergies, measure_ene);

   measurement->measure(rstep + nstep + 1);    
   correlation->measure(rstep + nstep + 1);  // TODO
   stopwatch.add("measurement");

   mqueue.finish();

   // Print remaining measurements
   measurement->flushMeasurements(rstep + nstep + 1);  // TODO
   correlation->flushCorrelations(gpuSim.cpuCorrelations, rstep + nstep + 1); 
   
   // Transfer GPU sample count back to Fortran for averaging
   if (FortranData::sc_nsamp_ptr != nullptr) {
       *FortranData::sc_nsamp_ptr = gpuSim.cpuCorrelations.sc_nsamp;
   }
      if (FortranData::sc_tidx_ptr != nullptr) {
         *FortranData::sc_tidx_ptr = gpuSim.cpuCorrelations.sc_tidx;
      }
   
   stopwatch.add("flush measurement");


   // Synchronize with device
   GPU_DEVICE_SYNCHRONIZE();
   stopwatch.add("final synchronize");
  
}

