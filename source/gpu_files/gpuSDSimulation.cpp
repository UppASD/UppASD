
#include "c_headers.hpp"
#include "c_helper.h"
#include "gpuDepondtIntegrator.hpp"
//#include "cudaGPUErrchk.hpp"
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

#include "gpu_wrappers.h"
#include "gpuCorrelations.hpp"

using ParallelizationHelper = GpuParallelizationHelper;

GpuSimulation::GpuSDSimulation::GpuSDSimulation() {
   // isInitiatedSD = false;
}

GpuSimulation::GpuSDSimulation::~GpuSDSimulation() {
}

// Printing simulation status
// Added copy to fortran line so that simulation status is printed correctly > Thomas Nystrand 14/09/09
void GpuSimulation::GpuSDSimulation::printMdStatus(std::size_t mstep, GpuSimulation& gpuSim) {
   if(gpuSim.SimParam.nstep > 20) {
      if(mstep % ((gpuSim.SimParam.rstep + gpuSim.SimParam.nstep) / 20) == 0) {
         gpuSim.copyToFortran();  // This is run so seldomly it has not impact on overall performance
         fortran_calc_simulation_status_variables(gpuSim.SimParam.mavg);
         std::printf("GPU: %3ld%% done. Mbar: %10.6f. U: %8.5f.\n",
                     mstep * 100 / (gpuSim.SimParam.rstep + gpuSim.SimParam.nstep),
                     *gpuSim.SimParam.mavg,
                     *gpuSim.SimParam.binderc);
      }
   } else {
      gpuSim.copyToFortran();
      fortran_calc_simulation_status_variables(gpuSim.SimParam.mavg);
      std::printf("GPU: Iteration %ld Mbar %13.6f\n", mstep, *gpuSim.SimParam.mavg);
   }
}

// Spin Dynamics initial phase
void GpuSimulation::GpuSDSimulation::SDiphase(GpuSimulation& gpuSim) {
   // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);
   std::printf("GpuSDSimulation: SD initial phase starting\n");

   // Initiated?
   if(!gpuSim.isInitiated) {
      std::fprintf(stderr, "GpuSimulation: not initiated!\n");
      return;
   }

   // Initiate default parallelization helper
   ParallelizationHelperInstance.initiate(gpuSim.SimParam.N, gpuSim.SimParam.M, gpuSim.SimParam.NH);
   // Timer
   StopwatchDeviceSync stopwatch(GlobalStopwatchPool::get("Gpu initial phase"),
                                 ParallelizationHelperInstance.getWorkStream());
   // Depontd integrator
   GpuDepondtIntegrator integrator;

   // Hamiltonian calculations
   GpuHamiltonianCalculations hamCalc;

   // Moment updater
   GpuMomentUpdater momUpdater(gpuSim.gpuLattice, gpuSim.SimParam.mompar, gpuSim.SimParam.initexc);


   // Initiate integrator and Hamiltonian
   if(!integrator.initiate(gpuSim.SimParam)) {
      std::fprintf(stderr, "GpuSDSimulation: integrator failed to initiate!\n");
      return;
   }
   if(!hamCalc.initiate(gpuSim.Flags, gpuSim.SimParam, gpuSim.gpuHamiltonian)) {
      std::fprintf(stderr, "GpuSDSimulation: Hamiltonian failed to initiate!\n");
      return;
   }

   // TEMPORARY PRINTING
   std::printf("\n");
   // std::printf("________DEBUG System Information:___________ \n");
   // std::printf("%zu\n", gpuSim.cpuHamiltonian.j_tensor.extent(0));
   // std::printf("%zu\n", gpuSim.cpuHamiltonian.j_tensor.extent(1));
   // std::printf("%zu\n", gpuSim.cpuHamiltonian.j_tensor.extent(2));
   // std::printf("%zu\n", gpuSim.cpuHamiltonian.j_tensor.extent(3));
   // std::printf("______________________________________\n");
   int mnn = gpuSim.cpuHamiltonian.j_tensor.extent(2);
   int l = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int N = gpuSim.SimParam.N;
   std::printf("_______________________________________________\n");
   Tensor<real, 1> ipTemp;
   ipTemp.AllocateHost(N);

   for(unsigned int k = 0; k < N; k++){
      ipTemp(k) = gpuSim.cpuLattice.ipTemp_array(k, 0);
   }
   // Initiate constants for integrator
   integrator.initiateConstants(gpuSim.SimParam, ipTemp);
   // Timing
   stopwatch.add("initiate");

   unsigned int steps;
   int ipnphase = gpuSim.SimParam.ipnphase; 
   for(unsigned int it = 0; it < ipnphase; it++){
       const real phaseDt = gpuSim.cpuLattice.ipdelta_t(it);
       const real phaseDamp = gpuSim.cpuLattice.iplambda1(it);
       gpuSim.SimParam.delta_t = phaseDt;
       gpuSim.SimParam.damping = phaseDamp;
       integrator.resetConstants(ipTemp, phaseDt, phaseDamp);
      steps = gpuSim.cpuLattice.ipnstep(it);

      // Time step loop
      for(std::size_t mstep = 1; mstep <= steps; mstep++) {
         // Print simulation status for each 5% of the simulation length
         printMdStatus(mstep, gpuSim); //-- Do we need it in initial phase?

         // Apply Hamiltonian to obtain effective field
         hamCalc.heisge(gpuSim.gpuLattice, gpuSim.gpuEnergies, false);
         stopwatch.add("hamiltonian");

         // Perform first step of SDE solver
         integrator.evolveFirst(gpuSim.gpuLattice);
         stopwatch.add("evolution");

         // Apply Hamiltonian to obtain effective field
         hamCalc.heisge(gpuSim.gpuLattice, gpuSim.gpuEnergies, false);
         stopwatch.add("hamiltonian");

         // Perform second (corrector) step of SDE solver
         integrator.evolveSecond(gpuSim.gpuLattice);
         stopwatch.add("evolution");

         // Update magnetic moments after time evolution step
         momUpdater.update();
         stopwatch.add("moments");

         // Check for error
         GPU_ERROR_T e = GPU_GET_LAST_ERROR();
         if(e != GPU_SUCCESS) {
            std::printf("Uncaught GPU error %d: %s\n", e, GPU_GET_ERROR_STRING(e));
            // Fatal GPU error: flush and _Exit. Do NOT device-reset then std::exit -
            // that runs the static GpuSimulation dtor, which cudaFrees on a torn-down
            // context and aborts with "driver shutting down". The OS reclaims the context.
            std::fflush(stdout);
            std::fflush(stderr);
            std::_Exit(EXIT_FAILURE);
         }
      }

      if(it < (ipnphase - 1)){
         for(unsigned int k = 0; k < N; k++){
            ipTemp(k) = gpuSim.cpuLattice.ipTemp_array(k, it + 1);
          }
         }
      }  
      
   // Synchronize with device
   GPU_STREAM_SYNC(ParallelizationHelperInstance.getWorkStream());
   // Explicitly export final initial-phase moments (emom/emom2/emomM/mmom/mmom2/mmomi)
   // so measurement phase can restart from this exact state.
   gpuSim.copyToFortran();
   stopwatch.add("final synchronize");
   ipTemp.FreeHost();
}

// Spin Dynamics measurement phase
void GpuSimulation::GpuSDSimulation::SDmphase(GpuSimulation& gpuSim) {
   // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);
   std::printf("GpuSDSimulation: SD measurement phase starting\n");

   // Initiated?
   if(!gpuSim.isInitiated) {
      std::fprintf(stderr, "GpuSimulation: not initiated!\n");
      return;
   }

   // initiateMatrices() copied the current Fortran state immediately before
   // gpuRunSimulation() called this routine.  A second full upload here would
   // repeat the immutable Hamiltonian/neighbour-list transfer as well as all
   // lattice arrays.

   // Make phase boundary explicit: start SDmphase with emom2 synchronized to emom.
   // This removes any dependence on historical emom2 content from prior phase bookkeeping.
   gpuSim.gpuLattice.emom2.copy_sync(gpuSim.gpuLattice.emom);
   // Only emom2 changed.  Exporting the complete lattice here used to add a
   // needless device-to-host transfer before the first measurement step.
   gpuSim.cpuLattice.emom2.copy_sync(gpuSim.gpuLattice.emom2);
   gpuSim.cpuLattice.emom2.copy_to(FortranData::emom2);

   // Timer
   StopwatchDeviceSync stopwatch(GlobalStopwatchPool::get("GPU measurement phase"),
                                 ParallelizationHelperInstance.getWorkStream());

   // Initiate default parallelization helper
    ParallelizationHelperInstance.initiate(gpuSim.SimParam.N, gpuSim.SimParam.M, gpuSim.SimParam.NH);
   // Depontd integrator
   GpuDepondtIntegrator integrator;

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

   // TEMPORARY PRINTING
   std::printf("\n");
   // std::printf("________DEBUG System Information:___________ \n");
   // std::printf("%zu\n", gpuSim.cpuHamiltonian.j_tensor.extent(0));
   // std::printf("%zu\n", gpuSim.cpuHamiltonian.j_tensor.extent(1));
   // std::printf("%zu\n", gpuSim.cpuHamiltonian.j_tensor.extent(2));
   // std::printf("%zu\n", gpuSim.cpuHamiltonian.j_tensor.extent(3));
   // std::printf("______________________________________\n");
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
      // Match the CPU convention: measure the state at the start of this step.
      measure_ene = ((gpuSim.Flags.do_ene > 0) && gpuSim.Flags.do_gpu_measurements &&
            (((mstep - 1) % gpuSim.SimParam.ene_step == 0) ||
             (gpuSim.Flags.do_cumu && ((mstep - 1) % gpuSim.SimParam.cumu_step == 0))));

      // Apply the Hamiltonian to the measured configuration.
      hamCalc.heisge(gpuSim.gpuLattice, gpuSim.gpuEnergies, measure_ene);
      stopwatch.add("hamiltonian");

      // Measure
      measurement->measure(mstep);
      correlation->measure(mstep);

      stopwatch.add("measurement");

      // Print simulation status for each 5% of the simulation length
      printMdStatus(mstep, gpuSim);

      // Perform first step of SDE solver
      integrator.evolveFirst(gpuSim.gpuLattice);
      stopwatch.add("evolution");

      // Apply the predictor field needed by the corrector without measuring it.
      hamCalc.heisge(gpuSim.gpuLattice, gpuSim.gpuEnergies, false);
      stopwatch.add("hamiltonian");
  

      // Perform second (corrector) step of SDE solver
      integrator.evolveSecond(gpuSim.gpuLattice);
      stopwatch.add("evolution");
      // Update magnetic moments after time evolution step
      momUpdater.update();
      stopwatch.add("moments");

      measurement->updateAC(mstep);

      // Check for error
      GPU_ERROR_T e = GPU_GET_LAST_ERROR();
      if(e != GPU_SUCCESS) {
         std::printf("Uncaught GPU error %d: %s\n", e, GPU_GET_ERROR_STRING(e));
         // Fatal GPU error: flush and _Exit. Do NOT device-reset then std::exit -
         // that runs the static GpuSimulation dtor, which cudaFrees on a torn-down
         // context and aborts with "driver shutting down". The OS reclaims the context.
         std::fflush(stdout);
         std::fflush(stderr);
         std::_Exit(EXIT_FAILURE);
      }


   }  // End loop over simulation steps

   // Final measure and print remaining measurements to file

   measure_ene = ((gpuSim.Flags.do_ene > 0 ) && (gpuSim.Flags.do_gpu_measurements));
   hamCalc.heisge(gpuSim.gpuLattice, gpuSim.gpuEnergies, measure_ene);

   measurement->measure(rstep + nstep + 1);    
   correlation->measure(rstep + nstep + 1);  // TODO
   stopwatch.add("measurement");

   // Drain the fast-copy path before retiring the queue: the final measure()
   // only queues a host callback on the copy stream, and that callback is what
   // pushes the measurement. Calling mqueue.finish() first would join the
   // consumer thread and destroy its mutex while that push is still pending.
   GPU_STREAM_SYNC(ParallelizationHelperInstance.getWorkStream());
   GPU_STREAM_SYNC(ParallelizationHelperInstance.getCopyStream());

   mqueue.finish();

   // Print remaining measurements
   measurement->flushMeasurements(rstep + nstep + 1);  // TODO
   correlation->flushCorrelations(gpuSim.cpuCorrelations, rstep + nstep + 1); 
   
   // Transfer GPU sample count back to Fortran for averaging. Only when the
   // correlations were actually computed on the GPU: with CPU correlations
   // (FortranCorrelation) the Fortran sampler already maintains sc%sc_tidx /
   // sc%sc_nsamp - and sc_tidx_ptr aliases sc%sc_tidx - so writing the GPU
   // counters (which stay 0 on this path) would clobber the CPU-accumulated
   // sample count and suppress all S(q,w) output.
   if (gpuSim.Flags.do_gpu_correlations) {
      if (FortranData::sc_nsamp_ptr != nullptr) {
          *FortranData::sc_nsamp_ptr = gpuSim.cpuCorrelations.sc_nsamp;
      }
      if (FortranData::sc_tidx_ptr != nullptr) {
         *FortranData::sc_tidx_ptr = gpuSim.cpuCorrelations.sc_tidx;
      }
   }
   
   stopwatch.add("flush measurement");


   // Finish computation and the fast-copy callback path before measurement
   // objects and their pinned buffers are released.
   GPU_STREAM_SYNC(ParallelizationHelperInstance.getWorkStream());
   GPU_STREAM_SYNC(ParallelizationHelperInstance.getCopyStream());

   // Explicitly export the final measurement-phase state, mirroring SDiphase.
   // Without this the Fortran arrays backing the restart file only hold
   // whatever printMdStatus last snapshotted: the slow-copy measurement path
   // happens to refresh them as a side effect, the fast-copy path does not.
   gpuSim.copyToFortran();
   stopwatch.add("final synchronize");
}

