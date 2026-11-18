
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

void GpuSimulation::GpuSDSimulation::printMdStatus_iphase(std::size_t mstep,GpuSimulation& gpuSim, std::size_t phaseSteps,
                                    unsigned int phaseIndex, real tempAvg, real tempMin, real tempMax) {

    if (phaseSteps == 0) {
        return;
    }

    const std::size_t printEvery = (phaseSteps > 20) ? phaseSteps / 20 : 1;
    if (mstep % printEvery != 0 && mstep != phaseSteps) {
        return;
    }

    gpuSim.copyToFortran();
    fortran_calc_simulation_status_variables(gpuSim.SimParam.mavg);

    const unsigned int progress = static_cast<unsigned int>(mstep * 100 / phaseSteps);

    std::printf(
        "GPU: IP phase %u/%u | Tavg: %.4f K "
        "(Tmin: %.4f, Tmax: %.4f) | "
        "%3u%% done (%zu/%zu steps). Mbar: %10.6f. U: %8.5f.\n",
        phaseIndex + 1,
        gpuSim.SimParam.ipnphase,
        static_cast<double>(tempAvg),
        static_cast<double>(tempMin),
        static_cast<double>(tempMax),
        progress,
        mstep,
        phaseSteps,
        static_cast<double>(*gpuSim.SimParam.mavg),
        static_cast<double>(*gpuSim.SimParam.binderc));
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

   // Timer
   StopwatchDeviceSync stopwatch = StopwatchDeviceSync(GlobalStopwatchPool::get("Gpu initial phase"));

   // Initiate default parallelization helper
   ParallelizationHelperInstance.initiate(gpuSim.SimParam.N, gpuSim.SimParam.M, gpuSim.SimParam.NH);
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

   std::printf("\n");

   int mnn = gpuSim.cpuHamiltonian.j_tensor.extent(2);
   int l = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int N = gpuSim.SimParam.N;
   std::printf("_______________________________________________\n");
   Tensor<real, 1> ipTemp;
   ipTemp.AllocateHost(N);


   if (FortranData::ipnphase == nullptr) {

      std::fprintf(stderr,
         "\nERROR: GPU Spin Dynamics initial phase is not configured correctly.\n"
         "The scalar input format using ip_temp and ip_nstep is not supported.\n"
         "\n"
         "Please use the conventional initial-phase input format:\n"
         "\n"
         "ip_mode S\n"
         "ip_nphase 1\n"
         "250000  0.01  1.0e-15  1.0\n"
         "\n"
         "The columns are:\n"
         "  Number of steps | Temperature | Time step | Damping\n"
         "\n"
         "The number of data rows must match ip_nphase.\n"
         "GPU initial phase will not be started.\n\n");

      return;
   }

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

      real tempMin = gpuSim.cpuLattice.ipTemp_array(0, it);
      real tempMax = tempMin;
      double tempSum = 0.0;

      for (std::size_t atom = 0; atom < gpuSim.SimParam.N; atom++) {

        const real T = gpuSim.cpuLattice.ipTemp_array(atom, it);

        tempSum += static_cast<double>(T);

        tempMin = (T < tempMin) ? T : tempMin;
        tempMax = (T > tempMax) ? T : tempMax;
      }
      const real tempAvg = static_cast<real>(tempSum / gpuSim.SimParam.N);

      std::printf(
        "\nGPU: Starting initial phase %u/%u\n"
        "    Temperature: %.4f K (min: %.4f, max: %.4f)\n"
        "    Steps: %u | dt: %.4e | damping: %.4f\n",
        it + 1,
        ipnphase,
        static_cast<double>(tempAvg),
        static_cast<double>(tempMin),
        static_cast<double>(tempMax),
        steps,
        static_cast<double>(phaseDt),
        static_cast<double>(phaseDamp));


      // Time step loop
      for(std::size_t mstep = 1; mstep <= steps; mstep++) {
         // Print simulation status for each 5% of the simulation length
         printMdStatus_iphase(mstep, gpuSim, steps, it, tempAvg, tempMin, tempMax);

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

         ASSERT_GPU_KERNEL("SDiphase");
      }

      if(it < (ipnphase - 1)){
         for(unsigned int k = 0; k < N; k++){
            ipTemp(k) = gpuSim.cpuLattice.ipTemp_array(k, it + 1);
          }
         }
      }  
      
   // Synchronize with device
   ASSERT_GPU(GPU_DEVICE_SYNCHRONIZE());
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
   GpuDepondtIntegrator integrator;

   // Hamiltonian calculations
   GpuHamiltonianCalculations hamCalc;

   // Moment updater
   GpuMomentUpdater momUpdater(gpuSim.gpuLattice, gpuSim.SimParam.mompar, gpuSim.SimParam.initexc);
   //Queue
   MeasurementQueue mqueue;
   // Measurement
   const auto measurement = MeasurementFactory::create(gpuSim.gpuLattice, gpuSim.cpuLattice, gpuSim.gpuEnergies, gpuSim.cpuProj, mqueue, gpuSim.Flags.do_jtensor);
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

   const bool energy_step = ((mstep - 1) % gpuSim.SimParam.ene_step == 0);
   const bool cumu_step = gpuSim.Flags.do_cumu && ((mstep - 1) % gpuSim.SimParam.cumu_step == 0);
   measure_ene = (gpuSim.Flags.do_ene > 0)&& (gpuSim.Flags.do_gpu_measurements) && (energy_step || cumu_step);

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

      // Apply Hamiltonian to obtain effective field
      hamCalc.heisge(gpuSim.gpuLattice, gpuSim.gpuEnergies, false);
      stopwatch.add("hamiltonian");
  

      // Perform second (corrector) step of SDE solver
      integrator.evolveSecond(gpuSim.gpuLattice);
      stopwatch.add("evolution");
      // Update magnetic moments after time evolution step
      momUpdater.update();
      stopwatch.add("moments");

      measurement->updateAC(mstep);

      ASSERT_GPU_KERNEL("SDmphase");

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
   ASSERT_GPU(GPU_DEVICE_SYNCHRONIZE());
   stopwatch.add("final synchronize");
}

