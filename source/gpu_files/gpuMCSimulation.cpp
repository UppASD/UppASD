#pragma once
#include "c_headers.hpp"
#include "c_helper.h"
//#include "cudaGPUErrchk.hpp"
#include "gpuHamiltonianCalculations.hpp"
#include "gpuMetropolis.hpp"
#include "gpuMetropolis_bruteforce.hpp"
#include "gpuParallelizationHelper.hpp"
#include "gpuSimulation.hpp"
#include "gpuStructures.hpp"
#include "fortranData.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "stopwatchPool.hpp"
#include "tensor.hpp"
#include "measurementFactory.hpp"

#include "gpu_wrappers.h"
#if defined(HIP_V)
#include <hip/hip_runtime.h>
#include <hiprand/hiprand.h>
#elif defined(CUDA_V)
#include <cuda.h>
#include <curand.h>
#endif
using ParallelizationHelper = GpuParallelizationHelper;

GpuSimulation::GpuMCSimulation::GpuMCSimulation() {
   // isInitiatedSD = false;
}

GpuSimulation::GpuMCSimulation::~GpuMCSimulation() {
}

void GpuSimulation::GpuMCSimulation::printMdStatus(std::size_t mstep, GpuSimulation& gpuSim) {
   if(gpuSim.SimParam.mcnstep > 20) {
      if(mstep % ((gpuSim.SimParam.rstep + gpuSim.SimParam.mcnstep) / 20) == 0) {
         gpuSim.copyToFortran();  // This is run so seldomly it has not impact on overall performance
         fortran_calc_simulation_status_variables(gpuSim.SimParam.mavg);
         std::printf("GPU: %3ld%% done. Mbar: %10.6f. U: %8.5f.\n",
                     mstep * 100 / (gpuSim.SimParam.rstep + gpuSim.SimParam.mcnstep),
                     *gpuSim.SimParam.mavg,
                     *gpuSim.SimParam.binderc);
      }
   } else {
      gpuSim.copyToFortran();
      fortran_calc_simulation_status_variables(gpuSim.SimParam.mavg);
      std::printf("GPU: Iteration %ld Mbar %13.6f\n", mstep, *gpuSim.SimParam.mavg);
   }
}

void GpuSimulation::GpuMCSimulation::printMdStatus_iphase(std::size_t mstep, GpuSimulation& gpuSim, int step) {
   if(step > 20) {
      if(mstep % ((step) / 20) == 0) {
         //cudaMC.mom_update(gpuSim.gpuLattice);

         gpuSim.copyToFortran();  // This is run so seldomly it has not impact on overall performance
         fortran_calc_simulation_status_variables(gpuSim.SimParam.mavg);
         std::printf("GPU: %3ld%% done. Mbar: %10.6f. U: %8.5f.\n",
                     mstep * 100 / (step),
                     *gpuSim.SimParam.mavg,
                     *gpuSim.SimParam.binderc);
      }
   } else {
      //cudaMC.mom_update(GpuSim.gpuLattice);
      gpuSim.copyToFortran();
      fortran_calc_simulation_status_variables(gpuSim.SimParam.mavg);
      std::printf("GPU: Iteration %ld Mbar %13.6f\n", mstep, *gpuSim.SimParam.mavg);
   }
}

// Monte Carlo initial phase (strict)
void GpuSimulation::GpuMCSimulation::MCiphase(GpuSimulation& gpuSim) {
   // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);
   std::printf("GpuMCSimulation: MC initial phase starting\n");

   // Initiated?
   if(!gpuSim.isInitiated) {
      std::fprintf(stderr, "GpuSimulation: not initiated!\n");
      return;
   }

   // Initiate default parallelization helper
    ParallelizationHelperInstance.initiate(gpuSim.SimParam.N, gpuSim.SimParam.M, gpuSim.SimParam.NH);

   // Timer
   StopwatchDeviceSync stopwatch = StopwatchDeviceSync(GlobalStopwatchPool::get("Gpu initial phase"));

    // Metropolis
   GpuMetropolis gpuMC;

   // Hamiltonian calculations
   GpuHamiltonianCalculations hamCalc;
   unsigned int num_subL = 0;
   //printf("HERE - 0\n");

   // Initiate MC and Hamiltonian
   num_subL = gpuMC.initiate(gpuSim.SimParam, gpuSim.cpuHamiltonian, gpuSim.cpuLattice);
   //printf("HERE - 1\n");
   if(!hamCalc.initiate(gpuSim.Flags, gpuSim.SimParam, gpuSim.gpuHamiltonian)) {
      std::fprintf(stderr, "GpuMCSimulation: Hamiltonian failed to initiate!\n");
      return;
   }

   int mnn = gpuSim.cpuHamiltonian.j_tensor.extent(2);
   int l = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int ipmcnphase = gpuSim.SimParam.ipmcnphase;
   real beta; 
   unsigned int mcs;
   // Timing
   stopwatch.add("initiate");
   printf("\n\\ninphase = %i\n\n", ipmcnphase);
   //printf("HERE - 2\n");

   for(unsigned int it = 0; it < ipmcnphase; it++){
   mcs = gpuSim.cpuLattice.ipmcnstep(it);
   beta = 1/(gpuSim.cpuLattice.ipTemp(it) * gpuSim.SimParam.k_bolt);
   // Apply Hamiltonian to obtain effective field
   hamCalc.heisge(gpuSim.gpuLattice);
   stopwatch.add("hamiltonian");
   //printf("HERE - 3\n");

   //printf("mcs = %i\n", mcs);   
   // Time step loop
   for(std::size_t mstep = 0; mstep <= mcs; mstep++) {
      printMdStatus_iphase(mstep, gpuSim, mcs);
      //printMdStatus(mstep, gpuSim);
         // Perform Metropolis sweep
         for(std::size_t sub = 0; sub < num_subL; sub++){
            gpuMC.MCrun(gpuSim.gpuLattice, beta, sub);
            stopwatch.add("montecarlo");
   //printf("HERE - 4\n");

             // Apply Hamiltonian to obtain effective field
            hamCalc.heisge(gpuSim.gpuLattice);
            stopwatch.add("hamiltonian");
         }

         // Check for error
         GPU_ERROR_T e = GPU_GET_LAST_ERROR();
         if(e != GPU_SUCCESS) {
            std::printf("Uncaught GPU error %d: %s\n", e, GPU_GET_ERROR_STRING(e));
            GPU_DEVICE_RESET();
            std::exit(EXIT_FAILURE);
         }
        // printf("mcs = %i\n", mstep);

   }  
   // End loop over simulation steps
   // Synchronize with device
   //cudaDeviceSynchronize();
   //printf("HERE - 3\n");
   //printf("HERE - 5\n");

   gpuMC.mom_update(gpuSim.gpuLattice);
   //printf("HERE - 6\n");

   GPU_DEVICE_SYNCHRONIZE();   
   }



   //cudaMC.release();
   stopwatch.add("final synchronize");
}

// Monte Carlo measurement phase (strict)
void GpuSimulation::GpuMCSimulation::MCmphase(GpuSimulation& gpuSim) {
   // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);
   std::printf("GpuMCSimulation: MC measurement phase starting\n");

   // Initiated?
   if(!gpuSim.isInitiated) {
      std::fprintf(stderr, "GpuSimulation: not initiated!\n");
      return;
   }
   // Initiate default parallelization helper
    ParallelizationHelperInstance.initiate(gpuSim.SimParam.N, gpuSim.SimParam.M, gpuSim.SimParam.NH);
   // Timer
   StopwatchDeviceSync stopwatch = StopwatchDeviceSync(GlobalStopwatchPool::get("Gpu measurement phase"));

    // MC
   GpuMetropolis gpuMC;

   // Hamiltonian calculations
   GpuHamiltonianCalculations hamCalc;
   unsigned int num_subL = 0;

 
   // Initiate MC and Hamiltonian
   num_subL = gpuMC.initiate(gpuSim.SimParam, gpuSim.cpuHamiltonian, gpuSim.cpuLattice);

   if(!hamCalc.initiate(gpuSim.Flags, gpuSim.SimParam, gpuSim.gpuHamiltonian)) {
      std::fprintf(stderr, "GpuMCSimulation: Hamiltonian failed to initiate!\n");
      return;
   }

 // Measurement
   const auto measurement = MeasurementFactory::create(gpuSim.gpuLattice, gpuSim.cpuLattice);
 
   int mnn = gpuSim.cpuHamiltonian.j_tensor.extent(2);
   int l = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   real beta = 1 / (gpuSim.SimParam.Temp * gpuSim.SimParam.k_bolt);

   // Timing
   stopwatch.add("initiate");

   size_t mcnstep = gpuSim.SimParam.mcnstep;

   // Time step loop
   for(std::size_t mstep = 1; mstep <= mcnstep; mstep++) {
      // Measure
      //printf("STEP = %i\n", mstep);
      measurement->measure(mstep);
      stopwatch.add("measurement");

      // Print simulation status for each 5% of the simulation length
      printMdStatus(mstep, gpuSim);

      // Apply Hamiltonian to obtain effective field
      //hamCalc.heisge(gpuSim.gpuLattice);
      //stopwatch.add("hamiltonian");
      for(unsigned int sub = 0; sub < num_subL; sub++){
         // Perform Metropolis sweep
         gpuMC.MCrun(gpuSim.gpuLattice, beta, sub);
         stopwatch.add("montecarlo");
         hamCalc.heisge(gpuSim.gpuLattice);
         stopwatch.add("hamiltonian");


      }



      // Check for error
      GPU_ERROR_T e = GPU_GET_LAST_ERROR();
      if(e != GPU_SUCCESS) {
         std::printf("Uncaught GPU error %d: %s\n", e, GPU_GET_ERROR_STRING(e));
         GPU_DEVICE_RESET();
         std::exit(EXIT_FAILURE);
      }

   }  // End loop over simulation steps

   // Final measure
   measurement->measure(mcnstep + 1);  // TODO
   stopwatch.add("measurement");

   // Print remaining measurements
   measurement->flushMeasurements(mcnstep + 1);  // TODO
   stopwatch.add("flush measurement");

   // Synchronize with device
   GPU_DEVICE_SYNCHRONIZE();
   //cudaMC.release();
   stopwatch.add("final synchronize");
}

// Monte Carlo initial phase (brute force)
void GpuSimulation::GpuMCSimulation::MCiphase_bf(GpuSimulation& gpuSim) {
   // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);
   std::printf("GpuMCSimulation: MC BF initial phase starting\n");

   // Initiated?
   if(!gpuSim.isInitiated) {
      std::fprintf(stderr, "GpuSimulation: not initiated!\n");
      return;
   }

   // Initiate default parallelization helper
    ParallelizationHelperInstance.initiate(gpuSim.SimParam.N, gpuSim.SimParam.M, gpuSim.SimParam.NH);

   // Timer
   StopwatchDeviceSync stopwatch = StopwatchDeviceSync(GlobalStopwatchPool::get("Gpu measurement phase"));

    // Metropolis
   GpuMetropolis_bruteforce gpuMC_bf;

   // Hamiltonian calculations
   GpuHamiltonianCalculations hamCalc;
   

   // Initiate MC and Hamiltonian
   //num_subL = cudaMC_bf.initiate(gpuSim.SimParam, gpuSim.cpuHamiltonian, gpuSim.cpuLattice);

   if(!hamCalc.initiate(gpuSim.Flags, gpuSim.SimParam, gpuSim.gpuHamiltonian)) {
      std::fprintf(stderr, "GpuMCSimulation: Hamiltonian failed to initiate!\n");
      return;
   }

   if(!gpuMC_bf.initiate(gpuSim.SimParam)) {
      std::fprintf(stderr, "GpuMCSimulation_bf: Hamiltonian failed to initiate!\n");
      return;
   }

   int mnn = gpuSim.cpuHamiltonian.j_tensor.extent(2);
   int l = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int ipmcnphase = gpuSim.SimParam.ipmcnphase;
   real beta; 
   unsigned int mcs;
   // Timing
   stopwatch.add("initiate");
   printf("\n\\ninphase = %i\n\n", ipmcnphase);

   for(unsigned int it = 0; it < ipmcnphase; it++){
   mcs = gpuSim.cpuLattice.ipmcnstep(it);
   beta = 1/(gpuSim.cpuLattice.ipTemp(it) * gpuSim.SimParam.k_bolt);
   // Apply Hamiltonian to obtain effective field
   hamCalc.heisge(gpuSim.gpuLattice);
   stopwatch.add("hamiltonian");

   //printf("mcs = %i\n", mcs);   
   // Time step loop
   for(std::size_t mstep = 0; mstep <= mcs; mstep++) {
      printMdStatus_iphase(mstep, gpuSim, mcs);
      //printMdStatus(mstep, gpuSim);
         // Perform Metropolis sweep

            gpuMC_bf.MCrun(gpuSim.gpuLattice, beta);
            stopwatch.add("montecarlo");
             // Apply Hamiltonian to obtain effective field
            hamCalc.heisge(gpuSim.gpuLattice);
            stopwatch.add("hamiltonian");
         

         // Check for error
      GPU_ERROR_T e = GPU_GET_LAST_ERROR();
      if(e != GPU_SUCCESS) {
         std::printf("Uncaught GPU error %d: %s\n", e, GPU_GET_ERROR_STRING(e));
         GPU_DEVICE_RESET();
         std::exit(EXIT_FAILURE);
      }
        // printf("mcs = %i\n", mstep);

   }  
   // End loop over simulation steps
   // Synchronize with device
   //cudaDeviceSynchronize();
   //printf("HERE - 3\n");
   gpuMC_bf.mom_update(gpuSim.gpuLattice);
   GPU_DEVICE_SYNCHRONIZE();   
   }



   //cudaMC.release();
   stopwatch.add("final synchronize");
}

// Monte Carlo measurement phase (brute force)
void GpuSimulation::GpuMCSimulation::MCmphase_bf(GpuSimulation& gpuSim) {
   // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);
   std::printf("GpuMCSimulation: MC BF measurement phase starting\n");

   // Initiated?
   if(!gpuSim.isInitiated) {
      std::fprintf(stderr, "GpuSimulation: not initiated!\n");
      return;
   }
   // Initiate default parallelization helper
   ParallelizationHelperInstance.initiate(gpuSim.SimParam.N, gpuSim.SimParam.M, gpuSim.SimParam.NH);
   // Timer
   StopwatchDeviceSync stopwatch = StopwatchDeviceSync(GlobalStopwatchPool::get("Cuda measurement phase"));

    // MC
   GpuMetropolis_bruteforce gpuMC_bf;

   // Hamiltonian calculations
   GpuHamiltonianCalculations hamCalc;
   unsigned int num_subL = 0;

 
   // Initiate MC and Hamiltonian
   if(!hamCalc.initiate(gpuSim.Flags, gpuSim.SimParam, gpuSim.gpuHamiltonian)) {
      std::fprintf(stderr, "GpuMCSimulation: Hamiltonian failed to initiate!\n");
      return;
   }

   if(!gpuMC_bf.initiate(gpuSim.SimParam)) {
      std::fprintf(stderr, "GpuMCSimulation_bf: Hamiltonian failed to initiate!\n");
      return;
   }

 // Measurement
const auto measurement = MeasurementFactory::create(gpuSim.gpuLattice, gpuSim.cpuLattice);
 
   int mnn = gpuSim.cpuHamiltonian.j_tensor.extent(2);
   int l = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   real beta = 1 / (gpuSim.SimParam.Temp * gpuSim.SimParam.k_bolt);



   // Timing
   stopwatch.add("initiate");

   size_t mcnstep = gpuSim.SimParam.mcnstep;

   // Time step loop
   for(std::size_t mstep = 1; mstep <= mcnstep; mstep++) {
      // Measure
      //printf("STEP = %i\n", mstep);
      measurement->measure(mstep);
      stopwatch.add("measurement");

      // Print simulation status for each 5% of the simulation length
      printMdStatus(mstep, gpuSim);

      // Apply Hamiltonian to obtain effective field
      //hamCalc.heisge(gpuSim.gpuLattice);
      //stopwatch.add("hamiltonian");

         // Perform Metropolis sweep
         gpuMC_bf.MCrun(gpuSim.gpuLattice, beta);
         stopwatch.add("montecarlo");
         hamCalc.heisge(gpuSim.gpuLattice);
         stopwatch.add("hamiltonian");


      



      // Check for error
      GPU_ERROR_T e = GPU_GET_LAST_ERROR();
      if(e != GPU_SUCCESS) {
         std::printf("Uncaught GPU error %d: %s\n", e, GPU_GET_ERROR_STRING(e));
         GPU_DEVICE_RESET();
         std::exit(EXIT_FAILURE);
      }

   }  // End loop over simulation steps

   // Final measure
   measurement->measure(mcnstep + 1);
   stopwatch.add("measurement");

   // Print remaining measurements
   measurement->flushMeasurements(mcnstep + 1);  // TODO
   stopwatch.add("flush measurement");

   // Synchronize with device
   GPU_DEVICE_SYNCHRONIZE();
   //cudaMC.release();
   stopwatch.add("final synchronize");
}
