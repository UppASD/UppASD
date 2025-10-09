
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

#include "gpu_wrappers.h"
#if defined(HIP_V)
#include <hip/hip_runtime.h>
#include <hiprand/hiprand.h>
#include "gpuCorrelations.hpp"
#elif defined(CUDA_V)
#include <cuda.h>
#include <curand.h>
#include "gpuCorrelations.cuh"
#endif

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

// Spin Dynamics measurement phase
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

   // TEMPORARY PRINTING
   std::printf("\n");
   std::printf("________DEBUG System Information:___________ \n");
   std::printf("%zu\n", gpuSim.cpuHamiltonian.j_tensor.extent(0));
   std::printf("%zu\n", gpuSim.cpuHamiltonian.j_tensor.extent(1));
   std::printf("%zu\n", gpuSim.cpuHamiltonian.j_tensor.extent(2));
   std::printf("%zu\n", gpuSim.cpuHamiltonian.j_tensor.extent(3));
   std::printf("______________________________________\n");
   int mnn = gpuSim.cpuHamiltonian.j_tensor.extent(2);
   int l = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int N = gpuSim.SimParam.N;
   std::printf("_______________________________________________\n");
   Tensor<real, 1> ipTemp;
   ipTemp.AllocateHost(N);

   for(unsigned int k = 0; k < N; k++){
      ipTemp(k) = gpuSim.cpuLattice.ipTemp_array(k, 1);
   }
   // Initiate constants for integrator
   integrator.initiateConstants(gpuSim.SimParam, ipTemp);
   // Timing
   stopwatch.add("initiate");

   unsigned int steps;
   int ipnphase = gpuSim.SimParam.ipnphase; 
   for(unsigned int it = 0; it < ipnphase; it++){
   steps = gpuSim.cpuLattice.ipnstep(it);

   // Time step loop
   for(std::size_t mstep = 1; mstep <= steps; mstep++) {
      // Print simulation status for each 5% of the simulation length
      // printMdStatus(mstep); -- Do we need it in initial phase?

      // Apply Hamiltonian to obtain effective field
      hamCalc.heisge(gpuSim.gpuLattice);
      stopwatch.add("hamiltonian");

      // Perform first step of SDE solver
      integrator.evolveFirst(gpuSim.gpuLattice);
      stopwatch.add("evolution");

      // Apply Hamiltonian to obtain effective field
      hamCalc.heisge(gpuSim.gpuLattice);
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
         GPU_DEVICE_RESET();
         std::exit(EXIT_FAILURE);
      }

      for(unsigned int k = 0; k < N; k++){
      ipTemp(k) = gpuSim.cpuLattice.ipTemp_array(k, it + 1);
      integrator.resetConstants(ipTemp);
   }

   }
   }  // End loop over simulation steps
   // Synchronize with device
   GPU_DEVICE_SYNCHRONIZE();
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
   // Measurement
   const auto measurement = MeasurementFactory::create(gpuSim.gpuLattice, gpuSim.cpuLattice);
   //Corrrelations
   const auto correlation = CorrelationFactory::create(gpuSim.gpuLattice, gpuSim.cpuLattice, gpuSim.Flags, gpuSim.SimParam, gpuSim.cpuCorrelations);

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
   std::printf("________DEBUG System Information:___________ \n");
   std::printf("%zu\n", gpuSim.cpuHamiltonian.j_tensor.extent(0));
   std::printf("%zu\n", gpuSim.cpuHamiltonian.j_tensor.extent(1));
   std::printf("%zu\n", gpuSim.cpuHamiltonian.j_tensor.extent(2));
   std::printf("%zu\n", gpuSim.cpuHamiltonian.j_tensor.extent(3));
   std::printf("______________________________________\n");
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

   // Time step loop
   for(std::size_t mstep = rstep + 1; mstep <= rstep + nstep; mstep++) {
      // Measure
      measurement->measure(mstep);
      stopwatch.add("measurement");

      // Print simulation status for each 5% of the simulation length
      printMdStatus(mstep, gpuSim);

      // Apply Hamiltonian to obtain effective field
      hamCalc.heisge(gpuSim.gpuLattice);
      stopwatch.add("hamiltonian");

      // Perform first step of SDE solver
      integrator.evolveFirst(gpuSim.gpuLattice);
      stopwatch.add("evolution");

      // Apply Hamiltonian to obtain effective field
      hamCalc.heisge(gpuSim.gpuLattice);
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
         GPU_DEVICE_RESET();
         std::exit(EXIT_FAILURE);
      }

   }  // End loop over simulation steps

   // Final measure
   measurement->measure(rstep + nstep + 1);  // TODO
   stopwatch.add("measurement");

   // Print remaining measurements
   measurement->flushMeasurements(rstep + nstep + 1);  // TODO
   stopwatch.add("flush measurement");

   // Synchronize with device
   GPU_DEVICE_SYNCHRONIZE();
   stopwatch.add("final synchronize");
}

