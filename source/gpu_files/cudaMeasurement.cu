#include <cstddef>
#include <cstdio>
#include <cstdlib>

using namespace std;

#include <cuda.h>
#include <pthread.h>

#include "c_helper.h"
#include "cudaMatrix.hpp"
#include "cudaMeasurement.hpp"
#include "cudaParallelizationHelper.hpp"
#include "fortranData.hpp"
#include "hostMatrix.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "stopwatchPool.hpp"

#ifdef NVPROF
#include <nvToolsExtCuda.h>
#endif

#include "measurementQueue.hpp"

// Constructor
CudaMeasurement::CudaMeasurement(const cudaMatrix<real, 3, 3> &p1, const cudaMatrix<real, 3, 3> &p2,
                                 const cudaMatrix<real, 2> &p3, bool p4, bool p5)
    : emomM(p1),
      emom(p2),
      mmom(p3),
      fastCopy(p4),
      alwaysCopy(p5),
      stopwatch(GlobalStopwatchPool::get("Cuda measurement")),
      parallel(CudaParallelizationHelper::def) {
#ifdef NVPROF
   nvtxNameOsThread(pthread_self(), "MAIN_THREAD");
#endif

   if(fastCopy) {
      // Initate temporary
      tmp_emomM.initiate(emomM);
      tmp_emom.initiate(emom);
      tmp_mmom.initiate(mmom);

      // Initiate pinned memory
      pinned_emomM = pinned_emom = pinned_mmom = nullptr;
      cudaHostAlloc(&pinned_emomM, emomM.data_size(), cudaHostAllocDefault);
      cudaHostAlloc(&pinned_emom, emom.data_size(), cudaHostAllocDefault);
      cudaHostAlloc(&pinned_mmom, mmom.data_size(), cudaHostAllocDefault);

      // Any out of memory?
      if(cudaPeekAtLastError() == cudaErrorMemoryAllocation) {
         // Flush error
         cudaGetLastError();

         // Free possible allocations
         tmp_emomM.free();
         tmp_emom.free();
         tmp_mmom.free();

         if(pinned_emomM != nullptr) {
            cudaFreeHost(pinned_emomM);
         }
         if(pinned_emom != nullptr) {
            cudaFreeHost(pinned_emom);
         }
         if(pinned_mmom != nullptr) {
            cudaFreeHost(pinned_mmom);
         }

         // Fall back to slow copy
         printf("Failed to allocate memory for fast copy in measurements!\n");
         fastCopy = false;
         return;
      }
   }
}

// Destructor
CudaMeasurement::~CudaMeasurement() {
   if(fastCopy) {
      cudaFreeHost(pinned_emomM);
      cudaFreeHost(pinned_emom);
      cudaFreeHost(pinned_mmom);
   }
}

// Callback
void CudaMeasurement::queue_callback(cudaStream_t, cudaError_t, void *data) {
#ifdef NVPROF
   nvtxRangePush("queue_callback");
#endif
   queue_callback_data *d = (queue_callback_data *)data;
   d->me->queueMeasurement(d->step);
   delete d;
#ifdef NVPROF
   nvtxRangePop();
#endif
}

// Callback method
void CudaMeasurement::queueMeasurement(size_t mstep) {
   measurementQueue.push(mstep, pinned_emomM, pinned_emom, pinned_mmom, mmom.size());
}

// Fast copy and measurement queueing (D -> D, D -> H (async), H -> H)
void CudaMeasurement::copyQueueFast(size_t mstep) {
   // Timing
   stopwatch.skip();

   // Streams
   cudaStream_t workStream = parallel.getWorkStream();
   cudaStream_t copyStream = parallel.getCopyStream();

   // Create new events
   CudaEventPool::Event &workDone = eventPool.get();
   CudaEventPool::Event &copyDone = eventPool.get();

   // The copying must wait for the work stream to finish
   cudaEventRecord(workDone.event(), workStream);
   cudaStreamWaitEvent(copyStream, workDone.event(), 0);

   // Async copy in copy stream (device -> temp. device)
   tmp_emomM.memcopy(emomM, copyStream);
   tmp_emom.memcopy(emom, copyStream);
   tmp_mmom.memcopy(mmom, copyStream);
   cudaEventRecord(copyDone.event(), copyStream);
   stopwatch.add("fast - D2D");

   // Then write to host in copy stream (asynchronously with work stream)
   tmp_emomM.writeAsync(pinned_emomM, copyStream);
   tmp_emom.writeAsync(pinned_emom, copyStream);
   tmp_mmom.writeAsync(pinned_mmom, copyStream);

   // Make the work stream wait out the copying
   cudaStreamWaitEvent(workStream, copyDone.event(), 0);
   copyDone.addDeactivateCallback(workStream);
   workDone.addDeactivateCallback(workStream);

   // Push to measure queue when done
   cudaStreamAddCallback(workStream, queue_callback, new queue_callback_data(this, mstep), 0);
}

// Slow copying (D -> H)
void CudaMeasurement::copyQueueSlow(size_t mstep) {
   // Timing
   stopwatch.skip();

   // Write directly to fortran
   // (this can't be done async, so it will block host until finished)
   emomM.write(FortranData::emomM);
   emom.write(FortranData::emom);
   mmom.write(FortranData::mmom);
   stopwatch.add("slow - D2H copy");

   // Queue measurement
   measurementQueue.push(mstep, FortranData::emomM, FortranData::emom, FortranData::mmom, mmom.size());
}

void CudaMeasurement::measure(size_t mstep) {
   // Copy required?
   bool copy = (alwaysCopy || fortran_do_measurements(mstep));

   if(copy) {
      // Copy and queue
      if(fastCopy) {
         copyQueueFast(mstep);
      } else {
         copyQueueSlow(mstep);
      }
   } else {
      // Push empty measurement
      measurementQueue.push(mstep);
   }
}

void CudaMeasurement::flushMeasurements(size_t mstep) {
   // Timing
   stopwatch.skip();

   // Wait out possible queue callbacks
   cudaStreamSynchronize(parallel.getWorkStream());

   // Flush internal queue
   measurementQueue.finish();

   // Print remaining measurements
   fortran_flush_measurements(mstep);
   stopwatch.add("measure");
}

