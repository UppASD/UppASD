#include <cuda.h>
#include <pthread.h>

#include "c_headers.hpp"
#include "c_helper.h"
#include "cudaMeasurement.hpp"
#include "cudaParallelizationHelper.hpp"
#include "fortranData.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "stopwatchPool.hpp"
#include "tensor.cuh"

#ifdef NVPROF
#include <nvToolsExtCuda.h>
#endif

#include "measurementQueue.hpp"

// Constructor
CudaMeasurement::CudaMeasurement(const CudaTensor<real, 3>& p1, const CudaTensor<real, 3>& p2,
                                 const CudaTensor<real, 2>& p3, Tensor<real, 3>& p4, Tensor<real, 3>& p5,
                                 Tensor<real, 2>& p6, bool p7, bool p8)
    : emomM(p1),
      emom(p2),
      mmom(p3),
      fortran_emomM(p4),
      fortran_emom(p5),
      fortran_mmom(p6),
      fastCopy(p7),
      alwaysCopy(p8),
      stopwatch(GlobalStopwatchPool::get("Cuda measurement")),
      parallel(CudaParallelizationHelper::def) {
#ifdef NVPROF
   nvtxNameOsThread(pthread_self(), "MAIN_THREAD");
#endif

   if(fastCopy) {
      // Initate temporary
      unsigned int N = emomM.extent(1);
      unsigned int M = emomM.extent(2);


      tmp_emomM.Allocate(3, N, M);
      tmp_emom.Allocate(3, N, M);
      tmp_mmom.Allocate(N, M);

      // Initiate pinned memory
      // pinned_emomM = pinned_emom = pinned_mmom = nullptr;

      pinned_emomM.AllocateHost(3, N, M);
      pinned_emom.AllocateHost(3, N, M);
      pinned_mmom.AllocateHost(N, M);
      // cudaHostAlloc(&pinned_emomM, emomM.bytes(), cudaHostAllocDefault);
      // cudaHostAlloc(&pinned_emom, emom.bytes(), cudaHostAllocDefault);
      // cudaHostAlloc(&pinned_mmom, mmom.bytes(), cudaHostAllocDefault);
   }
}

// Destructor
CudaMeasurement::~CudaMeasurement() {
   if(fastCopy) {
      tmp_emomM.Free();
      tmp_emom.Free();
      tmp_mmom.Free();

      pinned_emomM.FreeHost();
      pinned_emom.FreeHost();
      pinned_mmom.FreeHost();
   }
}

// Callback
void CudaMeasurement::queue_callback(cudaStream_t, cudaError_t, void* data) {
#ifdef NVPROF
   nvtxRangePush("queue_callback");
#endif
   queue_callback_data* d = (queue_callback_data*)data;
   d->me->queueMeasurement(d->step);
   delete d;
#ifdef NVPROF
   nvtxRangePop();
#endif
}

// Callback method
void CudaMeasurement::queueMeasurement(std::size_t mstep) {
   measurementQueue.push(mstep, pinned_emomM.data(), pinned_emom.data(), pinned_mmom.data(), mmom.size());
}

// Fast copy and measurement queueing (D -> D, D -> H (async), H -> H)
void CudaMeasurement::copyQueueFast(std::size_t mstep) {
   // Timing
   stopwatch.skip();

   // Streams
   cudaStream_t workStream = parallel.getWorkStream();
   cudaStream_t copyStream = parallel.getCopyStream();

   // Create new events
   CudaEventPool::Event& workDone = eventPool.get();
   CudaEventPool::Event& copyDone = eventPool.get();

   // The copying must wait for the work stream to finish
   cudaEventRecord(workDone.event(), workStream);
   cudaStreamWaitEvent(copyStream, workDone.event(), 0);

   // Async copy in copy stream (device -> temp. device)
   tmp_emomM.copy_async(emomM, copyStream);
   tmp_emom.copy_async(emom, copyStream);
   tmp_mmom.copy_async(mmom, copyStream);
   cudaEventRecord(copyDone.event(), copyStream);
   stopwatch.add("fast - D2D");

   // Then write to host in copy stream (asynchronously with work stream)

   pinned_emomM.copy_async(tmp_emomM, copyStream);
   pinned_emom.copy_async(tmp_emom, copyStream);
   pinned_mmom.copy_async(tmp_mmom, copyStream);

   // Make the work stream wait out the copying
   cudaStreamWaitEvent(workStream, copyDone.event(), 0);
   copyDone.addDeactivateCallback(workStream);
   workDone.addDeactivateCallback(workStream);

   // Push to measure queue when done
   cudaStreamAddCallback(workStream, queue_callback, new queue_callback_data(this, mstep), 0);
}

// Slow copying (D -> H)
void CudaMeasurement::copyQueueSlow(std::size_t mstep) {
   // Timing
   stopwatch.skip();

   // Write directly to fortran
   // (this can't be done async, so it will block host until finished)
   fortran_emomM.copy_sync(emomM);
   fortran_mmom.copy_sync(mmom);
   fortran_emom.copy_sync(emom);


   // emomM.write(FortranData::emomM);
   // emom.write(FortranData::emom);
   // mmom.write(FortranData::mmom);
   stopwatch.add("slow - D2H copy");

   // Queue measurement
   measurementQueue.push(mstep, fortran_emomM.data(), fortran_emom.data(), fortran_mmom.data(), mmom.size());
}

void CudaMeasurement::measure(std::size_t mstep) {
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

void CudaMeasurement::flushMeasurements(std::size_t mstep) {
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

