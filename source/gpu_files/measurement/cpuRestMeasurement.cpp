#include <pthread.h>
#include "cpuRestMeasurement.hpp"
#include "c_headers.hpp"
#include "c_helper.h"
#include "fortranData.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "stopwatchPool.hpp"
#include "tensor.hpp"
#include "gpuParallelizationHelper.hpp"
#if defined(HIP_V)
#include <hip/hip_runtime.h>
#elif defined(CUDA_V)
#include <cuda.h>
#endif
using ParallelizationHelper = GpuParallelizationHelper;

#ifdef NVPROF
#include <nvToolsExtCuda.h>
#endif

#include "measurementQueue.hpp"

/*   CpuMeasurement(const GpuTensor<real, 3>& emomM, const GpuTensor<real, 3>& emom,
                   const GpuTensor<real, 2>& mmom, const GpuTensor<real, 3>& beff, Tensor<real, 3>& f_emomM, Tensor<real, 3>& f_emom,
                   Tensor<real, 2>& f_mmom, Tensor<real, 3>& f_beff, bool fastCopy = DEFAULT_FAST_COPY,
                   bool alwaysCopy = false);*/
// Constructor
CpuRestMeasurement::CpuRestMeasurement(const GpuTensor<real, 3>& p_emomM, const GpuTensor<real, 3>& p_emom,
                                 const GpuTensor<real, 2>& p_mmom, const GpuTensor<real, 3>& p_beff, Tensor<real, 3>& p_femomM, 
                                 Tensor<real, 3>& p_femom, Tensor<real, 2>& p_fmmom, Tensor<real, 3>& p_fbeff, MeasurementQueue &mq, bool p7, bool p8)
    : emomM(p_emomM),
      emom(p_emom),
      mmom(p_mmom),
      beff(p_beff),
      fortran_emomM(p_femomM),
      fortran_emom(p_femom),
      fortran_mmom(p_fmmom),
      fortran_beff(p_fbeff),
      measurementQueue(mq),
      fastCopy(p7),
      alwaysCopy(p8),
      stopwatch(GlobalStopwatchPool::get("Fortran measurement")),
      parallel( ParallelizationHelperInstance) {
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
      tmp_beff.Allocate(3, N, M);


      // Initiate pinned memory
      // pinned_emomM = pinned_emom = pinned_mmom = nullptr;

      pinned_emomM.AllocateHost(3, N, M);
      pinned_emom.AllocateHost(3, N, M);
      pinned_mmom.AllocateHost(N, M);
      pinned_beff.AllocateHost(3, N, M);

      // cudaHostAlloc(&pinned_emomM, emomM.bytes(), cudaHostAllocDefault);
      // cudaHostAlloc(&pinned_emom, emom.bytes(), cudaHostAllocDefault);
      // cudaHostAlloc(&pinned_mmom, mmom.bytes(), cudaHostAllocDefault);
   }
}

// Destructor
CpuRestMeasurement::~CpuRestMeasurement() {
   if(fastCopy) {
      // An outstanding fast-copy callback still reads the pinned buffers and
      // this object, so drain the copy stream before releasing them.
      GPU_STREAM_SYNC(parallel.getCopyStream());
      tmp_emomM.Free();
      tmp_emom.Free();
      tmp_mmom.Free();
      tmp_beff.Free();

      pinned_emomM.FreeHost();
      pinned_emom.FreeHost();
      pinned_mmom.FreeHost();
      pinned_beff.FreeHost();
   }
}

// Callback
void CpuRestMeasurement::queue_callback(void* data) {
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
void CpuRestMeasurement::queueMeasurement(std::size_t mstep) {
   measurementQueue.push(mstep, pinned_emomM.data(), pinned_emom.data(), pinned_mmom.data(), pinned_beff.data(), mmom.size(),  MeasurementQueue::MeasurementType::Rest);
}

// Fast copy and measurement queueing (D -> D, D -> H (async), H -> H)
void CpuRestMeasurement::copyQueueFast(std::size_t mstep) {
   // Timing
   stopwatch.skip();

   // Streams
   GPU_STREAM_T workStream = parallel.getWorkStream();
   GPU_STREAM_T copyStream = parallel.getCopyStream();

   // Create new events
   GpuEventPool::Event& workDone = eventPool.get();
   GpuEventPool::Event& copyDone = eventPool.get();
   GpuEventPool::Event& d2hDone = eventPool.get();

   // The copying must wait for the work stream to finish
   GPU_EVENT_RECORD(workDone.event(), workStream);
   GPU_STREAM_WAIT_EVENT(copyStream, workDone.event(), 0);

   // Async copy in copy stream (device -> temp. device)
   tmp_emomM.copy_async(emomM, copyStream);
   tmp_emom.copy_async(emom, copyStream);
   tmp_mmom.copy_async(mmom, copyStream);
   tmp_beff.copy_async(beff, copyStream);
   GPU_EVENT_RECORD(copyDone.event(), copyStream);
   stopwatch.add("fast - D2D");

   // Then write to host in copy stream (asynchronously with work stream)

   pinned_emomM.copy_async(tmp_emomM, copyStream);
   pinned_emom.copy_async(tmp_emom, copyStream);
   pinned_mmom.copy_async(tmp_mmom, copyStream);
   pinned_beff.copy_async(tmp_beff, copyStream);
   GPU_EVENT_RECORD(d2hDone.event(), copyStream);

   // Release the work stream after the D2D staging copy completes.
   GPU_STREAM_WAIT_EVENT(workStream, copyDone.event(), 0);
   copyDone.addDeactivateCallback(workStream);
   workDone.addDeactivateCallback(workStream);

   // Deactivate d2hDone and queue the consumer after the D2H transfers.
   d2hDone.addDeactivateCallback(copyStream);
   GPU_STREAM_ADD_CALLBACK(copyStream, queue_callback, new queue_callback_data(this, mstep));
}

// Slow copying (D -> H)fortran
void CpuRestMeasurement::copyQueueSlow(std::size_t mstep) {
   // Timing
   stopwatch.skip();

   // Write directly to fortran
   // (this can't be done async, so it will block host until finished)
   fortran_emomM.copy_sync(emomM);
   fortran_mmom.copy_sync(mmom);
   fortran_emom.copy_sync(emom);
   fortran_beff.copy_sync(emom);


   // emomM.write(FortranData::emomM);
   // emom.write(FortranData::emom);
   // mmom.write(FortranData::mmom);
   stopwatch.add("slow - D2H copy");

   // Queue measurement
   measurementQueue.push(mstep, fortran_emomM.data(), fortran_emom.data(), fortran_mmom.data(), fortran_beff.data(), mmom.size(),  MeasurementQueue::MeasurementType::Rest);
}

void CpuRestMeasurement::measure(std::size_t mstep) {
   // Copy required?
   bool copy = (alwaysCopy || call_fortran_do_measurements(mstep));

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

void CpuRestMeasurement::flushMeasurements(std::size_t mstep) {
   // Timing
   stopwatch.skip();

   // Wait out possible queue callbacks. The fast-copy consumer callback is
   // queued on the copy stream (after D2H), so both streams must be drained
   // before the pending measurements are flushed.
   GPU_STREAM_SYNC(parallel.getWorkStream());
   GPU_STREAM_SYNC(parallel.getCopyStream());

   // Flush internal queue
   //measurementQueue.finish();

   // Print remaining measurements
   call_fortran_flush_measurements(mstep);
   stopwatch.add("measure");
}
