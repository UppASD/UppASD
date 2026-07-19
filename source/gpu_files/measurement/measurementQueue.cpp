/*
 * Measurement queue class
 *  Niklas Fejes 2013
 *
 * TODO allocation limit / handle out of memory
 */

#include "measurementQueue.hpp"

#include <pthread.h>
#include <sched.h>

#include <algorithm>
#include <queue>

#include "c_headers.hpp"
#include "c_helper.h"
#include "fortranData.hpp"
#include "real_type.h"

#ifdef NVPROF
#include <nvToolsExtCuda.h>
#endif

// Measurement class methods
MeasurementQueue::Measurement::Measurement(real* _emomM, real* _emom, real* _mmom, real* _beff, std::size_t NM,
                                           std::size_t _step, MeasurementType _type) {
   step = _step;
   type = _type;
   if(NM != 0) {
      emomM = new fortran_real[NM * 3];
      emom = new fortran_real[NM * 3];
      mmom = new fortran_real[NM];

      std::copy_n(_emomM, NM * 3, emomM);
      std::copy_n(_emom, NM * 3, emom);
      std::copy_n(_mmom, NM, mmom);

      if(type==MeasurementType::Rest){
         beff = new fortran_real[NM * 3];
         std::copy_n(_beff, NM * 3, beff);
      }
      else{
         beff = nullptr;
      }

   } else {
      emomM = emom = mmom = beff = nullptr;
   }
}

MeasurementQueue::Measurement::~Measurement() {
   if(emomM != nullptr) {
      delete[] emomM;
   }
   if(emom != nullptr) {
      delete[] emom;
   }
   if(mmom != nullptr) {
      delete[] mmom;
   }
   if(beff != nullptr) {
      delete[] beff;
   }
}

// MeasurementQueue class methods
void* MeasurementQueue::process_measurements(void* mqueue) {
   ((MeasurementQueue*)mqueue)->processMeasurements();
   pthread_exit(nullptr);
   return nullptr;
}

void MeasurementQueue::processMeasurements() {
#ifdef NVPROF
   nvtxNameOsThread(pthread_self(), "MEASURE_THREAD");
   nvtxRangePush("measuring");
#endif
   while(true) {
      // Lock the mutex
      pthread_mutex_lock(&mutex);

      // Anything to process?
      if(!q.empty()) {
         // Get the front measurement in the queue
         Measurement* m = q.front();
         q.pop();

         // Unlock the mutex (allow stack modification while working)
         pthread_mutex_unlock(&mutex);

         switch(m->type) {
            case MeasurementType::Moment:
               processMomentMeasurement(m);
               break;

            case MeasurementType::Rest:
               processRestMeasurement(m);
               break;
            case MeasurementType::Correlations:
               processCorrelationsMeasurement(m);
               break;
         }

         // Get data pointers
         /*const real* emomM = m->emomM;
         const real* emom = m->emom;
         const real* mmom = m->mmom;

         // Call with default data if no data needed
         if(emomM == nullptr) {
            emomM = FortranData::emomM;
         }
         if(emom == nullptr) {
            emom = FortranData::emom;
         }
         if(mmom == nullptr) {
            mmom = FortranData::mmom;
         }

         // Measure
         call_fortran_measure_moment(emomM, emom, mmom, m->step);
         */
         // Destroy measurement data
         delete m;
      }
      // Exit if the queue is empty, and the finish measurements flag is set
      else if(finishMeasurements) {
         // Unlock the mutex and exit
         pthread_mutex_unlock(&mutex);
         return;
      }
      // Else wait for condition change
      else {
#ifdef NVPROF
         nvtxRangePop();
#endif
         pthread_cond_wait(&cond, &mutex);
#ifdef NVPROF
         nvtxRangePush("measuring");
#endif
         pthread_mutex_unlock(&mutex);
      }
   }
}

void MeasurementQueue::processMomentMeasurement(Measurement* m)
{
   const fortran_real* emomM = (m->emomM) ? m->emomM : FortranData::emomM;
   const fortran_real* emom  = (m->emom)  ? m->emom  : FortranData::emom;
   const fortran_real* mmom  = (m->mmom)  ? m->mmom  : FortranData::mmom;

   call_fortran_measure_moment(emomM, emom, mmom, m->step);
}

void MeasurementQueue::processRestMeasurement(Measurement* m)
{
   const fortran_real* emomM = (m->emomM) ? m->emomM : FortranData::emomM;
   const fortran_real* emom  = (m->emom)  ? m->emom  : FortranData::emom;
   const fortran_real* mmom  = (m->mmom)  ? m->mmom  : FortranData::mmom;
   const fortran_real* beff  = (m->beff)  ? m->beff  : FortranData::beff;

   call_fortran_measure_rest(emomM, emom, mmom, beff, m->step);
}

void MeasurementQueue::processCorrelationsMeasurement(Measurement* m)
{
   const fortran_real* emomM = (m->emomM) ? m->emomM : FortranData::emomM;
   const fortran_real* emom  = (m->emom)  ? m->emom  : FortranData::emom;
   const fortran_real* mmom  = (m->mmom)  ? m->mmom  : FortranData::mmom;

   call_fortran_measure_correlations(emomM, emom, mmom, m->step);
}

void MeasurementQueue::startProcessThread() {
   if(!processThreadStarted) {
      // Flag started
      processThreadStarted = true;

      // Joinable thread
      pthread_attr_t attr;
      pthread_attr_init(&attr);
      pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

      // Set low priority for measurement thread
      struct sched_param sparam;
      sparam.sched_priority = sched_get_priority_min(SCHED_RR);
      pthread_attr_setschedpolicy(&attr, SCHED_RR);
      pthread_attr_setschedparam(&attr, &sparam);

      // Mutex and cond
      pthread_mutex_init(&mutex, nullptr);
      pthread_cond_init(&cond, nullptr);

      // Create
      pthread_create(&process_thread, &attr, process_measurements, this);
   }
}

void MeasurementQueue::finishProcessThread() {
   if(!finishMeasurements && processThreadStarted) {
      // Chage finish flag and send cond signal
      pthread_mutex_lock(&mutex);
      finishMeasurements = true;
      pthread_cond_signal(&cond);
      pthread_mutex_unlock(&mutex);

      // Join threads
      pthread_join(process_thread, nullptr);

      // Destroy mutex and cond
      pthread_mutex_destroy(&mutex);
      pthread_cond_destroy(&cond);

      // Flag that the thread is not running
      processThreadStarted = false;
   }
}

// Constructor / destructor
MeasurementQueue::MeasurementQueue() {
   finishMeasurements = false;
   processThreadStarted = false;
}

MeasurementQueue::~MeasurementQueue() {
   if(!finishMeasurements) {
      std::fprintf(stderr, "MeasurementQueue::finish() not called!\n");
   }
   finishProcessThread();
}

// Test if empty
bool MeasurementQueue::empty() {
   pthread_mutex_lock(&mutex);
   bool e = q.empty();
   pthread_mutex_unlock(&mutex);
   return e;
}

// Push a measurement with data to the queue
void MeasurementQueue::push(std::size_t mstep) {
   push(mstep, nullptr, nullptr, nullptr, nullptr, 0, MeasurementType::Moment);
}

void MeasurementQueue::push(std::size_t mstep, real* emomM,
                            real* emom,
                            real* mmom,
                            std::size_t NM)
{
   push(mstep, emomM, emom, mmom, nullptr, NM, MeasurementType::Moment);
}
void MeasurementQueue::push(std::size_t mstep,
                            real* emomM,
                            real* emom,
                            real* mmom,
                            real* beff,
                            std::size_t NM,
                            MeasurementType type) {
   // Finishing?
   if(finishMeasurements) {
      std::fprintf(stderr, "MeasurementQueue::push() called after finish()!\n");
      return;
   }

   // Initial push?
   if(!processThreadStarted) {
      // Push measurement to queue
      q.push(new Measurement(emomM, emom, mmom, beff, NM, mstep, type));

      // Start process thread
      startProcessThread();
   } else {
      // Lock mutex
      pthread_mutex_lock(&mutex);

      // Push measurement to queue
      q.push(new Measurement(emomM, emom, mmom, beff, NM, mstep, type));

      // Signal condition change if thread started
      if(processThreadStarted) {
         pthread_cond_signal(&cond);
      }

      // Unlock mutex
      pthread_mutex_unlock(&mutex);
   }
}

// Finish
void MeasurementQueue::finish() {
   if(finishMeasurements) {
      std::fprintf(stderr, "MeasurementQueue::finish() called multiple times!\n");
   }
   finishProcessThread();
}
