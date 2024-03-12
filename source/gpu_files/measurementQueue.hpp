/*
 * Measurement queue class
 *  Niklas Fejes 2013
 *
 */

#pragma once

#include <pthread.h>

#include <queue>

#include "real_type.h"

class MeasurementQueue {
   // Measurement class
   class Measurement {
      friend MeasurementQueue;
      real* emomM;
      real* emom;
      real* mmom;
      size_t step;

   public:
      // Constructor / destructor
      Measurement(real* emomM, real* emom, real* mmom, size_t NM, size_t step);
      ~Measurement();
   };

   // Private members
   volatile bool finishMeasurements;
   bool processThreadStarted;
   std::queue<Measurement*> q;
   pthread_mutex_t mutex;
   pthread_cond_t cond;
   pthread_t process_thread;

   // Helpers
   static void* process_measurements(void* mqueue);
   void processMeasurements();
   void startProcessThread();
   void finishProcessThread();

public:
   // Constructor / destructor
   MeasurementQueue();
   ~MeasurementQueue();

   // Test if empty
   bool empty();

   // Push a measurement with data to the queue
   void push(size_t mstep);
   void push(size_t mstep, real* emomM, real* emom, real* mmom, size_t NM);

   // Finish
   void finish();
};


