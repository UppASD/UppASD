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
      usd_int step;

   public:
      // Constructor / destructor
      Measurement(real* emomM, real* emom, real* mmom, usd_int NM, usd_int step);
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
   void push(usd_int mstep);
   void push(usd_int mstep, real* emomM, real* emom, real* mmom, usd_int NM);

   // Finish
   void finish();
};


