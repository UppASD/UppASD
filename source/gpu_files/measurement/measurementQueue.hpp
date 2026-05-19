

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
public:
   enum class MeasurementType {
      Moment,
      Rest,
      Correlations,
      All//TODO add all option on fortran side 
   };

   // Measurement class
   class Measurement {
      friend MeasurementQueue;
      real* emomM;
      real* emom;
      real* mmom;
      real* beff;
      std::size_t step;
      MeasurementType type; 

   public:
      // Constructor / destructor
      Measurement(real* emomM, real* emom, real* mmom, real* beff, std::size_t NM, std::size_t step, MeasurementType type);
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
   void processMomentMeasurement(Measurement* m);
   void processRestMeasurement(Measurement* m);
   void processCorrelationsMeasurement(Measurement* m);


public:
   // Constructor / destructor
   MeasurementQueue();
   ~MeasurementQueue();

   // Test if empty
   bool empty();

   // Push a measurement with data to the queue
   void push(std::size_t mstep);
   void push(std::size_t mstep, real* emomM, real* emom, real* mmom, std::size_t NM); //TODO: remove later
   void push(std::size_t mstep, real* emomM, real* emom, real* mmom, real* beff, std::size_t NM, MeasurementType type);

   // Finish
   void finish();
};


