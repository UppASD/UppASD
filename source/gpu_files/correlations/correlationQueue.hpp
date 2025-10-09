/*
 * Measurement queue class
 *  Niklas Fejes 2013
 *
 */

#pragma once

#include <pthread.h>

#include <queue>

#include "real_type.h"

class CorrelationQueue {
   // Measurement class
   class Correlation {
      friend CorrelationQueue;
      real* emomM;
      real* emom;
      real* mmom;
      std::size_t step;

   public:
      // Constructor / destructor
      Correlation(real* emomM, real* emom, real* mmom, std::size_t NM, std::size_t step);
      ~Correlation();
   };

   // Private members
   volatile bool finishCorrelations;
   bool processThreadStarted;
   std::queue<Correlation*> q;
   pthread_mutex_t mutex;
   pthread_cond_t cond;
   pthread_t process_thread;

   // Helpers
   static void* process_correlations(void* mqueue);
   void processCorrelations();
   void startProcessThread();
   void finishProcessThread();

public:
   // Constructor / destructor
   CorrelationQueue();
   ~CorrelationQueue();

   // Test if empty
   bool empty();

   // Push a measurement with data to the queue
   void push(std::size_t mstep);
   void push(std::size_t mstep, real* emomM, real* emom, real* mmom, std::size_t NM);

   // Finish
   void finish();
};


