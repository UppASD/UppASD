/*
 * Measurement queue class
 *  Niklas Fejes 2013
 *
 * TODO allocation limit / handle out of memory
 */

#include <cstdio>
#include <cstring>
#include <queue>
#include <pthread.h>
#include <sched.h>

#include "real_type.h"

#include "c_helper.h"

#include "measurementQueue.hpp"
#include "fortranData.hpp"

#ifdef NVPROF
#include <nvToolsExtCuda.h>
#endif

// Measurement class methods
MeasurementQueue::Measurement::Measurement(real * _emomM, real * _emom, real * _mmom,
		size_t NM, size_t _step) {
	step = _step;
	if (NM != 0) {
		emomM = new real[NM * 3];
		emom  = new real[NM * 3];
		mmom  = new real[NM * 1];
		memcpy(emomM, _emomM, NM * 3 * sizeof(real));
		memcpy(emom,  _emom,  NM * 3 * sizeof(real));
		memcpy(mmom,  _mmom,  NM * 1 * sizeof(real));
	}
	else {
		emomM = emom = mmom = NULL;
	}
}

MeasurementQueue::Measurement::~Measurement() {
	if (emomM != NULL) delete emomM;
	if (emom  != NULL) delete emom;
	if (mmom  != NULL) delete mmom;
}


// MeasurementQueue class methods
void * MeasurementQueue::process_measurements(void * mqueue) {
	((MeasurementQueue*)mqueue)->processMeasurements();
	pthread_exit(NULL);
	return NULL;
}

void MeasurementQueue::processMeasurements() {
#ifdef NVPROF
	nvtxNameOsThread(pthread_self(), "MEASURE_THREAD");
	nvtxRangePush("measuring");
#endif
	while (true) {
		// Lock the mutex
		pthread_mutex_lock(&mutex);

		// Anything to process?
		if (!q.empty()) {
			// Get the front measurement in the queue
			Measurement * m = q.front();
			q.pop();

			// Unlock the mutex (allow stack modification while working)
			pthread_mutex_unlock(&mutex);

			// Get data pointers
			const real * emomM = m->emomM;
			const real * emom  = m->emom;
			const real * mmom  = m->mmom;

			// Call with default data if no data needed
			if (emomM == NULL) emomM = FortranData::emomM;
			if (emom  == NULL) emom  = FortranData::emom;
			if (mmom  == NULL) mmom  = FortranData::mmom;

			// Measure
			fortran_measure_moment(emomM, emom, mmom, m->step);

			// Destroy measurement data
			delete m;
		}
		// Exit if the queue is empty, and the finish measurements flag is set
		else if (finishMeasurements) {
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

void MeasurementQueue::startProcessThread() {
	if (!processThreadStarted) {
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
		pthread_mutex_init(&mutex, NULL);
		pthread_cond_init(&cond, NULL);

		// Create
		pthread_create(&process_thread, &attr, process_measurements, this);
	}
}

void MeasurementQueue::finishProcessThread() {
	if (!finishMeasurements && processThreadStarted) {
		// Chage finish flag and send cond signal 
		pthread_mutex_lock(&mutex);
		finishMeasurements = true;
		pthread_cond_signal(&cond);
		pthread_mutex_unlock(&mutex);

		// Join threads
		pthread_join(process_thread, NULL);

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
	if (!finishMeasurements) fprintf(stderr, "MeasurementQueue::finish() not called!\n");
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
void MeasurementQueue::push(size_t mstep) { push(mstep, NULL, NULL, NULL, 0); }
void MeasurementQueue::push(size_t mstep, real * emomM, real * emom, real * mmom, size_t NM) {
	// Finishing?
	if (finishMeasurements) {
		fprintf(stderr, "MeasurementQueue::push() called after finish()!\n");
		return;
	}

	// Initial push?
	if (!processThreadStarted) {
		// Push measurement to queue
		q.push(new Measurement(emomM, emom, mmom, NM, mstep));

		// Start process thread
		startProcessThread();
	}
	else {
		// Lock mutex
		pthread_mutex_lock(&mutex);

		// Push measurement to queue
		q.push(new Measurement(emomM, emom, mmom, NM, mstep));

		// Signal condition change if thread started
		if (processThreadStarted)
			pthread_cond_signal(&cond);

		// Unlock mutex
		pthread_mutex_unlock(&mutex);
	}

}

// Finish
void MeasurementQueue::finish() {
	if (finishMeasurements) fprintf(stderr, "MeasurementQueue::finish() called multiple times!\n");
	finishProcessThread();
}

