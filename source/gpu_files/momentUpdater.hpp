#ifndef __MOMENT_UPDATER_HPP__
#define __MOMENT_UPDATER_HPP__

// Routines for updating magnetic moment after time evolution


#include "real_type.h"

#include "fortMatrix.hpp"
#include "stopwatch.hpp"
#include "stopwatchPool.hpp"

// Name      Type     Dimension  Intent  Description 
// emom      double   (3,N,M)    (out)   Current unit moment vector
// emom2     double   (3,N,M)    (in)    Final (or temporary) unit moment vector
// emomM     double   (3,N,M)    (out)   Current magnetic moment vector
// mmom      double   (N,M)      (out)   Magnitude of magnetic moments
// mmom0     double   (N,M)      (in)    Starting magnitude of magnetic moments
// mmom2     double   (N,M)      (out)   Temporary value of magnitude of magnetic moments
// mmomi     double   (N,M)      (out)   Inverse of magnitude of magnetic moments
// mompar:   int      1          (in)    Parametrization of magnetic moment magnitudes (0=no)
// initexc   char     1          (in)    Mode of excitation of initial magnetic moments (I=vacancies, R=two magnon Raman, F=no)


class MomentUpdater {
private:
	// Moments to update
	fortMatrix<real,2>   &mmom;
	fortMatrix<real,2>   &mmom0;
	fortMatrix<real,2>   &mmom2; 
	fortMatrix<real,3,3> &emom;
	fortMatrix<real,3,3> &emom2;
	fortMatrix<real,3,3> &emomM;
	fortMatrix<real,2>   &mmomi;

	// Parameters
	int mompar;
	char initexc;

	// Helpers
	void calcMoment();
	// Short
	// alt 0: mmom2 = mmom
	// alt 1: mmom2 = mmom0 * abs(emom_z)
	// alt 2: mmom2 = mmom0 * abs(emom_z^2)
	void copyMoment();
	// Short
	//   emom  = emom2
	//   mmom  = mmom2
	//   mmomi = 1/mmom2
	//   emomM = emom2 * mmom2

	// Timer
	Stopwatch &stopwatch;

public:

	// Constructor
	MomentUpdater(fortMatrix<real,2> &p1, fortMatrix<real,2> &p2, fortMatrix<real,2> &p3, 
		fortMatrix<real,3,3> &p4, fortMatrix<real,3,3> &p5, fortMatrix<real,3,3> &p6, 
		fortMatrix<real,2>  &p7, int p8, char p9) : 
			mmom(p1), mmom0(p2), mmom2(p3), emom(p4), emom2(p5), emomM(p6), mmomi(p7), 
			mompar(p8), initexc(p9), stopwatch(GlobalStopwatchPool::get("Moment")) {}

	// Updater
	void update();

};


#endif

