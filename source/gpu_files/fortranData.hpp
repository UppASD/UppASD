#pragma once

#include <cstddef>
#include "c_headers.hpp"

#include "real_type.h"

// NAME         TYPE     DIMENSION   DESCRIPTION 
//
// mompar       int       1          Parametrization of magnetic moment magnitudes (0=no)
// initexc      char      1          Mode of excitation of initial magnetic moments (I=vacancies, R=two magnon Raman, F=no)
//
// emom         real     (3,N,M)     Current unit moment vector
// emom2        real     (3,N,M)     Final (or temporary) unit moment vector
// emomM        real     (3,N,M)     Current magnetic moment vector
// mmom         real     (N,M)       Magnitude of magnetic moments
// mmom0        real     (N,M)       Starting magnitude of magnetic moments
// mmom2        real     (N,M)       Temporary value of magnitude of magnetic moments
// mmomi        real     (N,M)       Inverse of magnitude of magnetic moments
//
// mrod         real     (3,N,M)     Rotated magnetic moments
// btherm       real     (3,N,M)     Thermal stochastic field
// bloc         real     (3,N,M)     Local effective field
// bdup         real     (3,N,M)     Resulting effective field
//
// beff         real     (3,N,M)     Total effective field from application of Hamiltonian
// b2eff        real     (3,N,M)     Temporary storage of magnetic field
// btorque      real     (3,N,M)     Spin transfer torque
// emom         real     (3,N,M)     Current unit moment vector
// emom2        real     (3,N,M)     Final (or temporary) unit moment vector
// emomM        real     (3,N,M)     Current magnetic moment vector
// mmom         real     (N,M)       Magnitude of magnetic moments
// delta_t      real      1          Time step
// temperature  real     (N)         Temperature
//
// stt          char      1          Method to handle spin transfer torque 
// sb 			real 	 (N)	     Ratio between cubic and uniaxial anisotropy

class FortranData {
public:
	// Scalars
	static char * stt;
	static int *  SDEalgh;

	static unsigned int * rstep;
	static unsigned int * nstep;
	static unsigned int * Natom;
	static unsigned int * Mensemble;
	static unsigned int * max_no_neigh;

	static real * delta_t;
	static real * gamma;
	static real * k_bolt;
	static real * mub;
	static real * damping;

	static real * binderc;
	static real * mavg;

	static int *  mompar;
	static char * initexc;

	static unsigned int * do_dm;
	static unsigned int * max_no_dmneigh;

	static unsigned int * do_jtensor; // Information on weather the exchange coupling tensor should be used or not
	static unsigned int * do_aniso; // Information on weather the anisotropy should be used or not

	// Matrices / vectors
	static real *         ncoup;
	static unsigned int * nlist;
	static unsigned int * nlistsize;

	static real *         dmvect;
	static unsigned int * dmlist;
	static unsigned int * dmlistsize;
	
	static real * j_tensor;

	static real * kaniso;
	static real * eaniso;
	static unsigned int * taniso;
	static real * sb;

	static real * beff;
	static real * b2eff;
	static real * emomM;
	static real * emom;
	static real * emom2;
	static real * external_field;
	static real * mmom;
	static real * btorque;
	static real * temperature;
	static real * mmom0;
	static real * mmom2;
	static real * mmomi;               

	// Input
	static int * gpu_mode;
	static int * gpu_rng;
	static int * gpu_rng_seed;

	// Initiators
	static void setConstantPointers(char * p1, int * p2, unsigned int * p3,
		unsigned int * p4, unsigned int * p5, unsigned int * p6, unsigned int * p7, real * p8,
		real * p9, real * p10, real * p11, real * p12,
		real * p13, real * p14, int * p15, char * p16, unsigned int * p17, unsigned int * p18, unsigned int * p19, unsigned int * p20);

	static void setMatrixPointers(real * p1, unsigned int * p2, unsigned int * p3, real * p4,
		real * p5, real * p6, real * p7, real * p8, real * p9, 
		real * p10, real * p11, real * p12, real * p13, 
		real * p14, real * p15, real * p16, unsigned int * p17, unsigned int * p18, real * p19, real * p20, real * p21, unsigned int * p22, real * p23);

	static void setInputDataPointers(int * p1, int * p2, int * p3);

};



