#include <cstddef>
#include <cstring>
#include <cmath>

using namespace std;

#include "real_type.h"

#include "depondtIntegrator.hpp"

#include "randomnum.hpp"
#include "matrix.hpp"
#include "fortMatrix.hpp"
#include "stopwatch.hpp"
#include "stopwatchPool.hpp"

#include "thermfield.hpp"


// Constructor
DepondtIntegrator::DepondtIntegrator() :
		stopwatch(GlobalStopwatchPool::get("Depondt integrator")) {
	Natom = Mensemble = 0;
}

// Destructor
DepondtIntegrator::~DepondtIntegrator() {
	release();
}


// Initiator
bool DepondtIntegrator::initiate(size_t N, size_t M, char sttMode) {
	// Assert that we're not already initialized
	release();

	// Set sizes
	Natom     = N;
	Mensemble = M;

	// Param
	stt = sttMode;

	// Allocate "fake" fortran matrices
	size_t size = 3 * Natom * Mensemble;
	mrod  .set(new real[size], 3, Natom, Mensemble);
	blocal.set(new real[size], 3, Natom, Mensemble);
	bdup  .set(new real[size], 3, Natom, Mensemble);

	// All allocated?
	if (!mrod.has_data() || !blocal.has_data() || !bdup.has_data()) {
		release();
		return false;
	}

	// Initiate thermal field
	if (!tfield.initiate(N, M)) {
		release();
		return false;
	}		

	return true;
}


bool DepondtIntegrator::initiateConstants(real gamma_const, real k_bolt, 
		real mub, real damping_const, const hostMatrix<real,1> &temp, 
		real timestep_const) {
	// Set parameters
	gamma     = gamma_const;
	damping   = damping_const;
	timestep  = timestep_const;
	tfield.initiateConstants(temp, timestep, gamma, k_bolt, mub, damping);
	return true;
}


// Releaser
void DepondtIntegrator::release() {
	// Reset parameters
	Natom = Mensemble = 0;


	// Free data
	const real * data;
	data = mrod.get_data();   if (data != NULL) delete data;
	data = blocal.get_data(); if (data != NULL) delete data;
	data = bdup.get_data();   if (data != NULL) delete data;

	// Reset matrices
	mrod.set(NULL, 3, 1, 1);
	blocal.set(NULL, 3, 1, 1);
	bdup.set(NULL, 3, 1, 1);
}



// First step of Depond solver, calculates the stochastic field and rotates the
// magnetic moments according to the effective field 
void DepondtIntegrator::evolveFirst(
		const hostMatrix<real,3,3> &beff, 
		hostMatrix<real,3,3> &b2eff, 
		const hostMatrix<real,3,3> &btorque,
		hostMatrix<real,3,3> &emom, 
		hostMatrix<real,3,3> &emom2,
		hostMatrix<real,3,3> &emomM,
		const hostMatrix<real,2> &mmom) {
	// beff        - Total effective field from application of Hamiltonian
	// b2eff       - Temporary storage of magnetic field
	// btorque     - Spin transfer torque
	// emom        - Current unit moment vector
	// emom2       - Final (or temporary) unit moment vector
	// emomM       - Current magnetic moment vector
	// mmom        - Magnitude of magnetic moments

	// Dupont recipe J. Phys.: Condens. Matter 21 (2009) 336005

	// Timing
	stopwatch.skip();

	// Calculate stochastic field
	tfield.randomize(mmom);
	stopwatch.add("thermfield");

	// Construct local field
	const hostMatrix<real,3,3> &btherm = tfield.getField();
	#pragma omp parallel for collapse(2)
	for (size_t k = 0; k < Mensemble; k++) {
		for (size_t i = 0; i < Natom; i++) {
			blocal(0,i,k) = beff(0,i,k) + btherm(0,i,k);
			blocal(1,i,k) = beff(1,i,k) + btherm(1,i,k);
			blocal(2,i,k) = beff(2,i,k) + btherm(2,i,k);
		}
	}
	stopwatch.add("localfield");

	// Construct effective field (including damping term)
	buildbeff(emom, btorque);
	stopwatch.add("buildbeff");

	// Set up rotation matrices and perform rotations
	rotate(emom, timestep);
	stopwatch.add("rotate");

	// copy m(t) to emom2 and m(t+dt) to emom for heisge, save b(t)
	#pragma omp parallel for collapse(2)
	for (size_t k = 0; k < Mensemble; k++) {
		for (size_t i = 0; i < Natom; i++) {
			real m = mmom(i,k);
			for (int j = 0; j < 3; j++) {
				emom2(j,i,k) = emom(j,i,k);
				emomM(j,i,k) = mrod(j,i,k) * m;
				emom(j,i,k)  = mrod(j,i,k);
				b2eff(j,i,k) = bdup(j,i,k);
			}
		}
	}
	stopwatch.add("copy");
}

// Second step of Depond solver, calculates the corrected effective field from 
// the predicted effective fields. Rotates the moments in the corrected field
void DepondtIntegrator::evolveSecond(
		const hostMatrix<real,3,3> &beff, 
		const hostMatrix<real,3,3> &b2eff, 
		const hostMatrix<real,3,3> &btorque, 
		hostMatrix<real,3,3> &emom, 
		hostMatrix<real,3,3> &emom2) {
	// Timing
	stopwatch.skip();

	// Construct local field
	const hostMatrix<real,3,3> &btherm = tfield.getField();
	#pragma omp parallel for collapse(2)
	for (size_t k = 0; k < Mensemble; k++) {
		for (size_t i = 0; i < Natom; i++) {
			blocal(0,i,k) = beff(0,i,k) + btherm(0,i,k);
			blocal(1,i,k) = beff(1,i,k) + btherm(1,i,k);
			blocal(2,i,k) = beff(2,i,k) + btherm(2,i,k);
		}
	}
	stopwatch.add("localfield");

	// Construct effective field (including damping term)
	buildbeff(emom, btorque);
	stopwatch.add("buildbeff");

	// Corrected field
	#pragma omp parallel for collapse(2)
	for (size_t k = 0; k < Mensemble; k++) {
		for (size_t i = 0; i < Natom; i++) {
			for (int j = 0; j < 3; j++) {
				bdup(j,i,k) = 0.5*bdup(j,i,k) + 0.5*b2eff(j,i,k);
				emom(j,i,k) = emom2(j,i,k);
			}
		}
	}
	stopwatch.add("corrfield");

	// Final rotation
	rotate(emom, timestep);
	stopwatch.add("rotate");

	// Copy
	emom2.memcopy(mrod);
	stopwatch.add("copy");
}


// Performs a Rodrigues rotation of the magnetic moments in the 
// effective field.
bool DepondtIntegrator::rotate(const hostMatrix<real,3,3> &emom, real timestep) {
	// Initiated?
	if (Natom == 0) return false;

	// Rotate
	#pragma omp parallel for collapse(2)
	for (size_t k = 0; k < Mensemble; k++) {
		for (size_t i = 0; i < Natom; i++) {
			// Get effective field components and size
			real x = bdup(0,i,k);
			real y = bdup(1,i,k);
			real z = bdup(2,i,k);
			real norm = sqrt(x*x + y*y + z*z);

			// Normalize components
			x /= norm;
			y /= norm;
			z /= norm;

			// Get precession angle
			real angle = norm * timestep * gamma;

			// Calculate sin(angle) / cosine(angle)
			real cosv = cos(angle);
			real sinv = sin(angle);
			real u = 1 - cosv;

			// Calculate matrix
			real M[3][3];
			M[0][0]=x*x*u+  cosv; M[0][1]=y*x*u-z*sinv; M[0][2]=z*x*u+y*sinv;
			M[1][0]=x*y*u+z*sinv; M[1][1]=y*y*u+  cosv; M[1][2]=z*y*u-x*sinv;
			M[2][0]=x*z*u-y*sinv; M[2][1]=y*z*u+x*sinv; M[2][2]=z*z*u+  cosv;

			// Rotate
			real mx = emom(0,i,k), my = emom(1,i,k), mz = emom(2,i,k);
			mrod(0,i,k) = mx*M[0][0] + my*M[0][1] + mz*M[0][2];
			mrod(1,i,k) = mx*M[1][0] + my*M[1][1] + mz*M[1][2];
			mrod(2,i,k) = mx*M[2][0] + my*M[2][1] + mz*M[2][2];
//			real mx = emom(0,i,k), my = emom(1,i,k), mz = emom(2,i,k);
//			mrod(0,i,k) = mx*x*x*u+  cosv + my*y*x*u-z*sinv + mz*z*x*u+y*sinv;
//			mrod(1,i,k) = mx*x*y*u+z*sinv + my*y*y*u+  cosv + mz*z*y*u-x*sinv;
//			mrod(2,i,k) = mx*x*z*u-y*sinv + my*y*z*u+x*sinv + mz*z*z*u+  cosv;
		}
	}
	return true;
}

// Constructs the effective field (including damping term)
void DepondtIntegrator::buildbeff(const hostMatrix<real,3,3> &emom, const hostMatrix<real,3,3> &btorque) {

	#pragma omp parallel for collapse(2)
	for (size_t k = 0; k < Mensemble; k++) {
		for (size_t i = 0; i < Natom; i++) {
			bdup(0,i,k)=blocal(0,i,k)+damping*(emom(1,i,k)*blocal(2,i,k)-emom(2,i,k)*blocal(1,i,k));
			bdup(1,i,k)=blocal(1,i,k)+damping*(emom(2,i,k)*blocal(0,i,k)-emom(0,i,k)*blocal(2,i,k));
			bdup(2,i,k)=blocal(2,i,k)+damping*(emom(0,i,k)*blocal(1,i,k)-emom(1,i,k)*blocal(0,i,k));
		}
	}

	if (stt != 'N') {
		#pragma omp parallel for collapse(2)
		for (size_t k = 0; k < Natom; k++) {
			for (size_t i = 0; i < Mensemble; i++) {
				bdup(0,i,k) += btorque(0,i,k);
				bdup(1,i,k) += btorque(1,i,k);
				bdup(2,i,k) += btorque(2,i,k);
			}
		}
	}
}

