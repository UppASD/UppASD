#include <cstddef>
#include <cmath>

#include <cuda.h>
#include <curand.h>
#include <iostream>

using namespace std;

#include "real_type.h"
#include "printDebug.hpp"

#include "cudaDepondtIntegrator.hpp"

#include "cudaCommon.hpp"

#include "fortMatrix.hpp"
#include "cudaMatrix.hpp"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "stopwatchPool.hpp"


////////////////////////////////////////////////////////////////////////////////
// Parallelization helper classes
////////////////////////////////////////////////////////////////////////////////

class CudaDepondtIntegrator::BuildEffectiveField :
	public CudaParallelizationHelper::Atom {
private:
	real *       bdup;
	const real * blocal;
	const real * emom;
	real         damping;
public:
	BuildEffectiveField(real * p1, const real * p2, const real * p3,
			real p4) {
		bdup    = p1;
		blocal  = p2;
		emom    = p3;
		damping = p4;
	}

        __device__ void each(unsigned int atom) {
		real       * my_bdup = &bdup  [atom * 3];
		const real * my_bloc = &blocal[atom * 3];
		const real * my_emom = &emom  [atom * 3];
		my_bdup[0] = my_bloc[0] + damping*(my_emom[1]*my_bloc[2] - my_emom[2]*my_bloc[1]);
		my_bdup[1] = my_bloc[1] + damping*(my_emom[2]*my_bloc[0] - my_emom[0]*my_bloc[2]);
		my_bdup[2] = my_bloc[2] + damping*(my_emom[0]*my_bloc[1] - my_emom[1]*my_bloc[0]);
	}
};

class CudaDepondtIntegrator::Rotate :
	public CudaParallelizationHelper::Atom {
private:
	real *       mrod;
	const real * emom;
	const real * bdup;
	real         timestep;
	real         gamma;
	real 		 damping;
public:
	Rotate(real * p1, const real * p2, const real * p3,
			real p4, real p5, real p6) {
		mrod     = p1;
		emom     = p2;
		bdup     = p3;
		timestep = p4;
		gamma    = p5;
		damping  = p6;
	}

        __device__ void each(unsigned int atom) {
		real       * my_mrod = &mrod[atom * 3];
		const real * my_emom = &emom[atom * 3];
		const real * my_bdup = &bdup[atom * 3];

		// Get effective field components and size
		real x = my_bdup[0];
		real y = my_bdup[1];
		real z = my_bdup[2];
		real norm = sqrt(x*x + y*y + z*z);

		// Normalize components
		x /= norm;
		y /= norm;
		z /= norm;

		// Get precession angle
		real angle = norm * timestep * gamma;

		angle *= (real)1.0/((real)1.0+damping*damping);

		// Calculate sin and cos
		real cosv, sinv;
		sincos(angle, &sinv, &cosv);
		real u = 1 - cosv;

		// Calculate matrix
		real M[3][3];
		M[0][0]=x*x*u+  cosv; M[0][1]=x*y*u-z*sinv; M[0][2]=x*z*u+y*sinv;
		M[1][0]=x*y*u+z*sinv; M[1][1]=y*y*u+  cosv; M[1][2]=y*z*u-x*sinv;
		M[2][0]=x*z*u-y*sinv; M[2][1]=z*y*u+x*sinv; M[2][2]=z*z*u+  cosv;

		// Rotate
		real mx = my_emom[0], my = my_emom[1], mz = my_emom[2];
		my_mrod[0] = mx*M[0][0] + my*M[0][1] + mz*M[0][2];
		my_mrod[1] = mx*M[1][0] + my*M[1][1] + mz*M[1][2];
		my_mrod[2] = mx*M[2][0] + my*M[2][1] + mz*M[2][2];
// Alternative
//		real mx = my_emom[0], my = my_emom[1], mz = my_emom[2];
//		my_mrod[0] = mx*x*x*u+  cosv + my*y*x*u-z*sinv + mz*z*x*u+y*sinv;
//		my_mrod[1] = mx*x*y*u+z*sinv + my*y*y*u+  cosv + mz*z*y*u-x*sinv;
//		my_mrod[2] = mx*x*z*u-y*sinv + my*y*z*u+x*sinv + mz*z*z*u+  cosv;
	}
};



////////////////////////////////////////////////////////////////////////////////
// Class members
////////////////////////////////////////////////////////////////////////////////

// Constructor
CudaDepondtIntegrator::CudaDepondtIntegrator() :
		stopwatch(GlobalStopwatchPool::get("Cuda Depondt integrator")),
		parallel(CudaParallelizationHelper::def) {
}

// Destructor
CudaDepondtIntegrator::~CudaDepondtIntegrator() {
	release();
}


// Initiator
bool CudaDepondtIntegrator::initiate(size_t N, size_t M, char _stt, real _timestep, 
		curandRngType_t rng, unsigned long long seed) {
	// Assert that we're not already initialized
	release();

	// Param
	stt      = _stt;
	timestep = _timestep;

	// Initiate thermfield
	if (!thermfield.initiate(N, M, rng, seed)) {
		release();
		return false;
	}

	// Allocate device matrices
	mrod  .initiate(3, N, M);
	blocal.initiate(3, N, M);
	bdup  .initiate(3, N, M);

	// All initialized?
	if (cudaDeviceSynchronize() != cudaSuccess) {
		release();
		return false;
	}

	return true;
}

bool CudaDepondtIntegrator::initiateConstants(const fortMatrix<real,1> &temperature, 
		real timestep,
		real gamma_const, real k_bolt_const, 
		real mub_const, real damping_const) {
	// Set parameters
	gamma     = gamma_const;
	k_bolt    = k_bolt_const;
	mub       = mub_const;
	damping   = damping_const;
	dp_factor = (2.0*damping*k_bolt) / (gamma*mub*(1+damping*damping));

	// Initiate thermfield constants
	if (!thermfield.initiateConstants(temperature, timestep, gamma, k_bolt, mub, damping))
		return false;

	return true;
}

// Releaser
void CudaDepondtIntegrator::release() {
	mrod  .free();
	blocal.free();
	bdup  .free();
}


// First step of Depond solver, calculates the stochastic field and rotates the
// magnetic moments according to the effective field 
// Dupont recipe J. Phys.: Condens. Matter 21 (2009) 336005
void CudaDepondtIntegrator::evolveFirst(const cudaMatrix<real,3,3> &beff,
		cudaMatrix<real,3,3> &b2eff, 
		const cudaMatrix<real,3,3> &btorque,
		cudaMatrix<real,3,3> &emom, 
		cudaMatrix<real,3,3> &emom2,
		cudaMatrix<real,3,3> &emomM, 
		const cudaMatrix<real,2> &mmom) {

	// Timing
	stopwatch.skip();

	//_dpr;

	// Randomize stochastic field
	thermfield.randomize(mmom);
	stopwatch.add("thermfield");

	CudaCommon::Add cm(blocal, beff, thermfield.getField());
	
	//_dpr;
	// Construct local field
	parallel.cudaElementCall(CudaCommon::Add(blocal, beff, thermfield.getField()));
	stopwatch.add("localfield");

	//_dpr;
	// Construct effective field (including damping term)
	buildbeff(emom, btorque);
	stopwatch.add("buildbeff");

	//_dpr;
	// Set up rotation matrices and perform rotations
	rotate(emom, timestep);
	stopwatch.add("rotate");

	// copy m(t) to emom2 and m(t+dt) to emom for heisge, save b(t)
	parallel.cudaElementCall(CudaCommon::ScalarMult(emomM, mrod, mmom));
	emom2.swap(emom); // Previous emom will not be needed
	emom.swap(mrod);  // Previous mrod will not be needed
	b2eff.swap(bdup); // Previous bdup will not be needed
	stopwatch.add("copy");
}


// Second step of Depond solver, calculates the corrected effective field from 
// the predicted effective fields. Rotates the moments in the corrected field
void CudaDepondtIntegrator::evolveSecond(
		const cudaMatrix<real,3,3> &beff,
		const cudaMatrix<real,3,3> &b2eff, 
		const cudaMatrix<real,3,3> &btorque,
		cudaMatrix<real,3,3> &emom, 
		cudaMatrix<real,3,3> &emom2) {

	// Timing
	stopwatch.skip();

	// Construct local field
	parallel.cudaElementCall(CudaCommon::Add(blocal, beff, thermfield.getField()));
	stopwatch.add("localfield");

	// Construct effective field (including damping term)
	buildbeff(emom, btorque);
	stopwatch.add("buildbeff");

	// Corrected field
	parallel.cudaElementCall(CudaCommon::Avg(bdup, b2eff));
	emom.swap(emom2); // Vaild as emom2 wont be used after its coming swap
	stopwatch.add("corrfield");

	// Final rotation
	rotate(emom, timestep);
	stopwatch.add("rotate");

	// Swap
	emom2.swap(mrod); // Vaild as mrod wont be needed after this
	stopwatch.add("copy");
}


void CudaDepondtIntegrator::rotate(const cudaMatrix<real,3,3> &emom, real delta_t) {
	parallel.cudaAtomCall(Rotate(mrod, emom, bdup, timestep, gamma, damping));
}

// Constructs the effective field (including damping term)
void CudaDepondtIntegrator::buildbeff(const cudaMatrix<real,3,3> &emom,
	const cudaMatrix<real,3,3> &btorque) {
	parallel.cudaAtomCall(BuildEffectiveField(bdup, blocal, emom, damping));

	// TODO untested
	if (stt != 'N')
		parallel.cudaElementCall(CudaCommon::AddTo(bdup, btorque));
}


