#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <cuda.h>
#include <curand.h>

using namespace std;

#include "real_type.h"
#include "cudaGPUErrchk.hpp"

#include "cudaMdSimulation.hpp"

#include "cudaMatrix.hpp"
#include "fortMatrix.hpp"

#include "cudaParallelizationHelper.hpp"
#include "cudaDepondtIntegrator.hpp"
#include "cudaHamiltonianCalculations.hpp"
#include "cudaMomentUpdater.hpp"
#include "cudaMeasurement.hpp"

#include "fortranData.hpp"
#include "c_helper.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "stopwatchPool.hpp"

CudaMdSimulation::CudaMdSimulation() {
	isInitiated = false;
}
CudaMdSimulation::~CudaMdSimulation() {
	release();
}

void CudaMdSimulation::initiateConstants() {
	// Only heisge_jij allowed
	SDEalgh = *FortranData::SDEalgh;
	if(!(SDEalgh == 1 || SDEalgh == 4 || SDEalgh == 5 || SDEalgh == 11)) {
		fprintf(stderr, "Invalid SDEalgh!\n");
		exit(EXIT_FAILURE);
	}

	// Constants
	stt            = *FortranData::stt;
	rstep          = *FortranData::rstep;
	nstep          = *FortranData::nstep;
	Natom          = *FortranData::Natom;
	Mensemble      = *FortranData::Mensemble;
	max_no_neigh   = *FortranData::max_no_neigh;
	delta_t        = *FortranData::delta_t;
	gamma          = *FortranData::gamma;
	k_bolt         = *FortranData::k_bolt;
	mub            = *FortranData::mub;
	damping        = *FortranData::damping;
	binderc        =  FortranData::binderc;
	mavg           =  FortranData::mavg;
	mompar         = *FortranData::mompar;
	initexc        = *FortranData::initexc;
	do_dm          = static_cast<bool>(*FortranData::do_dm);
	do_jtensor     = static_cast<bool>(*FortranData::do_jtensor);
	max_no_dmneigh = *FortranData::max_no_dmneigh;


	// Thermfield
	switch (*FortranData::gpu_rng) {
	case 0: rngType = CURAND_RNG_PSEUDO_DEFAULT;  break;
	case 1: rngType = CURAND_RNG_PSEUDO_XORWOW;   break;
	case 2: rngType = CURAND_RNG_PSEUDO_MRG32K3A; break;
	case 3: rngType = CURAND_RNG_PSEUDO_MTGP32;   break;
	default:
		fprintf(stderr, "Unknown gpu_rng %d\n", *FortranData::gpu_rng);
		exit(EXIT_FAILURE);
		break;
	}
	randomSeed = (unsigned long long)*FortranData::gpu_rng_seed;
}

void CudaMdSimulation::initiate_fortran() {
	// Dimensions
	size_t N = Natom;
	size_t M = Mensemble;

	// Constants initiated?
	if (N == 0 || M == 0) {
		printf("MdSimulation: constants not initiated!\n");
		exit(EXIT_FAILURE);
	}

	// Inititate
	f_ncoup         .set(FortranData::ncoup         ,max_no_neigh,N);
	f_nlist         .set(FortranData::nlist         ,max_no_neigh,N);
	f_nlistsize     .set(FortranData::nlistsize     ,N);
	f_beff          .set(FortranData::beff          ,3,N,M);
	f_b2eff         .set(FortranData::b2eff         ,3,N,M);
	f_emomM         .set(FortranData::emomM         ,3,N,M);
	f_emom          .set(FortranData::emom          ,3,N,M);
	f_emom2         .set(FortranData::emom2         ,3,N,M);
	f_external_field.set(FortranData::external_field,3,N,M);
	f_mmom          .set(FortranData::mmom          ,N,M);
	f_btorque       .set(FortranData::btorque       ,3,N,M);
	f_temperature   .set(FortranData::temperature   ,N);
	f_mmom0         .set(FortranData::mmom0         ,N,M);
	f_mmom2         .set(FortranData::mmom2         ,N,M);
	f_mmomi         .set(FortranData::mmomi         ,N,M);
	f_dmvect        .set(FortranData::dmvect        ,3,max_no_dmneigh,N);
	f_dmlist        .set(FortranData::dmlist        ,max_no_dmneigh,N);
	f_dmlistsize    .set(FortranData::dmlistsize    ,N);
	
	if (*FortranData::do_jtensor == 1)
	{
		printf("\n CUDA: jTensor has been initialized \n");
		f_j_tensor  .set(FortranData::j_tensor      ,3,3,max_no_neigh,N);
	}
}

/*
static void printMemStat(const char * label) {
	size_t free;
	size_t total;
	CUresult result = cuMemGetInfo(&free, &total);
	printf("%s: free=%dk total=%dk (ret=%d)\n", label, free/1024, total/1024, result);
}
*/
bool CudaMdSimulation::initiateMatrices() {
	// Dimensions
	size_t N = Natom;
	size_t M = Mensemble;

	// Constants initiated?
	if (N == 0 || M == 0) {
		printf("CudaMdSimulation: constants not initiated!\n");
		exit(EXIT_FAILURE);
	}

	// Fortran initiate
	initiate_fortran();

	// Initiated?
	if (isInitiated) {
		printf("CudaMdSimulation: attempted to initiate already initiated CudaMdSimulation!\n");
		exit(EXIT_FAILURE);
	}

	// Inititate
	beff          .initiate(3,N,M);
	b2eff         .initiate(3,N,M);
	emomM         .initiate(3,N,M);
	emom          .initiate(3,N,M);
	emom2         .initiate(3,N,M);
	external_field.initiate(3,N,M);
	mmom          .initiate(N,M);
//	temperature   .initiate(N);
	mmom0         .initiate(N,M);
	mmom2         .initiate(N,M);
	mmomi         .initiate(N,M);

	if (FortranData::btorque)
		btorque.initiate(3,N,M);

	// Did we get the memory?
	if (!beff          .has_data() ||
	    !b2eff         .has_data() ||
	    !emomM         .has_data() ||
	    !emom          .has_data() ||
	    !emom2         .has_data() ||
	    !external_field.has_data() ||
	    !mmom          .has_data() ||
	    (!btorque.has_data() && (FortranData::btorque != NULL)) ||
//	    !temperature   .has_data() ||
	    !mmom0         .has_data() ||
	    !mmom2         .has_data() ||
	    !mmomi         .has_data()) {
		release();
		// Check for error
		const char * err = cudaGetErrorString(cudaGetLastError());
		fprintf(stderr, "CUDA: Failed to allocate memory: %s\n", err);
		return false;
	}

	// Flag that we're initiated
	isInitiated = true;

	// Initiate data
	copyFromFortran();
	return true;
}

void CudaMdSimulation::release() {
	isInitiated = false;

	beff          .free();
	b2eff         .free();
	emomM         .free();
	emom          .free();
	emom2         .free();
	external_field.free();
	mmom          .free();
	btorque       .free();
//	temperature   .free();
	mmom0         .free();
	mmom2         .free();
	mmomi         .free();
}

void CudaMdSimulation::copyFromFortran() {
	if (isInitiated) {
		beff          .read(FortranData::beff          );
		b2eff         .read(FortranData::b2eff         );
		emomM         .read(FortranData::emomM         );
		emom          .read(FortranData::emom          );
		emom2         .read(FortranData::emom2         );
		external_field.read(FortranData::external_field);
		mmom          .read(FortranData::mmom          );
//TODO		btorque       .read(FortranData::btorque       );
//		temperature   .read(FortranData::temperature   );
		mmom0         .read(FortranData::mmom0         );
		mmom2         .read(FortranData::mmom2         );
		mmomi         .read(FortranData::mmomi         );
	}
}

void CudaMdSimulation::copyToFortran() {
	if (isInitiated) {
		beff          .write(FortranData::beff          );
		b2eff         .write(FortranData::b2eff         );
		emomM         .write(FortranData::emomM         );
		emom          .write(FortranData::emom          );
		emom2         .write(FortranData::emom2         );
//		external_field.write(FortranData::external_field);
		mmom          .write(FortranData::mmom          );
//		btorque       .write(FortranData::btorque       );
//		temperature   .write(FortranData::temperature   );
		mmom0         .write(FortranData::mmom0         );
		mmom2         .write(FortranData::mmom2         );
		mmomi         .write(FortranData::mmomi         );
	}
}

void CudaMdSimulation::printConstants() {
	printf(
	"stt          : %c\n"
	"SDEalgh      : %d\n"
	"rstep        : %ld\n"
	"nstep        : %ld\n"
	"Natom        : %ld\n"
	"Mensemble    : %ld\n"
	"max_no_neigh : %ld\n"
	"delta_t      : %g\n"
	"gamma        : %g\n"
	"k_bolt       : %g\n"
	"mub          : %g\n"
	"damping      : %g\n"
	"binderc      : %g\n"
	"mavg         : %g\n",
	stt          ,
	SDEalgh      ,
	rstep        ,
	nstep        ,
	Natom        ,
	Mensemble    ,
	max_no_neigh ,
	delta_t      ,
	gamma        ,
	k_bolt       ,
	mub          ,
	damping      ,
	*binderc     ,
	*mavg        );
}

// Printing simulation status
// Added copy to fortran line so that simulation status is printed correctly > Thomas Nystrand 14/09/09
void CudaMdSimulation::printMdStatus(size_t mstep) {
	if (nstep > 20) {
		if (mstep % ((rstep + nstep) / 20) == 0) {
			copyToFortran(); // This is run so seldomly it has not impact on overall performance
			fortran_calc_simulation_status_variables(mavg);
			printf("CUDA: %3ld%% done. Mbar: %10.6f. U: %8.5f.\n",
				mstep * 100 / (rstep + nstep), *mavg, *binderc);
		}
	}
	else { 
		copyToFortran();
		fortran_calc_simulation_status_variables(mavg);
		printf("CUDA: Iteration %ld Mbar %13.6f\n", mstep, *mavg);
	}
}


// Spin Dynamics measurement phase
void CudaMdSimulation::measurementPhase() {
	// Unbuffered printf
	setbuf(stdout, NULL);
	setbuf(stderr, NULL);

	printf("CudaMdSimulation: md simulations starting\n");

	// Initiated?
	if (!isInitiated) {
		fprintf(stderr, "CudaMdSimulation: not initiated!\n");
		return;
	}


	// Timer
	StopwatchDeviceSync stopwatch = 
		StopwatchDeviceSync(GlobalStopwatchPool::get("Cuda measurement phase"));

	// Initiate default parallelization helper
	CudaParallelizationHelper::def.initiate(Natom, Mensemble);

	// Depontd integrator
	CudaDepondtIntegrator integrator;

	// Hamiltonian calculations
	CudaHamiltonianCalculations hamiltonian;

	// Moment updater
	CudaMomentUpdater momUpdater(mmom, mmom0, mmom2, emom, emom2, 
		emomM, mmomi, mompar, initexc);

	// Measurement
	CudaMeasurement measurement(emomM, emom, mmom);

	// Initiate integrator and Hamiltonian
	if (!integrator.initiate(Natom, Mensemble, stt, delta_t, rngType, randomSeed)) {
		fprintf(stderr, "CudaMdSimulation: integrator failed to initiate!\n");
		return;
	}
	if (!hamiltonian.initiate(f_ncoup, f_nlist, f_nlistsize, f_dmvect, f_dmlist, f_dmlistsize, do_dm, do_jtensor, f_j_tensor)) {
		fprintf(stderr, "CudaMdSimulation: Hamiltonian failed to initiate!\n");
		return;
	}

	// TEMPORARY PRINTING
	printf("\n");
	printf("________DEBUG System Information:___________ \n");
	printf("%zu\n", f_j_tensor.dimension_size(0));
	printf("%zu\n", f_j_tensor.dimension_size(1));
	printf("%zu\n", f_j_tensor.dimension_size(2));
	printf("%zu\n", f_j_tensor.dimension_size(3));
	printf("______________________________________\n");

	int mnn = f_j_tensor.dimension_size(2);
	int l = f_j_tensor.dimension_size(3);


	for (int l = 0 ; l < 1 ; l++  )
	{
		for (int k = 0 ; k < mnn ; k++  )
		{
			printf("__________\n");
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					unsigned int index = i + 3 * (j + 3 * (k + mnn * l));
					printf("%f\t", f_j_tensor[index]);
				}
        		printf("\n");
    		}
        	printf("\n");
		}
	}

	for(int i = 0 ; i < f_nlist.size(); i++)
	{
		//printf(" %d ", f_nlist.get_data()[i]);
		//real test = f_ncoup[i];
		//if (true) {
		//	printf("%f ",test);
		//}
	}
		

	


	printf("_______________________________________________\n");



	// Initiate constants for integrator
	integrator.initiateConstants(f_temperature, delta_t, gamma, k_bolt, mub, damping);

	// Debug
	//printConstants();

	// Timing
	stopwatch.add("initiate");

	// Time step loop
	for (size_t mstep = rstep + 1; mstep <= rstep + nstep; mstep++) {
		//export_mstep(mstep);

		// Measure
		measurement.measure(mstep);
		stopwatch.add("measurement");

		// Print simulation status for each 5% of the simulation length
		printMdStatus(mstep);

		// Apply Hamiltonian to obtain effective field
		hamiltonian.heisge(beff, emomM, external_field);
		stopwatch.add("hamiltonian");

		// Perform first step of SDE solver
		integrator.evolveFirst(beff, b2eff, btorque, emom, emom2, emomM, mmom);
		stopwatch.add("evolution");

		// Apply Hamiltonian to obtain effective field   
		hamiltonian.heisge(beff, emomM, external_field);
		stopwatch.add("hamiltonian");

		// Perform second (corrector) step of SDE solver
		integrator.evolveSecond(beff, b2eff, btorque, emom, emom2);
		stopwatch.add("evolution");

		// Update magnetic moments after time evolution step
		momUpdater.update();
		stopwatch.add("moments");

		// Check for error
		cudaError_t e = cudaGetLastError();
		if (e != cudaSuccess) {
			printf("Uncaught CUDA error %d: %s\n", e, cudaGetErrorString(e));
			cudaDeviceReset();
			exit(EXIT_FAILURE);
		}

	} // End loop over simulation steps

	// Final measure
	measurement.measure(rstep + nstep + 1);
	stopwatch.add("measurement");

	// Print remaining measurements 
	measurement.flushMeasurements(rstep + nstep + 1);
	stopwatch.add("flush measurement");

	// Synchronize with device
	cudaDeviceSynchronize();
	stopwatch.add("final synchronize");
}

