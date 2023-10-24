// Functions that help c call fortran
// by turning off C++ name_mangeling as well
// as following fortran name mangeling
// Niklas also defined short name wrappers to easy readability


#ifndef __C_HELPER_H__
#define __C_HELPER_H__

#include "real_type.h"

// Fortran naming
#if defined (__GNU__) && defined (__Intel__)
	#error "both __GNU__ and __Intel__ defined!"
#elif defined (__GNU__)
	#define FORTNAME(mod,func) __ ## mod ## _MOD_ ## func
#elif defined (__Intel__)
	#define FORTNAME(mod,func) mod ## _mp_ ## func ## _
#elif defined (__NVHPC__)
	#define FORTNAME(mod,func) mod ## _ ## func ## _
#else
	#error "FORTRAN compiler not specified! (compile with -D__GNU__ or -D__Intel__)"
#endif


// Don't use C++ object naming
#ifdef __cplusplus
extern "C" {
#endif


// Fortran definition
extern void FORTNAME(chelper,fortran_measure)(const size_t * mstep);
extern void FORTNAME(chelper,fortran_do_measurements)(const size_t * mstep, int * do_copy);
extern void FORTNAME(chelper,fortran_moment_update)();
extern void FORTNAME(chelper,fortran_flush_measurements)(const size_t * mstep);
extern void FORTNAME(chelper,cmdsim_initiate_constants)();
extern void FORTNAME(chelper,fortran_measure_moment)(const real * emomM, const real * emom, const real * mmom, const size_t * mstep);
extern void FORTNAME(chelper,fortran_calc_simulation_status_variables)(real * mavg);

// Short name wrappers
inline void fortran_measure(size_t mstep)            {FORTNAME(chelper,fortran_measure)(&mstep);}
inline void fortran_moment_update()                  {FORTNAME(chelper,fortran_moment_update)();}
inline void fortran_flush_measurements(size_t mstep) {FORTNAME(chelper,fortran_flush_measurements)(&mstep);}
inline void fortran_init_c_md_const()                {FORTNAME(chelper,cmdsim_initiate_constants)();}

// For the status variables
inline void fortran_calc_simulation_status_variables(real * mavg) {
         FORTNAME(chelper,fortran_calc_simulation_status_variables)(mavg);
}

// Checking if its time to copy data
inline int fortran_do_measurements(size_t mstep) {
	int do_copy = 0;
	FORTNAME(chelper,fortran_do_measurements)(&mstep, &do_copy);
	return do_copy;
}

inline void fortran_measure_moment(const real * emomM, const real * emom, const real * mmom, size_t mstep) {
	FORTNAME(chelper,fortran_measure_moment)(emomM, emom, mmom, &mstep);
}


// RNG
extern void FORTNAME(randomnumbers,fill_rngarray)(real * ranv, const size_t * arr);
inline void fill_rngarray(real * ranv, size_t len) {FORTNAME(randomnumbers,fill_rngarray)(ranv, &len);}


// Fortrans mstep variable
extern size_t FORTNAME(simulationdata,mstep);
inline void   export_mstep(size_t mstep) {FORTNAME(simulationdata,mstep) = mstep;}
inline size_t read_mstep()               {return FORTNAME(simulationdata,mstep);}


#ifdef __cplusplus
}
#endif

#endif

