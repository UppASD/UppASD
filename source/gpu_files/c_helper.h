// Functions that help c call fortran
// by turning off C++ name_mangeling as well
// as following fortran name mangeling
// Niklas also defined short name wrappers to easy readability

#ifndef __C_HELPER_H__
#define __C_HELPER_H__

#include "real_type.h"

// Don't use C++ object naming
#ifdef __cplusplus
extern "C" {
#endif

// Fortran definition
extern void fortran_print_measurables(const int* obs_step, const int* obs_buff, void* indxb_obs,
                                      const char* obs_name, const char* obs_label, const int* obs_dim,
                                      const void* obs_buffer, const int* mstep);
extern void fortran_measure(const int* mstep);
extern void fortran_do_measurements(const int* mstep, int* do_copy);
extern void fortran_moment_update();
extern void fortran_flush_measurements(const int* mstep);
extern void fortran_print_correlations();
extern void fortran_measure_correlations(const real* emomM, const real* emom, const real* mmom, const int* mstep);
extern void fortran_measure_rest(const real* emomM, const real* emom, const real* mmom, const real* beff, const int* mstep);
extern void fortran_measure_moment(const real* emomM, const real* emom, const real* mmom, const int* mstep);
extern void fortran_calc_simulation_status_variables(real* mavg);
extern void fortran_fill_rngarray(real* ranv, const int* len);
extern void fortran_set_mstep(const int* mstep);
extern void fortran_get_mstep(int* mstep);

// Short name wrappers
inline void call_fortran_print_measurables(int obs_step, int obs_buff, void* indxb_obs,
                              const char* obs_name, const char* obs_label, int obs_dim,
                              const void* obs_buffer, int mstep) {
   ::fortran_print_measurables(&obs_step, &obs_buff, indxb_obs, obs_name, obs_label, &obs_dim, obs_buffer, &mstep);
}


inline void call_fortran_measure(int mstep) {
   ::fortran_measure(&mstep);
}

inline void call_fortran_moment_update() {
   ::fortran_moment_update();
}

inline void call_fortran_flush_measurements(int mstep) {
   ::fortran_flush_measurements(&mstep);
}

inline void call_fortran_print_correlations() {
   ::fortran_print_correlations();
}

inline void call_fortran_measure_correlations(const real* emomM, const real* emom, const real* mmom, int mstep) {
   ::fortran_measure_correlations(emomM, emom, mmom, &mstep);
}

inline void call_fortran_measure_rest(const real* emomM, const real* emom, const real* mmom, const real* beff, int mstep) {
   ::fortran_measure_rest(emomM, emom, mmom, beff, &mstep);
}

// Checking if its time to copy data
inline int call_fortran_do_measurements(int mstep) {
   int do_copy = 0;
   ::fortran_do_measurements(&mstep, &do_copy);
   return do_copy;
}

inline void call_fortran_measure_moment(const real* emomM, const real* emom, const real* mmom, int mstep) {
   ::fortran_measure_moment(emomM, emom, mmom, &mstep);
}

// RNG
inline void fill_rngarray(real* ranv, int len) {
   ::fortran_fill_rngarray(ranv, &len);
}

inline void export_mstep(int mstep) {
   ::fortran_set_mstep(&mstep);
}

inline int read_mstep() {
   int mstep = 0;
   ::fortran_get_mstep(&mstep);
   return mstep;
}

#ifdef __cplusplus
}
#endif

#endif
