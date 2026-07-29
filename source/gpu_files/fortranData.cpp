#include "fortranData.hpp"

// Keep the out-of-class definitions and C ABI entry points in the Fortran
// storage precision declared by fortranData.hpp.
#define real fortran_real

#include "real_type.h"
#include <thrust/complex.h>


// Constants
char* FortranData::stt;
int* FortranData::SDEalgh;

unsigned int* FortranData::rstep;
unsigned int* FortranData::nstep;
unsigned int* FortranData::Natom;
unsigned int* FortranData::nHam;
unsigned int* FortranData::Mensemble;
unsigned int* FortranData::max_no_neigh;
unsigned int* FortranData::ipmcnphase;
unsigned int* FortranData::mcnstep;
unsigned int* FortranData::ipnphase;
unsigned int* FortranData::NA;
unsigned int* FortranData::Natom_full;
unsigned int* FortranData::N1 = nullptr;
unsigned int* FortranData::N2 = nullptr;
unsigned int* FortranData::N3 = nullptr;
char* FortranData::do_gpu_convolution = nullptr;
int* FortranData::gpu_dipole_mode = nullptr;
int* FortranData::gpu_dipole_surface = nullptr;
real* FortranData::gpu_dipole_alpha = nullptr;
real* FortranData::gpu_dipole_rcut = nullptr;
real* FortranData::gpu_dipole_tol = nullptr;
int* FortranData::gpu_dipole_mesh = nullptr;
char* FortranData::BC1 = nullptr;
char* FortranData::BC2 = nullptr;
char* FortranData::BC3 = nullptr;
real* FortranData::C1 = nullptr;
real* FortranData::C2 = nullptr;
real* FortranData::C3 = nullptr;
real* FortranData::Bas = nullptr;
real* FortranData::alat = nullptr;
int* FortranData::do_dip = nullptr;
unsigned int* FortranData::num_macro = nullptr;
unsigned int* FortranData::macro_block_x = nullptr;
unsigned int* FortranData::macro_block_y = nullptr;
unsigned int* FortranData::macro_block_z = nullptr;
unsigned int* FortranData::macro_cell_index = nullptr;
unsigned int* FortranData::macro_nlistsize = nullptr;
real* FortranData::macro_center = nullptr;
real* FortranData::macro_min_coord = nullptr;
real* FortranData::macro_max_coord = nullptr;
unsigned int* FortranData::pme_num_macro = nullptr;
unsigned int* FortranData::pme_macro_grid = nullptr;
unsigned int* FortranData::pme_cell_index = nullptr;
unsigned int* FortranData::pme_macro_nlistsize = nullptr;
real* FortranData::pme_macro_center = nullptr;
real* FortranData::pme_macro_min_coord = nullptr;
real* FortranData::pme_macro_max_coord = nullptr;
int* FortranData::adaptive_geometry_mode = nullptr;
unsigned int* FortranData::adaptive_atoms = nullptr;
unsigned int* FortranData::adaptive_blocks = nullptr;
unsigned int* FortranData::adaptive_basis = nullptr;
unsigned int* FortranData::adaptive_fft_channels = nullptr;
unsigned int* FortranData::adaptive_fft_grid_channels = nullptr;
unsigned int* FortranData::adaptive_dynamic_channels = nullptr;
unsigned int* FortranData::adaptive_ensembles = nullptr;
unsigned int* FortranData::adaptive_selector_criteria = nullptr;
int* FortranData::adaptive_repetition_shape = nullptr;
int* FortranData::adaptive_block_shape = nullptr;
int* FortranData::adaptive_block_grid = nullptr;
real* FortranData::adaptive_cell_vectors = nullptr;
real* FortranData::adaptive_block_vectors = nullptr;
int* FortranData::adaptive_atom_to_block = nullptr;
int* FortranData::adaptive_atom_to_basis = nullptr;
int* FortranData::adaptive_atom_to_dynamic_channel = nullptr;
int* FortranData::adaptive_atom_to_fft_channel = nullptr;
int* FortranData::adaptive_atom_to_fft_grid_index = nullptr;
int* FortranData::adaptive_basis_to_dynamic_channel = nullptr;
int* FortranData::adaptive_basis_to_fft_channel = nullptr;
int* FortranData::adaptive_block_atom_count = nullptr;
int* FortranData::adaptive_block_atom_offset = nullptr;
int* FortranData::adaptive_block_atoms = nullptr;
int* FortranData::adaptive_block_grid_coordinate = nullptr;
int* FortranData::adaptive_block_basis_population = nullptr;
int* FortranData::adaptive_block_fft_population = nullptr;
int* FortranData::adaptive_block_dynamic_population = nullptr;
real* FortranData::adaptive_block_center = nullptr;
real* FortranData::adaptive_block_volume = nullptr;
int* FortranData::adaptive_block_state = nullptr;
int* FortranData::adaptive_pending_state = nullptr;
unsigned int* FortranData::adaptive_state_age = nullptr;
unsigned int* FortranData::adaptive_transition_epoch = nullptr;
real* FortranData::adaptive_selector_scores = nullptr;
real* FortranData::adaptive_coarse_moment = nullptr;
real* FortranData::adaptive_coarse_direction = nullptr;
real* FortranData::adaptive_coarse_field = nullptr;
real* FortranData::adaptive_channel_moment_sum = nullptr;
real* FortranData::adaptive_atom_moment = nullptr;
int* FortranData::adaptive_projection_block = nullptr;
real* FortranData::adaptive_projection_weight = nullptr;
unsigned int* FortranData::adaptive_bonds = nullptr;
int* FortranData::adaptive_bond_atom = nullptr;
real* FortranData::adaptive_bond_matrix = nullptr;
unsigned int* FortranData::adaptive_selector_edges = nullptr;
int* FortranData::adaptive_selector_edge = nullptr;
real* FortranData::adaptive_inverse_block_transpose = nullptr;
real* FortranData::adaptive_exchange_stiffness = nullptr;
real* FortranData::adaptive_spiralization = nullptr;
int* FortranData::adaptive_anisotropy_axis_count = nullptr;
real* FortranData::adaptive_anisotropy_axis = nullptr;
real* FortranData::adaptive_anisotropy_k1 = nullptr;
real* FortranData::adaptive_anisotropy_k2 = nullptr;
real* FortranData::adaptive_normalization_floor = nullptr;
real* FortranData::adaptive_magnetic_moment_si = nullptr;
real* FortranData::adaptive_gamma_per_ts = nullptr;
real* FortranData::adaptive_damping = nullptr;
int* FortranData::adaptive_mask_mode = nullptr;
unsigned int* FortranData::adaptive_update_interval = nullptr;
real* FortranData::adaptive_refine_threshold = nullptr;
real* FortranData::adaptive_coarsen_threshold = nullptr;
unsigned int* FortranData::adaptive_minimum_dwell = nullptr;
unsigned int* FortranData::adaptive_buffer_dilation = nullptr;
int* FortranData::adaptive_reconstruction_scheme = nullptr;
real* FortranData::adaptive_cone_angle_rad = nullptr;

unsigned int* FortranData::nq;
unsigned int* FortranData::sc_step;
unsigned int* FortranData::sc_sep;
char* FortranData::do_sc;
char* FortranData::do_sc_proj;
char* FortranData::do_sc_projch;
unsigned int* FortranData::sc_max_nstep;
unsigned int* FortranData::sc_window_fun;
unsigned int* FortranData::nw;
unsigned int* FortranData::NT;
unsigned int* FortranData::Nchmax;

real* FortranData::r_mid;
real* FortranData::q;
real* FortranData::coord;
real* FortranData::w;
int* FortranData::atype;
int* FortranData::lattice_atype;
int* FortranData::achtype;
fortran_complex* FortranData::m_k = nullptr;
fortran_complex* FortranData::m_kw = nullptr;
fortran_complex* FortranData::m_kt = nullptr;
real* FortranData::deltat_corr = nullptr;
real* FortranData::scstep_arr = nullptr;
int* FortranData::sc_nsamp_ptr = nullptr;  // Pointer to GPU-computed sample count
int* FortranData::sc_tidx_ptr = nullptr;   // Pointer to GPU-computed time step index
fortran_complex* FortranData::m_k_proj = nullptr;
fortran_complex* FortranData::m_k_projch = nullptr;
fortran_complex* FortranData::m_kt_proj = nullptr;
fortran_complex* FortranData::m_kt_projch = nullptr;
fortran_complex* FortranData::m_kw_proj = nullptr;
fortran_complex* FortranData::m_kw_projch = nullptr;


real* FortranData::delta_t;
real* FortranData::gamma;
real* FortranData::k_bolt;
real* FortranData::mub;
real* FortranData::mry;
real* FortranData::damping;
real* FortranData::Temp;

real* FortranData::binderc;
real* FortranData::mavg;

int* FortranData::mompar;
char* FortranData::initexc;

unsigned int* FortranData::do_dm;
unsigned int* FortranData::do_jtensor;
unsigned int* FortranData::do_aniso;
unsigned int* FortranData::max_no_dmneigh;

char* FortranData::do_cuda_measurements;
char* FortranData::do_avrg;
char* FortranData::do_cumu;
char* FortranData::do_autocorr;
unsigned int* FortranData::do_ene;
char* FortranData::do_skyno;
char* FortranData::do_gpu_correlations;
char* FortranData::do_gpu_timings;
char* FortranData::real_time_measure;
char* FortranData::do_avrg_proj;
char* FortranData::do_avrg_projch;
char* FortranData::do_cumu_proj;
unsigned int* FortranData::do_ralloy;

// Matrices
unsigned int * FortranData::aHam;
real* FortranData::ncoup;
unsigned int* FortranData::nlist;
unsigned int* FortranData::nlistsize;
real* FortranData::dmvect;
unsigned int* FortranData::dmlist;
unsigned int* FortranData::dmlistsize;
real* FortranData::beff;
real* FortranData::b2eff;
real* FortranData::emomM;
real* FortranData::emom;
real* FortranData::emom2;
real* FortranData::external_field;
real* FortranData::mmom;
real* FortranData::btorque;
real* FortranData::temperature;
real* FortranData::ipTemp;
unsigned int* FortranData::ipmcnstep;
real* FortranData::ipTemp_array;
unsigned int* FortranData::ipnstep;
real* FortranData::ipdelta_t;
real* FortranData::iplambda1;
real* FortranData::mmom0;
real* FortranData::mmom2;
real* FortranData::mmomi;

real* FortranData::j_tensor;
real* FortranData::kaniso;
real* FortranData::eaniso;
unsigned int* FortranData::taniso;
real* FortranData::sb;

real* FortranData::dxyz_vec;
int* FortranData::dxyz_atom;
int* FortranData::dxyz_list;

// GPU stuff
int* FortranData::gpu_mode;
int* FortranData::gpu_rng;
int* FortranData::gpu_rng_seed;

real* FortranData::mavg_buff;
real* FortranData::mavg2_buff;
real* FortranData::mavg4_buff;
real* FortranData::eavg_buff;
real* FortranData::eavg2_buff;
real* FortranData::spinwait;
unsigned int* FortranData::spinwaittable; 
real* FortranData::autocorr_buff;
real* FortranData::indxb_ac; //TODO: get rid of once printed on C
int* FortranData::achem_ch;
int* FortranData::asite_ch;

unsigned int* FortranData::avrg_step;
unsigned int* FortranData::avrg_buff;
unsigned int* FortranData::cumu_step;
unsigned int* FortranData::cumu_buff;
unsigned int* FortranData::ene_step;
unsigned int* FortranData::ene_buff;
unsigned int* FortranData::skyno_step;
unsigned int* FortranData::skyno_buff;
unsigned int* FortranData::ac_step;
unsigned int* FortranData::ac_buff;
unsigned int* FortranData::nspinwait;

void FortranData::setFlagPointers(unsigned int* p_do_dm, unsigned int* p_do_jtensor, unsigned int* p_do_anisotropy,
                                  char* p_do_avrg, char* p_do_proj_avrg, char* p_do_projch_avrg, char* p_do_cumu, char* p_do_cumu_proj,
                                  unsigned int* p_plotenergy, char* p_do_autocorr, char* p_do_tottraj,
                                  unsigned int* p_ntraj, char* p_do_cuda_measurements, char* p_do_skyno, char* p_do_sc,
                                  char* p_do_gpu_correlations, char* p_do_gpu_timings, char* p_real_time_measure, char* p_do_sc_proj, char* p_do_sc_projch,
                                  unsigned int* p_do_ralloy){


   do_dm = p_do_dm;
   do_jtensor = p_do_jtensor;
   do_aniso = p_do_anisotropy;
   do_cuda_measurements = p_do_cuda_measurements;
   do_avrg = p_do_avrg;
   do_cumu = p_do_cumu;
   do_autocorr = p_do_autocorr;
   do_ene = p_plotenergy;
   do_skyno = p_do_skyno;
   do_gpu_correlations = p_do_gpu_correlations;
   do_gpu_timings = p_do_gpu_timings;
   do_sc = p_do_sc;
   real_time_measure = p_real_time_measure;
   do_sc_proj = p_do_sc_proj;
   do_sc_projch = p_do_sc_projch;
   do_avrg_proj = p_do_proj_avrg;
   do_avrg_projch = p_do_projch_avrg;
   do_cumu_proj = p_do_cumu_proj;
   do_ralloy = p_do_ralloy;
   
}

void FortranData::setConstantPointers(char* p_stt, int* p_SDEalgh, unsigned int* p_rstep, unsigned int* p_nstep,
                                      unsigned int* p_Natom, unsigned int* p_Mensemble, unsigned int* p_max_no_neigh, 
                                      real* p_delta_t, real* p_gamma, real* p_k_bolt, real* p_mub, real* p_mplambda1,
                                      real* p_binderc, real* p_mavg, int* p_mompar, char* p_initexc, unsigned int* p_max_no_dmneigh,
                                      unsigned int* p_nHam, real* p_Temp, unsigned int* p_ipmcnphase, unsigned int* p_mcnstep, unsigned int* p_ipnphase,
                                      unsigned int* p_avrg_step, unsigned int* p_avrg_buff, unsigned int* p_cumu_step, unsigned int* p_cumu_buff,
                                      unsigned int* p_ene_step, unsigned int* p_ene_buff,  unsigned int*p_tottraj_step, unsigned int*p_tottraj_buff,
                                      unsigned int* p_skyno_step, unsigned int* p_skyno_buff, unsigned int* p_nq, unsigned int* p_sc_window_fun, unsigned int* p_nw,
                                      unsigned int* p_sc_sep, unsigned int* p_sc_step, unsigned int* p_sc_max_nstep,
                                      unsigned int* p_nspinwait, unsigned int* p_ac_step, unsigned int* p_ac_buff,
                                      unsigned int* p_nt, unsigned int* p_Nchmax, real* p_mry, unsigned int* p_NA, unsigned int* p_Natom_full){

   stt = p_stt;
   SDEalgh = p_SDEalgh;
                                    
   rstep = p_rstep;
   nstep = p_nstep;
   Natom = p_Natom;
   nHam = p_nHam;
   Mensemble = p_Mensemble;
   max_no_neigh = p_max_no_neigh;
                                    
   delta_t = p_delta_t;
   gamma = p_gamma;
   k_bolt = p_k_bolt;
   mub = p_mub;
   mry = p_mry;
   damping = p_mplambda1;
                                    
   binderc = p_binderc;
   mavg = p_mavg;
                                    
   mompar = p_mompar;
   initexc = p_initexc;
                                    
                                      
   max_no_dmneigh = p_max_no_dmneigh;

   Temp = p_Temp;
   ipnphase = p_ipnphase;
   ipmcnphase = p_ipmcnphase;
   mcnstep = p_mcnstep;

   avrg_step = p_avrg_step;
   avrg_buff = p_avrg_buff;
   cumu_step = p_cumu_step;
   cumu_buff = p_cumu_buff;
   ene_step = p_ene_step;
   ene_buff = p_ene_buff;

   skyno_step = p_skyno_step;
   skyno_buff = p_skyno_buff;
   nq = p_nq;
   sc_window_fun = p_sc_window_fun;
   nw = p_nw;
   sc_sep = p_sc_sep;
   sc_step = p_sc_step; 
   sc_max_nstep = p_sc_max_nstep;
   nspinwait = p_nspinwait;
   ac_step = p_ac_step;
   ac_buff = p_ac_buff;
   NT = p_nt;
   Nchmax = p_Nchmax;
   NA = p_NA;
   Natom_full = p_Natom_full;
}

void FortranData::setGpuGeometryPointers(unsigned int* p_N1, unsigned int* p_N2, unsigned int* p_N3,
                                         unsigned int* p_NA, char* p_BC1, char* p_BC2, char* p_BC3,
                                         real* p_C1, real* p_C2, real* p_C3, real* p_Bas,
                                         real* p_alat, char* p_do_gpu_convolution) {
   N1 = p_N1;
   N2 = p_N2;
   N3 = p_N3;
   NA = p_NA;
   BC1 = p_BC1;
   BC2 = p_BC2;
   BC3 = p_BC3;
   C1 = p_C1;
   C2 = p_C2;
   C3 = p_C3;
   Bas = p_Bas;
   alat = p_alat;
   do_gpu_convolution = p_do_gpu_convolution;
}

void FortranData::setGpuDipolePointers(int* p_mode, int* p_surface, real* p_alpha, real* p_rcut,
                                        real* p_tol, int* p_mesh) {
   gpu_dipole_mode = p_mode;
   gpu_dipole_surface = p_surface;
   gpu_dipole_alpha = p_alpha;
   gpu_dipole_rcut = p_rcut;
   gpu_dipole_tol = p_tol;
   gpu_dipole_mesh = p_mesh;
}

void FortranData::setMacrocellPointers(int* p_do_dip, unsigned int* p_num_macro,
                                       unsigned int* p_block_x, unsigned int* p_block_y,
                                       unsigned int* p_block_z, unsigned int* p_cell_index,
                                       unsigned int* p_macro_nlistsize, real* p_macro_center,
                                       real* p_macro_min_coord, real* p_macro_max_coord) {
   do_dip = p_do_dip;
   num_macro = p_num_macro;
   macro_block_x = p_block_x;
   macro_block_y = p_block_y;
   macro_block_z = p_block_z;
   macro_cell_index = p_cell_index;
   macro_nlistsize = p_macro_nlistsize;
   macro_center = p_macro_center;
   macro_min_coord = p_macro_min_coord;
   macro_max_coord = p_macro_max_coord;
}

void FortranData::clearMacrocellPointers() {
   do_dip = nullptr;
   num_macro = nullptr;
   macro_block_x = nullptr;
   macro_block_y = nullptr;
   macro_block_z = nullptr;
   macro_cell_index = nullptr;
   macro_nlistsize = nullptr;
   macro_center = nullptr;
   macro_min_coord = nullptr;
   macro_max_coord = nullptr;
}

void FortranData::setPmeMacrocellPointers(unsigned int* p_num_macro, unsigned int* p_grid,
                                           unsigned int* p_cell_index, unsigned int* p_nlistsize,
                                           real* p_center, real* p_min_coord, real* p_max_coord) {
   pme_num_macro = p_num_macro;
   pme_macro_grid = p_grid;
   pme_cell_index = p_cell_index;
   pme_macro_nlistsize = p_nlistsize;
   pme_macro_center = p_center;
   pme_macro_min_coord = p_min_coord;
   pme_macro_max_coord = p_max_coord;
}

void FortranData::clearPmeMacrocellPointers() {
   pme_num_macro = nullptr;
   pme_macro_grid = nullptr;
   pme_cell_index = nullptr;
   pme_macro_nlistsize = nullptr;
   pme_macro_center = nullptr;
   pme_macro_min_coord = nullptr;
   pme_macro_max_coord = nullptr;
}

void FortranData::setAdaptivePointers(
   int* geometry_mode, unsigned int* atoms, unsigned int* blocks,
   unsigned int* basis, unsigned int* fft_channels,
   unsigned int* fft_grid_channels,
   unsigned int* dynamic_channels, unsigned int* ensembles,
   unsigned int* selector_criteria, int* repetition_shape, int* block_shape,
   int* block_grid, real* cell_vectors, real* block_vectors,
   int* atom_to_block, int* atom_to_basis,
   int* atom_to_dynamic_channel, int* atom_to_fft_channel,
   int* atom_to_fft_grid_index, int* basis_to_dynamic_channel,
   int* basis_to_fft_channel, int* block_atom_count, int* block_atom_offset,
   int* block_atoms, int* block_grid_coordinate, int* block_basis_population,
   int* block_fft_population, int* block_dynamic_population,
   real* block_center, real* block_volume, int* block_state,
   int* pending_state, unsigned int* state_age,
   unsigned int* transition_epoch, real* selector_scores,
   real* coarse_moment, real* coarse_direction, real* coarse_field,
   real* channel_moment_sum) {
   adaptive_geometry_mode = geometry_mode;
   adaptive_atoms = atoms;
   adaptive_blocks = blocks;
   adaptive_basis = basis;
   adaptive_fft_channels = fft_channels;
   adaptive_fft_grid_channels = fft_grid_channels;
   adaptive_dynamic_channels = dynamic_channels;
   adaptive_ensembles = ensembles;
   adaptive_selector_criteria = selector_criteria;
   adaptive_repetition_shape = repetition_shape;
   adaptive_block_shape = block_shape;
   adaptive_block_grid = block_grid;
   adaptive_cell_vectors = cell_vectors;
   adaptive_block_vectors = block_vectors;
   adaptive_atom_to_block = atom_to_block;
   adaptive_atom_to_basis = atom_to_basis;
   adaptive_atom_to_dynamic_channel = atom_to_dynamic_channel;
   adaptive_atom_to_fft_channel = atom_to_fft_channel;
   adaptive_atom_to_fft_grid_index = atom_to_fft_grid_index;
   adaptive_basis_to_dynamic_channel = basis_to_dynamic_channel;
   adaptive_basis_to_fft_channel = basis_to_fft_channel;
   adaptive_block_atom_count = block_atom_count;
   adaptive_block_atom_offset = block_atom_offset;
   adaptive_block_atoms = block_atoms;
   adaptive_block_grid_coordinate = block_grid_coordinate;
   adaptive_block_basis_population = block_basis_population;
   adaptive_block_fft_population = block_fft_population;
   adaptive_block_dynamic_population = block_dynamic_population;
   adaptive_block_center = block_center;
   adaptive_block_volume = block_volume;
   adaptive_block_state = block_state;
   adaptive_pending_state = pending_state;
   adaptive_state_age = state_age;
   adaptive_transition_epoch = transition_epoch;
   adaptive_selector_scores = selector_scores;
   adaptive_coarse_moment = coarse_moment;
   adaptive_coarse_direction = coarse_direction;
   adaptive_coarse_field = coarse_field;
   adaptive_channel_moment_sum = channel_moment_sum;
}

void FortranData::clearAdaptivePointers() {
   adaptive_geometry_mode = nullptr;
   adaptive_atoms = nullptr;
   adaptive_blocks = nullptr;
   adaptive_basis = nullptr;
   adaptive_fft_channels = nullptr;
   adaptive_fft_grid_channels = nullptr;
   adaptive_dynamic_channels = nullptr;
   adaptive_ensembles = nullptr;
   adaptive_selector_criteria = nullptr;
   adaptive_repetition_shape = nullptr;
   adaptive_block_shape = nullptr;
   adaptive_block_grid = nullptr;
   adaptive_cell_vectors = nullptr;
   adaptive_block_vectors = nullptr;
   adaptive_atom_to_block = nullptr;
   adaptive_atom_to_basis = nullptr;
   adaptive_atom_to_dynamic_channel = nullptr;
   adaptive_atom_to_fft_channel = nullptr;
   adaptive_atom_to_fft_grid_index = nullptr;
   adaptive_basis_to_dynamic_channel = nullptr;
   adaptive_basis_to_fft_channel = nullptr;
   adaptive_block_atom_count = nullptr;
   adaptive_block_atom_offset = nullptr;
   adaptive_block_atoms = nullptr;
   adaptive_block_grid_coordinate = nullptr;
   adaptive_block_basis_population = nullptr;
   adaptive_block_fft_population = nullptr;
   adaptive_block_dynamic_population = nullptr;
   adaptive_block_center = nullptr;
   adaptive_block_volume = nullptr;
   adaptive_block_state = nullptr;
   adaptive_pending_state = nullptr;
   adaptive_state_age = nullptr;
   adaptive_transition_epoch = nullptr;
   adaptive_selector_scores = nullptr;
   adaptive_coarse_moment = nullptr;
   adaptive_coarse_direction = nullptr;
   adaptive_coarse_field = nullptr;
   adaptive_channel_moment_sum = nullptr;
   adaptive_atom_moment = nullptr;
   adaptive_projection_block = nullptr;
   adaptive_projection_weight = nullptr;
   adaptive_bonds = nullptr;
   adaptive_bond_atom = nullptr;
   adaptive_bond_matrix = nullptr;
   adaptive_selector_edges = nullptr;
   adaptive_selector_edge = nullptr;
   adaptive_inverse_block_transpose = nullptr;
   adaptive_exchange_stiffness = nullptr;
   adaptive_spiralization = nullptr;
   adaptive_anisotropy_axis_count = nullptr;
   adaptive_anisotropy_axis = nullptr;
   adaptive_anisotropy_k1 = nullptr;
   adaptive_anisotropy_k2 = nullptr;
   adaptive_normalization_floor = nullptr;
   adaptive_magnetic_moment_si = nullptr;
   adaptive_gamma_per_ts = nullptr;
   adaptive_damping = nullptr;
   adaptive_mask_mode = nullptr;
   adaptive_update_interval = nullptr;
   adaptive_refine_threshold = nullptr;
   adaptive_coarsen_threshold = nullptr;
   adaptive_minimum_dwell = nullptr;
   adaptive_buffer_dilation = nullptr;
   adaptive_reconstruction_scheme = nullptr;
   adaptive_cone_angle_rad = nullptr;
}

void FortranData::setAdaptiveKernelPointers(
   real* atom_moment, int* projection_block, real* projection_weight,
   unsigned int* bonds, int* bond_atom, real* bond_matrix,
   unsigned int* selector_edges, int* selector_edge,
   real* inverse_block_transpose, real* exchange_stiffness,
   real* spiralization, int* anisotropy_axis_count,
   real* anisotropy_axis, real* anisotropy_k1, real* anisotropy_k2,
   real* normalization_floor, real* magnetic_moment_si,
   real* gamma_per_ts, real* damping, int* mask_mode,
   unsigned int* update_interval, real* refine_threshold,
   real* coarsen_threshold, unsigned int* minimum_dwell,
   unsigned int* buffer_dilation, int* reconstruction_scheme,
   real* cone_angle_rad) {
   adaptive_atom_moment = atom_moment;
   adaptive_projection_block = projection_block;
   adaptive_projection_weight = projection_weight;
   adaptive_bonds = bonds;
   adaptive_bond_atom = bond_atom;
   adaptive_bond_matrix = bond_matrix;
   adaptive_selector_edges = selector_edges;
   adaptive_selector_edge = selector_edge;
   adaptive_inverse_block_transpose = inverse_block_transpose;
   adaptive_exchange_stiffness = exchange_stiffness;
   adaptive_spiralization = spiralization;
   adaptive_anisotropy_axis_count = anisotropy_axis_count;
   adaptive_anisotropy_axis = anisotropy_axis;
   adaptive_anisotropy_k1 = anisotropy_k1;
   adaptive_anisotropy_k2 = anisotropy_k2;
   adaptive_normalization_floor = normalization_floor;
   adaptive_magnetic_moment_si = magnetic_moment_si;
   adaptive_gamma_per_ts = gamma_per_ts;
   adaptive_damping = damping;
   adaptive_mask_mode = mask_mode;
   adaptive_update_interval = update_interval;
   adaptive_refine_threshold = refine_threshold;
   adaptive_coarsen_threshold = coarsen_threshold;
   adaptive_minimum_dwell = minimum_dwell;
   adaptive_buffer_dilation = buffer_dilation;
   adaptive_reconstruction_scheme = reconstruction_scheme;
   adaptive_cone_angle_rad = cone_angle_rad;
}

void FortranData::setHamiltonianPointers(real* p_ncoup, unsigned int* p_nlist, unsigned int* p_nlistsize,
                                         real* p_dmvect, unsigned int* p_dmlist, unsigned int* p_dmlistsize,
                                         real* p_kaniso, real* p_eaniso, unsigned int* p_taniso, real* p_sb,
                                         real* p_j_tensor, unsigned int* p_aHam, 
                                         real* p_external_field, real* p_btorque, real* p_Temp_array, 
                                         real * p_ipTemp, unsigned int * p_ipmcnstep,
                                         real * p_ipTemp_array, unsigned int* p_ipnstep,
                                         real * p_ipdelta_t, real * p_iplambda1){

   ncoup = p_ncoup;
   nlist = p_nlist;
   nlistsize = p_nlistsize;
   external_field = p_external_field;
   btorque = p_btorque;
   temperature = p_Temp_array;
   dmvect = p_dmvect;
   dmlist = p_dmlist;
   dmlistsize = p_dmlistsize;
   j_tensor = p_j_tensor;
   kaniso = p_kaniso;
   eaniso = p_eaniso;
   taniso = p_taniso;
   sb = p_sb;
   aHam = p_aHam;
   ipTemp = p_ipTemp;
   ipmcnstep = p_ipmcnstep;
   ipTemp_array = p_ipTemp_array;
   ipnstep = p_ipnstep;
   ipdelta_t = p_ipdelta_t;
   iplambda1 = p_iplambda1;
}


void FortranData::setLatticePointers(real* p_beff, real* p_b2eff, real* p_emomM, real* p_emom, real* p_emom2, 
                                     real* p_mmom, real* p_mmom0, real* p_mmom2, real* p_mmomi,
                                     real* p_dxyz_vec, int* p_dxyz_atom, int* p_dxyz_list, int* p_atype){


   beff = p_beff;
   b2eff = p_b2eff;
   emomM = p_emomM;
   emom = p_emom;
   emom2 = p_emom2;
   mmom = p_mmom;
   mmom0 = p_mmom0;
   mmom2 = p_mmom2;
   mmomi = p_mmomi;

   dxyz_vec = p_dxyz_vec;

   dxyz_atom = p_dxyz_atom;
   dxyz_list = p_dxyz_list;
   lattice_atype = p_atype;
}

//TODO:binderc, autocorr_buff, spinwait
void FortranData::setMeasurablePointers(real* p_mavg_buff, real* p_mavg2_buff, real* p_mavg4_buff,
                                         real* p_mavg_buff_proj, real* p_mavg2_buff_proj, real* p_mavg4_buff_proj, 
                                         real* p_binderc, real* p_avrgmcum, real* p_avrgm2cum, real* p_avrgm4cum, 
                                         real* p_eavg_buff, real* p_eavg2_buff, 
                                         int* p_traj_step, int* p_traj_buff, int* p_traj_atom,
                                         real* p_mmomb, real* p_mmomb_traj, real* p_emomb, real* p_emomb_traj,
                                         unsigned int* p_spinwaitt, real* p_spinwait, real* p_autocorr_buff, real* p_indxb_ac,
                                         int* p_achem_ch, int* p_asite_ch){

   mavg_buff = p_mavg_buff;
   mavg2_buff = p_mavg2_buff;
   mavg4_buff = p_mavg4_buff;
   eavg_buff = p_eavg_buff;
   eavg2_buff = p_eavg2_buff;
   spinwaittable = p_spinwaitt;
   spinwait = p_spinwait;
   autocorr_buff = p_autocorr_buff;
   indxb_ac = p_indxb_ac;
   achem_ch = p_achem_ch;
   asite_ch = p_asite_ch;


}

void FortranData::setCorrelationPointers(real* p_q, real* p_r_mid, real* p_coord, real* p_w,  void* p_m_k, 
                                        void* p_m_kw, void* p_m_kt, real* p_deltat_corr, real* p_scstep_arr,
                                        int* p_sc_nsamp, int* p_sc_tidx, int* p_atype, int* p_achtype, void* p_m_k_proj, 
                                        void* p_m_k_projch, void* p_m_kt_proj, void* p_m_kt_projch, void* p_m_kw_proj, 
                                        void* p_m_kw_projch){


   q = p_q;
   r_mid = p_r_mid;
   coord = p_coord;
   w = p_w;
   m_k  = reinterpret_cast<fortran_complex*>(p_m_k);
   m_kw = reinterpret_cast<fortran_complex*>(p_m_kw);
   m_kt = reinterpret_cast<fortran_complex*>(p_m_kt);
   deltat_corr = p_deltat_corr;
   scstep_arr = p_scstep_arr;
   sc_nsamp_ptr = p_sc_nsamp;  // Store pointer to GPU sample count
   sc_tidx_ptr = p_sc_tidx;    // Store pointer to GPU time index
   atype = p_atype;
   achtype = p_achtype;
   m_k_proj = reinterpret_cast<fortran_complex*>(p_m_k_proj);
   m_k_projch = reinterpret_cast<fortran_complex*>(p_m_k_projch);
   m_kt_proj = reinterpret_cast<fortran_complex*>(p_m_kt_proj);
   m_kt_projch = reinterpret_cast<fortran_complex*>(p_m_kt_projch);
   m_kw_proj = reinterpret_cast<fortran_complex*>(p_m_kw_proj);
   m_kw_projch = reinterpret_cast<fortran_complex*>(p_m_kw_projch);

}
/*void FortranData::setConstantPointers(char* p1, int* p2, unsigned int* p3, unsigned int* p4, unsigned int* p5,
                                      unsigned int* p6, unsigned int* p7, real* p8, real* p9, real* p10,
                                      real* p11, real* p12, real* p13, real* p14, int* p15, char* p16,
                                      unsigned int* p17, unsigned int* p18, unsigned int* p19,
                                      unsigned int* p20, unsigned int* p21, real * p_Temp, unsigned int* p_ipmcnphase, unsigned int* p_mcnstep,
                                      unsigned int * p_ipnphase) {
   stt = p1;
   SDEalgh = p2;

   rstep = p3;
   nstep = p4;
   Natom = p5;
   nHam = p21;
   Mensemble = p6;
   max_no_neigh = p7;

   delta_t = p8;
   gamma = p9;
   k_bolt = p10;
   mub = p11;
   damping = p12;

   binderc = p13;
   mavg = p14;

   mompar = p15;
   initexc = p16;

   do_dm = p17;
   max_no_dmneigh = p18;
   do_jtensor = p19;
   do_aniso = p20;
   Temp = p_Temp;
   ipnphase = p_ipnphase;
   ipmcnphase = p_ipmcnphase;
   mcnstep = p_mcnstep;
}*/

/*void FortranData::setMatrixPointers(real* p1, unsigned int* p2, unsigned int* p3, real* p4, real* p5,
                                    real* p6, real* p7, real* p8, real* p9, real* p10, real* p11, real* p12,
                                    real* p13, real* p14, real* p15, real* p16, unsigned int* p17,
                                    unsigned int* p18, real* p19, real* p20, real* p21, unsigned int* p22,
                                    real* p23, unsigned int* p24, real * p_ipTem   m_k_proj = p_m_k_proj;
   m_k_projch = p_m_k_projch;p, unsigned int * p_ipmcnstep,
                                    real * p_ipTemp_array, unsigned int* p_ipnstep) {
   ncoup = p1;
   nlist = p2;
   nlistsize = p3;
   beff = p4;
   b2eff = p5;
   emomM = p6;
   emom = p7;
   emom2 = p8;
   external_field = p9;
   mmom = p10;
   btorque = p11;
   temperature = p12;
   mmom0 = p13;
   mmom2 = p14;
   mmomi = p15;
   dmvect = p16;
   dmlist = p17;
   dmlistsize = p18;
   j_tensor = p19;
   kaniso = p20;
   eaniso = p21;
   taniso = p22;
   sb = p23;
   aHam = p24;
   ipTemp = p_ipTemp;
   ipmcnstep = p_ipmcnstep;
   ipTemp_array = p_ipTemp_array;
   ipnstep = p_ipnstep;
}*/

void FortranData::setInputDataPointers(int* p1, int* p2, int* p3) {
   gpu_mode = p1;
   gpu_rng = p2;
   gpu_rng_seed = p3;
}

// Fortran helpers
extern "C" void fortrandata_setflags_(unsigned int* p_do_dm, unsigned int* p_do_jtensor, unsigned int* p_do_anisotropy, 
   char* p_do_avrg, char* p_do_proj_avrg, char* p_do_projch_avrg, char* p_do_cumu, char* p_do_cumu_proj,
   unsigned int* p_plotenergy, char* p_do_autocorr, char* p_do_tottraj,
   unsigned int* p_ntraj, char* p_do_cuda_measurements, char* p_do_skyno, char* p_do_sc, char* p_do_gpu_correlations,
   char* p_do_gpu_timings, char* p_real_time_measure, char* p_do_sc_proj, char* p_do_sc_projch, unsigned int* p_do_ralloy) {
FortranData::setFlagPointers(
   p_do_dm, p_do_jtensor, p_do_anisotropy, p_do_avrg, p_do_proj_avrg, p_do_projch_avrg, p_do_cumu, p_do_cumu_proj,  p_plotenergy, 
   p_do_autocorr, p_do_tottraj, p_ntraj, p_do_cuda_measurements, p_do_skyno, p_do_sc, p_do_gpu_correlations, p_do_gpu_timings,
   p_real_time_measure, p_do_sc_proj, p_do_sc_projch, p_do_ralloy);
}

extern "C" void fortrandata_setconstants_(char* p_stt, int* p_SDEalgh, unsigned int* p_rstep, unsigned int* p_nstep,
   unsigned int* p_Natom, unsigned int* p_Mensemble, unsigned int* p_max_no_neigh, 
   real* p_delta_t, real* p_gamma, real* p_k_bolt, real* p_mub, real* p_mplambda1,
   real* p_binderc, real* p_mavg, int* p_mompar, char* p_initexc, unsigned int* p_max_no_dmneigh,
   unsigned int* p_nHam, real* p_Temp, unsigned int* p_ipmcnphase, unsigned int* p_mcnstep, unsigned int* p_ipnphase,
   unsigned int* p_avrg_step, unsigned int* p_avrg_buff, unsigned int* p_cumu_step, unsigned int* p_cumu_buff,
   unsigned int* p_ene_step, unsigned int* p_ene_buff, unsigned int*p_tottraj_step, unsigned int*p_tottraj_buff,
   unsigned int* p_skyno_step, unsigned int* p_skyno_buff,  unsigned int* p_nq, unsigned int* p_sc_window_fun, unsigned int* p_nw,
   unsigned int* p_sc_sep, unsigned int* p_sc_step, unsigned int* p_sc_max_nstep,
   unsigned int* p_nspinwait, unsigned int* p_ac_step, unsigned int* p_ac_buff, unsigned int* p_nt, unsigned int* p_Nchmax,
   real* p_mry, unsigned int* p_NA, unsigned int* p_Natom_full) {
FortranData::setConstantPointers(
   p_stt, p_SDEalgh, p_rstep, p_nstep, p_Natom, p_Mensemble, p_max_no_neigh, p_delta_t, p_gamma, 
   p_k_bolt, p_mub, p_mplambda1, p_binderc, p_mavg, p_mompar, p_initexc, p_max_no_dmneigh, p_nHam, 
   p_Temp, p_ipmcnphase, p_mcnstep, p_ipnphase,
   p_avrg_step, p_avrg_buff, p_cumu_step, p_cumu_buff, p_ene_step, p_ene_buff, p_tottraj_step, p_tottraj_buff,
   p_skyno_step, p_skyno_buff, p_nq, p_sc_window_fun, p_nw, p_sc_sep, p_sc_step, p_sc_max_nstep,
   p_nspinwait, p_ac_step, p_ac_buff, p_nt, p_Nchmax, p_mry, p_NA, p_Natom_full);
}

extern "C" void fortrandata_setgpugeometry_(unsigned int* p_N1, unsigned int* p_N2, unsigned int* p_N3,
   unsigned int* p_NA, char* p_BC1, char* p_BC2, char* p_BC3,
   real* p_C1, real* p_C2, real* p_C3, real* p_Bas, real* p_alat, char* p_do_gpu_convolution) {
FortranData::setGpuGeometryPointers(p_N1, p_N2, p_N3, p_NA, p_BC1, p_BC2, p_BC3,
                                    p_C1, p_C2, p_C3, p_Bas, p_alat, p_do_gpu_convolution);
}

extern "C" void fortrandata_setgpudipole_(int* p_mode, int* p_surface, real* p_alpha, real* p_rcut,
                                            real* p_tol, int* p_mesh) {
   FortranData::setGpuDipolePointers(p_mode, p_surface, p_alpha, p_rcut, p_tol, p_mesh);
}

extern "C" void fortrandata_setmacrocell_(int* p_do_dip, unsigned int* p_num_macro,
   unsigned int* p_block_x, unsigned int* p_block_y, unsigned int* p_block_z,
   unsigned int* p_cell_index, unsigned int* p_macro_nlistsize, real* p_macro_center,
   real* p_macro_min_coord, real* p_macro_max_coord) {
   FortranData::setMacrocellPointers(p_do_dip, p_num_macro, p_block_x, p_block_y, p_block_z,
                                     p_cell_index, p_macro_nlistsize, p_macro_center,
                                     p_macro_min_coord, p_macro_max_coord);
}

extern "C" void fortrandata_clearmacell_() {
   FortranData::clearMacrocellPointers();
}

extern "C" void fortrandata_setpmemacrocell_(unsigned int* p_num_macro, unsigned int* p_macro_grid,
   unsigned int* p_cell_index, unsigned int* p_macro_nlistsize, real* p_macro_center,
   real* p_macro_min_coord, real* p_macro_max_coord) {
   FortranData::setPmeMacrocellPointers(p_num_macro, p_macro_grid, p_cell_index, p_macro_nlistsize,
                                         p_macro_center, p_macro_min_coord, p_macro_max_coord);
}

extern "C" void fortrandata_clearpmemacrocell_() {
   FortranData::clearPmeMacrocellPointers();
}

extern "C" void fortrandata_setadaptivetopology_(
   int* geometry_mode, unsigned int* atoms, unsigned int* blocks,
   unsigned int* basis, unsigned int* fft_channels,
   unsigned int* fft_grid_channels,
   unsigned int* dynamic_channels, unsigned int* ensembles,
   unsigned int* selector_criteria, int* repetition_shape, int* block_shape,
   int* block_grid, real* cell_vectors, real* block_vectors,
   int* atom_to_block, int* atom_to_basis,
   int* atom_to_dynamic_channel, int* atom_to_fft_channel,
   int* atom_to_fft_grid_index, int* basis_to_dynamic_channel,
   int* basis_to_fft_channel, int* block_atom_count, int* block_atom_offset,
   int* block_atoms, int* block_grid_coordinate, int* block_basis_population,
   int* block_fft_population, int* block_dynamic_population,
   real* block_center, real* block_volume, int* block_state,
   int* pending_state, unsigned int* state_age,
   unsigned int* transition_epoch, real* selector_scores,
   real* coarse_moment, real* coarse_direction, real* coarse_field,
   real* channel_moment_sum) {
   FortranData::setAdaptivePointers(
      geometry_mode, atoms, blocks, basis, fft_channels, fft_grid_channels,
      dynamic_channels, ensembles, selector_criteria, repetition_shape,
      block_shape, block_grid, cell_vectors, block_vectors, atom_to_block, atom_to_basis,
      atom_to_dynamic_channel, atom_to_fft_channel, atom_to_fft_grid_index,
      basis_to_dynamic_channel, basis_to_fft_channel, block_atom_count,
      block_atom_offset, block_atoms, block_grid_coordinate,
      block_basis_population, block_fft_population, block_dynamic_population,
      block_center, block_volume, block_state, pending_state, state_age,
      transition_epoch, selector_scores, coarse_moment, coarse_direction,
      coarse_field, channel_moment_sum);
}

extern "C" void fortrandata_clearadaptivetopology_() {
   FortranData::clearAdaptivePointers();
}

extern "C" void fortrandata_setadaptivekernels_(
   real* atom_moment, int* projection_block, real* projection_weight,
   unsigned int* bonds, int* bond_atom, real* bond_matrix,
   unsigned int* selector_edges, int* selector_edge,
   real* inverse_block_transpose, real* exchange_stiffness,
   real* spiralization, int* anisotropy_axis_count,
   real* anisotropy_axis, real* anisotropy_k1, real* anisotropy_k2,
   real* normalization_floor, real* magnetic_moment_si,
   real* gamma_per_ts, real* damping, int* mask_mode,
   unsigned int* update_interval, real* refine_threshold,
   real* coarsen_threshold, unsigned int* minimum_dwell,
   unsigned int* buffer_dilation, int* reconstruction_scheme,
   real* cone_angle_rad) {
   FortranData::setAdaptiveKernelPointers(
      atom_moment, projection_block, projection_weight, bonds, bond_atom,
      bond_matrix, selector_edges, selector_edge, inverse_block_transpose,
      exchange_stiffness, spiralization, anisotropy_axis_count,
      anisotropy_axis, anisotropy_k1, anisotropy_k2, normalization_floor,
      magnetic_moment_si, gamma_per_ts, damping, mask_mode, update_interval,
      refine_threshold, coarsen_threshold, minimum_dwell, buffer_dilation,
      reconstruction_scheme, cone_angle_rad);
}

extern "C" void fortrandata_sethamiltonian_(real* p_ncoup, unsigned int* p_nlist, unsigned int* p_nlistsize,
   real* p_dmvect, unsigned int* p_dmlist, unsigned int* p_dmlistsize,
   real* p_kaniso, real* p_eaniso, unsigned int* p_taniso, real* p_sb,
   real* p_j_tensor, unsigned int* p_aHam, 
   real* p_external_field, real* p_btorque, real* p_Temp_array, 
   real * p_ipTemp, unsigned int * p_ipmcnstep,
   real * p_ipTemp_array, unsigned int* p_ipnstep,
   real * p_ipdelta_t, real * p_iplambda1) {
FortranData::setHamiltonianPointers(
   p_ncoup, p_nlist, p_nlistsize, p_dmvect, p_dmlist, p_dmlistsize,  p_kaniso, p_eaniso, p_taniso, p_sb, 
   p_j_tensor, p_aHam, p_external_field, p_btorque, p_Temp_array, 
   p_ipTemp, p_ipmcnstep, p_ipTemp_array, p_ipnstep, p_ipdelta_t, p_iplambda1);
}

extern "C" void fortrandata_setlattice_(real* p_beff, real* p_b2eff, real* p_emomM, real* p_emom, real* p_emom2, 
   real* p_mmom, real* p_mmom0, real* p_mmom2, real* p_mmomi, real* p_dxyz_vec, int* p_dxyz_atom, int* p_dxyz_list,
   int* p_atype) {
FortranData::setLatticePointers(
   p_beff, p_b2eff, p_emomM, p_emom, p_emom2, p_mmom, p_mmom0, p_mmom2, p_mmomi, p_dxyz_vec, p_dxyz_atom, p_dxyz_list,
   p_atype);
}

extern "C" void fortrandata_setmeasurables_(real* p_mavg_buff, real* p_mavg2_buff, real* p_mavg4_buff,
   real* p_mavg_buff_proj, real* p_mavg2_buff_proj, real* p_mavg4_buff_proj, 
   real* p_binderc, real* p_avrgmcum, real* p_avrgm2cum, real* p_avrgm4cum, 
   real* p_eavg_buff, real* p_eavg2_buff, 
   int* p_traj_step, int* p_traj_buff, int* p_traj_atom,
   real* p_mmomb, real* p_mmomb_traj, real* p_emomb, real* p_emomb_traj,
   unsigned int* p_spinwaitt, real* p_spinwait, real* p_indxb_ac, real* p_autocorr_buff,
   int* p_achem_ch, int* p_asite_ch) {
FortranData::setMeasurablePointers(
   p_mavg_buff, p_mavg2_buff, p_mavg4_buff, p_mavg_buff_proj, p_mavg2_buff_proj, p_mavg4_buff_proj, 
   p_binderc, p_avrgmcum, p_avrgm2cum, p_avrgm4cum, p_eavg_buff, p_eavg2_buff, 
    p_traj_step, p_traj_buff, p_traj_atom,p_mmomb, p_mmomb_traj, p_emomb, p_emomb_traj,
   p_spinwaitt, p_spinwait, p_indxb_ac, p_autocorr_buff, p_achem_ch, p_asite_ch);
}

extern "C" void fortrandata_setcorrelations_(real* p_q, real* p_r_mid, real* p_coord, real* p_w, void* p_m_k, 
                                             void* p_m_kw, void* p_m_kt, real* p_deltat_corr, real* p_scstep_arr, 
                                             int* p_sc_nsamp, int* p_sc_tidx, int* p_atype, int* p_achtype, 
                                             void* p_m_k_proj, void* p_m_k_projch, void* p_m_kt_proj, void* p_m_kt_projch,
                                             void* p_m_kw_proj, void* p_m_kw_projch) {
FortranData::setCorrelationPointers(
   p_q, p_r_mid, p_coord, p_w,  p_m_k, p_m_kw, p_m_kt, p_deltat_corr, p_scstep_arr, p_sc_nsamp, p_sc_tidx,
   p_atype, p_achtype, p_m_k_proj, p_m_k_projch, p_m_kt_proj, p_m_kt_projch, p_m_kw_proj, p_m_kw_projch);
}
/*extern "C" void fortrandata_setconstants_(char* p1, int* p2, unsigned int* p3, unsigned int* p4,
                                          unsigned int* p5, unsigned int* p6, unsigned int* p7, real* p8,
                                          real* p9, real* p10, real* p11, real* p12, real* p13, real* p14,
                                          int* p15, char* p16, unsigned int* p17, unsigned int* p18,
                                          unsigned int* p19, unsigned int* p20, unsigned int* p21, real* p_Temp, unsigned int* p_ipmcnphase,
                                          unsigned int* p_mcnstep, unsigned int* p_ipnphase) {
   FortranData::setConstantPointers(
       p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p_Temp, p_ipmcnphase, p_mcnstep, p_ipnphase);
}*/

/*extern "C" void fortrandata_setmatrices_(real* p1, unsigned int* p2, unsigned int* p3, real* p4, real* p5,
                                         real* p6, real* p7, real* p8, real* p9, real* p10, real* p11,
                                         real* p12, real* p13, real* p14, real* p15, real* p16,
                                         unsigned int* p17, unsigned int* p18, real* p19, real* p20,
                                         real* p21, unsigned int* p22, real* p23, unsigned int* p24, 
                                         real* p_ipTemp, unsigned int* p_ipmcnstep, real* p_ipTemp_array, unsigned int* p_ipnstep) {
   FortranData::setMatrixPointers(p1,
                                  p2,
                                  p3,
                                  p4,
                                  p5,
                                  p6,
                                  p7,
                                  p8,
                                  p9,
                                  p10,
                                  p11,
                                  p12,
                                  p13,
                                  p14,
                                  p15,
                                  p16,
                                  p17,
                                  p18,
                                  p19,
                                  p20,
                                  p21,
                                  p22,
                                  p23,
                                  p24,
                                  p_ipTemp,
                                  p_ipmcnstep,
                                  p_ipTemp_array,
                                  p_ipnstep);
}*/

extern "C" void fortrandata_setinputdata_(int* p1, int* p2, int* p3) {
   FortranData::setInputDataPointers(p1, p2, p3);
}
