

#include "c_headers.hpp"
#include "c_helper.h"
//#include "cudaGPUErrchk.hpp"
#include "fortranData.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "stopwatchPool.hpp"
#include "gpuStructures.hpp"
#include "gpuSimulation.hpp"
#include "tensor.hpp"
#include "gpu_wrappers.h"
#include "gpuDepondtIntegrator.hpp"
#include "gpuAdaptiveRuntime.hpp"
#include "gpuLatticeConvolutionHamiltonian.hpp"
#include "gpuDipoleConvolution.hpp"
#include "gpuHamiltonianCalculations.hpp"
#include "measurement/gpuMeasurement.hpp"
#include "correlations/gpuCorrelations.hpp"
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#if defined(HIP_V)
#include <hiprand/hiprand.h>
#elif defined(CUDA_V)
#include <curand.h>
#endif

//#include "gpuSDSimulation.hpp"
GpuSimulation::GpuSimulation()
: maxThreads(256)
, maxBlocks(1024)
{
    isInitiated = false;
}

GpuSimulation::~GpuSimulation() {
    if(isInitiated && !isFreed) {release();}
}

void GpuSimulation::initiateConstants() {
   //printf("HERE - 0\n");

    // Read the Fortran input once per GPU-simulation initialization.  The
    // instrumentation sites then only test this cached C++ boolean.
    Stopwatch::setTimingMode(FortranData::do_gpu_timings ? *FortranData::do_gpu_timings : 'N');

    SimParam.SDEalgh = *FortranData::SDEalgh;
    if(!(SimParam.SDEalgh == 1 || SimParam.SDEalgh == 4 || SimParam.SDEalgh == 5 || SimParam.SDEalgh == 11)) {
        std::fprintf(stderr, "Invalid SDEalgh!\n");
        std::exit(EXIT_FAILURE);
    }

    Flags.do_dm = static_cast<bool>(*FortranData::do_dm);
    Flags.do_jtensor = static_cast<bool>(*FortranData::do_jtensor);
    Flags.do_aniso = static_cast<bool>(*FortranData::do_aniso);

    Flags.do_sc = *FortranData::do_sc;
    // Must match the Fortran allocation gate (chelper.f90 FortranData_Initiate,
    // 'do_gpu_correlations==Y'): a plain static_cast<bool>(char) is non-zero for
    // ANY letter including 'N', so the correlation buffers would be wired up even
    // when Fortran allocated nothing, then read through null pointers and crash.
    Flags.do_gpu_correlations = (*FortranData::do_gpu_correlations == 'Y');
    Flags.do_sc_proj = *FortranData::do_sc_proj;
    Flags.do_sc_projch = *FortranData::do_sc_projch;
    Flags.do_ene = *FortranData::do_ene;
    Flags.do_gpu_measurements = ((*FortranData::do_cuda_measurements == 'Y')? true : false);
    //Flags.do_avrg = static_cast<bool>(*FortranData::do_avrg);
    //Flags.do_cumu = static_cast<bool>(*FortranData::do_cumu);

    SimParam.N = *FortranData::Natom;
    SimParam.NH  = *FortranData::nHam;
    SimParam.M = *FortranData::Mensemble;
    SimParam.mnn = *FortranData::max_no_neigh;
    SimParam.mnndm = *FortranData::max_no_dmneigh;
    SimParam.N1 = FortranData::N1 ? *FortranData::N1 : 0;
    SimParam.N2 = FortranData::N2 ? *FortranData::N2 : 0;
    SimParam.N3 = FortranData::N3 ? *FortranData::N3 : 0;
    SimParam.NA = FortranData::NA ? *FortranData::NA : 0;
    SimParam.do_gpu_convolution = FortranData::do_gpu_convolution && *FortranData::do_gpu_convolution == 'Y';
    SimParam.BC1 = FortranData::BC1 ? *FortranData::BC1 : '0';
    SimParam.BC2 = FortranData::BC2 ? *FortranData::BC2 : '0';
    SimParam.BC3 = FortranData::BC3 ? *FortranData::BC3 : '0';
    geometryC1.set(FortranData::C1, 3);
    geometryC2.set(FortranData::C2, 3);
    geometryC3.set(FortranData::C3, 3);
    geometryBas.set(FortranData::Bas, 3 * SimParam.NA);
    SimParam.C1 = geometryC1.data();
    SimParam.C2 = geometryC2.data();
    SimParam.C3 = geometryC3.data();
    SimParam.Bas = geometryBas.data();
    SimParam.alat = FortranData::alat ? static_cast<double>(*FortranData::alat) : 0.0;
    SimParam.gpu_dipole_tol = FortranData::gpu_dipole_tol ?
       static_cast<double>(*FortranData::gpu_dipole_tol) : 1.0e-10;
    SimParam.ipmcnphase = *FortranData::ipmcnphase;
    SimParam.ipnphase = *FortranData::ipnphase;
    if(SimParam.ipnphase == 0) SimParam.ipnphase = 1;
    if(SimParam.ipmcnphase == 0) SimParam.ipmcnphase = 1;
    SimParam.mcnstep = *FortranData::mcnstep;
   // printf("cpp = %i, fortr = %i",  SimParam.ipnphase , *FortranData::ipnphase );
    // Constants 
    SimParam.stt = *FortranData::stt;
    SimParam.rstep = *FortranData::rstep;
    SimParam.nstep = *FortranData::nstep;

    SimParam.delta_t = *FortranData::delta_t;
    SimParam.gamma = *FortranData::gamma;
    SimParam.k_bolt = *FortranData::k_bolt;
    SimParam.mub = *FortranData::mub;
    SimParam.mry = *FortranData::mry;
    SimParam.Temp = *FortranData::Temp;
    SimParam.damping = *FortranData::damping;
    SimParam.mompar = *FortranData::mompar;
    SimParam.initexc = *FortranData::initexc;

    SimParam.binderc = FortranData::binderc;
    SimParam.mavg = FortranData::mavg;

    SimParam.sc_sep = *FortranData::sc_sep;
    SimParam.sc_step = *FortranData::sc_step;
    SimParam.sc_window_fun = *FortranData::sc_window_fun;
    SimParam.nq = *FortranData::nq;
    SimParam.nw = *FortranData::nw;
    SimParam.sc_max_nstep = *FortranData::sc_max_nstep;
    SimParam.NT = *FortranData::NT;
    SimParam.Nchmax= *FortranData::Nchmax;

    SimParam.ene_step = *FortranData::ene_step;
    SimParam.ene_buff = *FortranData::ene_buff;


   // SimParam.avrg_step = *FortranData::avrg_step;  
   // SimParam.avrg_buff = *FortranData::avrg_buff; 
    //SimParam.cumu_step = *FortranData::cumu_step; 
    //SimParam.cumu_buff = *FortranData::cumu_buff; 

 #if defined(CUDA_V)   
    switch(*FortranData::gpu_rng) {
    case 0: SimParam.rngType = CURAND_RNG_PSEUDO_DEFAULT; break;
    case 1: SimParam.rngType = CURAND_RNG_PSEUDO_XORWOW; break;
    case 2: SimParam.rngType = CURAND_RNG_PSEUDO_MRG32K3A; break;
    case 3: SimParam.rngType = CURAND_RNG_PSEUDO_MTGP32; break;
    default:
    std::fprintf(stderr, "Unknown gpu_rng %d\n", *FortranData::gpu_rng);
    std::exit(EXIT_FAILURE);
    break;
    }
 #elif defined(HIP_V)
switch(*FortranData::gpu_rng) {
    case 0: SimParam.rngType = HIPRAND_RNG_PSEUDO_DEFAULT; break;
    case 1: SimParam.rngType = HIPRAND_RNG_PSEUDO_XORWOW; break;
    case 2: SimParam.rngType = HIPRAND_RNG_PSEUDO_MRG32K3A; break;
    case 3: SimParam.rngType = HIPRAND_RNG_PSEUDO_MTGP32; break;
    default:
    std::fprintf(stderr, "Unknown gpu_rng %d\n", *FortranData::gpu_rng);
    std::exit(EXIT_FAILURE);
    break;
    }
#else
    #error "Either USE_CUDA or USE_HIP must be defined"
#endif

    SimParam.randomSeed = (unsigned long long)*FortranData::gpu_rng_seed;
   //printf("HERE - 1\n");

    
}

void GpuSimulation::initiate_fortran_cpu_matrices() {

    long int N = static_cast <long int>( SimParam.N);
    long int  NH = static_cast <long int>( SimParam.NH);
    long int  M = static_cast <long int>( SimParam.M);
   long int mnn = static_cast <long int>( SimParam.mnn);
    long int  mnndm = static_cast <long int>( SimParam.mnndm);
   long int ipnphase = static_cast <long int>( SimParam.ipnphase);
    long int ipmcnphase = static_cast <long int>( SimParam.ipmcnphase);
    long int nq = static_cast <long int>( SimParam.nq);
    long int nw= static_cast <long int>( SimParam.nw);
    long int NT= static_cast <long int>( SimParam.NT);
    long int Nchmax= static_cast <long int>( SimParam.Nchmax);


    // Constants initiated?
    if(N == 0 || M == 0 || NH == 0) {
        std::printf("GPUSimulation: constants not initiated!\n");
        std::exit(EXIT_FAILURE);
    }

    cpuHamiltonian.aHam.set(FortranData::aHam, N);                             
    cpuHamiltonian.ncoup.set(FortranData::ncoup, mnn, NH);           
    cpuHamiltonian.nlist.set(FortranData::nlist, mnn, N);
    cpuHamiltonian.nlistsize.set(FortranData::nlistsize, NH);

    if(Flags.do_dm != 0) {
        cpuHamiltonian.dmvect.set(FortranData::dmvect,  static_cast <long int>(3), mnndm, NH);    
        cpuHamiltonian.dmlist.set(FortranData::dmlist, mnndm, N);
        cpuHamiltonian.dmlistsize.set(FortranData::dmlistsize, NH);
        }
    if(Flags.do_jtensor != 0) {
        cpuHamiltonian.j_tensor.set(FortranData::j_tensor,  static_cast <long int>(3),  static_cast <long int>(3), mnn, NH);}
    if(Flags.do_aniso != 0) {
        cpuHamiltonian.kaniso.set(FortranData::kaniso,  static_cast <long int>(2), N);;
        cpuHamiltonian.eaniso.set(FortranData::eaniso,  static_cast <long int>(3), N);;
        cpuHamiltonian.taniso.set(FortranData::taniso, N);;
        cpuHamiltonian.sb.set(FortranData::sb, N);
    }
    cpuHamiltonian.extfield.set(FortranData::external_field,  static_cast <long int>(3), N, M);
    if(FortranData::gpu_dipole_mode && *FortranData::gpu_dipole_mode != 0 && FortranData::pme_num_macro &&
       *FortranData::pme_num_macro > 0 && FortranData::pme_cell_index && FortranData::pme_macro_nlistsize &&
       FortranData::pme_macro_center && FortranData::pme_macro_min_coord && FortranData::pme_macro_max_coord) {
        const std::size_t macroCount = *FortranData::pme_num_macro;
        std::vector<std::size_t> observedPopulation(macroCount, 0);
        for(long int atom = 0; atom < N; ++atom) {
            const unsigned int oneBasedCell = FortranData::pme_cell_index[atom];
            if(oneBasedCell == 0 || oneBasedCell > macroCount) {
                throw std::runtime_error("GPU macrocell map contains an out-of-range cell index");
            }
            ++observedPopulation[oneBasedCell - 1];
        }
        for(std::size_t cell = 0; cell < macroCount; ++cell) {
            if(observedPopulation[cell] != FortranData::pme_macro_nlistsize[cell]) {
                throw std::runtime_error("GPU macrocell map and macro_nlistsize disagree");
            }
        }
        cpuHamiltonian.macro_cell_index.set(FortranData::pme_cell_index, N);
        cpuHamiltonian.macro_nlistsize.set(FortranData::pme_macro_nlistsize, macroCount);
        cpuHamiltonian.macro_center.set(FortranData::pme_macro_center, 3, macroCount);
        cpuHamiltonian.macro_min_coord.set(FortranData::pme_macro_min_coord, 3, macroCount);
        cpuHamiltonian.macro_max_coord.set(FortranData::pme_macro_max_coord, 3, macroCount);
    }
    cpuLattice.beff.set(FortranData::beff,  static_cast <long int>(3), N, M);
    cpuLattice.b2eff.set(FortranData::b2eff,  static_cast <long int>(3), N, M);
    cpuLattice.emomM.set(FortranData::emomM,  static_cast <long int>(3), N, M);
    cpuLattice.emom.set(FortranData::emom,  static_cast <long int>(3), N, M);
    cpuLattice.emom2.set(FortranData::emom2,  static_cast <long int>(3), N, M);
    cpuLattice.mmom.set(FortranData::mmom, N, M);
    cpuLattice.mmom0.set(FortranData::mmom0, N, M);
    cpuLattice.mmom2.set(FortranData::mmom2, N, M);
    cpuLattice.mmomi.set(FortranData::mmomi, N, M);
    // btorque is allocated only for the spin-transfer-torque path.  In fp64
    // Tensor::set merely retained a null optional pointer; fp32 staging must
    // not dereference it while converting Fortran doubles to device storage.
    if(FortranData::btorque)
        cpuLattice.btorque.set(FortranData::btorque,  static_cast <long int>(3), N, M);
    cpuLattice.temperature.set(FortranData::temperature, N);
    cpuLattice.ipTemp.set(FortranData::ipTemp, ipmcnphase);
    cpuLattice.ipmcnstep.set(FortranData::ipmcnstep, ipmcnphase);
    cpuLattice.ipTemp_array.set(FortranData::ipTemp_array, N, ipnphase);
    cpuLattice.ipnstep.set(FortranData::ipnstep, ipnphase);
    cpuLattice.ipdelta_t.set(FortranData::ipdelta_t, ipnphase);
    cpuLattice.iplambda1.set(FortranData::iplambda1, ipnphase);
    if(Flags.do_gpu_correlations){
        cpuCorrelations.r_mid.set(FortranData::r_mid, static_cast <long int>(3));
        cpuCorrelations.q.set(FortranData::q, static_cast <long int>(3), nq);
        cpuCorrelations.w.set(FortranData::w, nw);
        cpuCorrelations.coord.set(FortranData::coord, static_cast <long int>(3), N);

        // The Fortran side (chelper.f90:FortranData_Initiate) only allocates the
        // static buffer (m_k) for do_sc 'C'/'Y' and the dynamic buffers
        // (m_kt/m_kw) for do_sc 'Q'/'Y'.  The complex set() overload copies from
        // its source (cpu_complex != fortran_complex), so calling it on a buffer
        // that was never allocated reads through a null/unassociated Fortran
        // pointer and segfaults.  Gate each set() on the mode that allocated it.
        const bool sc_static  = (Flags.do_sc == 'C') || (Flags.do_sc == 'Y');
        const bool sc_dynamic = (Flags.do_sc == 'Q') || (Flags.do_sc == 'Y');
        if(sc_static)
            cpuCorrelations.m_k.set(FortranData::m_k, static_cast <long int>(3), nq);
        if(sc_dynamic){
            cpuCorrelations.m_kt.set(FortranData::m_kt, static_cast <long int>(3), static_cast <long int>(nq), static_cast <long int>(SimParam.sc_max_nstep));
            cpuCorrelations.m_kw.set(FortranData::m_kw, static_cast <long int>(3), nq, nw);
            cpuCorrelations.deltat_corr.set(FortranData::deltat_corr, static_cast <long int>(SimParam.sc_max_nstep + 1));
            cpuCorrelations.scstep_arr.set(FortranData::scstep_arr, static_cast <long int>(SimParam.sc_max_nstep + 1));
        }
        // Projected buffers follow the same static/dynamic split (chelper.f90):
        // m_k_proj for 'C'/'Y', m_kt_proj/m_kw_proj for 'Q'/'T'/'Y'.
        if((Flags.do_sc_proj=='C')||(Flags.do_sc_proj=='Q')||(Flags.do_sc_proj=='Y')||(Flags.do_sc_proj=='T')){
        cpuCorrelations.atype.set(FortranData::atype, N);
        if((Flags.do_sc_proj=='C')||(Flags.do_sc_proj=='Y'))
            cpuCorrelations.m_k_proj.set(FortranData::m_k_proj, static_cast <long int>(3), NT, nq);
        if((Flags.do_sc_proj=='Q')||(Flags.do_sc_proj=='T')||(Flags.do_sc_proj=='Y')){
            cpuCorrelations.m_kt_proj.set(FortranData::m_kt_proj, static_cast <long int>(3), NT, nq, static_cast <long int>(SimParam.sc_max_nstep));
            cpuCorrelations.m_kw_proj.set(FortranData::m_kw_proj, static_cast <long int>(3), NT, nq, nw);
        }
        }
        if((Flags.do_sc_projch=='C')||(Flags.do_sc_projch=='Q')||(Flags.do_sc_projch=='Y')||(Flags.do_sc_projch=='T')){
        cpuCorrelations.achtype.set(FortranData::achtype, N);
        if((Flags.do_sc_projch=='C')||(Flags.do_sc_projch=='Y'))
            cpuCorrelations.m_k_projch.set(FortranData::m_k_projch, static_cast <long int>(3), Nchmax, nq);
        if((Flags.do_sc_projch=='Q')||(Flags.do_sc_projch=='T')||(Flags.do_sc_projch=='Y')){
            cpuCorrelations.m_kt_projch.set(FortranData::m_kt_projch, static_cast <long int>(3), Nchmax, nq, static_cast <long int>(SimParam.sc_max_nstep));
            cpuCorrelations.m_kw_projch.set(FortranData::m_kw_projch, static_cast <long int>(3), Nchmax, nq, nw);
        }
        }

    }
  // printf("HERE - 2\n");

  //  if(FortranData::ipnstep == nullptr)printf("ITS EMPTY\n");

   /* if (Flags.do_mphase_now != 0){
        if (Flags.do_avrg !=0){
            cpuMeasurables.mavg_buff.set(FortranData::mavrg_buff, SimParam.avrg_buff);
            //cpuMeasurables.eavg_buff.set(FortranData::eavrg_buff, SimParam.avrg_buff);
        }
        if (Flags.do_cumu !=0){
            cpuMeasurables.mcumu_buff.set(FortranData::mcumu_buff, SimParam.avrg_buff);
            //cpuMeasurables.ecumu_buff.set(FortranData::ecumu_buff, SimParam.cumu_buff); 
        }
    }*/
        

    

}

// ----------------------------------------------------------------------------
// M1: upfront device-memory budget. Each helper below mirrors the Allocate
// calls in initiateMatrices; the integrator/thermfield/convolution/correlations
// expose their own estimateBytes() next to their allocations. The release()-time
// self-check (estimate vs TensorMemoryTracker peak) guards against drift.
// ----------------------------------------------------------------------------
namespace {

std::string formatBytes(std::size_t bytes) {
   char buf[64];
   const double b = static_cast<double>(bytes);
   if(bytes >= 1000000000ULL)   std::snprintf(buf, sizeof(buf), "%.3f GB", b / 1e9);
   else if(bytes >= 1000000ULL) std::snprintf(buf, sizeof(buf), "%.1f MB", b / 1e6);
   else if(bytes >= 1000ULL)    std::snprintf(buf, sizeof(buf), "%.1f kB", b / 1e3);
   else                         std::snprintf(buf, sizeof(buf), "%zu B", bytes);
   return std::string(buf);
}

struct BudgetLine { const char* name; std::size_t bytes; };

std::size_t hamiltonianBytes(const Flag& F, const SimulationParameters& P) {
   const std::size_t N = P.N, NH = P.NH, mnn = P.mnn, mnndm = P.mnndm, M = P.M;
   std::size_t b = N * sizeof(unsigned int);                 // aHam(N)
   if(F.do_jtensor) {
      b += NH * 3 * 3 * mnn * sizeof(real);                  // j_tensor(NH,3,3,mnn)
      b += mnn * N * sizeof(unsigned int);                   // nlist(mnn,N)
      b += NH * sizeof(unsigned int);                        // nlistsize(NH)
   } else {
      b += NH * mnn * sizeof(real);                          // ncoup(NH,mnn)
      b += N * mnn * sizeof(unsigned int);                   // nlist(N,mnn)
      b += NH * sizeof(unsigned int);                        // nlistsize(NH)
      if(F.do_dm) {
         b += 3 * mnndm * NH * sizeof(real);                 // dmvect(3,mnndm,NH)
         b += mnndm * N * sizeof(unsigned int);              // dmlist(mnndm,N)
         b += NH * sizeof(unsigned int);                     // dmlistsize(NH)
      }
   }
   if(F.do_aniso) {
      b += 2 * N * sizeof(real);                             // kaniso(2,N)
      b += 3 * N * sizeof(real);                             // eaniso(3,N)
      b += N * sizeof(unsigned int);                         // taniso(N)
      b += N * sizeof(real);                                 // sb(N)
   }
   b += 3 * N * M * sizeof(real);                            // extfield(3,N,M)
   // The map and geometry are created with the base matrices.  The macro
   // moment tensor itself belongs to the phase-local dipole operator and is
   // accounted together with its FFT allocations below.
   if(FortranData::gpu_dipole_mode && *FortranData::gpu_dipole_mode != 0 && FortranData::pme_num_macro) {
      const std::size_t macro = *FortranData::pme_num_macro;
      b += N * sizeof(unsigned int);                          // cell_index(N)
      b += macro * sizeof(unsigned int);                      // macro_nlistsize(Nmacro)
      b += 9 * macro * sizeof(real);                          // centre/min/max geometry
   }
   return b;
}

std::size_t latticeBytes(const Flag& F, const SimulationParameters& P, bool is_mc) {
   const std::size_t nm3 = 3 * P.N * P.M, nm = P.N * P.M;
   std::size_t b = 5 * nm3 * sizeof(real);                   // beff,b2eff,emomM,emom,emom2
   b += 2 * nm * sizeof(real);                               // mmom,mmomi
   // M2: eneff aliases beff for SD without aniso; own buffer only for MC or aniso.
   if(is_mc || F.do_aniso != 0) b += nm3 * sizeof(real);     // eneff
   // M2: mmom0/mmom2 only for MC or mompar!=0.
   if(is_mc || P.mompar != 0)   b += 2 * nm * sizeof(real);  // mmom0, mmom2
   if(FortranData::btorque) b += nm3 * sizeof(real);         // btorque(3,N,M)
   return b;
}

std::size_t energiesBytes(const Flag& F, const SimulationParameters& P) {
   return (F.do_ene > 0 ? P.M * 7 : static_cast<std::size_t>(1)) * sizeof(real); // energyM
}

// Mirror of GpuHamiltonianCalculations::canUseLatticeConvolution parameter gates:
// convolution only activates for a reduced Hamiltonian (NH==NA) on a full
// Bravais grid, otherwise the run falls back to neighbour lists (no FFT arrays).
bool willUseConvolution(const SimulationParameters& P) {
   if(!P.do_gpu_convolution) return false;
   if(P.N1 == 0 || P.N2 == 0 || P.N3 == 0 || P.NA == 0) return false;
   if(P.N != P.N1 * P.N2 * P.N3 * P.NA) return false;
   if(P.NH != P.NA) return false;
   return true;
}

std::size_t cv6DipoleBytes(const SimulationParameters& P) {
   if(!FortranData::gpu_dipole_mode || *FortranData::gpu_dipole_mode == 0) return 0;
   const int mode = *FortranData::gpu_dipole_mode;
   if(mode != 1 && mode != 2) throw std::runtime_error("unknown GPU dipole mode before device allocation");
#if defined(HIP_V)
   // CUDA fp32 OPEN_FFT has its finite-oracle, production-E2E, and memcheck
   // evidence.  HIP has none on this project host, so keep that backend
   // explicitly closed until it satisfies the identical acceptance matrix.
   if(mode == 2 && sizeof(real) != sizeof(double)) {
      throw std::runtime_error("OPEN_FFT fp32 is enabled for accepted CUDA only; HIP fp32 requires its own oracle and sanitizer acceptance");
   }
#endif
   if(mode == 2 && P.gpu_dipole_tol != 1.0e-10) {
      throw std::runtime_error("OPEN_FFT rejects Ewald tolerance overrides before device allocation");
   }
   if(!FortranData::pme_num_macro || *FortranData::pme_num_macro == 0 || !FortranData::pme_macro_grid ||
      !FortranData::pme_macro_center || P.NA == 0 || P.M == 0) {
      throw std::runtime_error("GPU dipole mode was requested without a complete GPU macrocell layout");
   }
   GpuDipoleDescriptorInput input{};
   input.atomistic_grid = {P.N1, P.N2, P.N3};
   input.macro_grid = {static_cast<std::size_t>(FortranData::pme_macro_grid[0]),
                       static_cast<std::size_t>(FortranData::pme_macro_grid[1]),
                       static_cast<std::size_t>(FortranData::pme_macro_grid[2])};
   input.basis = static_cast<unsigned int>(P.NA);
   input.ensembles = static_cast<unsigned int>(P.M);
   input.boundaries = {{P.BC1, P.BC2, P.BC3}};
   input.c1 = P.C1;
   input.c2 = P.C2;
   input.c3 = P.C3;
   input.basis_offsets = P.Bas;
   // The preflight owns no device geometry; the same descriptor factory only
   // needs its shape/count contract here. Runtime supplies the device centres.
   input.macro_centers = nullptr;
   input.macro_count = *FortranData::pme_num_macro;
   input.alat = P.alat;
   input.tolerance = P.gpu_dipole_tol;
   GpuDipoleConvolutionDescriptor descriptor{};
   const bool descriptor_ok = mode == 1 ? makeEwald3dFftDipoleDescriptor(input, descriptor) :
                                          makeOpenFftDipoleDescriptor(input, descriptor);
   if(!descriptor_ok) {
      throw std::runtime_error(mode == 1 ? "GPU EWALD3D_FFT descriptor is invalid before device allocation" :
                                          "GPU OPEN_FFT descriptor is invalid before device allocation; it cannot fall back to periodic FFT");
   }
   const std::size_t convolution = GpuDipoleConvolution::estimateBytes(descriptor);
   const std::size_t macroCount = static_cast<std::size_t>(*FortranData::pme_num_macro);
   if(macroCount > std::numeric_limits<std::size_t>::max() / 3 / P.M ||
      3 * macroCount * P.M > std::numeric_limits<std::size_t>::max() / sizeof(real))
      throw std::runtime_error("GPU dipole macro-moment memory estimate overflowed");
   const std::size_t macroMoments = 3 * macroCount * P.M * sizeof(real);
   if(convolution == 0 || macroMoments > std::numeric_limits<std::size_t>::max() - convolution)
      throw std::runtime_error(mode == 1 ? "GPU EWALD3D_FFT memory estimate overflowed or failed" :
                                           "GPU OPEN_FFT memory estimate overflowed or failed");
   return convolution + macroMoments;
}

GpuAdaptiveTopologyInput adaptiveTopologyInput() {
   GpuAdaptiveTopologyInput t;
   if(!FortranData::adaptive_geometry_mode) return t;
   t.geometryMode = *FortranData::adaptive_geometry_mode;
   t.atoms = *FortranData::adaptive_atoms;
   t.blocks = *FortranData::adaptive_blocks;
   t.basis = *FortranData::adaptive_basis;
   t.fftChannelsPerBlock = *FortranData::adaptive_fft_channels;
   t.fftGridChannels = *FortranData::adaptive_fft_grid_channels;
   t.dynamicChannels = *FortranData::adaptive_dynamic_channels;
   t.ensembles = *FortranData::adaptive_ensembles;
   t.repetitionShape = FortranData::adaptive_repetition_shape;
   t.blockShape = FortranData::adaptive_block_shape;
   t.blockGrid = FortranData::adaptive_block_grid;
   t.cellVectors = FortranData::adaptive_cell_vectors;
   t.blockVectors = FortranData::adaptive_block_vectors;
   t.atomToBlock = FortranData::adaptive_atom_to_block;
   t.atomToBasis = FortranData::adaptive_atom_to_basis;
   t.atomToDynamicChannel = FortranData::adaptive_atom_to_dynamic_channel;
   t.atomToFftChannel = FortranData::adaptive_atom_to_fft_channel;
   t.atomToFftGridIndex = FortranData::adaptive_atom_to_fft_grid_index;
   t.basisToDynamicChannel = FortranData::adaptive_basis_to_dynamic_channel;
   t.basisToFftChannel = FortranData::adaptive_basis_to_fft_channel;
   t.blockAtomCount = FortranData::adaptive_block_atom_count;
   t.blockAtomOffset = FortranData::adaptive_block_atom_offset;
   t.blockAtoms = FortranData::adaptive_block_atoms;
   t.blockGridCoordinate = FortranData::adaptive_block_grid_coordinate;
   t.blockBasisPopulation = FortranData::adaptive_block_basis_population;
   t.blockFftChannelPopulation = FortranData::adaptive_block_fft_population;
   t.blockDynamicChannelPopulation = FortranData::adaptive_block_dynamic_population;
   t.blockCenter = FortranData::adaptive_block_center;
   t.blockVolume = FortranData::adaptive_block_volume;
   return t;
}

GpuAdaptiveRuntimeInput adaptiveRuntimeInput() {
   GpuAdaptiveRuntimeInput r;
   if(!FortranData::adaptive_geometry_mode) return r;
   r.blockState = FortranData::adaptive_block_state;
   r.pendingState = FortranData::adaptive_pending_state;
   r.stateAge = FortranData::adaptive_state_age;
   r.transitionEpoch = FortranData::adaptive_transition_epoch;
   r.selectorCriteria = *FortranData::adaptive_selector_criteria;
   r.selectorScores = FortranData::adaptive_selector_scores;
   r.coarseMoment = FortranData::adaptive_coarse_moment;
   r.coarseDirection = FortranData::adaptive_coarse_direction;
   r.coarseField = FortranData::adaptive_coarse_field;
   r.channelMomentSum = FortranData::adaptive_channel_moment_sum;
   r.kernels.atomMoment = FortranData::adaptive_atom_moment;
   r.kernels.atomAnisotropyAxisCount = FortranData::adaptive_atom_anisotropy_axis_count;
   r.kernels.atomAnisotropyAxis = FortranData::adaptive_atom_anisotropy_axis;
   r.kernels.atomAnisotropyK1 = FortranData::adaptive_atom_anisotropy_k1;
   r.kernels.atomAnisotropyK2 = FortranData::adaptive_atom_anisotropy_k2;
   r.kernels.projectionBlock = FortranData::adaptive_projection_block;
   r.kernels.projectionWeight = FortranData::adaptive_projection_weight;
   r.kernels.bonds = FortranData::adaptive_bonds ? *FortranData::adaptive_bonds : 0;
   r.kernels.bondAtom = FortranData::adaptive_bond_atom;
   r.kernels.bondMatrix = FortranData::adaptive_bond_matrix;
   r.kernels.selectorEdges =
      FortranData::adaptive_selector_edges ? *FortranData::adaptive_selector_edges : 0;
   r.kernels.selectorEdge = FortranData::adaptive_selector_edge;
   r.kernels.inverseBlockTranspose = FortranData::adaptive_inverse_block_transpose;
   r.kernels.exchangeStiffness = FortranData::adaptive_exchange_stiffness;
   r.kernels.spiralization = FortranData::adaptive_spiralization;
   r.kernels.anisotropyAxisCount = FortranData::adaptive_anisotropy_axis_count;
   r.kernels.anisotropyAxis = FortranData::adaptive_anisotropy_axis;
   r.kernels.anisotropyK1 = FortranData::adaptive_anisotropy_k1;
   r.kernels.anisotropyK2 = FortranData::adaptive_anisotropy_k2;
   if(FortranData::adaptive_normalization_floor)
      r.kernels.normalizationFloor = *FortranData::adaptive_normalization_floor;
   if(FortranData::adaptive_magnetic_moment_si)
      r.kernels.magneticMomentSi = *FortranData::adaptive_magnetic_moment_si;
   if(FortranData::adaptive_gamma_per_ts)
      r.kernels.gammaPerTs = *FortranData::adaptive_gamma_per_ts;
   if(FortranData::adaptive_damping)
      r.kernels.damping = *FortranData::adaptive_damping;
   return r;
}

std::size_t adaptiveRuntimeBytes(const SimulationParameters& parameters) {
   if(!FortranData::adaptive_geometry_mode) return 0;
   if(!FortranData::adaptive_atoms || !FortranData::adaptive_blocks ||
      !FortranData::adaptive_basis || !FortranData::adaptive_fft_channels ||
      !FortranData::adaptive_fft_grid_channels ||
      !FortranData::adaptive_dynamic_channels || !FortranData::adaptive_ensembles ||
      !FortranData::adaptive_selector_criteria) {
      throw std::runtime_error("GPU adaptive topology staging is missing scalar counts");
   }
   const auto topology = adaptiveTopologyInput();
   const auto runtime = adaptiveRuntimeInput();
   std::string diagnostic;
   if(!GpuAdaptiveRuntime::validate(topology, runtime, parameters.N, parameters.M, diagnostic))
      throw std::runtime_error(diagnostic + " (adaptive preflight; no device allocation attempted)");
   return GpuAdaptiveRuntime::estimateBytes(topology, runtime);
}

// Sum every device allocation the run will make, compare to free VRAM, and abort
// with a table before the first Allocate if it will not fit. Returns the total.
std::size_t computeAndCheckDeviceBudget(const Flag& F, const SimulationParameters& P, bool is_mc) {
   std::vector<BudgetLine> lines = {
      {"Hamiltonian (Jij/DM/aniso) + ext field", hamiltonianBytes(F, P)},
      {"Lattice streaming arrays",               latticeBytes(F, P, is_mc)},
      {"Depondt integrator + thermfield",        GpuDepondtIntegrator::estimateBytes(P)},
      {"Energy buffers",                         energiesBytes(F, P)},
      {"GPU measurement phase",                  GpuMeasurement::estimateBytes(P.N, P.M)},
      {"FFT convolution",                        willUseConvolution(P) ? GpuLatticeConvolutionHamiltonian::estimateBytes(P) : static_cast<std::size_t>(0)},
      {"CV6 dipole FFT, macro moments + staging", cv6DipoleBytes(P)},
      {"Correlations",                           (FortranData::do_gpu_correlations && *FortranData::do_gpu_correlations == 'Y') ? GpuCorrelations::estimateBytes(F, P) : static_cast<std::size_t>(0)},
      {"Adaptive CG topology + runtime",          adaptiveRuntimeBytes(P)},
   };
   std::size_t total = 0;
   for(const auto& l : lines) {
      if(l.bytes > std::numeric_limits<std::size_t>::max() - total) {
         throw std::runtime_error("GPU device-memory budget overflow");
      }
      total += l.bytes;
   }
   if(FortranData::gpu_dipole_mode && *FortranData::gpu_dipole_mode != 0) {
      std::printf("Gpu: planned device tensor inventory (full-run preflight scope):\n");
      for(const auto& line : lines)
         std::printf("Gpu:   %-40s %zu B\n", line.name, line.bytes);
      std::printf("Gpu:   %-40s %zu B\n", "TOTAL", total);
   }

   std::size_t freeB = 0, totalB = 0;
   const bool haveInfo = (GPU_MEM_GET_INFO(&freeB, &totalB) == GPU_SUCCESS);
   // Test-only deterministic preflight seam.  A real free-memory query is
   // inherently machine-dependent, so acceptance coverage can cap it and
   // prove that the rejection happens before the normal allocation phase.
   // It is intentionally opt-in and has no effect on normal evaluations.
   bool effectiveHaveInfo = haveInfo;
   if(const char* testFree = std::getenv("UPPASD_GPU_TEST_FREE_BYTES")) {
      char* end = nullptr;
      const unsigned long long parsed = std::strtoull(testFree, &end, 10);
      if(end != testFree && *end == '\0' && parsed > 0) {
         freeB = static_cast<std::size_t>(parsed);
         totalB = std::max(totalB, freeB);
         effectiveHaveInfo = true;
      }
   }
   const char* skip = std::getenv("UPPASD_GPU_SKIP_BUDGET");
   const bool bypass = (skip != nullptr && skip[0] != '\0');
   const bool overBudget = effectiveHaveInfo && (static_cast<double>(total) > 0.90 * static_cast<double>(freeB));

   if(overBudget && !bypass) {
      std::sort(lines.begin(), lines.end(),
                [](const BudgetLine& a, const BudgetLine& b) { return a.bytes > b.bytes; });
      std::fflush(stdout);
      std::fprintf(stderr, "\n========================================================================\n");
      std::fprintf(stderr, " GPU MEMORY BUDGET: this run will not fit and would abort mid-allocation\n");
      std::fprintf(stderr, "------------------------------------------------------------------------\n");
      std::fprintf(stderr, " Projected device use: %s   (approximate estimate)\n", formatBytes(total).c_str());
      std::fprintf(stderr, " Free / total on GPU:  %.3f GB / %.3f GB\n",
                   static_cast<double>(freeB) / 1e9, static_cast<double>(totalB) / 1e9);
      std::fprintf(stderr, "------------------------------------------------------------------------\n");
      std::fprintf(stderr, " Top consumers:\n");
      for(const auto& l : lines) {
         if(l.bytes == 0) continue;
         std::fprintf(stderr, "   %-40s %s\n", l.name, formatBytes(l.bytes).c_str());
      }
      std::fprintf(stderr, "------------------------------------------------------------------------\n");
      std::fprintf(stderr, " To make it fit:\n");
      std::fprintf(stderr, "   * Rebuild in single precision to roughly halve real-valued arrays.\n");
      if(P.M > 1)
         std::fprintf(stderr, "   * Reduce Mensemble (now %zu); device memory scales with ensembles.\n", P.M);
      std::fprintf(stderr, "   * Reduce the system size (Natom now %zu).\n", P.N);
      std::fprintf(stderr, "   * Set UPPASD_GPU_SKIP_BUDGET=1 to bypass this check (estimate is\n");
      std::fprintf(stderr, "     approximate; the run may still run out).\n");
      std::fprintf(stderr, "========================================================================\n");
      std::fflush(stderr);
      std::exit(EXIT_FAILURE);
   }

   if(effectiveHaveInfo) {
      std::printf("Gpu: projected device use %s of %.3f GB free (%.3f GB total)%s\n",
                  formatBytes(total).c_str(),
                  static_cast<double>(freeB) / 1e9, static_cast<double>(totalB) / 1e9,
                  (overBudget && bypass) ? " [OVER BUDGET - bypassed via UPPASD_GPU_SKIP_BUDGET]" : "");
   } else {
      std::printf("Gpu: projected device use %s (free-memory query unavailable)\n",
                  formatBytes(total).c_str());
   }
   return total;
}

} // namespace

bool GpuSimulation::initiateMatrices(int is_mc) {
   runIsMC = (is_mc != 0);
   if(runIsMC && FortranData::gpu_dipole_mode && *FortranData::gpu_dipole_mode != 0) {
      throw std::runtime_error("GPU EWALD3D_FFT is not available with GPU Monte Carlo");
   }
   // Dimensions
   printf("Initiate matrices GPU -1 (is_mc=%d)\n", is_mc);
    long int N = static_cast <long int>( SimParam.N);
    long int NH = static_cast <long int>(SimParam.NH);
    long int M = static_cast <long int>( SimParam.M);
    long int mnn = static_cast <long int>( SimParam.mnn);
    long int mnndm = static_cast <long int>( SimParam.mnndm);

   // Constants initiated?
   if(N == 0 || M == 0 || NH == 0) {
      std::printf("GpuSimulation: constants not initiated!\n");
      std::exit(EXIT_FAILURE);
   }

   // initiate corresponding cpu matrices from Fortran
   initiate_fortran_cpu_matrices();

   // Initiated?
   if(isInitiated) {
      std::printf("GpuSimulation: attempted to initiate already initiated GpuSimulation!\n");
      std::exit(EXIT_FAILURE);
   }

   // Allocate
   // M1: project the full-run device footprint and abort here (before the first
   // Allocate) if it will not fit. Stored for the release()-time self-check.
   try {
      estimatedDeviceBytes = computeAndCheckDeviceBudget(Flags, SimParam, runIsMC);
   } catch(...) {
      FortranData::clearAdaptivePointers();
      throw;
   }



    gpuHamiltonian.aHam.Allocate(N);  
     if(Flags.do_jtensor != 0) {
        std::printf("\n GPU: jTensor has been initialized \n");
        gpuHamiltonian.j_tensor.Allocate(NH, static_cast <long int>(3), static_cast <long int>(3), mnn);
        gpuHamiltonian.nlist.Allocate(mnn, N);
        gpuHamiltonian.nlistsize.Allocate(NH);

    }
    else {
        
        gpuHamiltonian.ncoup.Allocate(NH, mnn);            
        gpuHamiltonian.nlist.Allocate(N, mnn);
        gpuHamiltonian.nlistsize.Allocate(NH);



        if(Flags.do_dm != 0) {
            gpuHamiltonian.dmvect.Allocate( static_cast <long int>(3), mnndm, NH);     
            gpuHamiltonian.dmlist.Allocate(mnndm, N);
            gpuHamiltonian.dmlistsize.Allocate(NH);

        }
    }


    if(Flags.do_aniso != 0) {
        gpuHamiltonian.kaniso.Allocate( static_cast <long int>(2), N);;
        gpuHamiltonian.eaniso.Allocate( static_cast <long int>(3), N);;
        gpuHamiltonian.taniso.Allocate(N);;
        gpuHamiltonian.sb.Allocate(N);

    }
    gpuHamiltonian.extfield.Allocate( static_cast <long int>(3), N, M);
    if(!cpuHamiltonian.macro_cell_index.empty()) {
        gpuHamiltonian.macro_cell_index.Allocate(N);
        gpuHamiltonian.macro_nlistsize.Allocate(cpuHamiltonian.macro_nlistsize.size());
        gpuHamiltonian.macro_center.Allocate(3, cpuHamiltonian.macro_nlistsize.size());
        gpuHamiltonian.macro_min_coord.Allocate(3, cpuHamiltonian.macro_nlistsize.size());
        gpuHamiltonian.macro_max_coord.Allocate(3, cpuHamiltonian.macro_nlistsize.size());
    }
    gpuLattice.beff.Allocate( static_cast <long int>(3), N, M);
    gpuLattice.b2eff.Allocate( static_cast <long int>(3), N, M);
    // M2: eneff is written by Heisge but read only by MC; for SD without aniso it
    // equals beff bit-for-bit, so alias it (no separate buffer) instead of allocating.
    eneffAliased = (!runIsMC && Flags.do_aniso == 0);
    if(eneffAliased) {
        gpuLattice.eneff = gpuLattice.beff;   // shallow: shares beff's device buffer
    } else {
        gpuLattice.eneff.Allocate( static_cast <long int>(3), N, M);
    }
    //gpuLattice.energy.Allocate(1);
    gpuLattice.emomM.Allocate( static_cast <long int>(3), N, M);
    gpuLattice.emom.Allocate( static_cast <long int>(3), N, M);
    gpuLattice.emom2.Allocate( static_cast <long int>(3), N, M);
    gpuLattice.mmom.Allocate(N, M);
    // M2: mmom0/mmom2 are used only by MC (moms) and the mompar!=0 MomentUpdater.
    const bool needMomExtra = runIsMC || (SimParam.mompar != 0);
    if(needMomExtra) {
        gpuLattice.mmom0.Allocate(N, M);
        gpuLattice.mmom2.Allocate(N, M);
    }
    gpuLattice.mmomi.Allocate(N, M);
    //gpuLattice.ipTemp_array.Allocate(N);

    if(!eneffAliased) gpuLattice.eneff.zeros();

    if(Flags.do_ene > 0) {
        //gpuEnergies.totalM.Allocate(M);
        //gpuEnergies.exchM.Allocate(M);
        //gpuEnergies.aniM.Allocate(M);
        //gpuEnergies.dmM.Allocate(M);
        //gpuEnergies.extM.Allocate(M);
        //gpuEnergies.tensorM.Allocate(M);
        gpuEnergies.energyM.Allocate(M, static_cast <long int>(7));

    }
    else{ 
        //gpuEnergies.totalM.Allocate(1);
        //gpuEnergies.exchM.Allocate(1);
        //gpuEnergies.aniM.Allocate(1);
        //gpuEnergies.dmM.Allocate(1);
        //gpuEnergies.extM.Allocate(1);
        //gpuEnergies.tensorM.Allocate(1);
        gpuEnergies.energyM.Allocate(1, 1);
    }

        //gpuEnergies.totalM.zeros();
       // gpuEnergies.exchM.zeros();
       // gpuEnergies.aniM.zeros();
        //gpuEnergies.dmM.zeros();
        //gpuEnergies.extM.zeros();
        //gpuEnergies.tensorM.zeros();
        gpuEnergies.energyM.zeros();



    //gpuLattice.temperature.initiate(N); //is initiated if we run SD or MC simulation inside corresponding classes where they are requires
    if(FortranData::btorque) {gpuLattice.btorque.Allocate(static_cast <long int>(3), N, M);}
    if(FortranData::adaptive_geometry_mode) {
        try {
            adaptiveMaskEnabled =
               FortranData::adaptive_mask_mode && *FortranData::adaptive_mask_mode != 0;
            adaptiveUpdateInterval = FortranData::adaptive_update_interval ?
               *FortranData::adaptive_update_interval : 1;
            if(FortranData::adaptive_refine_threshold)
               adaptiveSelectorPolicy.refineThreshold =
                  static_cast<real>(*FortranData::adaptive_refine_threshold);
            if(FortranData::adaptive_coarsen_threshold)
               adaptiveSelectorPolicy.coarsenThreshold =
                  static_cast<real>(*FortranData::adaptive_coarsen_threshold);
            if(FortranData::adaptive_minimum_dwell)
               adaptiveSelectorPolicy.minimumDwellUpdates =
                  *FortranData::adaptive_minimum_dwell;
            if(FortranData::adaptive_buffer_dilation)
               adaptiveSelectorPolicy.bufferDilationBlocks =
                  *FortranData::adaptive_buffer_dilation;
            if(FortranData::adaptive_reconstruction_scheme &&
               *FortranData::adaptive_reconstruction_scheme == 2)
               adaptiveReconstructionPolicy.scheme =
                  GpuAdaptiveReconstruction::ConstrainedCone;
            if(FortranData::adaptive_cone_angle_rad)
               adaptiveReconstructionPolicy.coneAngleRadians =
                  static_cast<real>(*FortranData::adaptive_cone_angle_rad);
            adaptiveDiagnostics = FortranData::adaptive_diagnostics ?
               *FortranData::adaptive_diagnostics : 0;
            adaptiveEnergyJumpLimitJ = FortranData::adaptive_energy_jump_limit_j ?
               static_cast<double>(*FortranData::adaptive_energy_jump_limit_j) : 0.0;
            gpuAdaptiveRuntime.initialize(adaptiveTopologyInput(), adaptiveRuntimeInput(),
                                          SimParam.N, SimParam.M);
            const auto work = gpuAdaptiveRuntime.downloadWorkSnapshot();
            std::printf("Gpu: AdaptiveCG initial active_atoms=%zu active_blocks=%zu "
                        "interface_atoms=%zu device_bytes=%zu\n",
                        work.activeAtoms.size(), work.activeBlocks.size(),
                        work.interfaceAtoms.size(),
                        gpuAdaptiveRuntime.allocatedBytes());
            FortranData::clearAdaptivePointers();
        } catch(...) {
            FortranData::clearAdaptivePointers();
            throw;
        }
    }
    //e.Allocate( static_cast <long int>(3), N, M);} 

   /* if (Flags.do_mphase_now != 0){
        if (Flags.do_avrg !=0){
            gpuMeasurables.mavg_buff.Allocate(N, M);
            //gpuMeasurables.eavg_buff.Allocate(N, M);
        }
        if (Flags.do_cumu !=0){
            gpuMeasurables.mcumu_buff.Allocate(N, M);
            //gpuMeasurables.ecumu_buff.Allocate(N, M); 
        }
    }*/

   // Did we get the memory?
    if(gpuHasNoData()){
      FortranData::clearAdaptivePointers();
      release();
      // Check for error
      const char* err = GPU_GET_ERROR_STRING(GPU_GET_LAST_ERROR());
      std::fprintf(stderr, "Gpu: Failed to allocate memory: %s\n", err);
      return false;
   }

   // Flag that we're initiated
   isInitiated = true;
   isFreed = false;
   // Initiate data
   copyFromFortran();
   // Post-init report (V3 part 2): what initiateMatrices actually placed on device.
   {
      std::size_t freeB = 0, totalB = 0;
      const int64_t allocated = TensorMemoryTracker::current_device();
      if(GPU_MEM_GET_INFO(&freeB, &totalB) == GPU_SUCCESS)
         std::printf("Gpu: device base inventory after init: %lld B; %.3f GB free / %.3f GB total\n",
                     static_cast<long long>(allocated),
                     static_cast<double>(freeB) / 1e9, static_cast<double>(totalB) / 1e9);
   }
   return true;
}

bool GpuSimulation::gpuHasNoData(){
    bool check = (  gpuHamiltonian.aHam.empty() ||                             
                    (gpuHamiltonian.ncoup.empty() && (FortranData::j_tensor == nullptr))||            
                    gpuHamiltonian.nlist.empty() || 
                    gpuHamiltonian.nlistsize.empty() ||
                    (gpuHamiltonian.dmvect.empty() && (FortranData::dmvect != nullptr)) ||
                    (gpuHamiltonian.dmlist.empty() && (FortranData::dmlist != nullptr)) ||    
                    (gpuHamiltonian.dmlistsize.empty() && (FortranData::dmlistsize != nullptr)) ||    
                    (gpuHamiltonian.j_tensor.empty() && (FortranData::j_tensor != nullptr)) ||     
                    (gpuHamiltonian.kaniso.empty() && (FortranData::kaniso != nullptr)) ||     
                    (gpuHamiltonian.eaniso.empty() && (FortranData::eaniso != nullptr)) ||     
                    (gpuHamiltonian.taniso.empty() && (FortranData::taniso != nullptr)) ||     
                    (gpuHamiltonian.sb.empty() && (FortranData::sb != nullptr)) ||     
                    gpuHamiltonian.extfield.empty() || 
                    gpuLattice.beff.empty() || 
                    gpuLattice.b2eff.empty() ||
                    gpuLattice.eneff.empty() ||
                    //gpuLattice.energy.empty() ||
                    gpuLattice.emomM.empty() || 
                    gpuLattice.emom.empty() || 
                    gpuLattice.emom2.empty() || 
                    gpuLattice.mmom.empty() || 
                    (gpuLattice.mmom0.empty() && (runIsMC || SimParam.mompar != 0)) || 
                    (gpuLattice.mmom2.empty() && (runIsMC || SimParam.mompar != 0)) || 
                    gpuLattice.mmomi.empty() ||
                   // gpuLattice.ipTemp.empty()||
                    (gpuLattice.btorque.empty()&& (FortranData::btorque != nullptr)));
    if(FortranData::adaptive_geometry_mode) check = check || !gpuAdaptiveRuntime.ready();
    //TODO: add measurables
    return check;
}

void GpuSimulation::release() {
    if(isInitiated && !isFreed)
    isInitiated = false;
    isFreed = true;
    gpuHamiltonian.neighbourListsPrepared = false;
    if(gpuAdaptiveRuntime.ready()) {
      const auto work = gpuAdaptiveRuntime.downloadWorkSnapshot();
      const auto& phase = gpuAdaptiveRuntime.phaseMetrics();
      std::printf("Gpu: AdaptiveCG final active_atoms=%zu active_blocks=%zu "
                  "interface_atoms=%zu phases_ms(atom=%.3f coarse=%.3f "
                  "interface=%.3f selector=%.3f compaction=%.3f integration=%.3f)\n",
                  work.activeAtoms.size(), work.activeBlocks.size(),
                  work.interfaceAtoms.size(), phase.atomisticMilliseconds,
                  phase.coarseMilliseconds, phase.interfaceMilliseconds,
                  phase.selectorMilliseconds, phase.compactionMilliseconds,
                  phase.integrationMilliseconds);
      if(adaptiveDiagnostics > 0) {
         const auto diagnostic =
            gpuAdaptiveRuntime.diagnosticSnapshot(gpuLattice.emom.data());
         std::printf("Gpu: AdaptiveCG resolved diagnostics=%d energy_jump_limit_j=%.16e "
                     "fft_source_mapping=basis-resolved-to-single-dynamical-channel\n",
                     adaptiveDiagnostics, adaptiveEnergyJumpLimitJ);
         std::printf("Gpu: AdaptiveCG final_state values=");
         for(std::size_t block = 0; block < diagnostic.blockState.size(); ++block)
            std::printf("%s%d", block == 0 ? "" : ",",
                        diagnostic.blockState[block]);
         std::printf("\n");
         const auto countState = [&](int state) {
            return static_cast<unsigned long long>(std::count(
               diagnostic.blockState.begin(), diagnostic.blockState.end(), state));
         };
         unsigned long long accepted = 0;
         for(const auto epoch : diagnostic.transitionEpoch) accepted += epoch;
         std::printf("Gpu: AdaptiveCG resolution_counts fine=%llu interface=%llu "
                     "coarse=%llu accepted_transitions=%llu rejected_transitions=0\n",
                     countState(2), countState(1), countState(0), accepted);
         std::printf("Gpu: AdaptiveCG last_energy_j atomistic_bilinear=%.16e "
                     "atomistic_onsite=%.16e "
                     "coarse_exchange=%.16e coarse_spiralization=%.16e "
                     "coarse_anisotropy=%.16e coarse_external=%.16e "
                     "coarse_dipole=%.16e total=%.16e\n",
                     diagnostic.energy.atomisticBilinearJ,
                     diagnostic.energy.atomisticOnsiteJ,
                     diagnostic.energy.coarseExchangeJ,
                     diagnostic.energy.coarseSpiralizationJ,
                     diagnostic.energy.coarseAnisotropyJ,
                     diagnostic.energy.coarseExternalJ,
                     diagnostic.energy.dipoleJ, diagnostic.energy.totalJ);
         std::printf("Gpu: AdaptiveCG last_field_checksums_t atom_sum=%.16e "
                     "atom_norm2=%.16e coarse_sum=%.16e coarse_norm2=%.16e\n",
                     diagnostic.atomFieldSumT, diagnostic.atomFieldNorm2T2,
                     diagnostic.coarseFieldSumT,
                     diagnostic.coarseFieldNorm2T2);
         std::printf("Gpu: AdaptiveCG trajectory_checksums direction_sum=%.16e "
                     "direction_norm2=%.16e\n",
                     diagnostic.directionSum, diagnostic.directionNorm2);
         std::printf("Gpu: AdaptiveCG phase_ms field_atom=%.6f field_coarse=%.6f "
                     "interface=%.6f selector=%.6f compaction=%.6f fft=%.6f "
                     "integration=%.6f device_bytes=%zu\n",
                     phase.atomisticMilliseconds, phase.coarseMilliseconds,
                     phase.interfaceMilliseconds, phase.selectorMilliseconds,
                     phase.compactionMilliseconds, phase.fftMilliseconds,
                     phase.integrationMilliseconds,
                     gpuAdaptiveRuntime.allocatedBytes());
      }
      gpuAdaptiveRuntime.release();
    }
    gpuHamiltonian.aHam.Free();
         
    gpuHamiltonian.nlist.Free();  
    gpuHamiltonian.nlistsize.Free(); 

    gpuEnergies.energyM.Free();
    //gpuEnergies.totalM.Free();
    //gpuEnergies.extM.Free();
    //gpuEnergies.tensorM.Free();
    //gpuEnergies.exchM.Free();
    //gpuEnergies.dmM.Free();
    //gpuEnergies.aniM.Free();

    if(Flags.do_jtensor != 0) {
        gpuHamiltonian.j_tensor.Free();
    }
    else {
        gpuHamiltonian.ncoup.Free();

        if(Flags.do_dm != 0) { 
            gpuHamiltonian.dmvect.Free();  
            gpuHamiltonian.dmlist.Free();     
            gpuHamiltonian.dmlistsize.Free();

        }
    }

    
    if(Flags.do_aniso != 0) {
        gpuHamiltonian.kaniso.Free();    
        gpuHamiltonian.eaniso.Free();     
        gpuHamiltonian.taniso.Free(); 
        gpuHamiltonian.sb.Free();

    }   
     
    gpuHamiltonian.extfield.Free();
    gpuHamiltonian.macro_cell_index.Free();
    gpuHamiltonian.macro_nlistsize.Free();
    gpuHamiltonian.macro_center.Free();
    gpuHamiltonian.macro_min_coord.Free();
    gpuHamiltonian.macro_max_coord.Free();
    gpuLattice.beff.Free();  
    gpuLattice.b2eff.Free();   
    if(eneffAliased) { gpuLattice.eneff = GpuTensor<real, 3>{}; eneffAliased = false; }
    else            { gpuLattice.eneff.Free(); }
    //gpuLattice.energy.Free();
    gpuLattice.emomM.Free();  
    gpuLattice.emom.Free();  
    gpuLattice.emom2.Free();   
    gpuLattice.mmom.Free();  
    if(!gpuLattice.mmom0.empty()) gpuLattice.mmom0.Free();
    if(!gpuLattice.mmom2.empty()) gpuLattice.mmom2.Free();
    gpuLattice.mmomi.Free();
   // gpuLattice.ipTemp.Free();
     if(FortranData::btorque) {gpuLattice.btorque.Free();  }

    

   // gpuMeasurables.mavg_buff.Free();  
   // gpuMeasurables.mcumu_buff.Free();  
  
    if(FortranData::gpu_dipole_mode && *FortranData::gpu_dipole_mode == 2) {
      std::printf("Gpu: OPEN_FFT release device inventory: %lld B.\n",
                  static_cast<long long>(TensorMemoryTracker::current_device()));
    }
    TensorMemoryTracker::printResults();
    // The tracker peak covers every phase-specific allocation made over the
    // process lifetime; print its scope explicitly rather than comparing it
    // to an unrelated single-component diagnostic.
    if(estimatedDeviceBytes > 0) {
      const int64_t peak = TensorMemoryTracker::peak_device();
      if(peak > 0) {
         const double est = static_cast<double>(estimatedDeviceBytes);
         const double pk  = static_cast<double>(peak);
         std::printf("Gpu: full-run planned device inventory %zu B vs TensorMemoryTracker process peak %lld B (estimate %+.1f%%).\n",
                     estimatedDeviceBytes, static_cast<long long>(peak), 100.0 * (est - pk) / pk);
      }
    }
    // TensorMemoryTracker::saveToFile();
    TensorDataMovementTracker::printResults();
}

void GpuSimulation::updateAdaptiveBlockState(const int* blockState, std::size_t count) {
   gpuAdaptiveRuntime.updateBlockState(blockState, count);
   const auto& metrics = gpuAdaptiveRuntime.compactionMetrics();
   std::printf("Gpu: adaptive selector compaction sync #%llu: %zu block bytes, "
               "%.3f ms wall / %.3f ms host wait / %.3f ms device cumulative.\n",
               static_cast<unsigned long long>(metrics.hostSynchronizations),
               count * sizeof(int), metrics.elapsedMilliseconds,
               metrics.hostWaitMilliseconds, metrics.deviceMilliseconds);
}

void GpuSimulation::advanceAdaptiveStep(
   std::size_t step, GpuHamiltonianCalculations* hamiltonian) {
   if(!gpuAdaptiveRuntime.ready() || !gpuAdaptiveRuntime.kernelsReady())
      throw std::logic_error("GPU adaptive production step requires a ready CG-10 runtime");
   GpuAdaptiveFftEvaluator fftEvaluator{};
   if(hamiltonian && hamiltonian->hasAdaptiveFftDipole()) {
      fftEvaluator = [this, hamiltonian](const real* direction) {
         const auto started = std::chrono::steady_clock::now();
         auto view = hamiltonian->evaluateAdaptiveFftDipole(
            direction, gpuLattice.mmom,
            gpuAdaptiveRuntime.deviceTopology(),
            gpuAdaptiveRuntime.deviceRuntime());
         ASSERT_GPU(GPU_STREAM_SYNC(
            ParallelizationHelperInstance.getWorkStream()));
         const auto stopped = std::chrono::steady_clock::now();
         gpuAdaptiveRuntime.recordFftMilliseconds(
            std::chrono::duration<double, std::milli>(
               stopped - started).count());
         return view;
      };
   }
   gpuAdaptiveRuntime.integrateHeun(
      static_cast<real>(SimParam.delta_t), gpuLattice.emom.data(),
      nullptr, nullptr, fftEvaluator);
   gpuAdaptiveRuntime.synchronizeAtomicState(
      gpuLattice.emom.data(), adaptiveReconstructionPolicy);
   if(adaptiveMaskEnabled && adaptiveUpdateInterval > 0 &&
      step % adaptiveUpdateInterval == 0) {
      gpuAdaptiveRuntime.restrictMoments(gpuLattice.emom.data());
      gpuAdaptiveRuntime.evaluateSelectorScores(gpuLattice.emom.data());
      gpuAdaptiveRuntime.proposeSelectorState(adaptiveSelectorPolicy);
      gpuAdaptiveRuntime.publishProposedState(
         gpuLattice.emom.data(), adaptiveReconstructionPolicy, true);
      gpuAdaptiveRuntime.synchronizeAtomicState(
         gpuLattice.emom.data(), adaptiveReconstructionPolicy);
   }
   ASSERT_GPU(GPU_STREAM_SYNC(gpuAdaptiveRuntime.stream()));
}

void GpuSimulation::copyFromFortran() {
   if(isInitiated) {
   //printf("HERE - 5\n");

    gpuHamiltonian.neighbourListsPrepared = false;

    gpuHamiltonian.aHam.copy_sync(cpuHamiltonian.aHam);
    if(Flags.do_jtensor != 0) {
        cpuHamiltonian.j_tensor.tensor_exchange_copy_to(gpuHamiltonian.j_tensor);
        gpuHamiltonian.nlist.copy_sync(cpuHamiltonian.nlist);  
        gpuHamiltonian.nlistsize.copy_sync(cpuHamiltonian.nlistsize);         
    }
    else {
        cpuHamiltonian.ncoup.transposed_copy_to(gpuHamiltonian.ncoup);
        cpuHamiltonian.nlist.transposed_copy_to(gpuHamiltonian.nlist);
        gpuHamiltonian.nlistsize.copy_sync(cpuHamiltonian.nlistsize);
        if(Flags.do_dm != 0) {
            gpuHamiltonian.dmvect.copy_sync(cpuHamiltonian.dmvect);  
            gpuHamiltonian.dmlist.copy_sync(cpuHamiltonian.dmlist);     
            gpuHamiltonian.dmlistsize.copy_sync(cpuHamiltonian.dmlistsize);  }
    }

    if(Flags.do_aniso != 0) {
        gpuHamiltonian.kaniso.copy_sync(cpuHamiltonian.kaniso);    
        gpuHamiltonian.eaniso.copy_sync(cpuHamiltonian.eaniso);     
        gpuHamiltonian.taniso.copy_sync(cpuHamiltonian.taniso);  
        gpuHamiltonian.sb.copy_sync(cpuHamiltonian.sb);  }     
    gpuHamiltonian.extfield.copy_sync(cpuHamiltonian.extfield);
    if(!gpuHamiltonian.macro_cell_index.empty()) {
        gpuHamiltonian.macro_cell_index.copy_sync(cpuHamiltonian.macro_cell_index);
        gpuHamiltonian.macro_nlistsize.copy_sync(cpuHamiltonian.macro_nlistsize);
        gpuHamiltonian.macro_center.copy_sync(cpuHamiltonian.macro_center);
        gpuHamiltonian.macro_min_coord.copy_sync(cpuHamiltonian.macro_min_coord);
        gpuHamiltonian.macro_max_coord.copy_sync(cpuHamiltonian.macro_max_coord);
    }
   // printf("HERE - 6\n");
    gpuLattice.beff.copy_sync(cpuLattice.beff);  
    gpuLattice.b2eff.copy_sync(cpuLattice.b2eff);   
    gpuLattice.emomM.copy_sync(cpuLattice.emomM);
    //printf("HERE - 7\n");  
    gpuLattice.emom.copy_sync(cpuLattice.emom);  
    gpuLattice.emom2.copy_sync(cpuLattice.emom2);   
    gpuLattice.mmom.copy_sync(cpuLattice.mmom);    
    //printf("HERE - 8\n");  
    if(!gpuLattice.mmom0.empty()) gpuLattice.mmom0.copy_sync(cpuLattice.mmom0);
    if(!gpuLattice.mmom2.empty()) gpuLattice.mmom2.copy_sync(cpuLattice.mmom2);
    gpuLattice.mmomi.copy_sync(cpuLattice.mmomi);  
    //printf("HERE - 9\n");
   // gpuLattice.ipTemp.copy_sync(cpuLattice.ipTemp);
   // printf("HERE - 10\n");
    if(FortranData::btorque) {gpuLattice.btorque.copy_sync(cpuLattice.btorque); } 

   // gpuMeasurables.mavg_buff.copy_sync(cpuMeasurables.mavg_buff);  
   // gpuMeasurables.mcumu_buff.copy_sync(cpuMeasurables.mcumu_buff);
   }
}
void GpuSimulation::copyToFortran() {
   if(isInitiated) {
   // printf("COPIED\n");
    cpuLattice.beff.copy_sync(gpuLattice.beff);  
    cpuLattice.b2eff.copy_sync(gpuLattice.b2eff);   
    cpuLattice.emomM.copy_sync(gpuLattice.emomM);  
    cpuLattice.emom.copy_sync(gpuLattice.emom);  
    cpuLattice.emom2.copy_sync(gpuLattice.emom2);   
    cpuLattice.mmom.copy_sync(gpuLattice.mmom);  
    if(!gpuLattice.mmom0.empty()) cpuLattice.mmom0.copy_sync(gpuLattice.mmom0);
    if(!gpuLattice.mmom2.empty()) cpuLattice.mmom2.copy_sync(gpuLattice.mmom2);
    cpuLattice.mmomi.copy_sync(gpuLattice.mmomi);
    // cpuLattice owns precision-converted staging in SINGLE_PREC builds.
    // Convert results back into the double-precision Fortran buffers.
    cpuLattice.beff.copy_to(FortranData::beff);
    cpuLattice.b2eff.copy_to(FortranData::b2eff);
    cpuLattice.emomM.copy_to(FortranData::emomM);
    cpuLattice.emom.copy_to(FortranData::emom);
    cpuLattice.emom2.copy_to(FortranData::emom2);
    cpuLattice.mmom.copy_to(FortranData::mmom);
    if(!gpuLattice.mmom0.empty()) cpuLattice.mmom0.copy_to(FortranData::mmom0);
    if(!gpuLattice.mmom2.empty()) cpuLattice.mmom2.copy_to(FortranData::mmom2);
    cpuLattice.mmomi.copy_to(FortranData::mmomi);
    if(FortranData::btorque) {
        cpuLattice.btorque.copy_to(FortranData::btorque);
    }
   // gpuMeasurables.mavg_buff.copy_sync(cpuMeasurables.mavg_buff);
    //gpuMeasurables.mcumu_buff.copy_sync(cpuMeasurables.mcumu_buff);
   }
}


void GpuSimulation::gpuRunSimulation(const int whichsim, const int whichphase, const char bf){
printf("current type %i\n", whichsim);
    if(whichsim == 0){
        GpuSDSimulation GpuSD;
        if(whichphase == 0) {
            GpuSD.SDiphase(*this);
            copyToFortran();
        }
        else if(whichphase == 1) {
            // initiateMatrices() has just created this GPU session and copied the
            // current Fortran state.  Do not upload it again here: in particular,
            // the full neighbour list is immutable during a phase and can be very
            // large for long-ranged Hamiltonians.
            GpuSD.SDmphase(*this);
        }
        else {printf("Wrong phase! 0 - initial, 1 - measurement");}
    }
    else if(whichsim == 1){
         GpuMCSimulation GpuMC;//TODO

        if(whichphase == 0) {
            if(bf == 'Y') GpuMC.MCiphase_bf(*this);
            else GpuMC.MCiphase(*this);
            copyToFortran();
        }
        else if(whichphase == 1) {
            if(bf == 'Y') GpuMC.MCmphase_bf(*this);
            else GpuMC.MCmphase(*this);
        }
        else {printf("Wrong phase! 0 - initial, 1 - measurement");}
    }
    else {printf("Wrong simulation type! 0 - SD, 1 - MC; current type %i\n", whichsim);}
    //release();
}
