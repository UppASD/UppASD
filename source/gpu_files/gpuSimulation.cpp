

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
#include "gpuLatticeConvolutionHamiltonian.hpp"
#include "correlations/gpuCorrelations.hpp"
#include <algorithm>
#include <cstdlib>
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

    SimParam.SDEalgh = *FortranData::SDEalgh;
    if(!(SimParam.SDEalgh == 1 || SimParam.SDEalgh == 4 || SimParam.SDEalgh == 5 || SimParam.SDEalgh == 11)) {
        std::fprintf(stderr, "Invalid SDEalgh!\n");
        std::exit(EXIT_FAILURE);
    }

    Flags.do_dm = static_cast<bool>(*FortranData::do_dm);
    Flags.do_jtensor = static_cast<bool>(*FortranData::do_jtensor);
    Flags.do_aniso = static_cast<bool>(*FortranData::do_aniso);

    Flags.do_sc = *FortranData::do_sc;
    Flags.do_gpu_correlations = static_cast<bool>(*FortranData::do_gpu_correlations);
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
    cpuLattice.beff.set(FortranData::beff,  static_cast <long int>(3), N, M);
    cpuLattice.b2eff.set(FortranData::b2eff,  static_cast <long int>(3), N, M);
    cpuLattice.emomM.set(FortranData::emomM,  static_cast <long int>(3), N, M);
    cpuLattice.emom.set(FortranData::emom,  static_cast <long int>(3), N, M);
    cpuLattice.emom2.set(FortranData::emom2,  static_cast <long int>(3), N, M);
    cpuLattice.mmom.set(FortranData::mmom, N, M);
    cpuLattice.mmom0.set(FortranData::mmom0, N, M);
    cpuLattice.mmom2.set(FortranData::mmom2, N, M);
    cpuLattice.mmomi.set(FortranData::mmomi, N, M);
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

        cpuCorrelations.m_k.set(FortranData::m_k, static_cast <long int>(3), nq);
        cpuCorrelations.m_kt.set(FortranData::m_kt, static_cast <long int>(3), static_cast <long int>(nq), static_cast <long int>(SimParam.sc_max_nstep));
        cpuCorrelations.m_kw.set(FortranData::m_kw, static_cast <long int>(3), nq, nw);
        cpuCorrelations.deltat_corr.set(FortranData::deltat_corr, static_cast <long int>(SimParam.sc_max_nstep + 1));
        cpuCorrelations.scstep_arr.set(FortranData::scstep_arr, static_cast <long int>(SimParam.sc_max_nstep + 1));
        if((Flags.do_sc_proj=='C')||(Flags.do_sc_proj=='Q')||(Flags.do_sc_proj=='Y')){          
        cpuCorrelations.atype.set(FortranData::atype, N);
        cpuCorrelations.m_k_proj.set(FortranData::m_k_proj, static_cast <long int>(3), NT, nq);
        cpuCorrelations.m_kt_proj.set(FortranData::m_kt_proj, static_cast <long int>(3), NT, nq, static_cast <long int>(SimParam.sc_max_nstep));
        cpuCorrelations.m_kw_proj.set(FortranData::m_kw_proj, static_cast <long int>(3), NT, nq, nw);
        }
        if((Flags.do_sc_projch=='C')||(Flags.do_sc_projch=='Q')||(Flags.do_sc_projch=='Y')){   
        cpuCorrelations.achtype.set(FortranData::achtype, N);
        cpuCorrelations.m_k_projch.set(FortranData::m_k_projch, static_cast <long int>(3), Nchmax, nq);
        cpuCorrelations.m_kt_projch.set(FortranData::m_kt_projch, static_cast <long int>(3), Nchmax, nq, static_cast <long int>(SimParam.sc_max_nstep));
        cpuCorrelations.m_kw_projch.set(FortranData::m_kw_projch, static_cast <long int>(3), Nchmax, nq, nw);
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
   return (F.do_ene > 0 ? P.M * 6 : static_cast<std::size_t>(1)) * sizeof(real); // energyM
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

// Sum every device allocation the run will make, compare to free VRAM, and abort
// with a table before the first Allocate if it will not fit. Returns the total.
std::size_t computeAndCheckDeviceBudget(const Flag& F, const SimulationParameters& P, bool is_mc) {
   std::vector<BudgetLine> lines = {
      {"Hamiltonian (Jij/DM/aniso) + ext field", hamiltonianBytes(F, P)},
      {"Lattice streaming arrays",               latticeBytes(F, P, is_mc)},
      {"Depondt integrator + thermfield",        GpuDepondtIntegrator::estimateBytes(P)},
      {"Energy buffers",                         energiesBytes(F, P)},
      {"FFT convolution",                        willUseConvolution(P) ? GpuLatticeConvolutionHamiltonian::estimateBytes(P) : static_cast<std::size_t>(0)},
      {"Correlations",                           (FortranData::do_gpu_correlations && *FortranData::do_gpu_correlations == 'Y') ? GpuCorrelations::estimateBytes(F, P) : static_cast<std::size_t>(0)},
   };
   std::size_t total = 0;
   for(const auto& l : lines) total += l.bytes;

   std::size_t freeB = 0, totalB = 0;
   const bool haveInfo = (GPU_MEM_GET_INFO(&freeB, &totalB) == GPU_SUCCESS);
   const char* skip = std::getenv("UPPASD_GPU_SKIP_BUDGET");
   const bool bypass = (skip != nullptr && skip[0] != '\0');
   const bool overBudget = haveInfo && (static_cast<double>(total) > 0.90 * static_cast<double>(freeB));

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
      std::fprintf(stderr, "     approximate and excludes FFT workspace; the run may still run out).\n");
      std::fprintf(stderr, "========================================================================\n");
      std::fflush(stderr);
      std::exit(EXIT_FAILURE);
   }

   if(haveInfo) {
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
   estimatedDeviceBytes = computeAndCheckDeviceBudget(Flags, SimParam, runIsMC);



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
        gpuEnergies.energyM.Allocate(M, static_cast <long int>(6));

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
      const int64_t allocated = TensorMemoryTracker::peak_device();
      if(GPU_MEM_GET_INFO(&freeB, &totalB) == GPU_SUCCESS)
         std::printf("Gpu: device allocated %.3f GB after init; %.3f GB free / %.3f GB total\n",
                     static_cast<double>(allocated) / 1e9,
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
    //TODO: add measurables
    return check;
}

void GpuSimulation::release() {
    if(isInitiated && !isFreed)
    isInitiated = false;
    isFreed = true;
    gpuHamiltonian.neighbourListsPrepared = false;
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
  
    TensorMemoryTracker::printResults();
    // M1: 5% self-check - upfront estimate vs the peak the tracker actually saw.
    if(estimatedDeviceBytes > 0) {
       const int64_t peak = TensorMemoryTracker::peak_device();
       if(peak > 0) {
          const double est = static_cast<double>(estimatedDeviceBytes);
          const double pk  = static_cast<double>(peak);
          std::printf("Gpu: device estimate %.3f GB vs measured peak %.3f GB (estimate %+.1f%%)\n",
                      est / 1e9, pk / 1e9, 100.0 * (est - pk) / pk);
       }
    }
    // TensorMemoryTracker::saveToFile();
    TensorDataMovementTracker::printResults();
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
            // Ensure measurement phase starts from the latest moments stored on the Fortran side.
            // This is important when initial and measurement phases are executed in separate GPU sessions.
            copyFromFortran();
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
