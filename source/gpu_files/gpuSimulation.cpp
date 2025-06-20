

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
#if defined(HIP_V)
#include <hiprand/hiprand.h>
#elif defined(CUDA_V)
#include <curand.h>
#endif

//#include "gpuSDSimulation.hpp"
GpuSimulation::GpuSimulation() {
    isInitiated = false;
}

GpuSimulation::~GpuSimulation() {
    if(isInitiated && !isFreed) {release();}
}

void GpuSimulation::initiateConstants() {
    SimParam.SDEalgh = *FortranData::SDEalgh;
    if(!(SimParam.SDEalgh == 1 || SimParam.SDEalgh == 4 || SimParam.SDEalgh == 5 || SimParam.SDEalgh == 11)) {
        std::fprintf(stderr, "Invalid SDEalgh!\n");
        std::exit(EXIT_FAILURE);
    }
    Flags.do_dm = static_cast<bool>(*FortranData::do_dm);
    Flags.do_jtensor = static_cast<bool>(*FortranData::do_jtensor);
    Flags.do_aniso = static_cast<bool>(*FortranData::do_aniso);
    //Flags.do_avrg = static_cast<bool>(*FortranData::do_avrg);
    //Flags.do_cumu = static_cast<bool>(*FortranData::do_cumu);

    SimParam.N = *FortranData::Natom;
    SimParam.NH  = *FortranData::nHam;
    SimParam.M = *FortranData::Mensemble;
    SimParam.mnn = *FortranData::max_no_neigh;
    SimParam.mnndm = *FortranData::max_no_dmneigh;
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
    SimParam.Temp = *FortranData::Temp;
    SimParam.damping = *FortranData::damping;
    SimParam.mompar = *FortranData::mompar;
    SimParam.initexc = *FortranData::initexc;

    SimParam.binderc = FortranData::binderc;
    SimParam.mavg = FortranData::mavg;

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

bool GpuSimulation::initiateMatrices() {
   // Dimensions
    long int N = static_cast <long int>( SimParam.N);
    long int NH = static_cast <long int>(SimParam.NH);
    long int M = static_cast <long int>( SimParam.M);
    long int mnn = static_cast <long int>( SimParam.mnn);
    long int mnndm = static_cast <long int>( SimParam.mnndm);
    long int ipnphase = static_cast <long int>( SimParam.ipnphase);
    long int ipmcnphase = static_cast <long int>( SimParam.ipmcnphase);

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
    gpuHamiltonian.aHam.Allocate(N);  
     if(Flags.do_jtensor != 0) {
        std::printf("\n GPU: jTensor has been initialized \n");
        gpuHamiltonian.j_tensor.Allocate( static_cast <long int>(3),  static_cast <long int>(3), mnn, NH);        
         gpuHamiltonian.nlist.Allocate(mnn, N);
         gpuHamiltonian.nlistsize.Allocate(NH);}
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
    gpuLattice.eneff.Allocate( static_cast <long int>(3), N, M);
    gpuLattice.emomM.Allocate( static_cast <long int>(3), N, M);
    gpuLattice.emom.Allocate( static_cast <long int>(3), N, M);
    gpuLattice.emom2.Allocate( static_cast <long int>(3), N, M);
    gpuLattice.mmom.Allocate(N, M);
    gpuLattice.mmom0.Allocate(N, M);
    gpuLattice.mmom2.Allocate(N, M);
    gpuLattice.mmomi.Allocate(N, M);
    //gpuLattice.ipTemp_array.Allocate(N);

    gpuLattice.eneff.zeros();


   
    //gpuLattice.temperature.initiate(N); //is initiated if we run SD or MC simulation inside corresponding classes where they are requires
    if(FortranData::btorque) {gpuLattice.btorque.Allocate( static_cast <long int>(3), N, M);} 

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
                    gpuLattice.emomM.empty() || 
                    gpuLattice.emom.empty() || 
                    gpuLattice.emom2.empty() || 
                    gpuLattice.mmom.empty() || 
                    gpuLattice.mmom0.empty() || 
                    gpuLattice.mmom2.empty() || 
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
    gpuHamiltonian.aHam.Free();                         
         
    gpuHamiltonian.nlist.Free();  
    gpuHamiltonian.nlistsize.Free(); 

  if(Flags.do_dm != 0) { gpuHamiltonian.dmvect.Free();  
    gpuHamiltonian.dmlist.Free();     
    gpuHamiltonian.dmlistsize.Free();  }
     if(Flags.do_jtensor != 0) {gpuHamiltonian.j_tensor.Free();}
     else {gpuHamiltonian.ncoup.Free();}

    
    if(Flags.do_aniso != 0) {
        gpuHamiltonian.kaniso.Free();    
    gpuHamiltonian.eaniso.Free();     
    gpuHamiltonian.taniso.Free(); 
    gpuHamiltonian.sb.Free();}   
     
   gpuHamiltonian.extfield.Free();  
    gpuLattice.beff.Free();  
    gpuLattice.b2eff.Free();   
    gpuLattice.eneff.Free();   
    gpuLattice.emomM.Free();  
    gpuLattice.emom.Free();  
    gpuLattice.emom2.Free();   
    gpuLattice.mmom.Free();  
    gpuLattice.mmom0.Free();  
    gpuLattice.mmom2.Free();  
    gpuLattice.mmomi.Free();
   // gpuLattice.ipTemp.Free();
     if(FortranData::btorque) {gpuLattice.btorque.Free();  }

    

   // gpuMeasurables.mavg_buff.Free();  
   // gpuMeasurables.mcumu_buff.Free();  
  

}

void GpuSimulation::copyFromFortran() {
   if(isInitiated) {
   //printf("HERE - 5\n");

    gpuHamiltonian.aHam.copy_sync(cpuHamiltonian.aHam);     
    if(Flags.do_jtensor != 0) {
        gpuHamiltonian.j_tensor.copy_sync(cpuHamiltonian.j_tensor);            
        gpuHamiltonian.nlist.copy_sync(cpuHamiltonian.nlist);  
        gpuHamiltonian.nlistsize.copy_sync(cpuHamiltonian.nlistsize);         
    }     
    else { 
        cpuHamiltonian.ncoup.transpose();
        cpuHamiltonian.nlist.transpose();                      
        gpuHamiltonian.ncoup.copy_sync(cpuHamiltonian.ncoup);            
        gpuHamiltonian.nlist.copy_sync(cpuHamiltonian.nlist);  
        gpuHamiltonian.nlistsize.copy_sync(cpuHamiltonian.nlistsize); 
        cpuHamiltonian.ncoup.transpose();
        cpuHamiltonian.nlist.transpose();
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
    gpuLattice.mmom0.copy_sync(cpuLattice.mmom0);  
    gpuLattice.mmom2.copy_sync(cpuLattice.mmom2);  
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
    cpuLattice.mmom0.copy_sync(gpuLattice.mmom0);  
    cpuLattice.mmom2.copy_sync(gpuLattice.mmom2);  
    cpuLattice.mmomi.copy_sync(gpuLattice.mmomi);  
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
            GpuSD.SDmphase(*this);
        }
        else {printf("Wrong phase! 0 - initial, 1 - measurement");}
    }
    else if(whichsim == 1){
      /*   CudaMCSimulation CudaMC;//TODO

        if(whichphase == 0) {
            if(bf == 'Y') CudaMC.MCiphase_bf(*this);
            else CudaMC.MCiphase(*this);
            copyToFortran();
        }
        else if(whichphase == 1) {
            if(bf == 'Y') CudaMC.MCmphase_bf(*this);
            else CudaMC.MCmphase(*this);
        }
        else {printf("Wrong phase! 0 - initial, 1 - measurement");}*/
    }
    else {printf("Wrong simulation type! 0 - SD, 1 - MC; current type %i\n", whichsim);}
    //release();
}
