#include "mdSimulation.hpp"

#include "c_headers.hpp"
#include "c_helper.h"
#include "depondtIntegrator.hpp"
#include "fortranData.hpp"
#include "hamiltonianCalculations.hpp"
#include "momentUpdater.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchPool.hpp"

MdSimulation::~MdSimulation() {
   freeOwn();
}

void MdSimulation::initiateConstants() {
   // Only heisge_jij allowed
   SDEalgh = *FortranData::SDEalgh;
   if(!(SDEalgh == 1 || SDEalgh == 4 || SDEalgh == 5 || SDEalgh == 11)) {
      std::fprintf(stderr, "Invalid SDEalgh!\n");
      std::exit(EXIT_FAILURE);
   }

   // Constants
   stt = *FortranData::stt;
   rstep = *FortranData::rstep;
   nstep = *FortranData::nstep;
   Natom = *FortranData::Natom;
   Mensemble = *FortranData::Mensemble;
   max_no_neigh = *FortranData::max_no_neigh;
   delta_t = *FortranData::delta_t;
   gamma = *FortranData::gamma;
   k_bolt = *FortranData::k_bolt;
   mub = *FortranData::mub;
   damping = *FortranData::damping;
   binderc = FortranData::binderc;
   mavg = FortranData::mavg;
   mompar = *FortranData::mompar;
   initexc = *FortranData::initexc;
   do_dm = *FortranData::do_dm;
   max_no_dmneigh = *FortranData::max_no_dmneigh;
}

void MdSimulation::initiateFortran() {
   // Dimensions
   usd_int N = Natom;
   usd_int M = Mensemble;

   // Constants initiated?
   if(N == 0 || M == 0) {
      std::printf("MdSimulation: constants not initiated!\n");
      std::exit(EXIT_FAILURE);
   }

   // Inititate
   ncoup.set(FortranData::ncoup, max_no_neigh, N);
   nlist.set(FortranData::nlist, max_no_neigh, N);
   nlistsize.set(FortranData::nlistsize, N);
   dmvect.set(FortranData::dmvect, 3, max_no_dmneigh, N);
   dmlist.set(FortranData::dmlist, max_no_dmneigh, N);
   dmlistsize.set(FortranData::dmlistsize, N);
   beff.set(FortranData::beff, 3, N, M);
   b2eff.set(FortranData::b2eff, 3, N, M);
   emomM.set(FortranData::emomM, 3, N, M);
   emom.set(FortranData::emom, 3, N, M);
   emom2.set(FortranData::emom2, 3, N, M);
   external_field.set(FortranData::external_field, 3, N, M);
   mmom.set(FortranData::mmom, N, M);
   btorque.set(FortranData::btorque, 3, N, M);
   temperature.set(FortranData::temperature, N);
   mmom0.set(FortranData::mmom0, N, M);
   mmom2.set(FortranData::mmom2, N, M);
   mmomi.set(FortranData::mmomi, N, M);

   isOwnData = false;
}

void MdSimulation::printConstants() {
   std::printf(
       "stt             : %c\n"
       "SDEalgh         : %d\n"
       "rstep           : %ld\n"
       "nstep           : %ld\n"
       "Natom           : %ld\n"
       "Mensemble       : %ld\n"
       "max_no_neigh    : %ld\n"
       "delta_t         : %g\n"
       "gamma           : %g\n"
       "k_bolt          : %g\n"
       "mub             : %g\n"
       "damping         : %g\n"
       "binderc         : %g\n"
       "mavg            : %g\n"
       "do_dm           : %d\n"
       "max_no_dmneigh  : %ld\n",
       stt,
       SDEalgh,
       rstep,
       nstep,
       Natom,
       Mensemble,
       max_no_neigh,
       delta_t,
       gamma,
       k_bolt,
       mub,
       damping,
       *binderc,
       *mavg,
       do_dm,
       max_no_dmneigh);
}

// Spin Dynamics measurement phase
void MdSimulation::measurementPhase() {
   // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);

   std::printf("C/C++: md simulations starting\n");

   // Timer
   // Stopwatch stopwatch;
   Stopwatch& stopwatch = GlobalStopwatchPool::get("Measurement phase");

   // Depontd integrator
   DepondtIntegrator integrator;

   // Hamiltonian calculations
   HamiltonianCalculations ham(ncoup, nlist, nlistsize, do_dm, dmvect, dmlist, dmlistsize);

   // Moment updater
   MomentUpdater momUpdater(mmom, mmom0, mmom2, emom, emom2, emomM, mmomi, mompar, initexc);

   // Attempt initiation (C/C++ memory allocation)
   if(!integrator.initiate(Natom, Mensemble, stt)) {
      return;
   }

   // Initiate constants for integrator
   integrator.initiateConstants(gamma, k_bolt, mub, damping, temperature, delta_t);

   // Debug
   // fortran_init_c_md_const(); printConstants();

   stopwatch.add("initiate");
   // Time step loop
   for(usd_int mstep = rstep; mstep < rstep + nstep; mstep++) {
      // export_mstep(mstep);

      // Measure averages and trajectories (through fortran call)
      if(isOwnData) {
         copyToFortran();
         stopwatch.add("F/C copy");
      }
      fortran_measure(mstep);  // through c_helper.h to chelper.f90

      // Print simulation status for each 5% of the simulation length given that
      // length is larger than 20, otherwise print each step
      if(nstep > 20) {
         if((mstep + 1) % ((rstep + nstep) / 20) == 0) {
            // Debug
            // fortran_init_c_md_const(); printConstants();

            // Update mavrg and binderc if necessary
            fortran_calc_simulation_status_variables(mavg);
            // Print status
            std::printf("C/C++: %2ld%% done. Mbar: %10.6f. U: %8.5f.\n",
                        (mstep + 1) * 100 / (rstep + nstep),
                        *mavg,
                        *binderc);
         }
      } else {
         // Update mavrg and binderc if necessary
         fortran_calc_simulation_status_variables(mavg);
         std::printf("C/C++: Iteration %ld Mbar %13.6f\n", mstep, *mavg);
      }
      stopwatch.add("measurement");

      // Apply Hamiltonian to obtain effective field
      ham.heisge_jij(beff, emomM, emom, external_field);

      stopwatch.add("hamiltonian");

      // Perform first step of SDE solver
      integrator.evolveFirst(beff, b2eff, btorque, emom, emom2, emomM, mmom);

      stopwatch.add("evolution");

      // Apply Hamiltonian to obtain effective field
      ham.heisge_jij(beff, emomM, emom, external_field);

      stopwatch.add("hamiltonian");

      // Perform second (corrector) step of SDE solver
      integrator.evolveSecond(beff, b2eff, btorque, emom, emom2);

      stopwatch.add("evolution");

      // Update magnetic moments after time evolution step
      //		fortran_moment_update();
      momUpdater.update();

      stopwatch.add("moments");
   }  // End loop over simulation steps

   if(isOwnData) {
      copyToFortran();
      stopwatch.add("F/C copy");
   }

   // Measure averages and trajectories
   fortran_measure(rstep + nstep);

   // Print remaining measurements
   fortran_flush_measurements(rstep + nstep);
   stopwatch.add("measurement");

   //	std::printf("C/C++: md simulations done!\n");
}

// Safe copy (allows nullptr pointer)
static inline void* scopy(void* p1, void* p2, usd_int s) {
   //	std::printf("memcpy(%10p, %10p, %ld);\n", p1, p2, s);
   return (p1 && p2) ? memcpy(p1, p2, s) : p1;
}

void MdSimulation::copyFromFortran() {
   usd_int N = Natom;
   usd_int M = Mensemble;

   if(!isOwnData) {
      return;
   }

   // Copy data
   scopy(ncoup.get_data(), FortranData::ncoup, ncoup.data_size());
   scopy(nlist.get_data(), FortranData::nlist, nlist.data_size());
   scopy(nlistsize.get_data(), FortranData::nlistsize, nlistsize.data_size());
   scopy(beff.get_data(), FortranData::beff, beff.data_size());
   scopy(b2eff.get_data(), FortranData::b2eff, b2eff.data_size());
   scopy(emomM.get_data(), FortranData::emomM, emomM.data_size());
   scopy(emom.get_data(), FortranData::emom, emom.data_size());
   scopy(emom2.get_data(), FortranData::emom2, emom2.data_size());
   scopy(external_field.get_data(), FortranData::external_field, external_field.data_size());
   scopy(mmom.get_data(), FortranData::mmom, mmom.data_size());
   scopy(btorque.get_data(), FortranData::btorque, btorque.data_size());
   scopy(temperature.get_data(), FortranData::temperature, temperature.data_size());
   scopy(mmom0.get_data(), FortranData::mmom0, mmom0.data_size());
   scopy(mmom2.get_data(), FortranData::mmom2, mmom2.data_size());
   scopy(mmomi.get_data(), FortranData::mmomi, mmomi.data_size());
}

void MdSimulation::copyToFortran() {
   usd_int N = Natom;
   usd_int M = Mensemble;

   if(!isOwnData) {
      return;
   }

   // Copy data
   scopy(FortranData::ncoup, ncoup.get_data(), ncoup.data_size());
   scopy(FortranData::nlist, nlist.get_data(), nlist.data_size());
   scopy(FortranData::nlistsize, nlistsize.get_data(), nlistsize.data_size());
   scopy(FortranData::beff, beff.get_data(), beff.data_size());
   scopy(FortranData::b2eff, b2eff.get_data(), b2eff.data_size());
   scopy(FortranData::emomM, emomM.get_data(), emomM.data_size());
   scopy(FortranData::emom, emom.get_data(), emom.data_size());
   scopy(FortranData::emom2, emom2.get_data(), emom2.data_size());
   scopy(FortranData::external_field, external_field.get_data(), external_field.data_size());
   scopy(FortranData::mmom, mmom.get_data(), mmom.data_size());
   scopy(FortranData::btorque, btorque.get_data(), btorque.data_size());
   scopy(FortranData::temperature, temperature.get_data(), temperature.data_size());
   scopy(FortranData::mmom0, mmom0.get_data(), mmom0.data_size());
   scopy(FortranData::mmom2, mmom2.get_data(), mmom2.data_size());
   scopy(FortranData::mmomi, mmomi.get_data(), mmomi.data_size());
}

/////////////////////////////////////////////////////////
// This does not work right now and is not used either //
// Kept by Niklas for later use?                       //
/////////////////////////////////////////////////////////
// TODO: figure out what to do with this
void MdSimulation::initiateOwn() {
   // Dimensions
   usd_int N = Natom;
   usd_int M = Mensemble;

   // Constants initiated?
   if(N == 0 || M == 0) {
      std::printf("MdSimulation: constants not initiated!\n");
      std::exit(EXIT_FAILURE);
   }

   // Inititate
   ncoup.set(new real[max_no_neigh * N], max_no_neigh, N);
   nlist.set(new usd_int[max_no_neigh * N], max_no_neigh, N);
   nlistsize.set(new usd_int[N], N);
   beff.set(new real[3 * N * M], 3, N, M);
   b2eff.set(new real[3 * N * M], 3, N, M);
   emomM.set(new real[3 * N * M], 3, N, M);
   emom.set(new real[3 * N * M], 3, N, M);
   emom2.set(new real[3 * N * M], 3, N, M);
   external_field.set(new real[3 * N * M], 3, N, M);
   mmom.set(new real[N * M], N, M);
   btorque.set(new real[3 * N * M], 3, N, M);
   temperature.set(new real[N], N);
   mmom0.set(new real[N * M], N, M);
   mmom2.set(new real[N * M], N, M);
   mmomi.set(new real[N * M], N, M);

   isOwnData = true;

   // Initiate data
   copyFromFortran();
}

void MdSimulation::freeOwn() {
   if(isOwnData) {
      // Delete
      delete ncoup.get_data();
      delete nlist.get_data();
      delete nlistsize.get_data();
      delete beff.get_data();
      delete b2eff.get_data();
      delete emomM.get_data();
      delete emom.get_data();
      delete emom2.get_data();
      delete external_field.get_data();
      delete mmom.get_data();
      delete btorque.get_data();
      delete temperature.get_data();
      delete mmom0.get_data();
      delete mmom2.get_data();
      delete mmomi.get_data();

      // Reset
      ncoup.set(nullptr, 0, 0, 0);
      nlist.set(nullptr, 0, 0, 0);
      nlistsize.set(nullptr, 0, 0, 0);
      beff.set(nullptr, 0, 0, 0);
      b2eff.set(nullptr, 0, 0, 0);
      emomM.set(nullptr, 0, 0, 0);
      emom.set(nullptr, 0, 0, 0);
      emom2.set(nullptr, 0, 0, 0);
      external_field.set(nullptr, 0, 0, 0);
      btorque.set(nullptr, 0, 0, 0);
      mmom0.set(nullptr, 0, 0);
      mmom2.set(nullptr, 0, 0);
      mmomi.set(nullptr, 0, 0);
      mmom.set(nullptr, 0, 0);
      temperature.set(nullptr, 0);
   }
}

