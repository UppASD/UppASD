#include "momentUpdater.hpp"

#include <cmath>

#include "fortMatrix.hpp"
#include "real_type.h"
#include "stopwatch.hpp"

// Wrapper routine for updating the magnetic moments
void MomentUpdater::update() {
   // Timing
   stopwatch.skip();

   // Calculate
   calcMoment();
   char calc_str[32];
   sprintf(calc_str, "calculate-%d", mompar);
   stopwatch.add(calc_str);

   // Copy
   copyMoment();
   stopwatch.add("copy");

   // Short - calc
   //   alt 0: mmom2 = mmom
   //   alt 1: mmom2 = mmom0 * abs(emom_z)
   //   alt 2: mmom2 = mmom0 * abs(emom_z^2)

   // Short - copy
   //   emom  = emom2
   //   mmom  = mmom2
   //   mmomi = 1/mmom2
   //   emomM = emom2 * mmom2
}

// Transfer moments from emom2 to emom and emomM
void MomentUpdater::copyMoment() {
   // Tolerance
   real momtol = (initexc != 'I') ? -1.0 : 0.000001;

   // Dimensions
   std::size_t M = mmom.dimension_size(1);  // Number of ensembles
   std::size_t N = mmom.dimension_size(0);  // Number of atoms

// Transfer moments
#pragma omp parallel for collapse(2)
   for(std::size_t j = 0; j < M; j++) {
      for(std::size_t i = 0; i < N; i++) {
         real m2 = mmom2(i, j);
         emom(0, i, j) = emom2(0, i, j);
         emom(1, i, j) = emom2(1, i, j);
         emom(2, i, j) = emom2(2, i, j);
         emomM(0, i, j) = emom2(0, i, j) * m2;
         emomM(1, i, j) = emom2(1, i, j) * m2;
         emomM(2, i, j) = emom2(2, i, j) * m2;
         mmom(i, j) = m2;
         mmomi(i, j) = (m2 < momtol) ? (1.0) : (1.0 / m2);
      }
   }
}

// Update the magnitude of the magnetic moment if wanted.
void MomentUpdater::calcMoment() {
   // Dimensions
   std::size_t M = mmom.dimension_size(1);  // Number of ensembles
   std::size_t N = mmom.dimension_size(0);  // Number of atoms

   switch(mompar) {
      case 1:  // M = M0 * mz
#pragma omp parallel for collapse(2)
         for(std::size_t j = 0; j < M; j++) {
            for(std::size_t i = 0; i < N; i++) {
               real m = mmom0(i, j) * std::abs(emom2(2, i, j));
               mmom2(i, j) = (m > 1e-4) ? m : 1e-4;
            }
         }
         break;
      case 2:  // M = M0 * mz^2
#pragma omp parallel for collapse(2)
         for(std::size_t j = 0; j < M; j++) {
            for(std::size_t i = 0; i < N; i++) {
               real mz = emom2(2, i, j);
               real m = mmom0(i, j) * mz * mz;
               mmom2(i, j) = (m > 1e-4) ? m : 1e-4;
            }
         }
         break;
      case 3:  // M=parametrized from Nature Nano. (Pt wire)
         printf("mompar == 3 Not implemented!");
         /*
                         for (int j = 0; j < M; j++) {
                                 for (int i = 0; i < N; i++) {
                                         mmom2(i,j) = max(ptnanowire(abs(emom2(2,i,j)),mmom0(i,j)), 1e-4);
                                 }
                         }
         */
         break;
      default:  // Normal case:
         mmom2.memcopy(mmom);
         break;
   }
}

