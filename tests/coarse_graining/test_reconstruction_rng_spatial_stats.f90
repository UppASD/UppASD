program test_reconstruction_rng_spatial_stats

   ! RCG-06D (F-20) RECONSTRUCTION-RNG-SPATIAL-STATS raw-draw generator.
   !
   ! Prints CSV rows of (global_seed,channel,ensemble,epoch,block,u1,u2) for
   ! a sweep of block indices at several fixed (global_seed,channel,ensemble,
   ! epoch) tuples, using the literal production seed formula
   ! (deterministic_reconstruction_seed) and RNG step (next_uniform) from
   ! AdaptiveHybridSolver -- not a reimplementation that could silently
   ! drift, matching the gpuAtomicDouble.hpp/ENERGY-FP32-ACCUM precedent.
   ! run_reconstruction_rng_spatial_stats.py computes the actual spatial-
   ! correlation statistics from this raw output; this program only emits
   ! numbers, deliberately NOT linked against asdlib beyond this one module
   ! and its own dependencies, so it stays valid evidence for the RNG
   ! defect class independent of later production source changes elsewhere.

   use, intrinsic :: iso_fortran_env, only : int64, output_unit
   use Parameters, only : dblprec
   use AdaptiveHybridSolver, only : deterministic_reconstruction_seed, next_uniform

   implicit none

   integer, parameter :: nblocks = 4096
   integer(int64), parameter :: global_seeds(3) = &
      (/ 1_int64, 99173_int64, 20260812_int64 /)
   integer, parameter :: channels(2) = (/ 1, 2 /)
   integer, parameter :: ensembles(2) = (/ 1, 2 /)
   integer, parameter :: epochs(2) = (/ 0, 7 /)

   integer :: si, ci, ei, pi, block
   integer(int64) :: seed, rng
   real(dblprec) :: u1, u2

   write(output_unit,'(a)') '# RCG-06D reconstruction RNG spatial statistics raw draws'
   write(output_unit,'(a)') '# global_seed channel ensemble epoch block u1 u2'
   do si = 1, size(global_seeds)
      do ci = 1, size(channels)
         do ei = 1, size(ensembles)
            do pi = 1, size(epochs)
               do block = 1, nblocks
                  seed = deterministic_reconstruction_seed(global_seeds(si),block, &
                     channels(ci),ensembles(ei),epochs(pi))
                  rng = seed
                  u1 = next_uniform(rng)
                  u2 = next_uniform(rng)
                  write(output_unit,'(i0,1x,i0,1x,i0,1x,i0,1x,i0,1x,es24.16,1x,es24.16)') &
                     global_seeds(si),channels(ci),ensembles(ei),epochs(pi),block,u1,u2
               end do
            end do
         end do
      end do
   end do

end program test_reconstruction_rng_spatial_stats
