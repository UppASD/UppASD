!-------------------------------------------------------------------------------
!> MODULE: Evolution
!> @brief
!> Driver module for integrating the stochastic LLG-equations
!! @todo See if unit lenght magnetic vectors emom can be removed
!> @author
!> Anders Bergman, Jonathan Chico, Johan Hellsvik et. al
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
!-------------------------------------------------------------------------------
module Evolution
   use Parameters
   use Profiling
   use MidPoint
   use SphMidPoint
   use Heun_single
   use Heun_proper
   use Depondt
   use LLGI
   use RandomNumbers, only : rannum, ranv
   use Constants, only : gama
   use Fielddata, only : thermal_field
   use optimizationroutines, only : modeulermpf_constl,smodeulermpt_constl,invalidationCheck

   implicit none

   private

   public :: evolve_first, evolve_second, evolve_imp

contains

   !-------------------------------------------------------------------------------
   !> SUBROUTINE: evolve_first
   !> @brief
   !> First step of solver, calculates the predictor. Only step for SDEalgh=2,3
   !-------------------------------------------------------------------------------
   subroutine evolve_first(Natom, Mensemble, Landeg, llg, SDEalgh, bn, &
         lambda1_array, lambda2_array, NA, compensate_drift, delta_t, relaxtime, Temp_array, temprescale, &
         beff, b2eff, thermal_field,beff2, btorque, field1, field2, emom, emom2, emomM, mmom, mmomi, stt,&
         do_site_damping,nlist,nlistsize,constellationsUnitVec,constellationsUnitVec2,constellationsMag, &
         constellations,unitCellType,OPT_flag,cos_thr,max_no_constellations,do_she,she_btorque)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      integer, intent(in) :: llg !< Type of equation of motion (1=LLG)
      integer, intent(in) :: SDEalgh !< Solver for equations of motion (1-5)
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0d0)
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: compensate_drift !< Correct for drift in RNG
      real(dblprec), intent(in) :: delta_t !< Time step
      real(dblprec), intent(in) :: relaxtime !< Relaxation time for inertial regime
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(Natom), intent(in) :: lambda2_array !< Additional damping parameter (not used for llg=1)
      real(dblprec), dimension(Natom), intent(in) :: Temp_array !< Temperature
      real(dblprec), intent(inout) :: temprescale !< Temperature rescaling from QHB
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff2 !< External field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: btorque !< Spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: she_btorque !< Spin Hall effect transfer torque
      real(dblprec), dimension(3), intent(in) :: field1 !< Average internal effective field
      real(dblprec), dimension(3), intent(in) :: field2 !< Average external effective field
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmomi !< Inverse of magnitude of magnetic moments
      character(len=1), intent(in) :: STT    !< Treat spin transfer torque?
      character(len=1), intent(in) :: do_she !< Treat spin hall effect torque
      character(len=1), intent(in) :: do_site_damping

      ! Optimization used variables
      integer, dimension(:,:), intent(in)            :: nlist
      integer, dimension(:), intent(in)              :: nlistsize
      real(dblprec), dimension(:,:,:), intent(in)    :: constellationsUnitVec
      real(dblprec), dimension(:,:,:), intent(inout) :: constellationsUnitVec2
      real(dblprec), dimension(:,:), intent(in)      :: constellationsMag
      real(dblprec), dimension(:,:,:), intent(inout) :: constellations
      integer, dimension(:,:),intent(inout)          :: unitCellType
      real(dblprec),intent(in)                       :: cos_thr
      integer, intent(in)                            :: max_no_constellations
      logical, intent(in)                            :: OPT_flag

      real(dblprec) :: lambdatol = 1d-12

      integer :: ij,i,j

      thermal_field=0.0d0
      ! Move optimization cores
      if (OPT_flag) then
         call smodeulermpt_constl(max_no_constellations,Mensemble,Landeg,bn,lambda1_array, &
            beff2(:,1:max_no_constellations,:),constellationsUnitVec,constellationsUnitVec2,&
            constellations,constellationsMag,delta_t)
      end if

      ! Single point Heun
      if(SDEalgh==2) then
         call evolve2(Natom, Mensemble, Landeg, llg, bn, lambda1_array, lambda2_array, &
            beff, beff2, field1, field2, mmomi, emom2, compensate_drift, delta_t, Temp_array,temprescale,thermal_field)

         ! Single point Heun
      elseif(SDEalgh==3) then
         call evolve3(Natom, Mensemble, Landeg, llg, lambda1_array, beff, &
            field1, field2, mmomi, emom2, compensate_drift, delta_t, Temp_array,temprescale, thermal_field)

         ! Predictor-corrector Heun
      elseif(SDEalgh==4) then
         !JohanH-19Feb
         !Heun with effective field updated after
         !predictor step. In evolve4 calls are done
         !to heun4p (predictor step), to heisge (update field)
         !and heun4f (final step)
         call evolve4p(Natom, Mensemble, Landeg, llg, lambda1_array, compensate_drift, &
            beff, b2eff, mmomi, emom2, emom, emomM, mmom,delta_t,Temp_array,temprescale,thermal_field)

         ! Mentink's midpoint solver
      elseif(SDEalgh==1) then
         ! random numbers the same for calculation et and emom2
         ! original subroutine evolve.f90 split
         ! - calculation random numbers
         ! - heun step t gives et (temporarily, actually euler)
         ! - calculation of effective field with et
         ! - heun step f gives emom2 (final)
         ! - normalization in final step
         ! If the damping is negliably small ( < lambdatol ), there is
         ! no need to generate random numbers
         call timing(0,'Evolution     ','OF')
         call timing(0,'RNG           ','ON')

         if (do_site_damping=='Y') then
            call rannum(Natom, Mensemble, NA,  llg, lambda1_array, lambda2_array, &
               compensate_drift, bn, field1, field2, mmomi, Temp_array,temprescale)
         else
            if(minval(lambda1_array) .gt. lambdatol .or. minval(lambda2_array) .gt. lambdatol) then
               call rannum(Natom, Mensemble, NA,  llg, lambda1_array, lambda2_array, &
                  compensate_drift, bn, field1, field2, mmomi, Temp_array,temprescale)
            else
               ranv=0_dblprec
            end if
         endif

         call timing(0,'RNG           ','OF')
         call timing(0,'Evolution     ','ON')

         call smodeulermpt(Natom, Mensemble, Landeg,bn, lambda1_array, beff, &
            emom, emom2, emomM, mmom, delta_t,thermal_field )

         ! Mentink's Semi-implicit midpoint solver with fixed point iteration
      elseif(SDEalgh==6) then
         ! random numbers the same for calculation et and emom2
         ! original subroutine evolve.f90 split
         ! - calculation random numbers
         ! - heun step t gives et (temporarily, actually euler)
         ! - calculation of effective field with et
         ! - heun step f gives emom2 (final)
         ! - normalization in final step
         ! If the damping is negliably small ( < lambdatol ), there is
         ! no need to generate random numbers
         if (do_site_damping=='Y') then
            call rannum(Natom, Mensemble, NA,  llg, lambda1_array, lambda2_array, &
               compensate_drift, bn, field1, field2, mmomi, Temp_array,temprescale)
         else
            if(minval(lambda1_array) .gt. lambdatol .or. minval(lambda2_array) .gt. lambdatol) then
               call rannum(Natom, Mensemble, NA,  llg, lambda1_array, lambda2_array, &
                  compensate_drift, bn, field1, field2, mmomi, Temp_array,temprescale)
            else
               ranv=0_dblprec
            end if
         endif

         write(*,*) 'calls sibt'
         call sibt(Natom, Mensemble, Landeg,bn, lambda1_array, beff, &
            emom, emom2, emomM, mmom, delta_t,thermal_field )

         ! Semi-implicit spherical midpoint solver with fixed point iteration
      elseif(SDEalgh==7) then
         ! random numbers the same for calculation et and emom2
         ! original subroutine evolve.f90 split
         ! - calculation random numbers
         ! - heun step t gives et (temporarily, actually euler)
         ! - calculation of effective field with et
         ! - heun step f gives emom2 (final)
         ! - normalization in final step
         ! If the damping is negliably small ( < lambdatol ), there is
         ! no need to generate random numbers
         if (do_site_damping=='Y') then
            call rannum(Natom, Mensemble, NA,  llg, lambda1_array, lambda2_array, &
               compensate_drift, bn, field1, field2, mmomi, Temp_array,temprescale)
         else
            if(minval(lambda1_array) .gt. lambdatol .or. minval(lambda2_array) .gt. lambdatol) then
               call rannum(Natom, Mensemble, NA,  llg, lambda1_array, lambda2_array, &
                  compensate_drift, bn, field1, field2, mmomi, Temp_array,temprescale)
            else
               ranv=0_dblprec
            end if
         endif

         call sispht(Natom, Mensemble, Landeg,bn, lambda1_array, beff, &
            emom, emom2, emomM, mmom, delta_t,thermal_field )

         ! The Depondt solver
      elseif(SDEalgh==5) then
         call depondt_evolve_first(Natom, Mensemble, lambda1_array, beff, b2eff, &
            btorque, emom, emom2, emomM, mmom, delta_t, &
            Temp_array, temprescale, stt,thermal_field,do_she,she_btorque)

         ! Predictor-corrector solver for the LLGI equation
      elseif(SDEalgh==11) then
         call llgi_evolve_first(Natom, Mensemble, lambda1_array, beff, b2eff, &
            emom, emom2, emomM, mmom, delta_t, relaxtime,Temp_array,temprescale,thermal_field)
      end if

      if(OPT_flag) then
         !$omp parallel do default(shared) private(ij,i,j)
         do ij=1,Natom*Mensemble
            i=mod(ij-1,Natom)+1
            j=int((ij-1)/Natom)+1
            ! Invalidation check: Check whether an atom, after randomization and time stepping, has deviated from the threshold value
            call invalidationCheck(i, j, nlist, nlistsize, constellationsUnitVec, unitCellType, cos_thr, emom2)
         end do
         !$omp end parallel do
      end if

   end subroutine evolve_first

   !-------------------------------------------------------------------------------
   !> SUBROUTINE: evolve_second
   !> @brief
   !> Second step of solver, calculates the corrected solution. Only used for SDEalgh=1,5
   !-------------------------------------------------------------------------------
   subroutine evolve_second(Natom, Mensemble, Landeg, llg, SDEalgh, bn, lambda1_array, &
         delta_t, relaxtime, beff, beff2,b2eff, btorque, emom, emom2,stt, &
         nlist,nlistsize,constellationsUnitVec,constellationsUnitVec2,constellationsMag, &
         constellations,unitCellType,OPT_flag,cos_thr,max_no_constellations,do_she,she_btorque)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      integer, intent(in) :: llg !< Type of equation of motion (1=LLG)
      integer, intent(in) :: SDEalgh !< Solver for equations of motion (1-5)
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0d0)
      real(dblprec), intent(in) :: delta_t !< Time step
      real(dblprec), intent(in) :: relaxtime !< Relaxation time for inertial regime
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff    !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff2 !< External field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: b2eff   !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: btorque !< Spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: she_btorque !< SHE generated spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom2  !< Final (or temporary) unit moment vector
      character(len=1), intent(in) :: STT    !< Treat spin transfer torque?
      character(len=1), intent(in) :: do_she !< Treat the spin hall effect spin transfer torque?

      ! Optimization used variables
      integer, dimension(:,:), intent(in)            :: nlist
      integer, dimension(:), intent(in)              :: nlistsize
      real(dblprec), dimension(:,:,:), intent(inout) :: constellationsUnitVec
      real(dblprec), dimension(:,:,:), intent(inout) :: constellationsUnitVec2
      real(dblprec), dimension(:,:), intent(in)      :: constellationsMag
      real(dblprec), dimension(:,:,:), intent(inout) :: constellations
      integer, dimension(:,:),intent(inout)          :: unitCellType
      real(dblprec),intent(in)                       :: cos_thr
      integer, intent(in)                            :: max_no_constellations
      logical, intent(in)                            :: OPT_flag
      integer                                        :: i,j,k,ij

      if (OPT_flag) then
         call modeulermpf_constl(max_no_constellations,Mensemble,Landeg,bn,lambda1_array,&
            beff2(:,1:max_no_constellations,:), constellationsUnitVec,&
            constellationsUnitVec2,delta_t)

         do k=1,Mensemble
            do i=1,3
               constellations(i,:,k) = constellationsUnitVec(i,:,k)*constellationsMag(:,k)
            end do
         end do
      end if

      ! Predictor-corrector Heun
      if(SDEalgh==4) then
         call evolve4f(Natom, Mensemble, Landeg, lambda1_array, llg, beff,  b2eff, emom2, emom,delta_t)

         ! Mentink's midpoint solver
      elseif(SDEalgh==1) then
         call modeulermpf(Natom, Mensemble, Landeg, bn, lambda1_array, beff, emom, emom2, delta_t)

         ! Mentink's Semi-implicit midpoint solver with fixed point iteration
      elseif(SDEalgh==6) then
         call sibf(Natom, Mensemble, Landeg, bn, lambda1_array, beff, emom, emom2, delta_t)

         ! Semi-implicit spherical midpoint solver
      elseif(SDEalgh==7) then
         call sisphf(Natom, Mensemble, Landeg, bn, lambda1_array, beff, emom, emom2, delta_t)

         ! The Depondt solver
      elseif(SDEalgh==5) then
         call depondt_evolve_second(Natom, Mensemble, lambda1_array, beff, b2eff, &
            btorque, emom, emom2, delta_t, stt,do_she,she_btorque)

         ! Predictor-corrector solver for the LLGI equation
      elseif(SDEalgh==11) then
         call llgi_evolve_second(Natom, Mensemble, lambda1_array, beff, b2eff, &
            emom, emom2, delta_t,relaxtime)

      end if

      if(OPT_flag) then
         !$omp parallel do default(shared) private(ij,i,j)
         do ij=1,Natom*Mensemble
            i=mod(ij-1,Natom)+1
            j=int((ij-1)/Natom)+1
            ! Invalidation check: Check whether an atom, after randomization and time stepping, has deviated from the threshold value
            call invalidationCheck(i, j, nlist, nlistsize, constellationsUnitVec2, unitCellType, cos_thr, emom2)
         end do
         !$omp end parallel do
      end if

   end subroutine evolve_second


   !-------------------------------------------------------------------------------
   !> SUBROUTINE: evolve_imp
   !> @brief
   !> For fully implicit methods SDEalgh=21, 22
   !-------------------------------------------------------------------------------
   subroutine evolve_imp(Natom, Mensemble, Landeg, llg, SDEalgh, bn, &
         lambda1_array, lambda2_array, NA, compensate_drift, delta_t, relaxtime, Temp_array, temprescale, &
         beff, b2eff, thermal_field,beff2, btorque, field1, field2, emom, emom2, emomM, mmom, mmomi, stt,&
         do_site_damping,nlist,nlistsize,constellationsUnitVec,constellationsUnitVec2,constellationsMag, &
         constellations,unitCellType,OPT_flag,cos_thr,max_no_constellations,do_she,she_btorque,converged)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      integer, intent(in) :: llg !< Type of equation of motion (1=LLG)
      integer, intent(in) :: SDEalgh !< Solver for equations of motion (1-5)
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0d0)
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: compensate_drift !< Correct for drift in RNG
      real(dblprec), intent(in) :: delta_t !< Time step
      real(dblprec), intent(in) :: relaxtime !< Relaxation time for inertial regime
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(Natom), intent(in) :: lambda2_array !< Additional damping parameter (not used for llg=1)
      real(dblprec), dimension(Natom), intent(in) :: Temp_array !< Temperature
      real(dblprec), intent(inout) :: temprescale !< Temperature rescaling from QHB
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff2 !< External field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: btorque !< Spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: she_btorque !< Spin Hall effect transfer torque
      real(dblprec), dimension(3), intent(in) :: field1 !< Average internal effective field
      real(dblprec), dimension(3), intent(in) :: field2 !< Average external effective field
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmomi !< Inverse of magnitude of magnetic moments
      character(len=1), intent(in) :: STT    !< Treat spin transfer torque?
      character(len=1), intent(in) :: do_she !< Treat spin hall effect torque
      character(len=1), intent(in) :: do_site_damping
      logical, intent(inout) :: converged !< Fix-point iteration converged or not?

      ! Optimization used variables
      integer, dimension(:,:), intent(in)            :: nlist
      integer, dimension(:), intent(in)              :: nlistsize
      real(dblprec), dimension(:,:,:), intent(in)    :: constellationsUnitVec
      real(dblprec), dimension(:,:,:), intent(inout) :: constellationsUnitVec2
      real(dblprec), dimension(:,:), intent(in)      :: constellationsMag
      real(dblprec), dimension(:,:,:), intent(inout) :: constellations
      integer, dimension(:,:),intent(inout)          :: unitCellType
      real(dblprec),intent(in)                       :: cos_thr
      integer, intent(in)                            :: max_no_constellations
      logical, intent(in)                            :: OPT_flag

      real(dblprec) :: lambdatol = 1e-7

      integer :: ij,i,j

      thermal_field=0.0d0
      ! Move optimization cores
      if (OPT_flag) then
         call smodeulermpt_constl(max_no_constellations,Mensemble,Landeg,bn,lambda1_array, &
            beff2(:,1:max_no_constellations,:),constellationsUnitVec,constellationsUnitVec2,&
            constellations,constellationsMag,delta_t)
      end if

      ! Implicit midpoint solver with fixed point iteration
      ! This call enacts only one iteration.
      ! Calls in order to recalculate the fields are done from 0sd.
      if(SDEalgh==21) then
         ! If the damping is negliably small ( < lambdatol ), there is
         ! no need to generate random numbers
         if (do_site_damping=='Y') then
            call rannum(Natom, Mensemble, NA,  llg, lambda1_array, lambda2_array, &
               compensate_drift, bn, field1, field2, mmomi, Temp_array,temprescale)
         else
            if(minval(lambda1_array) .gt. lambdatol .or. minval(lambda2_array) .gt. lambdatol) then
               call rannum(Natom, Mensemble, NA,  llg, lambda1_array, lambda2_array, &
                  compensate_drift, bn, field1, field2, mmomi, Temp_array, temprescale)
            else
               ranv=0_dblprec
            end if
         endif

         write(*,*) 'calls imp_fp'
         call imp_fp(Natom, Mensemble, Landeg, bn, lambda1_array, beff, &
            emom, emom2, emomM, mmom, delta_t, thermal_field, SDEalgh, converged)

         ! Implicit spherical midpoint solver with fixed point iteration
         ! This call enacts only one iteration.
         ! Calls in order to recalculate the fields are done from 0sd.
      elseif(SDEalgh==22) then
         ! If the damping is negliably small ( < lambdatol ), there is
         ! no need to generate random numbers
         if (do_site_damping=='Y') then
            call rannum(Natom, Mensemble, NA,  llg, lambda1_array, lambda2_array, &
               compensate_drift, bn, field1, field2, mmomi, Temp_array,temprescale)
         else
            if(minval(lambda1_array) .gt. lambdatol .or. minval(lambda2_array) .gt. lambdatol) then
               call rannum(Natom, Mensemble, NA,  llg, lambda1_array, lambda2_array, &
                  compensate_drift, bn, field1, field2, mmomi, Temp_array,temprescale)
            else
               ranv=0_dblprec
            end if
         endif

         ! Temporary fix!
         do i=1,Natom
            do j=1,Mensemble
               b2eff(1:3,i,j) = beff(1:3,i,j)
            end do
         end do

         if(OPT_flag) then
            !$omp parallel do default(shared) private(ij,i,j)
            do ij=1,Natom*Mensemble
               i=mod(ij-1,Natom)+1
               j=int((ij-1)/Natom)+1
               ! Invalidation check: Check whether an atom, after randomization and time stepping, has deviated from the threshold value
               call invalidationCheck(i, j, nlist, nlistsize, constellationsUnitVec, unitCellType, cos_thr, emom2)
            end do
            !$omp end parallel do
         end if

      end if

   end subroutine evolve_imp


end module Evolution
