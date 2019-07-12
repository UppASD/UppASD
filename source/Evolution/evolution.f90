!-------------------------------------------------------------------------------
!  MODULE: Evolution
!> @brief
!> Driver module for integrating the stochastic LLG-equations
!> @todo See if unit lenght magnetic vectors emom can be removed
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module Evolution

   use LLGI
   use Depondt
   use MidPoint
   use Profiling
   use Constants,             only : gama
   use Fielddata,             only : thermal_field
   use InputData,             only : perp
   use Parameters
   use SphMidPoint
   use Heun_single
   use Heun_proper
   use RandomNumbers,         only : rannum, ranv
   use optimizationroutines,  only : modeulermpf_constl,smodeulermpt_constl,invalidationCheck
   use InputData, only : perp

   implicit none

   private

   public :: evolve_first, evolve_second

contains

   !----------------------------------------------------------------------------
   !  SUBROUTINE: evolve_first
   !> @brief
   !> First step of solver, calculates the predictor. Only step for SDEalgh=2,3
   !----------------------------------------------------------------------------
   subroutine evolve_first(Natom,Mensemble,Landeg,llg,SDEalgh,bn,lambda1_array,     &
      lambda2_array,NA,compensate_drift,delta_t,relaxtime,Temp_array,temprescale,   &
      beff,b2eff,thermal_field,beff2,btorque,field1,field2,emom,emom2,emomM,mmom,   &
      mmomi,stt,do_site_damping,nlist,nlistsize,constellationsUnitVec,              &
      constellationsUnitVec2,constellationsMag,constellations,unitCellType,OPT_flag,&
      cos_thr,max_no_constellations,do_she,she_btorque,Nred,red_atom_list,          &
      do_fixed_mom,do_sot,sot_btorque)

      implicit none

      integer, intent(in) :: NA        !< Number of atoms in one cell
      integer, intent(in) :: llg       !< Type of equation of motion (1=LLG)
      integer, intent(in) :: Nred      !< Number of moments that can be updated
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: SDEalgh   !< Solver for equations of motion (1-5)
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: compensate_drift !< Correct for drift in RNG
      real(dblprec), intent(in) :: bn        !< Scaling factor for LLG equation (parameter=1.0_dblprec)
      real(dblprec), intent(in) :: delta_t   !< Time step
      real(dblprec), intent(in) :: relaxtime !< Relaxation time for inertial regime
      character(len=1), intent(in) :: STT    !< Treat spin transfer torque?
      character(len=1), intent(in) :: do_she !< Treat spin hall effect torque
      character(len=1), intent(in) :: do_sot !< Treat the general SOT model
      character(len=1), intent(in) :: do_fixed_mom 	!< Do Fixed moment calculation (Y/N)
      character(len=1), intent(in) :: do_site_damping	!< Flag for site dependent damping in measurement phase
      integer, dimension(Nred), intent(in) :: red_atom_list !< Reduced list containing atoms allowed to evolve in a fixed moment calculation
      real(dblprec), dimension(3), intent(in) :: field1 !< Average internal effective field
      real(dblprec), dimension(3), intent(in) :: field2 !< Average external effective field
      real(dblprec), dimension(Natom), intent(in) :: Landeg          !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array   !< Damping parameter
      real(dblprec), dimension(Natom), intent(in) :: lambda2_array   !< Additional damping parameter (not used for llg=1)
      real(dblprec), dimension(Natom), intent(in) :: Temp_array      !< Temperature
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom  !< Magnitude of magnetic moments
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmomi !< Inverse of magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff2        !< External field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: btorque      !< Spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: she_btorque  !< Spin Hall effect transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: sot_btorque  !< Spin orbit torque
      ! .. In/out variables
      real(dblprec), intent(inout) :: temprescale !< Temperature rescaling from QHB
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: beff   !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field !< Thermal field
      ! .. Output variables
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM  !< Current magnetic moment vector

      ! Optimization used variables
      integer, intent(in)                            :: max_no_constellations  !< Max number of constellations for all ensembles
      real(dblprec), intent(in)                      :: cos_thr   !< Cosine similarity threshold
      logical, intent(in)                            :: OPT_flag  !< Optimization flag
      integer, dimension(:), intent(in)              :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(:,:), intent(in)            :: nlist  !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(:,:), intent(in)      :: constellationsMag  !< Magnitude of magnetic moments
      real(dblprec), dimension(:,:,:), intent(in)    :: constellationsUnitVec !< Normalized constellation matrix It is used for efficient cosinus comparisons in evolve step
      integer, dimension(:,:), intent(inout)         :: unitCellType !< Array of constellation id and classification (core, boundary, or noise) per atom
      real(dblprec), dimension(:,:,:), intent(inout) :: constellationsUnitVec2
      real(dblprec), dimension(:,:,:), intent(inout) :: constellations  !< Saved fixed unit cell configurations, these represent a configuration present in the domain with same configuration of unit cells in the neighborhood

      ! .. Local variables
      integer :: ij,i,j,k
      real(dblprec) :: lambdatol = 1d-12

      thermal_field=0.0_dblprec

      if(perp=='Y') then
         !$omp parallel do default(shared) private(k,i) collapse(2)
         do k=1,Mensemble
            do i=1,Natom
               beff(:,i,k)=beff(:,i,k)-emom(:,i,k)*sum(beff(:,i,k)*emom(:,i,k))
            end do
         end do
         !$omp end parallel do
      end if
      !------------------------------------------------------------------------------
      ! Move optimization cores
      !------------------------------------------------------------------------------
      if (OPT_flag) then
         call smodeulermpt_constl(max_no_constellations,Mensemble,Landeg,bn,        &
            lambda1_array,beff2(:,1:max_no_constellations,:),constellationsUnitVec, &
            constellationsUnitVec2,constellations,constellationsMag,delta_t)
      end if
      !------------------------------------------------------------------------------
      ! Mentink's midpoint solver
      !------------------------------------------------------------------------------
      if(SDEalgh==1) then
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
         ! If one is using site dependent damping generate the random numbers
         if (do_site_damping=='Y') then
            call rannum(Natom,Mensemble,NA,llg,lambda1_array,lambda2_array,         &
               compensate_drift,bn,field1,field2,mmomi,Temp_array,temprescale)
            else
               if(minval(lambda1_array).gt.lambdatol.or.minval(lambda2_array).gt.lambdatol) then
                  call rannum(Natom,Mensemble,NA,llg,lambda1_array,lambda2_array,   &
                  compensate_drift,bn,field1,field2,mmomi,Temp_array,temprescale)
               else
                  ranv=0_dblprec
               end if
            endif

         call timing(0,'RNG           ','OF')
         call timing(0,'Evolution     ','ON')

         call smodeulermpt(Natom,Mensemble,Landeg,bn,lambda1_array,beff,emom,emom2, &
            emomM,mmom,delta_t,thermal_field,STT,do_she,do_sot,btorque,she_btorque, &
            sot_btorque,Nred,red_atom_list)
      !------------------------------------------------------------------------------
      ! Single point Heun
      !------------------------------------------------------------------------------
      elseif(SDEalgh==2) then
         call evolve2(Natom,Mensemble,Landeg,llg,bn,lambda1_array,lambda2_array,    &
            beff,beff2,field1,field2,mmomi,emom2,compensate_drift,delta_t,          &
            Temp_array,temprescale,thermal_field,Nred,red_atom_list)
      !------------------------------------------------------------------------------
      ! Single point Heun
      !------------------------------------------------------------------------------
      elseif(SDEalgh==3) then
         call evolve3(Natom,Mensemble,Landeg,llg,lambda1_array,beff,field1,field2,  &
            mmomi,emom2,compensate_drift,delta_t,Temp_array,temprescale,            &
            thermal_field,Nred,red_atom_list)
      !------------------------------------------------------------------------------
      ! Predictor-corrector Heun
      !------------------------------------------------------------------------------
      elseif(SDEalgh==4) then
         !JohanH-19Feb
         !Heun with effective field updated after
         !predictor step. In evolve4 calls are done
         !to heun4p (predictor step), to heisge (update field)
         !and heun4f (final step)
         call evolve4p(Natom,Mensemble,Landeg,llg,lambda1_array,compensate_drift,   &
            beff,b2eff,mmomi,emom2,emom,emomM,mmom,delta_t,Temp_array,temprescale,  &
            thermal_field,STT,do_she,do_sot,btorque,she_btorque,sot_btorque,Nred,   &
            red_atom_list)
      !------------------------------------------------------------------------------
      ! The Depondt solver
      !------------------------------------------------------------------------------
      elseif(SDEalgh==5) then
         ! Selective updater only some moments are updated
         call depondt_evolve_first(Natom,Nred,Mensemble,lambda1_array,beff,b2eff,   &
            btorque,emom,emom2,emomM,mmom,delta_t,Temp_array,temprescale,stt,       &
            thermal_field,do_she,she_btorque,do_sot,sot_btorque,red_atom_list)
      !------------------------------------------------------------------------------
      ! Mentink's Semi-implicit midpoint solver with fixed point iteration
      !------------------------------------------------------------------------------
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
         call sibt(Natom,Mensemble,Landeg,bn,lambda1_array,beff,emom,emom2,emomM,   &
            mmom,delta_t,thermal_field)
      !------------------------------------------------------------------------------
      ! Semi-implicit spherical midpoint solver with fixed point iteration
      !------------------------------------------------------------------------------
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
               ranv=0.0_dblprec
            end if
         endif

         call sispht(Natom,Mensemble,Landeg,bn,lambda1_array,beff,emom,emom2,emomM, &
            mmom,delta_t,thermal_field )
      !------------------------------------------------------------------------------
      ! Predictor-corrector solver for the LLGI equation
      !------------------------------------------------------------------------------
      elseif(SDEalgh==11) then
         call llgi_evolve_first(Natom,Mensemble,lambda1_array,beff,b2eff,emom,emom2,&
            emomM,mmom,delta_t,relaxtime,Temp_array,temprescale,thermal_field,Nred, &
            red_atom_list)
      end if
      !------------------------------------------------------------------------------
      ! Optimization region
      !------------------------------------------------------------------------------
      if(OPT_flag) then
         !$omp parallel do default(shared) private(ij,i,j)
         do ij=1,Natom*Mensemble
            i=mod(ij-1,Natom)+1
            j=int((ij-1)/Natom)+1
            ! Invalidation check: Check whether an atom, after randomization and time stepping, has deviated from the threshold value
            call invalidationCheck(i,j,nlist,nlistsize,constellationsUnitVec,       &
               unitCellType, cos_thr, emom2)
         end do
         !$omp end parallel do
      end if

   end subroutine evolve_first

   !----------------------------------------------------------------------------
   !  SUBROUTINE: evolve_second
   !> @brief
   !> Second step of solver, calculates the corrected solution. Only used for SDEalgh=1,5
   !----------------------------------------------------------------------------
   subroutine evolve_second(Natom,Mensemble,Landeg,llg,SDEalgh,bn,lambda1_array,    &
      delta_t,relaxtime,beff,beff2,b2eff,btorque,emom,emom2,stt,nlist,nlistsize,    &
      constellationsUnitVec,constellationsUnitVec2,constellationsMag,constellations,&
      unitCellType,OPT_flag,cos_thr,max_no_constellations,do_she,she_btorque,Nred,  &
      red_atom_list,do_fixed_mom,do_sot,sot_btorque)

      implicit none

      integer, intent(in) :: llg       !< Type of equation of motion (1=LLG)
      integer, intent(in) :: Nred      !< Number of moments that can be updated
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: SDEalgh   !< Solver for equations of motion (1-5)
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in) :: bn        !< Scaling factor for LLG equation (parameter=1.0_dblprec)
      real(dblprec), intent(in) :: delta_t   !< Time step
      real(dblprec), intent(in) :: relaxtime !< Relaxation time for inertial regime
      character(len=1), intent(in) :: STT    !< Treat spin transfer torque?
      character(len=1), intent(in) :: do_she !< Treat the spin hall effect spin transfer torque?
      character(len=1), intent(in) :: do_sot !< Treat the general SOT model
      character(len=1), intent(in) :: do_fixed_mom !< Do Fixed moment calculation (Y/N)
      integer, dimension(Nred), intent(in) :: red_atom_list !< Reduced list containing atoms allowed to evolve in a fixed moment calculation
      real(dblprec), dimension(Natom), intent(in) :: Landeg !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array            !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: beff		!< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff2        !< External field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: b2eff        !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: btorque      !< Spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: she_btorque  !< SHE generated spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: sot_btorque  !< Spin orbit torque
      ! .. Output variables
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom        !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom2       !< Final (or temporary) unit moment vector

      ! Optimization used variables
      integer, intent(inout)                         :: max_no_constellations  !< Max number of constellations for all ensembles
      real(dblprec), intent(in)                      :: cos_thr   !< Cosine similarity threshold
      logical, intent(in)                            :: OPT_flag  !< Optimization flag
      integer, dimension(:), intent(in)              :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(:,:), intent(in)            :: nlist  !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(:,:), intent(in)      :: constellationsMag  !< Magnitude of magnetic moments
      integer, dimension(:,:), intent(inout)         :: unitCellType !< Array of constellation id and classification (core, boundary, or noise) per atom
      real(dblprec), dimension(:,:,:), intent(inout) :: constellationsUnitVec !< Normalized constellation matrix It is used for efficient cosinus comparisons in evolve step
      real(dblprec), dimension(:,:,:), intent(inout) :: constellationsUnitVec2
      real(dblprec), dimension(:,:,:), intent(inout) :: constellations  !< Saved fixed unit cell configurations, these represent a configuration present in the domain with same configuration of unit cells in the neighborhood
      ! .. Local variables
      integer                                        :: i,j,k,ij

      if(perp=='Y') then
         !$omp parallel do default(shared) private(k,i) collapse(2)
         do k=1,Mensemble
            do i=1,Natom
               beff(:,i,k)=beff(:,i,k)-emom(:,i,k)*sum(beff(:,i,k)*emom(:,i,k))
            end do
         end do
         !$omp end parallel do
      end if

      !------------------------------------------------------------------------------
      ! Optimization region
      !------------------------------------------------------------------------------
      if (OPT_flag) then
         call modeulermpf_constl(max_no_constellations,Mensemble,Landeg,bn,         &
            lambda1_array,beff2(:,1:max_no_constellations,:),constellationsUnitVec, &
            constellationsUnitVec2,delta_t)

         do k=1,Mensemble
            do i=1,3
               constellations(i,:,k) = constellationsUnitVec(i,:,k)*constellationsMag(:,k)
            end do
         end do
      end if
      !------------------------------------------------------------------------------
      ! Mentink's midpoint solver
      !------------------------------------------------------------------------------
      if(SDEalgh==1) then
         call modeulermpf(Natom,Mensemble,Landeg,bn,lambda1_array,beff, emom,emom2, &
            delta_t,STT,do_she,do_sot,btorque,she_btorque,sot_btorque,Nred,red_atom_list)
      !------------------------------------------------------------------------------
      ! Predictor-corrector Heun
      !------------------------------------------------------------------------------
      elseif(SDEalgh==4) then
         call evolve4f(Natom,Mensemble,Landeg,lambda1_array,llg,beff, b2eff,emom2,  &
            emom,delta_t,STT,do_she,do_sot,btorque,she_btorque,sot_btorque,Nred,    &
            red_atom_list)
      !------------------------------------------------------------------------------
      ! The Depondt solver
      !------------------------------------------------------------------------------
      elseif(SDEalgh==5) then
         ! Selective moment update only some moments are evolved
         call depondt_evolve_second(Natom,Nred,Mensemble,lambda1_array,beff,b2eff,  &
            btorque,emom,emom2,delta_t,stt,do_she,she_btorque,do_sot,sot_btorque,   &
            red_atom_list)
      !------------------------------------------------------------------------------
      ! Mentink's Semi-implicit midpoint solver with fixed point iteration
      !------------------------------------------------------------------------------
      elseif(SDEalgh==6) then
         call sibf(Natom,Mensemble,Landeg,bn,lambda1_array,beff,emom,emom2,delta_t)
      !------------------------------------------------------------------------------
      ! Semi-implicit spherical midpoint solver
      !------------------------------------------------------------------------------
      elseif(SDEalgh==7) then
         call sisphf(Natom,Mensemble,Landeg,bn,lambda1_array,beff,emom,emom2,delta_t)
      !------------------------------------------------------------------------------
      ! Predictor-corrector solver for the LLGI equation
      !------------------------------------------------------------------------------
      elseif(SDEalgh==11) then
         call llgi_evolve_second(Natom,Mensemble,lambda1_array,beff,b2eff,emom,     &
            emom2,delta_t,relaxtime,Nred,red_atom_list)
      end if
      !------------------------------------------------------------------------------
      ! Optimization region
      !------------------------------------------------------------------------------
      if(OPT_flag) then
         !$omp parallel do default(shared) private(ij,i,j)
         do ij=1,Natom*Mensemble
            i=mod(ij-1,Natom)+1
            j=int((ij-1)/Natom)+1
            ! Invalidation check: Check whether an atom, after randomization and
            ! time stepping, has deviated from the threshold value
            call invalidationCheck(i,j,nlist,nlistsize,constellationsUnitVec2,      &
               unitCellType,cos_thr,emom2)
         end do
         !$omp end parallel do
      end if

   end subroutine evolve_second

end module Evolution
