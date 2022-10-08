!> Calculate effective interatomic force field by applying the derivative of the Hamiltonian
!> \details The effective field, \f$\mathbf{F}_i\f$, on an atom \f$\textit{i}\f$, is calculated from
!> \f$ \mathbf{F}_i=-\frac{\partial \mathbf{H}}{\partial \mathbf{u}_i}
module LatticeHamiltonianActions
   use Parameters
   use Profiling
   use HamiltonianData, only :  ham
   use LatticeHamiltonianData
   use inputdata, only : do_hoc_debug

contains


   !> Calculate effective field by applying the derivative of the Hamiltonian
   subroutine effective_latticefield(Natom, Mensemble, start_atom, stop_atom, &
      do_ll, do_ml, do_mml, mode, &
      uvec, emomM, latt_external_field, latt_time_external_field, &
      eeff, eeff1, eeff2, eeff3)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: start_atom !< Atom to start loop for
      integer, intent(in) :: stop_atom !< Atom to end loop for

      integer, intent(in) :: do_ll   !< Add LL term to lattice Hamiltonian (0/1)
      integer, intent(in) :: do_ml   !< Add ML term to lattice Hamiltonian (0/1)
      integer, intent(in) :: do_mml   !< Add MML term to lattice Hamiltonian (0/1)
      character(len=1) :: mode   !< Simulation mode (S=SD, M=MC, H=MC Heat Bath, P=LD, R=SLD)

      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: uvec  !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: latt_external_field !< External electric field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: latt_time_external_field !< External time-dependent electric field

      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: eeff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: eeff1 !< Internal effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: eeff2 !< External field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: eeff3 !< Internal effective field from mixed spin-lattice Hamiltonians

      !.. Local scalars
      integer :: i, k
      real(dblprec), dimension(3) :: tfield, eeff_l, eeff_m

      !$omp parallel do default(shared) schedule(static) private(i, k, eeff_l, eeff_m, tfield) collapse(2)
      do i=start_atom, stop_atom
         do k=1, Mensemble

            eeff_l=0.0_dblprec
            eeff_m=0.0_dblprec

            ! SLKTODO change to use eeff1 and eeff3 as arguments to avoid copying?

            ! LL term
            if(do_ll==1) call ll_field(i, k, eeff_l)

            ! MML term
            if(do_ml==1 .and. mode=='R') call ml_efield(i, k, eeff_m)

            ! MML term
            if(do_mml==1 .and. mode=='R') call mml_efield(i, k, eeff_m)
            if(do_mml==2 .and. mode=='R') call mml_efield_diag(i, k, eeff_m)

            eeff1(1:3,i,k) = eeff_l
            eeff2(1:3,i,k) = latt_external_field(1:3,i,k) + latt_time_external_field(1:3,i,k)
            eeff3(1:3,i,k) = eeff_m
            eeff(1:3,i,k)  = eeff1(1:3,i,k) + eeff2(1:3,i,k) + eeff3(1:3,i,k)

         end do
      end do
      !$omp end parallel do


10004 format (i8,30es16.8)


   contains


      !---------------ll_field---------------!
      subroutine ll_field(i, k, field)
         !
         !.. Implicit declarations
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         !
         integer :: j

         !LL term
         do j=1,ll_listsize(i)
            field(1) = field(1) - ( &
               ll_tens(1,1,j,i) * uvec(1,ll_list(1,j,i),k) + &
               ll_tens(2,1,j,i) * uvec(2,ll_list(1,j,i),k) + &
               ll_tens(3,1,j,i) * uvec(3,ll_list(1,j,i),k) )
            field(2) = field(2) - ( &
               ll_tens(1,2,j,i) * uvec(1,ll_list(1,j,i),k) + &
               ll_tens(2,2,j,i) * uvec(2,ll_list(1,j,i),k) + &
               ll_tens(3,2,j,i) * uvec(3,ll_list(1,j,i),k) )
            field(3) = field(3) - ( &
               ll_tens(1,3,j,i) * uvec(1,ll_list(1,j,i),k) + &
               ll_tens(2,3,j,i) * uvec(2,ll_list(1,j,i),k) + &
               ll_tens(3,3,j,i) * uvec(3,ll_list(1,j,i),k) )
         end do
         !write(301,'(a,3f10.6)') 'll-field hc   ', field(1:3)
         !field=0_dblprec

         !do j=1,ll_listsize(i)
         !   do ja=1,3
         !      do ia=1,3
         !         write(303,*) j, ja, ia
         !         field(ja) = field(ja) - ll_tens(ia,ja,j,i) &
         !              * uvec(ia,ll_list(1,j,i),k)
         !      end do
         !   end do
         !end do
         !write(302,'(a,3f10.6)') 'll-field loop ', field(1:3)

         ! The factor from double counting is canceled with
         ! the factor 1/2 for Phi_ij in the LL Hamiltonian
         ! Compare for the Heisenberg Hamiltonian and magnetic fields
         !field(1:3) = field(1:3) * 2.0_dblprec 

      end subroutine ll_field


      !---------------ml_field---------------!
      subroutine ml_efield(i, k, field)
         !
         !.. Implicit declarations
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         !
         integer :: j

         ! ML term
         !! SLKTODO Change sign convention of lattice Hamiltonian? Would be logical
         !! in particular with regard to the mixed Hamiltonian terms
         do j=1,lm_listsize(i)
            field(1) = field(1) - ( &
               lm_tens(1,1,j,i) * emomM(1,lm_list(1,j,i),k) + &
               lm_tens(2,1,j,i) * emomM(2,lm_list(1,j,i),k) + &
               lm_tens(3,1,j,i) * emomM(3,lm_list(1,j,i),k) )
            field(2) = field(2) - ( &
               lm_tens(1,2,j,i) * emomM(1,lm_list(1,j,i),k) + &
               lm_tens(2,2,j,i) * emomM(2,lm_list(1,j,i),k) + &
               lm_tens(3,2,j,i) * emomM(3,lm_list(1,j,i),k) )
            field(3) = field(3) - ( &
               lm_tens(1,3,j,i) * emomM(1,lm_list(1,j,i),k) + &
               lm_tens(2,3,j,i) * emomM(2,lm_list(1,j,i),k) + &
               lm_tens(3,3,j,i) * emomM(3,lm_list(1,j,i),k) )
         end do

      end subroutine ml_efield


      !---------------mml_efield---------------!
      subroutine mml_efield(i, k, field)
         !
         !.. Implicit declarations
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         !
         integer :: j
         integer :: ja, ka

         real(dblprec) :: mdot,hdot

         ! MML term
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
!!$omp simd
#endif
         do j=1,lmm_listsize(i)

            !dbg_count = 0
            !if(do_hoc_debug==1) then
            !   write(231,*) 'i ', i, 'j ', j
            !end if
            do ka=1,3
               !mdot=0.0_dblprec
               hdot=0.5_dblprec*emomM(ka,lmm_list(2,j,i),k)
               do ja=1,3
                  mdot=emomM(ja,lmm_list(1,j,i),k) * hdot
                  !do ia=1,3
                     field(:) = field(:) - lmm_tens(:,ja,ka,j,i) * mdot
                        !* emomM(ja,lmm_list(1,j,i),k) * emomM(ka,lmm_list(2,j,i),k) * 0.5_dblprec
                     !if(do_hoc_debug==1) then
                     !   dbg_field = 0.0_dblprec
                     !   dbg_count = dbg_count + 1
                     !   dbg_field(ia) = -lmm_tens(ia,ja,ka,j,i) * emomM(ja,lmm_list(1,j,i),k) * emomM(ka,lmm_list(2,j,i),k)
                     !   write(231,'(a,i4,a,3f10.6)') ' dbg_count ', dbg_count, ' dbg_mml_efield ', dbg_field(1:3)
                     !end if
                  !end do
               end do
            end do

         end do

      end subroutine mml_efield


      !---------------mml_efield---------------!
      subroutine mml_efield_diag(i, k, field)
         !
         !.. Implicit declarations
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         !
         integer :: j

         real(dblprec), dimension(3) :: lfield
         real(dblprec) :: mdot

         ! MML term
         lfield=0.0_dblprec
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
!!$omp simd
#endif
         do j=1,lmm_listsize(i)

            mdot=emomM(1,lmm_list(1,j,i),k) * emomM(1,lmm_list(2,j,i),k)  + emomM(2,lmm_list(1,j,i),k) * emomM(2,lmm_list(2,j,i),k)  + emomM(3,lmm_list(1,j,i),k) * emomM(3,lmm_list(2,j,i),k)
            lfield(:) = lfield(:) + lmm_tens_diag(:,j,i) * mdot
            !lfield(:) = lfield(:) + lmm_tens(:,1,1,j,i) * mdot
            !field(:) = field(:) - 0.5_dblprec * lmm_tens(:,1,1,j,i) * mdot

         end do
         field(:) = field(:) - 0.5_dblprec*lfield(:)

      end subroutine mml_efield_diag



   end subroutine effective_latticefield


   !> Calculate effective magnetic field from mixed spin-lattice Hamiltonians
   subroutine effective_bmixedfield(Natom, Mensemble, start_atom, stop_atom, &
      do_ml, do_mml, uvec, emomM, beff, beff3)

      use Constants

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: start_atom !< Atom to start loop for
      integer, intent(in) :: stop_atom !< Atom to end loop for

      integer, intent(in) :: do_ml   !< Add ML term to lattice Hamiltonian (0/1)
      integer, intent(in) :: do_mml   !< Add MML term to lattice Hamiltonian (0/1)
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: uvec  !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff3 !< Internal effective magnetic field from 
      !mixed spin-lattice Hamiltonians
      !.. Local scalars
      integer :: i, k
      real(dblprec), dimension(3) :: beff_l
      real(dblprec) :: fc

      ! Factors for mRy energy conversion
      ! Here fc is used, the double counting is treated 
      ! for each coupling separately
      fc = mry/mub
      !fc2 = 2.0_dblprec*mry/mub

      !$omp parallel do default(shared) schedule(static) private(i, k, beff_l) collapse(2)
      do i=start_atom, stop_atom
         do k=1, Mensemble

            beff_l=0.0_dblprec

            ! SLDTODO Change to use beff3 as argument to avoid copying?

            ! ML term
            if(do_ml==1) call ml_bfield(i, k, beff_l)

            ! MML term
            if(do_mml==1) call mml_bfield(i, k, beff_l)
            if(do_mml==2) call mml_bfield_diag(i, k, beff_l)


            beff3(1:3,i,k)= fc * beff_l
            !beff3(1:3,i,k)= beff_l
            beff(1:3,i,k) = beff(1:3,i,k)+beff3(1:3,i,k)

         end do
      end do
      !$omp end parallel do


10004 format (i8,30es16.8)


   contains


      !---------------ml_bfield---------------!
      subroutine ml_bfield(i, k, field)
         !
         !.. Implicit declarations
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         !
         integer :: j
         integer :: ia, ja


         ! ML term
         do j=1,ml_listsize(i)

            do ja=1,3
               do ia=1,3
                  field(ja) = field(ja) - ml_tens(ia,ja,j,i) &
                     * uvec(ia,ml_list(1,j,i),k)
               end do
            end do

         end do

      end subroutine ml_bfield

      !---------------mml_bfield---------------!
      subroutine mml_bfield(i, k, field)
         !
         !.. Implicit declarations
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         !
         integer :: j
         integer :: ia, ka

         real(dblprec) :: mdot

         ! MML term
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
!!$omp simd
#endif
         do j=1,mml_listsize(i)

            !dbg_count = 0
            !if(do_hoc_debug==1) then
            !   write(232,*) 'i ', i, 'j ', j
            !end if

            do ka=1,3
                  do ia=1,3
                  mdot=emomM(ia,mml_list(1,j,i),k) * uvec(ka,mml_list(2,j,i),k)
!                 do ja=1,3
                     field(:) = field(:) - mml_tens(ia,:,ka,j,i)  * mdot
                        !* emomM(ia,mml_list(1,j,i),k) * uvec(ka,mml_list(2,j,i),k)
                     !if(do_hoc_debug==1) then
                     !   dbg_field = 0.0_dblprec
                     !   dbg_count = dbg_count + 1
                     !   dbg_field(ja) = -mml_tens(ia,ja,ka,j,i) &
                     !     * emomM(ia,mml_list(1,j,i),k) * uvec(ka,mml_list(2,j,i),k)
                     !   write(232,'(a,i4,a,3f10.6)') ' dbg_count ', dbg_count, ' dbg_mml_bfield ', dbg_field(1:3)
                     !end if
!                 end do
               end do
            end do

         end do

         ! Factor of two to account for double counting in Hamiltonian
         !field(1:3) = field(1:3) * 2.0_dblprec

      end subroutine mml_bfield

      !---------------mml_bfield---------------!
      subroutine mml_bfield_diag(i, k, field)
         !
         !.. Implicit declarations
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         !
         integer :: j

         real(dblprec) :: audot

         ! MML term
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
!!$omp simd
#endif
         do j=1,mml_listsize(i)

            !audot=mml_tens(1,1,1,j,i)*uvec(1,mml_list(2,j,i),k) +  mml_tens(1,1,2,j,i)*uvec(2,mml_list(2,j,i),k) +  mml_tens(1,1,3,j,i)*uvec(3,mml_list(2,j,i),k)
            audot=mml_tens_diag(1,j,i)*uvec(1,mml_list(2,j,i),k) +  mml_tens_diag(2,j,i)*uvec(2,mml_list(2,j,i),k) +  mml_tens_diag(3,j,i)*uvec(3,mml_list(2,j,i),k)
            field(:) = field(:) - audot * emomM(:,mml_list(1,j,i),k)

         end do



         ! Factor of two to account for double counting in Hamiltonian
         !field(1:3) = field(1:3) * 2.0_dblprec

      end subroutine mml_bfield_diag


   end subroutine effective_bmixedfield


   !> Calculate lattice Hamiltonian and spin-lattice Hamiltonian energies
   subroutine calc_lattenergies(Natom, Mensemble, start_ensemb, stop_ensemb, start_atom, stop_atom, &
                                !subroutine calc_lattenergies(Natom, Mensemble, start_atom, stop_atom, &
      do_ll, do_ml, do_mml, mode, &
      uvec, emomM, latt_external_field, latt_time_external_field, &
      ll_energy, ml_energy, mml_energy,  &
      ldpot_energy, sdpot_energy, sldpot_energy, totpot_energy, sld_single_energy, &
      mm_energy0, ammom_inp, aemom_inp, NA)
      !ldpot_energy, sdpot_energy, sldpot_energy, totpot_energy, sld_single_energy)
      !ldpot_energy, sdpot_energy, sldpot_energy, totpot_energy)

      use Constants

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: start_ensemb !< Ensemble to start loop for
      integer, intent(in) :: stop_ensemb !< Ensemble to stop loop for
      integer, intent(in) :: start_atom !< Atom to start loop for
      integer, intent(in) :: stop_atom !< Atom to end loop for

      integer, intent(in) :: do_ll   !< Add LL term to lattice Hamiltonian (0/1)
      integer, intent(in) :: do_ml   !< Add ML term to lattice Hamiltonian (0/1)
      integer, intent(in) :: do_mml   !< Add MML term to lattice Hamiltonian (0/1)
      character(len=1) :: mode   !< Simulation mode (S=SD, M=MC, H=MC Heat Bath, P=LD, R=SLD, B=SLD MC)

      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: uvec  !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: latt_external_field !< External electric field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: latt_time_external_field !< External time-dependent electric field

      real(dblprec), dimension(Mensemble), intent(out) :: ll_energy        !< Total LL energy
      real(dblprec), dimension(Mensemble), intent(out) :: ml_energy        !< Total ML energy
      real(dblprec), dimension(Mensemble), intent(out) :: mml_energy       !< Total MML energy

      !real(dblprec), dimension(Mensemble), intent(out) :: mm_energy        !< Total MM energy     ***************************
      real(dblprec), dimension(Mensemble), intent(out) :: ldpot_energy     !< LD potential energy
      real(dblprec), dimension(Mensemble), intent(in)  :: sdpot_energy     !< SD potential energy
      real(dblprec), dimension(Mensemble), intent(out) :: sldpot_energy    !< SLD potential energy (without pure LD or SD potential energies)
      real(dblprec), dimension(Mensemble), intent(out) :: totpot_energy    !< Total potential energy: LD + SD + SLD. No kinetic energy!
      real(dblprec), dimension(Mensemble), intent(out) :: sld_single_energy !< Trial potential energy: LD + Heisenberg + SLD. No kinetic energy!
      real(dblprec), dimension(Mensemble), intent(out) :: mm_energy0        !< Total MM ground state energy     ***************************

      real(dblprec), intent(in), dimension(:,:,:) :: ammom_inp  !< Magnetic moment magnitudes from input (for alloys)
      real(dblprec), intent(in), allocatable, dimension(:,:,:,:) :: aemom_inp  !< Magnetic moment directions from input (for alloys)
      integer, intent(in) :: NA  !< Number of atoms in one cell

      !.. Local scalars
      integer :: i, k
      real(dblprec) :: tmp_energy, tmp_ll
      real(dblprec) :: tmp_ml, tmp_mml,  tmp_mm, tmp_mm0

      real(dblprec), dimension(Mensemble) :: mm_energy

      ll_energy = 0.0_dblprec
      ml_energy = 0.0_dblprec
      mml_energy = 0.0_dblprec

      ldpot_energy = 0.0_dblprec
      sldpot_energy = 0.0_dblprec
      totpot_energy = 0.0_dblprec

      mm_energy = 0.0_dblprec
      mm_energy0 = 0.0_dblprec
      sld_single_energy = 0.0_dblprec

      tmp_energy = 0.0_dblprec

      !SLDTODO Add calculation of the energy from coupling to external electric field

      !SLDTODO Use OMP reduce statements to sum up the contributions to the energies
      !SLDTODO But, right now this subroutine is also used for the Monte Carlo trial moves
      !SLDTODO for which reduce statements are not relevant.
      do k=start_ensemb, stop_ensemb

         tmp_ll=0.0_dblprec;tmp_ml=0.0_dblprec;tmp_mml=0.0_dblprec;tmp_mm=0.0_dblprec;tmp_mm0=0.0_dblprec

!$omp parallel do default(shared) schedule(static) private(i, tmp_energy)  &
!$omp& reduction(+:tmp_ll,tmp_ml,tmp_mml,tmp_mm,tmp_mm0)
         do i=start_atom, stop_atom

            !do k=1, Mensemble

            ! LL term
            if(do_ll==1) then
               call calc_ll_energy(i, k, tmp_energy)
               tmp_ll=tmp_ll+tmp_energy
               !ll_energy(k) = ll_energy(k) + tmp_energy
            end if

            ! ML term
            if(do_ml==1 .and. (mode=='R' .or. mode=='B') ) then
               call calc_ml_energy(i, k, tmp_energy)
               tmp_ml = tmp_ml + tmp_energy
               !ml_energy(k) = ml_energy(k) + tmp_energy
            end if

            ! MML term
            if(do_mml>=1 .and. (mode=='R' .or. mode=='B') ) then
               call calc_mml_energy(i, k, tmp_energy)
               tmp_mml = tmp_mml + tmp_energy
               !mml_energy(k) = mml_energy(k) + tmp_energy
            end if


            ! MM term ground state
            if(mode=='R'.and..not.mm_energy0_calc) then
               call calc_mm_energy0(i, k, tmp_energy, ammom_inp, aemom_inp, NA)
               tmp_mm0 = tmp_mm0 + tmp_energy
            end if

         end do
!$omp end parallel do

         ll_energy(k) = tmp_ll
         ml_energy(k) = tmp_ml
         mml_energy(k) = tmp_mml
         mm_energy(k) = tmp_mm
         mm_energy0(k) = tmp_mm0
         
         ! For energy measurements

         ldpot_energy(k) = ll_energy(k) 
         sldpot_energy(k) = ml_energy(k) + mml_energy(k) 

         ! Substracts the ground state Heisenberg energy
         totpot_energy(k) = sdpot_energy(k) + ldpot_energy(k) + sldpot_energy(k) - mm_energy0(k)


         ! For Monte Carlo trial moves
         sld_single_energy(k) = mm_energy(k) + ldpot_energy(k) + sldpot_energy(k)

      end do


10004 format (i8,30es16.8)


   contains


      !---------------ll_energy---------------!
      subroutine calc_ll_energy(i, k, energy)
         !
         !.. Implicit declarations
         implicit none

         integer, intent(in) :: i !< Current atom
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), intent(out) :: energy !< LL energy
         !
         integer :: j

         ! LL term
         energy=0.0_dblprec
         do j=1,ll_listsize(i)
            energy = energy + &
               ll_tens(1,1,j,i) * uvec(1,ll_list(1,j,i),k) * uvec(1,i,k) + &
               ll_tens(2,1,j,i) * uvec(2,ll_list(1,j,i),k) * uvec(1,i,k) + &
               ll_tens(3,1,j,i) * uvec(3,ll_list(1,j,i),k) * uvec(1,i,k) + &
               
               ll_tens(1,2,j,i) * uvec(1,ll_list(1,j,i),k) * uvec(2,i,k) + &
               ll_tens(2,2,j,i) * uvec(2,ll_list(1,j,i),k) * uvec(2,i,k) + &
               ll_tens(3,2,j,i) * uvec(3,ll_list(1,j,i),k) * uvec(2,i,k) + &
               
               ll_tens(1,3,j,i) * uvec(1,ll_list(1,j,i),k) * uvec(3,i,k) + &
               ll_tens(2,3,j,i) * uvec(2,ll_list(1,j,i),k) * uvec(3,i,k) + &
               ll_tens(3,3,j,i) * uvec(3,ll_list(1,j,i),k) * uvec(3,i,k)
         end do

         ! Here double counting correction is necessary, c.f. ll_field
         energy = energy * 0.5_dblprec

      end subroutine calc_ll_energy


      !---------------ml_energy---------------!
      subroutine calc_ml_energy(i, k, energy)
         !
         !.. Implicit declarations
         implicit none

         integer, intent(in) :: i !< Current atom
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), intent(out) :: energy !< ML energy
         !
         integer :: j

         ! ML term
         energy=0.0_dblprec
         do j=1,ml_listsize(i)
            energy = energy + &
               ml_tens(1,1,j,i) * uvec(1,ml_list(1,j,i),k) * emomM(1,i,k) + &
               ml_tens(2,1,j,i) * uvec(2,ml_list(1,j,i),k) * emomM(1,i,k) + &
               ml_tens(3,1,j,i) * uvec(3,ml_list(1,j,i),k) * emomM(1,i,k) + &
               
               ml_tens(1,2,j,i) * uvec(1,ml_list(1,j,i),k) * emomM(2,i,k) + &
               ml_tens(2,2,j,i) * uvec(2,ml_list(1,j,i),k) * emomM(2,i,k) + &
               ml_tens(3,2,j,i) * uvec(3,ml_list(1,j,i),k) * emomM(2,i,k) + &
               
               ml_tens(1,3,j,i) * uvec(1,ml_list(1,j,i),k) * emomM(3,i,k) + &
               ml_tens(2,3,j,i) * uvec(2,ml_list(1,j,i),k) * emomM(3,i,k) + &
               ml_tens(3,3,j,i) * uvec(3,ml_list(1,j,i),k) * emomM(3,i,k)
         end do

         ! single counting in displacement and in spin
         !energy = energy*1.0_dblprec

      end subroutine calc_ml_energy


      !---------------mml_energy---------------!
      subroutine calc_mml_energy(i, k, energy)
         !
         !.. Implicit declarations
         implicit none

         integer, intent(in) :: i !< Current atom
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), intent(out) :: energy !< MML energy

         !
         integer :: j
         real(dblprec) :: fenergy

         ! MML energy
         energy=0.0_dblprec
         do j=1,mml_listsize(i)

            fenergy =  &
               mml_tens(1,1,1,j,i) * emomM(1,mml_list(1,j,i),k) * uvec(1,mml_list(2,j,i),k) * emomM(1,i,k) + &
               mml_tens(2,1,1,j,i) * emomM(2,mml_list(1,j,i),k) * uvec(1,mml_list(2,j,i),k) * emomM(1,i,k) + &
               mml_tens(3,1,1,j,i) * emomM(3,mml_list(1,j,i),k) * uvec(1,mml_list(2,j,i),k) * emomM(1,i,k) + &
               mml_tens(1,1,2,j,i) * emomM(1,mml_list(1,j,i),k) * uvec(2,mml_list(2,j,i),k) * emomM(1,i,k) + &
               mml_tens(2,1,2,j,i) * emomM(2,mml_list(1,j,i),k) * uvec(2,mml_list(2,j,i),k) * emomM(1,i,k) + &
               mml_tens(3,1,2,j,i) * emomM(3,mml_list(1,j,i),k) * uvec(2,mml_list(2,j,i),k) * emomM(1,i,k) + &
               mml_tens(1,1,3,j,i) * emomM(1,mml_list(1,j,i),k) * uvec(3,mml_list(2,j,i),k) * emomM(1,i,k) + &
               mml_tens(2,1,3,j,i) * emomM(2,mml_list(1,j,i),k) * uvec(3,mml_list(2,j,i),k) * emomM(1,i,k) + &
               mml_tens(3,1,3,j,i) * emomM(3,mml_list(1,j,i),k) * uvec(3,mml_list(2,j,i),k) * emomM(1,i,k) + &
               
               mml_tens(1,2,1,j,i) * emomM(1,mml_list(1,j,i),k) * uvec(1,mml_list(2,j,i),k) * emomM(2,i,k) + &
               mml_tens(2,2,1,j,i) * emomM(2,mml_list(1,j,i),k) * uvec(1,mml_list(2,j,i),k) * emomM(2,i,k) + &
               mml_tens(3,2,1,j,i) * emomM(3,mml_list(1,j,i),k) * uvec(1,mml_list(2,j,i),k) * emomM(2,i,k) + &
               mml_tens(1,2,2,j,i) * emomM(1,mml_list(1,j,i),k) * uvec(2,mml_list(2,j,i),k) * emomM(2,i,k) + &
               mml_tens(2,2,2,j,i) * emomM(2,mml_list(1,j,i),k) * uvec(2,mml_list(2,j,i),k) * emomM(2,i,k) + &
               mml_tens(3,2,2,j,i) * emomM(3,mml_list(1,j,i),k) * uvec(2,mml_list(2,j,i),k) * emomM(2,i,k) + &
               mml_tens(1,2,3,j,i) * emomM(1,mml_list(1,j,i),k) * uvec(3,mml_list(2,j,i),k) * emomM(2,i,k) + &
               mml_tens(2,2,3,j,i) * emomM(2,mml_list(1,j,i),k) * uvec(3,mml_list(2,j,i),k) * emomM(2,i,k) + &
               mml_tens(3,2,3,j,i) * emomM(3,mml_list(1,j,i),k) * uvec(3,mml_list(2,j,i),k) * emomM(2,i,k) + &
               
               mml_tens(1,3,1,j,i) * emomM(1,mml_list(1,j,i),k) * uvec(1,mml_list(2,j,i),k) * emomM(3,i,k) + &
               mml_tens(2,3,1,j,i) * emomM(2,mml_list(1,j,i),k) * uvec(1,mml_list(2,j,i),k) * emomM(3,i,k) + &
               mml_tens(3,3,1,j,i) * emomM(3,mml_list(1,j,i),k) * uvec(1,mml_list(2,j,i),k) * emomM(3,i,k) + &
               mml_tens(1,3,2,j,i) * emomM(1,mml_list(1,j,i),k) * uvec(2,mml_list(2,j,i),k) * emomM(3,i,k) + &
               mml_tens(2,3,2,j,i) * emomM(2,mml_list(1,j,i),k) * uvec(2,mml_list(2,j,i),k) * emomM(3,i,k) + &
               mml_tens(3,3,2,j,i) * emomM(3,mml_list(1,j,i),k) * uvec(2,mml_list(2,j,i),k) * emomM(3,i,k) + &
               mml_tens(1,3,3,j,i) * emomM(1,mml_list(1,j,i),k) * uvec(3,mml_list(2,j,i),k) * emomM(3,i,k) + &
               mml_tens(2,3,3,j,i) * emomM(2,mml_list(1,j,i),k) * uvec(3,mml_list(2,j,i),k) * emomM(3,i,k) + &
               mml_tens(3,3,3,j,i) * emomM(3,mml_list(1,j,i),k) * uvec(3,mml_list(2,j,i),k) * emomM(3,i,k) 

            ! Original version
            energy = energy + fenergy
            ! Dual version
            !energy = energy + fenergy * 0.5_dblprec + benergy * 0.5_dblprec

         end do

         ! Here double counting correction is necessary, c.f. mml_efield
         ! and mml_bfield
         energy = energy * 0.5_dblprec

      end subroutine calc_mml_energy


      !--------------- mm_energy ---------------!
      !> Heisenberg 
      subroutine calc_mm_energy(i, k, energy)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), intent(out) :: energy !< MM energy

         real(dblprec), dimension(3) :: field !< Effective field

         integer :: j, ih

         field = 0.0_dblprec
         energy = 0.0_dblprec
         ih=ham%aHam(i)
         do j=1,ham%nlistsize(ih)
            field = field + ham%ncoup(j,ih,1)*emomM(:,ham%nlist(j,i),k)
         end do
         energy = -field(1) * emomM(1,i,k) - field(2) * emomM(2,i,k) - field(3) * emomM(3,i,k)
         !write(*,'(a,3f10.6,a,f10.6)') 'field ', field(1:3), ' energy ', energy

      end subroutine calc_mm_energy


      !--------------- mm_energy0 ---------------!
      !> Heisenberg ground state energy
      subroutine calc_mm_energy0(i, k, energy, ammom_inp, aemom_inp, NA)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), intent(out) :: energy !< MM energy
         real(dblprec), intent(in), dimension(:,:,:) :: ammom_inp  !< Magnetic moment magnitudes from input (for alloys)
         real(dblprec), intent(in), allocatable, dimension(:,:,:,:) :: aemom_inp  !< Magnetic moment directions from input (for alloys)
         integer, intent(in) :: NA  !< Number of atoms in one cell

         real(dblprec), dimension(3) :: field !< Effective field
         real(dblprec) :: fcinv, fc

         integer :: j, ih
         integer :: iatom, jatom

         real(dblprec), dimension(3) :: emomI, emomJ

         ! Factor for energy scale
         fcinv = mub/mry
         fc = mry/mub

         field = 0.0_dblprec
         energy = 0.0_dblprec
         ih=ham%aHam(i)
         do j=1,ham%nlistsize(ih)
            jatom = mod(ham%nlist(j,i)-1,NA) + 1
            !write(*,*) 'ham%nlist(j,i), jatom ', ham%nlist(j,i), jatom
            emomJ(1:3) = aemom_inp(1:3,jatom,1,1) * ammom_inp(jatom,1,1)
            field(1:3) = field(1:3) + ham%ncoup(j,ih,1) * emomJ(1:3)
         end do
         iatom = mod(i-1,NA) + 1
         !write(*,*) 'i, iatom ', i, iatom
         emomI(1:3) = aemom_inp(1:3,iatom,1,1) * ammom_inp(iatom,1,1)
         energy = -field(1) * emomI(1) - field(2) * emomI(2) - field(3) * emomI(3)

         !double counting in spins
         energy = energy * 0.5_dblprec

         !Conversion to mRyd
         energy = energy * fcinv

         !energy = -field(1) * emomM(1,i,k) - field(2) * emomM(2,i,k) - field(3) * emomM(3,i,k)
         !write(*,'(a,3f10.6,a,f10.6)') 'field ', field(1:3), ' energy ', energy

      end subroutine calc_mm_energy0


   end subroutine calc_lattenergies


end module LatticeHamiltonianActions
