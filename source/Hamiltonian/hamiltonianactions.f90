!-------------------------------------------------------------------------------
! MODULE: HamiltonianActions
!> @brief
!> Calculate effective field by applying the derivative of the Hamiltonian
!> @details The effective field, \f$\mathbf{B}_i\f$, on an atom \f$\textit{i}\f$, is calculated from
!> \f$ \mathbf{B}_i=-\frac{\partial \mathbf{H}}{\partial \mathbf{m}_i},\f$ where primarily the part of
!> the Hamiltonian, \f$\mathbf{H}\f$, which represents interatomic exchange interactions,
!> \f$\mathbf{H}_\mathrm{ex}\f$, are considered. For this we use the classical Heisenberg Hamiltonian,
!> \f$ \mathbf{H}_\mathrm{ex}=-\frac{1}{2}\sum_{i\neq j}J_{ij}\mathbf{m}_i\cdot\mathbf{m}_j,\f$ where
!> \f$i\f$ and \f$j\f$ are atomic indices and \f$J_{ij}\f$ is the strength of the exchange interaction,
!> which is calculated from first principles theory.
!> @author Anders Bergman, Lars Bergqvist, Johan Hellsvik, Nikos Ntallis
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module HamiltonianActions

   use Profiling
   use Parameters
   use HamiltonianData
   use InputData, only : ham_inp

   implicit none

   interface effective_field
      module procedure effective_field_bare, effective_field_full
   end interface

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: effective_field
   !> @brief
   !> Calculate effective field by applying the derivative of the Hamiltonian
   !> @author Anders Bergman, Lars Bergqvist, Johan Hellsvik
   !> @todo Check consistency of terms wrt the input parameters, especially the anisotropies
   !> @todo Replace moment unit vectors emom with full length vectors emomM
   !ham%dm_vect !< DM vector \f$H_{DM}=\sum{D_{ij}\dot(m_i \times m_j)}\f$
   !> @todo Check the sign of the dipolar field
   !----------------------------------------------------------------------------
   ! Aliases with less arguments
   subroutine effective_field_bare(measure_energy)
      use InputData, only : Natom, Mensemble, NA, N1, N2, N3
      use MomentData, only : emomM, mmom
      use FieldData, only : external_field, time_external_field, beff, beff1, beff2
      use macrocells, only : Num_macro, cell_index, emomM_macro, macro_nlistsize
      implicit none

      real(dblprec) :: energy
      logical, intent(in), optional :: measure_energy


      call effective_field_full(Natom,Mensemble,1,Natom,   &
         emomM,mmom,external_field,time_external_field,beff,beff1,beff2,energy,       &
         Num_macro,cell_index,emomM_macro,    &
         macro_nlistsize,NA,N1,N2,N3,measure_energy)

   end subroutine effective_field_bare
      !
   subroutine effective_field_full(Natom,Mensemble,start_atom,stop_atom,   &
      emomM,mmom,external_field,time_external_field,beff,beff1,beff2,energy,         &
      Num_macro,cell_index,emomM_macro,    &
      macro_nlistsize,NA,N1,N2,N3,measure_energy)
      !
      use Constants, only : mry,mub
      use DipoleManager, only : dipole_field_calculation
      use HamiltonianTargetOrder, only : target_order_range
#ifdef _OPENMP
      use omp_lib, only : omp_get_num_threads, omp_get_thread_num
#endif
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: NA
      integer, intent(in) :: N1
      integer, intent(in) :: N2
      integer, intent(in) :: N3
      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Mensemble    !< Number of ensembles
      integer, intent(in) :: start_atom   !< Atom to start loop for
      integer, intent(in) :: stop_atom    !< Atom to end loop for
      integer, intent(in) :: Num_macro    !< Number of macrocells in the system
      integer, dimension(Natom), intent(in) :: cell_index            !< Macrocell index for each atom
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize   !< Number of atoms per macrocell
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: external_field  !< External magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: time_external_field !< External time-dependent magnetic field
      ! .. Output Variables
      real(dblprec), intent(out) :: energy !< Total energy
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff  !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff1 !< Internal effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff2 !< External field from application of Hamiltonian
      logical, intent(in), optional :: measure_energy !< Request the canonical energy reduction

      !.. Local scalars
      integer :: i,k,q,thread_id,nthreads,q_start,q_stop
      logical :: calculate_energy,ordered_targets,weighted_targets
      real(dblprec) :: atom_energy

      !.. Executable statements
      calculate_energy=.true.
      if (present(measure_energy)) calculate_energy=measure_energy

      ! Initialization of the energy
      energy=0.0_dblprec
      ! Initialization if the effective field
      beff=0.0_dblprec
      ! Wrapper for the calculation of the dipole-dipole interaction field
      ! The field is stored in the bfield array which then is passed to the main loop
      ! This is inefficient for the brute-force methods, but it the best way to ensure
      ! that the FFT approaches can be used in an appropriate way
      if (ham_inp%do_dip>0) then
         call timing(0,'Hamiltonian   ','OF')
         call timing(0,'Dipolar Int.  ','ON')
         call dipole_field_calculation(NA,N1,N2,N3,Natom,ham_inp%do_dip,Num_macro,          &
            Mensemble,stop_atom,start_atom,cell_index,macro_nlistsize,emomM,        &
            emomM_macro,ham%Qdip,ham%Qdip_macro,energy,beff,calculate_energy)
         call timing(0,'Dipolar Int.  ','OF')
         call timing(0,'Hamiltonian   ','ON')
      endif
      ! The target order is a permutation of physical atom IDs.  Partial
      ! callers retain the historical natural loop because their range is not
      ! generally a contiguous interval in the permutation.
      ordered_targets=.false.
      if (allocated(ham%target_order)) then
         ordered_targets=(size(ham%target_order)==Natom .and. start_atom==1 .and. stop_atom==Natom)
      endif
      weighted_targets=.false.
      if (ordered_targets .and. allocated(ham%target_work_prefix)) then
         weighted_targets=ham%target_order_weighted
      endif

      if (ordered_targets) then
         if (weighted_targets) then
            ! One contiguous target-order interval per thread, cut by the
            ! persistent cumulative neighbour-work array.
            if (calculate_energy) then
               !$omp parallel default(shared) private(i,k,q,thread_id,nthreads,q_start,q_stop,atom_energy) reduction(+:energy)
#ifdef _OPENMP
               thread_id=omp_get_thread_num()
               nthreads=omp_get_num_threads()
#else
               thread_id=0
               nthreads=1
#endif
               call target_order_range(ham%target_work_prefix,ham%target_total_work, &
                  nthreads,thread_id,q_start,q_stop)
               do q=q_start,q_stop
                  i=ham%target_order(q)
                  do k=1,Mensemble
                     call effective_field_atom(i,k,Natom,Mensemble,emomM,mmom,external_field, &
                        time_external_field,beff(:,i,k),beff1(:,i,k),beff2(:,i,k),atom_energy,.true.)
                     energy=energy+atom_energy
                  end do
               end do
               !$omp end parallel
            else
               !$omp parallel default(shared) private(i,k,q,thread_id,nthreads,q_start,q_stop,atom_energy)
#ifdef _OPENMP
               thread_id=omp_get_thread_num()
               nthreads=omp_get_num_threads()
#else
               thread_id=0
               nthreads=1
#endif
               call target_order_range(ham%target_work_prefix,ham%target_total_work, &
                  nthreads,thread_id,q_start,q_stop)
               do q=q_start,q_stop
                  i=ham%target_order(q)
                  do k=1,Mensemble
                     call effective_field_atom(i,k,Natom,Mensemble,emomM,mmom,external_field, &
                        time_external_field,beff(:,i,k),beff1(:,i,k),beff2(:,i,k),atom_energy,.false.)
                  end do
               end do
               !$omp end parallel
            endif
         else if (calculate_energy) then
            !$omp parallel do default(shared) schedule(static) private(i,k,q,atom_energy) reduction(+:energy)
            do q=1,Natom
               i=ham%target_order(q)
               do k=1,Mensemble
                  call effective_field_atom(i,k,Natom,Mensemble,emomM,mmom,external_field, &
                     time_external_field,beff(:,i,k),beff1(:,i,k),beff2(:,i,k),atom_energy,.true.)
                  energy=energy+atom_energy
               end do
            end do
            !$omp end parallel do
         else
            !$omp parallel do default(shared) schedule(static) private(i,k,q,atom_energy)
            do q=1,Natom
               i=ham%target_order(q)
               do k=1,Mensemble
                  call effective_field_atom(i,k,Natom,Mensemble,emomM,mmom,external_field, &
                     time_external_field,beff(:,i,k),beff1(:,i,k),beff2(:,i,k),atom_energy,.false.)
               end do
            end do
            !$omp end parallel do
         endif
      else if (calculate_energy) then
         !$omp parallel do default(shared) schedule(static) private(i,k,atom_energy) collapse(2) reduction(+:energy)
         do k=1, Mensemble
            do i=start_atom, stop_atom
               call effective_field_atom(i,k,Natom,Mensemble,emomM,mmom,external_field, &
                  time_external_field,beff(:,i,k),beff1(:,i,k),beff2(:,i,k),atom_energy,.true.)
               energy=energy+atom_energy
            end do
         end do
         !$omp end parallel do
      else
         !$omp parallel do default(shared) schedule(static) private(i,k,atom_energy) collapse(2)
         do k=1, Mensemble
            do i=start_atom, stop_atom
               call effective_field_atom(i,k,Natom,Mensemble,emomM,mmom,external_field, &
                  time_external_field,beff(:,i,k),beff1(:,i,k),beff2(:,i,k),atom_energy,.false.)
            end do
         end do
         !$omp end parallel do
      endif
      if (calculate_energy) energy = energy * mub / mry

   end subroutine effective_field_full

   ! Shared canonical field assembly for the energy-enabled and field-only paths.
   subroutine effective_field_atom(i,k,Natom,Mensemble,emomM,mmom,external_field, &
      time_external_field,beff,beff1,beff2,atom_energy,calculate_energy)
      implicit none

      integer, intent(in) :: i,k,Natom,Mensemble
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: external_field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: time_external_field
      real(dblprec), dimension(3), intent(inout) :: beff
      real(dblprec), dimension(3), intent(out) :: beff1
      real(dblprec), dimension(3), intent(out) :: beff2
      real(dblprec), intent(out) :: atom_energy
      logical, intent(in) :: calculate_energy

      real(dblprec), dimension(3) :: tfield, beff_s, beff_q

      beff_s=0.0_dblprec
      beff_q=0.0_dblprec

      if(ham_inp%do_jtensor/=1) then
         ! Heisenberg exchange term
         if(ham_inp%exc_inter=='N') then
            call heisenberg_field(i, k, beff_s,Natom,Mensemble,emomM)
         else
            call heisenberg_rescaling_field(i, k, beff_s,Natom,Mensemble,mmom,emomM)
         endif
      else
         call tensor_field(i, k, beff_s,Natom,Mensemble,emomM)
      end if

      ! Dzyaloshinskii-Moriya term
      if(ham_inp%do_dm==1) call dzyaloshinskii_moriya_field(i, k, beff_s,Natom,Mensemble,emomM)

      ! Symmetric anisotropic term
      if(ham_inp%do_sa==1) call symmetric_anisotropic_field(i, k, beff_s,Natom,Mensemble,emomM)

      ! Pseudo-Dipolar term
      if(ham_inp%do_pd==1) call pseudo_dipolar_field(i, k, beff_s,Natom,Mensemble,emomM)

      ! BIQDM term
      if(ham_inp%do_biqdm==1) call dzyaloshinskii_moriya_bq_field(i, k, beff_q,Natom,Mensemble,emomM)

      ! Biquadratic exchange term
      if(ham_inp%do_bq==1) call biquadratic_field(i, k, beff_q,Natom,Mensemble,emomM)

      ! Four-spin ring exchange term
      if(ham_inp%do_ring==1) call ring_field(i, k, beff_s,Natom,Mensemble,emomM)

      ! Scalar chiral term
      if(ham_inp%do_chir==1) call chirality_field(i, k, beff_s,Natom,Mensemble,emomM)

      ! Anisotropy
      if (ham_inp%do_anisotropy==1) then
         if (ham%taniso(i)==1) then
            call uniaxial_anisotropy_field(i, k, beff_s,Natom,Mensemble,ham_inp%mult_axis,emomM)
         elseif (ham%taniso(i)==2) then
            call cubic_anisotropy_field(i, k, beff_s,Natom,Mensemble,ham_inp%mult_axis,emomM)
         elseif (ham%taniso(i)==7)then
            call uniaxial_anisotropy_field(i, k, beff_s,Natom,Mensemble,ham_inp%mult_axis,emomM)
            tfield=0.0_dblprec
            call cubic_anisotropy_field(i, k, tfield,Natom,Mensemble,ham_inp%mult_axis,emomM)
            beff_q=beff_q+tfield*ham%sb(i)
         endif
      endif

      beff1=beff_s
      beff2=beff_q+external_field(1:3,i,k)+time_external_field(1:3,i,k)
      ! Here the dipole contribution is added since beff != 0 in that case.
      beff = beff + beff1 + beff2

      atom_energy=0.0_dblprec
      if (calculate_energy) then
         tfield=0.50_dblprec*(beff_s+2.0_dblprec*beff_q+2.0_dblprec*external_field(1:3,i,k)+time_external_field(1:3,i,k))
         atom_energy = -emomM(1,i,k)*tfield(1)-emomM(2,i,k)*tfield(2)-emomM(3,i,k)*tfield(3)
      endif

   end subroutine effective_field_atom

      !---------------heisenberg_field---------------!
      !> Heisenberg
      subroutine heisenberg_field(i, k, field,Natom,Mensemble,emomM)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         integer, intent(in) :: Natom        !< Number of atoms in system
         integer, intent(in) :: Mensemble    !< Number of ensembles
         real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field

         integer :: j, ih, x, n_neigh
         real(dblprec) :: bx, by, bz, coup

         ih=ham%aHam(i)
         n_neigh=ham%nlistsize(ih)

         ! Keep the gather formulation, but make the three scalar reductions
         ! explicit. This avoids constructing and repeatedly updating a
         ! three-element array expression in the neighbour loop.
         bx=field(1)
         by=field(2)
         bz=field(3)
         do j=1,n_neigh
            x=ham%nlist(j,i)
            coup=ham%ncoup(j,ih,1)
            bx=bx+coup*emomM(1,x,k)
            by=by+coup*emomM(2,x,k)
            bz=bz+coup*emomM(3,x,k)
         end do
         field(1)=bx
         field(2)=by
         field(3)=bz
      end subroutine heisenberg_field

      !---------------heisenberg_rescaling_field---------------!
      !> Heisenberg
      subroutine heisenberg_rescaling_field(i,k,field,Natom,Mensemble,mmom,emomM)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         integer, intent(in) :: Natom        !< Number of atoms in system
         integer, intent(in) :: Mensemble    !< Number of ensembles
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment
         real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector

         integer :: j, ih
         real(dblprec) :: excscale

         ih=ham%aHam(i)
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422) && __INTEL_COMPILER < 1800
!         !$omp simd private(excscale) reduction(+:field)
#endif
         do j=1,ham%nlistsize(ih)
            excscale=abs(sum(emomM(:,ham%nlist(j,i),k)*emomM(:,i,k)))/(mmom(ham%nlist(j,i),k)*mmom(i,k))
            field = field + ((excscale*ham%ncoup(j,ih,1)+(1.0_dblprec-excscale)*ham%ncoupD(j,ih,1)))*emomM(:,ham%nlist(j,i),k)
         end do
      end subroutine heisenberg_rescaling_field

      !---------------tensor_heisenberg_field---------------!
      !> @brief Calculates heisenberg, DM and anisotropy through one tensor
      !>
      !> @Date 09/15/2014 - Thomas Nystrand
      !> - Added 0.5*m_i*I part
      subroutine tensor_field(i, k, field,Natom,Mensemble,emomM)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         integer, intent(in) :: Natom        !< Number of atoms in system
         integer, intent(in) :: Mensemble    !< Number of ensembles
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
         !
         integer :: j ! Neighbourlist index
         integer :: x ! Exchange index
         integer :: ih ! Hamiltonian index


         !Exchange term
         ih=ham%aHam(i)
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422) && __INTEL_COMPILER < 1800
!         !$omp simd reduction(+:field)
#endif
         do j=1,ham%nlistsize(ih)
            x = ham%nlist(j,i);
            ! Matrix tensor multiplication: f = f+0.5*I*m_i+0.5*m_i*I
            ! Faster then: field = field + 0.5*MATMUL(ham%j_tens(:,:,j,i),emomM(:,x,k)) + 0.5*MATMUL(emomM(:,x,k),ham%j_tens(:,:,j,i))
            !field = field + ham%ncoup(j,i)*emomM(:,x,k)
            !field = field + MATMUL(ham%j_tens(:,:,j,i),emomM(:,x,k)) 
            field = field + ham%j_tens(:,1,j,i)*emomM(1,x,k) + ham%j_tens(:,2,j,i)*emomM(2,x,k) + ham%j_tens(:,3,j,i)*emomM(3,x,k)
            !!! field(1) = field(1) &
            !!!    + 0.50_dblprec*(        &
            !!!    + ham%j_tens(1,1,j,ih)*emomM(1,x,k) + ham%j_tens(1,2,j,ih)*emomM(2,x,k) + ham%j_tens(1,3,j,ih)*emomM(3,x,k) &
            !!!    + emomM(1,x,k)*ham%j_tens(1,1,j,ih) + emomM(2,x,k)*ham%j_tens(2,1,j,ih) + emomM(3,x,k)*ham%j_tens(3,1,j,ih) &
            !!!    )
            !!! field(2) = field(2) &
            !!!    + 0.50_dblprec*(        &
            !!!    + ham%j_tens(2,1,j,ih)*emomM(1,x,k) + ham%j_tens(2,2,j,ih)*emomM(2,x,k) + ham%j_tens(2,3,j,ih)*emomM(3,x,k) &
            !!!    + emomM(1,x,k)*ham%j_tens(1,2,j,ih) + emomM(2,x,k)*ham%j_tens(2,2,j,ih) + emomM(3,x,k)*ham%j_tens(3,2,j,ih) &
            !!!    )
            !!! field(3) = field(3) &
            !!!    + 0.50_dblprec*(        &
            !!!    + ham%j_tens(3,1,j,ih)*emomM(1,x,k) + ham%j_tens(3,2,j,ih)*emomM(2,x,k) + ham%j_tens(3,3,j,ih)*emomM(3,x,k) &
            !!!    + emomM(1,x,k)*ham%j_tens(1,3,j,ih) + emomM(2,x,k)*ham%j_tens(2,3,j,ih) + emomM(3,x,k)*ham%j_tens(3,3,j,ih) &
            !!!    )
         end do
      end subroutine tensor_field

      !---------------dzyaloshinskii_moriya_field---------------!
      !> DM-field
      subroutine dzyaloshinskii_moriya_field(i, k, field,Natom,Mensemble,emomM)
         !
         !.. Implicit declarations
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         integer, intent(in) :: Natom        !< Number of atoms in system
         integer, intent(in) :: Mensemble    !< Number of ensembles
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
         !
         integer :: j, ih

         ! Dzyaloshinskii_moriya term
         ih=ham%aHam(i)
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422) && __INTEL_COMPILER < 1800
!         !$omp simd reduction(+:field)
#endif
         do j=1,ham%dmlistsize(ih)
            field(1) = field(1) - ham%dm_vect(3,j,ih)*emomM(2,ham%dmlist(j,i),k) +&
               ham%dm_vect(2,j,ih)*emomM(3,ham%dmlist(j,i),k)
            field(2) = field(2) - ham%dm_vect(1,j,ih)*emomM(3,ham%dmlist(j,i),k) +&
               ham%dm_vect(3,j,ih)*emomM(1,ham%dmlist(j,i),k)
            field(3) = field(3) - ham%dm_vect(2,j,ih)*emomM(1,ham%dmlist(j,i),k) +&
               ham%dm_vect(1,j,ih)*emomM(2,ham%dmlist(j,i),k)
         end do

      end subroutine dzyaloshinskii_moriya_field

      !---------------symmetric_anisotropic_field---------------!
      !> SA-field
      subroutine symmetric_anisotropic_field(i, k, field,Natom,Mensemble,emomM)
         !
         !.. Implicit declarations
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         integer, intent(in) :: Natom        !< Number of atoms in system
         integer, intent(in) :: Mensemble    !< Number of ensembles
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
         !
         integer :: j, ih

         ! Symmetric anisotropic term
         ih=ham%aHam(i)
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422) && __INTEL_COMPILER < 1800
!         !$omp simd reduction(+:field)
#endif
         do j=1,ham%salistsize(ih)
            field(1) = field(1) + ham%sa_vect(3,j,ih)*emomM(2,ham%salist(j,i),k) +&
               ham%sa_vect(2,j,ih)*emomM(3,ham%salist(j,i),k)
            field(2) = field(2) + ham%sa_vect(1,j,ih)*emomM(3,ham%salist(j,i),k) +&
               ham%sa_vect(3,j,ih)*emomM(1,ham%salist(j,i),k)
            field(3) = field(3) + ham%sa_vect(2,j,ih)*emomM(1,ham%salist(j,i),k) +&
               ham%sa_vect(1,j,ih)*emomM(2,ham%salist(j,i),k)
         end do

      end subroutine symmetric_anisotropic_field

      !---------------chirality_field---------------!
      !> CHIR-field
      subroutine chirality_field(i, k, field,Natom,Mensemble,emomM)
         !
         !.. Implicit declarations
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         integer, intent(in) :: Natom        !< Number of atoms in system
         integer, intent(in) :: Mensemble    !< Number of ensembles
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
         real(dblprec), dimension(3) :: tmp_field
         !
         integer :: j, ih, ip1, ip2, im1, im2

         ! Chirality  term
         tmp_field=field
         ih=ham%aHam(i)
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422) && __INTEL_COMPILER < 1800
         !$omp simd reduction(+:field)
#endif
         do j=1,ham%chirlistsize(ih)

            im1=ham%chirlist(2,j,i)
            ip1=ham%chirlist(1,j,i)
            im2=ham%chirlist(2,j,im1)
            if(im2==0) im2=im1
            ip2=ham%chirlist(1,j,ip1)
            if(ip2==0) ip2=ip1
            !print '(2x,a,2i4,5x,5i4)','->  ', i,j,im2,im1,i,ip1,ip2
            field(1) = field(1)  &
               - ham%chir_coup(j,ih)*emomM(2,ip1,k)*emomM(3,im1,k) + ham%chir_coup(j,ih)*emomM(3,ip1,k)*emomM(2,im1,k) &
               - ham%chir_coup(j,ih)*emomM(2,im2,k)*emomM(3,im1,k) + ham%chir_coup(j,ih)*emomM(3,im2,k)*emomM(2,im1,k) &
               - ham%chir_coup(j,ih)*emomM(2,ip1,k)*emomM(3,ip2,k) + ham%chir_coup(j,ih)*emomM(3,ip1,k)*emomM(2,ip2,k)

            field(2) = field(2)  &
               - ham%chir_coup(j,ih)*emomM(3,ip1,k)*emomM(1,im1,k) + ham%chir_coup(j,ih)*emomM(1,ip1,k)*emomM(3,im1,k) &
               - ham%chir_coup(j,ih)*emomM(3,im2,k)*emomM(1,im1,k) + ham%chir_coup(j,ih)*emomM(1,im2,k)*emomM(3,im1,k) &
               - ham%chir_coup(j,ih)*emomM(3,ip1,k)*emomM(1,ip2,k) + ham%chir_coup(j,ih)*emomM(1,ip1,k)*emomM(3,ip2,k)

            field(3) = field(3)  &
               - ham%chir_coup(j,ih)*emomM(1,ip1,k)*emomM(2,im1,k) + ham%chir_coup(j,ih)*emomM(2,ip1,k)*emomM(1,im1,k) &
               - ham%chir_coup(j,ih)*emomM(1,im2,k)*emomM(2,im1,k) + ham%chir_coup(j,ih)*emomM(2,im2,k)*emomM(1,im1,k) &
               - ham%chir_coup(j,ih)*emomM(1,ip1,k)*emomM(2,ip2,k) + ham%chir_coup(j,ih)*emomM(2,ip1,k)*emomM(1,ip2,k)

            !!! field(1) = field(1) - ham%chir_coup(j,ih)*emomM(2,ham%chirlist(1,j,i),k)*emomM(3,ham%chirlist(2,j,i),k) &
            !!!                     + ham%chir_coup(j,ih)*emomM(3,ham%chirlist(1,j,i),k)*emomM(2,ham%chirlist(2,j,i),k) &
            !!!   - ham%chir_coup(j,ih)*emomM(2,ham%chirlist(2,j,ham%chirlist(2,j,i)),k)*emomM(3,ham%chirlist(2,j,i),k)  &
            !!!   + ham%chir_coup(j,ih)*emomM(3,ham%chirlist(2,j,ham%chirlist(2,j,i)),k)*emomM(2,ham%chirlist(2,j,i),k)  &
            !!!   - ham%chir_coup(j,ih)*emomM(2,ham%chirlist(1,j,i),k)*emomM(3,ham%chirlist(1,j,ham%chirlist(1,j,i)),k)  &
            !!!   + ham%chir_coup(j,ih)*emomM(3,ham%chirlist(1,j,i),k)*emomM(2,ham%chirlist(1,j,ham%chirlist(1,j,i)),k)

            !!! field(2) = field(2) - ham%chir_coup(j,ih)*emomM(3,ham%chirlist(1,j,i),k)*emomM(1,ham%chirlist(2,j,i),k) &
            !!!                     + ham%chir_coup(j,ih)*emomM(1,ham%chirlist(1,j,i),k)*emomM(3,ham%chirlist(2,j,i),k) &
            !!!   - ham%chir_coup(j,ih)*emomM(3,ham%chirlist(2,j,ham%chirlist(2,j,i)),k)*emomM(1,ham%chirlist(2,j,i),k)  &
            !!!   + ham%chir_coup(j,ih)*emomM(1,ham%chirlist(2,j,ham%chirlist(2,j,i)),k)*emomM(3,ham%chirlist(2,j,i),k)  &
            !!!   - ham%chir_coup(j,ih)*emomM(3,ham%chirlist(1,j,i),k)*emomM(1,ham%chirlist(1,j,ham%chirlist(1,j,i)),k)  &
            !!!   + ham%chir_coup(j,ih)*emomM(1,ham%chirlist(1,j,i),k)*emomM(3,ham%chirlist(1,j,ham%chirlist(1,j,i)),k)

            !!! field(3) = field(3) - ham%chir_coup(j,ih)*emomM(1,ham%chirlist(1,j,i),k)*emomM(2,ham%chirlist(2,j,i),k) &
            !!!                     + ham%chir_coup(j,ih)*emomM(2,ham%chirlist(1,j,i),k)*emomM(1,ham%chirlist(2,j,i),k) &
            !!!   - ham%chir_coup(j,ih)*emomM(1,ham%chirlist(2,j,ham%chirlist(2,j,i)),k)*emomM(2,ham%chirlist(2,j,i),k)  &
            !!!   + ham%chir_coup(j,ih)*emomM(2,ham%chirlist(2,j,ham%chirlist(2,j,i)),k)*emomM(1,ham%chirlist(2,j,i),k)  &
            !!!   - ham%chir_coup(j,ih)*emomM(1,ham%chirlist(1,j,i),k)*emomM(2,ham%chirlist(1,j,ham%chirlist(1,j,i)),k)  &
            !!!   + ham%chir_coup(j,ih)*emomM(2,ham%chirlist(1,j,i),k)*emomM(1,ham%chirlist(1,j,ham%chirlist(1,j,i)),k)

         !print '(f12.6)',ham%chir_coup(j,ih)
         !print '(a,i4,3f12.6)','moment i  : ', i,emomM(:,i,k)
         !print '(a,i4,3f12.6)','moment j  : ', ham%chirlist(1,j,i),emomM(:,ham%chirlist(1,j,i),k)
         !print '(a,i4,3f12.6)','moment k  : ', ham%chirlist(2,j,i),emomM(:,ham%chirlist(2,j,i),k)
         end do
         !print '(a,3f12.6)','chir field:         ',field-tmp_field

      end subroutine chirality_field

      !---------------pseudo_dipolar_field---------------!
      !> PD-field
      subroutine pseudo_dipolar_field(i, k, field,Natom,Mensemble,emomM)
         !
         !.. Implicit declarations
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         integer, intent(in) :: Natom        !< Number of atoms in system
         integer, intent(in) :: Mensemble    !< Number of ensembles
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
         !
         integer :: j, ih

         ! Pseudo-Dipolar term
         ih=ham%aHam(i)
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422) && __INTEL_COMPILER < 1800
!         !$omp simd reduction(+:field)
#endif
        ! do j=1,ham%pdlistsize(ih)
        !    field(1) = field(1) + ham%pd_vect(1,j,ih)*emomM(1,ham%pdlist(j,i),k) +&
        !       ham%pd_vect(4,j,ih)*emomM(2,ham%pdlist(j,i),k) +&
        !       ham%pd_vect(5,j,ih)*emomM(3,ham%pdlist(j,i),k)
        !    field(2) = field(2) + ham%pd_vect(4,j,ih)*emomM(1,ham%pdlist(j,i),k) +&
        !       ham%pd_vect(2,j,ih)*emomM(2,ham%pdlist(j,i),k) +&
        !       ham%pd_vect(6,j,ih)*emomM(3,ham%pdlist(j,i),k)
        !    field(3) = field(3) + ham%pd_vect(5,j,ih)*emomM(1,ham%pdlist(j,i),k) +&
        !       ham%pd_vect(6,j,ih)*emomM(2,ham%pdlist(j,i),k) +&
        !       ham%pd_vect(3,j,ih)*emomM(3,ham%pdlist(j,i),k)
        ! end do
          do j=1,ham%pdlistsize(ih)
            field(1) = field(1) + ham%pd_vect(1,j,ih)*emomM(1,ham%pdlist(j,i),k) +&
               ham%pd_vect(2,j,ih)*emomM(2,ham%pdlist(j,i),k) +&
               ham%pd_vect(3,j,ih)*emomM(3,ham%pdlist(j,i),k)
            field(2) = field(2) + ham%pd_vect(4,j,ih)*emomM(1,ham%pdlist(j,i),k) +&
               ham%pd_vect(5,j,ih)*emomM(2,ham%pdlist(j,i),k) +&
               ham%pd_vect(6,j,ih)*emomM(3,ham%pdlist(j,i),k)
            field(3) = field(3) + ham%pd_vect(7,j,ih)*emomM(1,ham%pdlist(j,i),k) +&
               ham%pd_vect(8,j,ih)*emomM(2,ham%pdlist(j,i),k) +&
               ham%pd_vect(9,j,ih)*emomM(3,ham%pdlist(j,i),k)
         end do
      end subroutine pseudo_dipolar_field

      !---------------dzyaloshinskii_moriya_bq_field---------------!
      !> DM-BQ field
      subroutine dzyaloshinskii_moriya_bq_field(i, k, field,Natom,Mensemble,emomM)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         integer, intent(in) :: Natom        !< Number of atoms in system
         integer, intent(in) :: Mensemble    !< Number of ensembles
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
         !
         integer :: j, ih
         real(dblprec), dimension(3) :: dot !< Work array

         ! BIQDM term
         ih=ham%aHam(i)
         do j=1,ham%biqdmlistsize(ih)
            dot(1) = emomM(1,i,k)*emomM(2,ham%biqdmlist(j,i),k)-&
               emomM(2,i,k)*emomM(1,ham%biqdmlist(j,i),k)
            dot(2) = emomM(2,i,k)*emomM(3,ham%biqdmlist(j,i),k)-&
               emomM(3,i,k)*emomM(2,ham%biqdmlist(j,i),k)
            dot(3) = emomM(3,i,k)*emomM(1,ham%biqdmlist(j,i),k)-&
               emomM(1,i,k)*emomM(3,ham%biqdmlist(j,i),k)
            field(1) = field(1) + 2.0_dblprec*ham%biqdm_vect(1,j,ih)*(&
               dot(1)*emomM(3,ham%biqdmlist(j,i),k)-&
               dot(2)*emomM(2,ham%biqdmlist(j,i),k))
            field(2) = field(2) + 2.0_dblprec*ham%biqdm_vect(1,j,ih)*(&
               dot(2)*emomM(1,ham%biqdmlist(j,i),k)-&
               dot(3)*emomM(3,ham%biqdmlist(j,i),k))
            field(3) = field(3) + 2.0_dblprec*ham%biqdm_vect(1,j,ih)*(&
               dot(3)*emomM(2,ham%biqdmlist(j,i),k)-&
               dot(1)*emomM(1,ham%biqdmlist(j,i),k))
         end do
      end subroutine dzyaloshinskii_moriya_bq_field


      !---------------biquadratic_field---------------!
      !> BQ-field
      subroutine biquadratic_field(i, k, field,Natom,Mensemble,emomM)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         integer, intent(in) :: Natom        !< Number of atoms in system
         integer, intent(in) :: Mensemble    !< Number of ensembles
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
         !
         integer :: j, ih
         real(dblprec) :: dot

         ! Biquadratic exchange term
         ih=ham%aHam(i)
         do j=1,ham%bqlistsize(ih)
            dot=emomM(1,ham%bqlist(j,i),k)*emomM(1,i,k)+&
               emomM(2,ham%bqlist(j,i),k)*emomM(2,i,k)+&
               emomM(3,ham%bqlist(j,i),k)*emomM(3,i,k)
            field = field + 2.0_dblprec*ham%j_bq(j,ih)*dot*emomM(1:3,ham%bqlist(j,i),k)
         end do
      end subroutine biquadratic_field
      
	  !---------------------ring_field--------------------!
      subroutine ring_field(i, k, field,Natom,Mensemble,emomM)
	  !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: i !< Atom to calculate effective field for
      integer, intent(in) :: k !< Current ensemble
      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Mensemble    !< Number of ensembles
      real(dblprec), dimension(3), intent(inout) :: field !< Effective field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      !
      integer :: j
      real(dblprec) :: dotkl,dotkj,dotjl
      real(dblprec), dimension(3) :: tmpfield

      tmpfield=field
       ! Four-spin ring exchange term 
          do j=1,ham%ringlistsize(i)
             dotkl=emomM(1,ham%ringlist(i,j,2),k)*emomM(1,ham%ringlist(i,j,3),k)+&
                  emomM(2,ham%ringlist(i,j,2),k)*emomM(2,ham%ringlist(i,j,3),k)+&
                  emomM(3,ham%ringlist(i,j,2),k)*emomM(3,ham%ringlist(i,j,3),k)

             dotkj=emomM(1,ham%ringlist(i,j,2),k)*emomM(1,ham%ringlist(i,j,1),k)+&
                  emomM(2,ham%ringlist(i,j,2),k)*emomM(2,ham%ringlist(i,j,1),k)+&
                  emomM(3,ham%ringlist(i,j,2),k)*emomM(3,ham%ringlist(i,j,1),k)

             dotjl=emomM(1,ham%ringlist(i,j,1),k)*emomM(1,ham%ringlist(i,j,3),k)+&
                  emomM(2,ham%ringlist(i,j,1),k)*emomM(2,ham%ringlist(i,j,3),k)+&
                  emomM(3,ham%ringlist(i,j,1),k)*emomM(3,ham%ringlist(i,j,3),k)

             field = field - ham%j_ring(i,j)*dotkl*emomM(1:3,ham%ringlist(i,j,1),k)-&
             ham%j_ring(i,j)*dotkj*emomM(1:3,ham%ringlist(i,j,3),k)+&
             ham%j_ring(i,j)*dotjl*emomM(1:3,ham%ringlist(i,j,2),k)

          end do
    end subroutine ring_field
      
      

      !---------------uniaxial_anisotropy---------------------------------------
      !> @brief Field from the uniaxial anisotropy
      !-------------------------------------------------------------------------
      subroutine uniaxial_anisotropy_field(i, k, field,Natom,Mensemble,mult_axis,emomM)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         integer, intent(in) :: Natom        !< Number of atoms in system
         integer, intent(in) :: Mensemble    !< Number of ensembles
         character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
         !
         real(dblprec) :: tt1,tt2,tt3
         real(dblprec) :: tt1_d,tt2_d,tt3_d

         ! Uniaxial anisotropy
         ! cos(theta)
         tt1=emomM(1,i,k)*ham%eaniso(1,i)+emomM(2,i,k)*ham%eaniso(2,i)+emomM(3,i,k)*ham%eaniso(3,i)
         ! k1 + 2*k2*sin^2(theta) = k1 + 2*k2*(1-cos^2(theta))
         tt2=ham%kaniso(1,i)+2.0_dblprec*ham%kaniso(2,i)*(1.0_dblprec-tt1*tt1)
         ! 2 * cos(theta)* [k1 + 2*k2*sin^2(theta)]
         tt3= 2.0_dblprec*tt1*tt2

         if (mult_axis=='Y') then
            ! Uniaxial anisotropy
            ! cos(theta)
            tt1_d=emomM(1,i,k)*ham%eaniso_diff(1,i)+emomM(2,i,k)*ham%eaniso_diff(2,i)+emomM(3,i,k)*ham%eaniso_diff(3,i)
            ! k1 + 2*k2*sin^2(theta) = k1 + 2*k2*(1-cos^2(theta))
            tt2_d=ham%kaniso_diff(1,i)+2.0_dblprec*ham%kaniso_diff(2,i)*(1.0_dblprec-tt1_d*tt1_d)
            ! 2 * cos(theta)* [k1 + 2*k2*sin^2(theta)]
            tt3_d= 2.0_dblprec*tt1_d*tt2_d

            field  = field - tt3*ham%eaniso(1:3,i)-tt3_d*ham%eaniso_diff(1:3,i)
         else

            field  = field - tt3*ham%eaniso(1:3,i)
         endif

      end subroutine uniaxial_anisotropy_field

      !---------------cubic_anisotropy_field---------------!
      !> Cubic anisotropy
      subroutine cubic_anisotropy_field(i, k, field,Natom,Mensemble,mult_axis,emomM)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         integer, intent(in) :: Natom        !< Number of atoms in system
         integer, intent(in) :: Mensemble    !< Number of ensembles
         character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
         !

         field(1) = field(1)  &
            + 2.0_dblprec*ham%kaniso(1,i)*emomM(1,i,k)*(emomM(2,i,k)**2+emomM(3,i,k)**2) &
            + 2.0_dblprec*ham%kaniso(2,i)*emomM(1,i,k)*(emomM(2,i,k)**2*emomM(3,i,k)**2)
         field(2) = field(2)  &
            + 2.0_dblprec*ham%kaniso(1,i)*emomM(2,i,k)*(emomM(3,i,k)**2+emomM(1,i,k)**2) &
            + 2.0_dblprec*ham%kaniso(2,i)*emomM(2,i,k)*(emomM(3,i,k)**2*emomM(1,i,k)**2)
         field(3) = field(3)  &
            + 2.0_dblprec*ham%kaniso(1,i)*emomM(3,i,k)*(emomM(1,i,k)**2+emomM(2,i,k)**2) &
            + 2.0_dblprec*ham%kaniso(2,i)*emomM(3,i,k)*(emomM(1,i,k)**2*emomM(2,i,k)**2)

         if (mult_axis=='Y') then

            field(1) = field(1)  &
               + 2.0_dblprec*ham%kaniso_diff(1,i)*emomM(1,i,k)*(emomM(2,i,k)**2+emomM(3,i,k)**2) &
               + 2.0_dblprec*ham%kaniso_diff(2,i)*emomM(1,i,k)*emomM(2,i,k)**2*emomM(3,i,k)**2
            field(2) = field(2)  &
               + 2.0_dblprec*ham%kaniso_diff(1,i)*emomM(2,i,k)*(emomM(3,i,k)**2+emomM(1,i,k)**2) &
               + 2.0_dblprec*ham%kaniso_diff(2,i)*emomM(2,i,k)*emomM(3,i,k)**2*emomM(1,i,k)**2
            field(3) = field(3)  &
               + 2.0_dblprec*ham%kaniso_diff(1,i)*emomM(3,i,k)*(emomM(1,i,k)**2+emomM(2,i,k)**2) &
               + 2.0_dblprec*ham%kaniso_diff(2,i)*emomM(3,i,k)*emomM(1,i,k)**2*emomM(2,i,k)**2

         endif

      end subroutine cubic_anisotropy_field

end module HamiltonianActions
