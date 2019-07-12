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
!> @author Anders Bergman, Lars Bergqvist, Johan Hellsvik
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module HamiltonianActions_lite

   use Profiling
   use Parameters
   use HamiltonianData

   implicit none

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
   subroutine effective_field_lite(Natom,Mensemble,start_atom,stop_atom,do_jtensor,      &
      do_anisotropy,exc_inter,do_dm,do_pd,do_biqdm,do_bq,do_chir,do_dip,emomM,mmom, &
      external_field,time_external_field,beff,beff1,beff2,mult_axis,energy,NA,N1,N2,N3)
      !
      use Constants, only : mry,mub
      use DipoleManager, only : dipole_field_calculation
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
      integer, intent(in) :: do_jtensor   !< Use SKKR style exchange tensor (0=off, 1=on, 2=with biquadratic exchange)
      integer, intent(in) :: do_anisotropy !< Add on-site anisotropy terms to Hamiltonian
      integer, intent(in) :: do_dm        !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_pd        !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
      integer, intent(in) :: do_biqdm     !< Add Biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_bq        !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      integer, intent(in) :: do_chir      !< Add scalar chirality exchane (CHIR) term to Hamiltonian (0/1)
      integer, intent(in) :: do_dip       !< Calculate dipole-dipole contribution (0=Off, 1= Brute Force, 2=macrocell)
      !integer, intent(in) :: Num_macro    !< Number of macrocells in the system
      character(len=1), intent(in) :: exc_inter !< Interpolation of Jij (Y/N)
      character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
      !integer, dimension(Natom), intent(in) :: cell_index            !< Macrocell index for each atom
      !integer, dimension(Num_macro), intent(in) :: macro_nlistsize   !< Number of atoms per macrocell
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      !real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: external_field  !< External magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: time_external_field !< External time-dependent magnetic field
      ! .. Output Variables
      real(dblprec), intent(out) :: energy !< Total energy
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff  !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff1 !< Internal effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff2 !< External field from application of Hamiltonian

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!! +++ Optimization Routines Variables +++ !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!integer, intent(in) :: max_no_constellations !< The maximum (global) length of the constellation matrix
      !!!logical, intent(in) :: OPT_flag
      !!!integer, dimension(:), intent(in) :: maxNoConstl
      !!!integer, dimension(:,:), intent(in) :: unitCellType !< Array of constellation id and classification (core, boundary, or noise) per atom
      !!!integer, dimension(:,:,:), intent(in) :: constellationsNeighType
      !!!real(dblprec), dimension(:,:,:), intent(in) :: constlNCoup
      !!!real(dblprec), dimension(:,:,:), intent(in) :: constellations
      !!!real(dblprec), dimension(3,max_no_constellations,Mensemble) :: beff1_constellations
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!! +++ End Region +++ !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !.. Local scalars
      integer :: i,k
      real(dblprec), dimension(3) :: tfield, beff_s, beff_q, beff_m

      ! Initialization of the energy
      energy=0.0_dblprec
      ! Initialization if the effective field
      beff=0.0_dblprec
      
      !!reduction(+:energy)
      !!$omp parallel do default(shared) schedule(static) private(i,k,beff_s,beff_q,tfield,beff_m) collapse(2) reduction(+:energy)
      do k=1, Mensemble
         do i=start_atom, stop_atom

            beff_s=0.0_dblprec
            beff_q=0.0_dblprec
            beff_m=0.0_dblprec

             call heisenberg_field_lite(i, k, beff_s,Natom,Mensemble,emomM)
            !!! if(do_jtensor/=1) then
            !!!    ! Heisenberg exchange term
            !!!    if(exc_inter=='N') then
            !!!       call heisenberg_field(i, k, beff_s,Natom,Mensemble,OPT_flag,      &
            !!!          beff1_constellations,unitCellType,emomM,max_no_constellations)
            !!!    else
            !!!       call heisenberg_rescaling_field(i, k, beff_s,Natom,Mensemble,mmom,emomM)
            !!!    endif
            !!!    ! Dzyaloshinskii-Moriya term
            !!!    if(do_dm==1) call dzyaloshinskii_moriya_field(i, k, beff_s,Natom,Mensemble,emomM)
            !!! else
            !!!    call tensor_field(i, k, beff_s,Natom,Mensemble,emomM)
            !!! end if

            beff1(1:3,i,k)= beff_s
            beff2(1:3,i,k)= beff_q+external_field(1:3,i,k)!+time_external_field(1:3,i,k)
            beff(1:3,i,k) = beff(1:3,i,k)+ beff1(1:3,i,k)+beff2(1:3,i,k)

            tfield=0.50_dblprec*(beff_s+2.0_dblprec*beff_q+external_field(1:3,i,k))!+time_external_field(1:3,i,k))
            energy=energy - emomM(1,i,k)*tfield(1)-emomM(2,i,k)*tfield(2)-emomM(3,i,k)*tfield(3)
         end do
      end do
      !!$omp end parallel do
      energy = energy * mub / mry

   end subroutine effective_field_lite

   subroutine effective_field_extralite(Natom,Mensemble,start_atom,stop_atom,emomM,mmom,energy,beff)
      !
      use Constants, only : mry,mub
      !use DipoleManager, only : dipole_field_calculation
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Mensemble    !< Number of ensembles
      integer, intent(in) :: start_atom   !< Atom to start loop for
      integer, intent(in) :: stop_atom    !< Atom to end loop for
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment
      !real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: external_field  !< External magnetic field
      !real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: time_external_field !< External time-dependent magnetic field
      ! .. Output Variables
      real(dblprec), intent(out) :: energy !< Total energy
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff  !< Total effective field from application of Hamiltonian

      !.. Local scalars
      integer :: i,k
      real(dblprec), dimension(3) :: tfield, beff_s, beff_q, beff_m

      ! Initialization of the energy
      energy=0.0_dblprec
      ! Initialization if the effective field
      beff=0.0_dblprec
      
      !!reduction(+:energy)
      !!$omp parallel do default(shared) schedule(static) private(i,k,beff_s,beff_q,tfield,beff_m) collapse(2) reduction(+:energy)
      do k=1, Mensemble
         do i=start_atom, stop_atom

            beff_s=0.0_dblprec
            beff_q=0.0_dblprec
            beff_m=0.0_dblprec

             call heisenberg_field_lite(i, k, beff_s,Natom,Mensemble,emomM)
            !!! if(do_jtensor/=1) then
            !!!    ! Heisenberg exchange term
            !!!    if(exc_inter=='N') then
            !!!       call heisenberg_field(i, k, beff_s,Natom,Mensemble,OPT_flag,      &
            !!!          beff1_constellations,unitCellType,emomM,max_no_constellations)
            !!!    else
            !!!       call heisenberg_rescaling_field(i, k, beff_s,Natom,Mensemble,mmom,emomM)
            !!!    endif
            !!!    ! Dzyaloshinskii-Moriya term
            !!!    if(do_dm==1) call dzyaloshinskii_moriya_field(i, k, beff_s,Natom,Mensemble,emomM)
            !!! else
            !!!    call tensor_field(i, k, beff_s,Natom,Mensemble,emomM)
            !!! end if

            beff(1:3,i,k) = beff(1:3,i,k)+beff_s+beff_q

            tfield=0.50_dblprec*(beff_s+2.0_dblprec*beff_q)!+external_field(1:3,i,k))!+time_external_field(1:3,i,k))
            energy=energy - emomM(1,i,k)*tfield(1)-emomM(2,i,k)*tfield(2)-emomM(3,i,k)*tfield(3)
         end do
      end do
      !!$omp end parallel do
      energy = energy * mub / mry

   end subroutine effective_field_extralite
   !contains

      !---------------heisenberg_field---------------!
      !> Heisenberg
      subroutine heisenberg_field_lite(i, k, field,Natom,Mensemble,emomM)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         integer, intent(in) :: Natom        !< Number of atoms in system
         integer, intent(in) :: Mensemble    !< Number of ensembles
         real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field

         integer :: j, ih
    
         ih=ham%aHam(i)

#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422) && __INTEL_COMPILER < 1800
!            !$omp simd reduction(+:field)
#endif
            do j=1,ham%nlistsize(ih)
               !           !DIR$ vector always aligned
               field = field + ham%ncoup(j,ih,1)*emomM(:,ham%nlist(j,i),k)
            end do
      end subroutine heisenberg_field_lite

end module HamiltonianActions_lite
