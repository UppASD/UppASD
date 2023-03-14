!-------------------------------------------------------------------------------
! MODULE: HamiltonianActions_cpu
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
module HamiltonianActions_cpu

   use Profiling
   use Parameters
   use HamiltonianData
   use InputData, only : ham_inp

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
   subroutine effective_field_cpu(Natom,Mensemble,start_atom,stop_atom,   &
      emomM,mmom,external_field,time_external_field,beff,beff1,beff2,OPT_flag,      &
      max_no_constellations,maxNoConstl,unitCellType,constlNCoup,constellations,    &
      constellationsNeighType,energy,Num_macro,cell_index,emomM_macro,    &
      macro_nlistsize,NA,N1,N2,N3)
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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! +++ Optimization Routines Variables +++ !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer, intent(in) :: max_no_constellations !< The maximum (global) length of the constellation matrix
      logical, intent(in) :: OPT_flag
      integer, dimension(:), intent(in) :: maxNoConstl
      integer, dimension(:,:), intent(in) :: unitCellType !< Array of constellation id and classification (core, boundary, or noise) per atom
      integer, dimension(:,:,:), intent(in) :: constellationsNeighType
      real(dblprec), dimension(:,:,:), intent(in) :: constlNCoup
      real(dblprec), dimension(:,:,:), intent(in) :: constellations
      real(dblprec), dimension(3,max_no_constellations,Mensemble) :: beff1_constellations
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! +++ End Region +++ !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !.. Local scalars
      integer :: i,k
      real(dblprec), dimension(3) :: tfield, beff_s, beff_q, beff_m

      !.. Executable statements
      if(OPT_flag) call pre_optimize_cpu(Natom,Mensemble,max_no_constellations,         &
         maxNoConstl,constellations,constlNCoup,beff1_constellations,               &
         constellationsNeighType)

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
            emomM_macro,ham%Qdip,ham%Qdip_macro,energy,beff)
         call timing(0,'Dipolar Int.  ','OF')
         call timing(0,'Hamiltonian   ','ON')
      endif
      !reduction(+:energy)
      !$omp parallel do default(shared) schedule(static) private(i,k,beff_s,beff_q,tfield,beff_m) collapse(2) reduction(+:energy)
      do k=1, Mensemble
         do i=start_atom, stop_atom

            beff_s=0.0_dblprec
            beff_q=0.0_dblprec
            beff_m=0.0_dblprec

            if(ham_inp%do_jtensor/=1) then
               ! Heisenberg exchange term
               call heisenberg_field_cpu(i, k, beff_s,Natom,Mensemble,OPT_flag,      &
                    beff1_constellations,unitCellType,emomM,max_no_constellations)
            else
               call tensor_field_cpu(i, k, beff_s,Natom,Mensemble,emomM)
            end if
            beff1(1:3,i,k)= beff_s
            beff2(1:3,i,k)= beff_q+external_field(1:3,i,k)+time_external_field(1:3,i,k)
            beff(1:3,i,k) = beff(1:3,i,k)+ beff1(1:3,i,k)+beff2(1:3,i,k)

            tfield=0.50_dblprec*(beff_s+2.0_dblprec*beff_q+2.0_dblprec*external_field(1:3,i,k)+time_external_field(1:3,i,k))
            energy=energy - emomM(1,i,k)*tfield(1)-emomM(2,i,k)*tfield(2)-emomM(3,i,k)*tfield(3)
         end do
      end do
      !$omp end parallel do
      energy = energy * mub / mry

   end subroutine effective_field_cpu

   !contains

      !---------------pre_optimize---------------!
      !> Preoptimization
      subroutine pre_optimize_cpu(Natom,Mensemble,max_no_constellations,maxNoConstl,    &
         constellations,constlNCoup,beff1_constellations,constellationsNeighType)

         implicit none

         integer, intent(in) :: Natom        !< Number of atoms in system
         integer, intent(in) :: Mensemble    !< Number of ensembles
         integer, intent(in) :: max_no_constellations ! The maximum (global) length of the constellation matrix
         integer, dimension(:), intent(in) :: maxNoConstl
         integer, dimension(:,:,:), intent(in) :: constellationsNeighType
         real(dblprec), dimension(:,:,:), intent(in) :: constellations
         real(dblprec), dimension(:,:,:), intent(in) :: constlNCoup
         real(dblprec), dimension(3,max_no_constellations,Mensemble), intent(inout) :: beff1_constellations

         integer :: k,i,j

         beff1_constellations(:,:,:) = 0.0_dblprec
         do k=1,Mensemble
            ! Compute the Heissenberg exchange term for the (representative) unit cells in each constellation
            do i=1,maxNoConstl(k)
               ! Computation of the exchange term, going to ham%max_no_neigh since the exchange factor,
               ! governed by constlNCoup, is zero for iterations beyond the neighbourhood of atom i
               do j=1,ham%max_no_neigh
                  beff1_constellations(:,i,k) = beff1_constellations(:,i,k) + &
                     constlNCoup(j,i,k)*constellations(:,constellationsNeighType(j,i,k),k)
               end do
            end do
         end do

      end subroutine pre_optimize_cpu

      !---------------heisenberg_field---------------!
      !> Heisenberg
      subroutine heisenberg_field_cpu(i, k, field,Natom,Mensemble,OPT_flag,             &
         beff1_constellations,unitCellType,emomM,max_no_constellations)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         integer, intent(in) :: Natom        !< Number of atoms in system
         integer, intent(in) :: Mensemble    !< Number of ensembles
         logical, intent(in) :: OPT_flag
         integer, intent(in) :: max_no_constellations ! The maximum (global) length of the constellation matrix
         integer, dimension(:,:), intent(in) :: unitCellType ! Array of constellation id and classification (core, boundary, or noise) per atom
         real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
         real(dblprec), dimension(3,max_no_constellations,Mensemble), intent(inout) :: beff1_constellations
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field

         integer :: j, ih

         ! WARNING LSF is not working for ASD ham%ncoup set to first configuration only ASK FAN
         !Exchange term, run optimization or use the standard scheme
         if(OPT_flag) then
            do j=1,ham%nlistsize(i)
               if(unitCellType(i,k).ge.1) then ! If in a constellation
                  field = field+beff1_constellations(:,unitCellType(i,k),k)
               endif
            end do
         else
            ih=ham%aHam(i)
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422) && __INTEL_COMPILER < 1800
!            !$omp simd reduction(+:field)
#endif
            do j=1,ham%nlistsize(ih)
               !           !DIR$ vector always aligned
               field = field + ham%ncoup(j,ih,1)*emomM(:,ham%nlist(j,i),k)
            end do
         endif
      end subroutine heisenberg_field_cpu

      !---------------tensor_heisenberg_field---------------!
      !> @brief Calculates heisenberg, DM and anisotropy through one tensor
      !>
      !> @Date 09/15/2014 - Thomas Nystrand
      !> - Added 0.5*m_i*I part
      subroutine tensor_field_cpu(i, k, field,Natom,Mensemble,emomM)
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
            field = field + ham%j_tens(:,1,j,i)*emomM(1,x,k) + ham%j_tens(:,2,j,i)*emomM(2,x,k) + ham%j_tens(:,3,j,i)*emomM(3,x,k)
         end do
      end subroutine tensor_field_cpu

end module HamiltonianActions_cpu
