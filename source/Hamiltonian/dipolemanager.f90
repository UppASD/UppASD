!-------------------------------------------------------------------------------
! MODULE: DipoleManager
!> @brief Wrapper module containing the needed subroutines for the dipole-dipole
!> interaction.
!> @details This module handles the caller routines for allocation/deallocation,
!> intialization and actual calculation of the dipole-dipole interactions for the
!> different methods implemented.
!> @author Jonathan Chico
!-------------------------------------------------------------------------------
module DipoleManager

   use Profiling
   use Parameters
#ifdef USE_FFTW
   use fftdipole_fftw
#endif
#ifdef USE_MKL_FFT
   use fftdipole_mkl
#endif

   implicit none

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: dipole_setup
   !> @brief Wrapper routine for the setup of the dipole-dipole interaction.
   !> @details This routine takes care of the setup of the dipole-dipole geometrical
   !> tensor for the brute force methods. For the FFT methods it takes care of the
   !> calculation of the geometrical part of the dipole-dipole interaction
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine dipole_setup(NA,N1,N2,N3,Natom,do_dip,Num_macro,Mensemble,            &
      max_num_atom_macro_cell,macro_nlistsize,macro_atom_nlist,alat,C1,C2,C3,Bas,   &
      coord,simid,read_dipole,print_dip_tensor,qdip_files,Qdip,Qdip_macro)

      use HamiltonianData, only: allocate_dipole,allocate_macro_dipole
      use DipoleCommon, only : setup_qdip,setup_macro_qdip,read_qdip

      implicit none

      integer, intent(in) :: NA
      integer, intent(in) :: N1
      integer, intent(in) :: N2
      integer, intent(in) :: N3
      integer, intent(in) :: Natom
      integer, intent(in) :: do_dip
      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, intent(in) :: Mensemble
      integer, intent(in) :: max_num_atom_macro_cell  !< Maximum number of atoms in  a macrocell
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize   !< Number of atoms per macrocell
      integer, dimension(Num_macro,max_num_atom_macro_cell), intent(in) :: macro_atom_nlist  !< List containing the information of which atoms are in a given macrocell
      real(dblprec), intent(in) :: alat
      real(dblprec), dimension(3), intent(in) :: C1
      real(dblprec), dimension(3), intent(in) :: C2
      real(dblprec), dimension(3), intent(in) :: C3
      real(dblprec), dimension(3,NA), intent(in) :: Bas
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      character(len=8), intent(in) :: simid
      character(len=1), intent(in) :: read_dipole !< Flag to read the dipole-dipole tensor
      character(len=1), intent(in) :: print_dip_tensor
      character(len=30), intent(in) :: qdip_files !< Input file that contains the dipole-dipole tensor
      real(dblprec), dimension(:,:,:,:), allocatable, intent(inout) :: Qdip
      real(dblprec), dimension(:,:,:,:), allocatable, intent(inout):: Qdip_macro !< Matrix for the macro spin dipole-dipole interaction

      !  Set up dipole-tensor for the brute force calculation
      if(do_dip==1) then
         write(*,'(2x,a)',advance='no') "Set up dipole-dipole matrix"
         call allocate_dipole(Natom,1)
         if (read_dipole=='Y') then
            call read_qdip(Natom,Qdip,qdip_files)
         else
            call setup_qdip(Natom, coord, alat, Qdip,simid,print_dip_tensor)
         endif
         write(*,'(a)') '  done'
      ! Setup dipole-tensor for the macro-cell calculation
      else if (do_dip==2) then
         write(*,'(2x,a)',advance='no') "Set up macro dipole-dipole matrix"
         call allocate_macro_dipole(Num_macro,1)
         if (read_dipole=='Y') then
            call read_qdip(Num_macro,Qdip_macro,qdip_files)
         else
            call setup_macro_qdip(Natom,Num_macro,max_num_atom_macro_cell,          &
               macro_nlistsize,macro_atom_nlist,coord,alat,Qdip_macro,simid,        &
               print_dip_tensor)
         endif
         write(*,'(a)') '  done'
      ! Setup the calculation of the dipole-dipole FFT calculation, also calculate
      ! the FFT of the geometrical part of the dipole interaction
      else if (do_dip==3) then
#ifdef USE_FFTW
         write(*,'(2x,a)',advance='no') "Set up FFT dipole-dipole matrix"
         call geometrical_DFFT_wrapper(NA,N1,N2,N3,Natom,Mensemble,alat,C1,C2,C3,Bas)
         write(*,'(a)') '  done'
#endif
#ifdef USE_MKL_FFT
         write(*,'(2x,a)',advance='no') "Set up MKL FFT dipole-dipole matrix"
         call geometrical_DFFT_wrapper_MKL(NA,N1,N2,N3,Natom,Mensemble,alat,C1,C2,  &
            C3,Bas)
         write(*,'(a)') '  done'
#endif
#if (! defined USE_FFTW && ! defined USE_MKL_FFT)
         write(*,'(2x,a)') 'UppASD was not compiled with the FFTW library or MKL support, needed for this mode. Shutting down'
         stop
#endif
      end if

   end subroutine dipole_setup

   !----------------------------------------------------------------------------
   ! SUBROUTINE: dipole_field_calculation
   !> @brief Wrapper for the calculation of the field steming from the dipole-dipole interaction
   !> @details Calculation of the dipole-dipole field for the different approaches,
   !> for the brute force terms it calculates the double loop of the interaction.
   !> For the FFT methods it performs the FFT of the magnetic moments, calculates
   !> the field in reciprocal space and then it transforms it back to real space.
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine dipole_field_calculation(NA,N1,N2,N3,Natom,do_dip,Num_macro,          &
      Mensemble,stop_atom,start_atom,cell_index,macro_nlistsize,emomM,emomM_macro,  &
      Qdip,Qdip_macro,energy,bfield)

      implicit none

      integer, intent(in) :: NA
      integer, intent(in) :: N1
      integer, intent(in) :: N2
      integer, intent(in) :: N3
      integer, intent(in) :: Natom
      integer, intent(in) :: do_dip
      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, intent(in) :: Mensemble
      integer, intent(in) :: stop_atom
      integer, intent(in) :: start_atom
      integer, dimension(Natom), intent(in) :: cell_index            !< Macrocell index for each atom
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize   !< Number of atoms per macrocell
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
      real(dblprec), dimension(3,3,Natom,Natom), intent(in) :: Qdip
      real(dblprec), dimension(3,3,Num_macro,Num_macro), intent(in):: Qdip_macro !< Matrix for the macro spin dipole-dipole interaction
      real(dblprec), intent(inout) :: energy
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: bfield

      integer :: ii,kk

      ! Brute force calculation of the dipole-dipole field
      if (do_dip==1) then
         !$omp parallel do default(shared) schedule(static) private(ii,kk) collapse(2) reduction(+:energy)
         do kk=1, Mensemble
            do ii=start_atom, stop_atom
               ! Dipolar term
               call dipolar_field(ii,kk,bfield(:,ii,kk),Natom,Mensemble,emomM,Qdip)
               energy=energy-0.50_dblprec*(bfield(1,ii,kk)*emomM(1,ii,kk)+          &
                  bfield(2,ii,kk)*emomM(2,ii,kk)+bfield(3,ii,kk)*emomM(3,ii,kk))
            enddo
         enddo
         !$omp end parallel do
      ! Macro cell calculation of the dipole-dipole field
      else if (do_dip==2) then

         !$omp parallel do default(shared) schedule(static) private(ii,kk,energy) collapse(2)
         do kk=1, Mensemble
            do ii=start_atom, stop_atom
               ! Macro cell dipolar field
               call macro_dipolar_field(ii,kk,bfield(:,ii,kk),Natom,Num_macro,      &
                  Mensemble,cell_index,emomM_macro,Qdip_macro)
               ! Calculate the contribution to the energy comming from the macrocell treatment of the
               ! dipole-dipole interaction
               call calc_macro_energy(ii,kk,bfield(:,ii,kk),energy,Natom,Num_macro,Mensemble, &
                  cell_index,emomM_macro,macro_nlistsize)
            enddo
         enddo
         !$omp end parallel do
      ! Calculation of the dipole-dipole field making use of the FFT and the convolution theorem
      else if (do_dip==3) then
#ifdef USE_FFTW
         call calculate_fft_dipolar_field(NA,N1,N2,N3,Natom,Mensemble,emomM,bfield)
#endif
#ifdef USE_MKL_FFT
         call calculate_fft_dipolar_field_MKL(NA,N1,N2,N3,Natom,Mensemble,emomM,bfield)
#endif
      endif

   end subroutine dipole_field_calculation

   !-------------------------------------------------------------------------
   ! SUBROUTINE: calc_macro_energy
   !
   !> @brief Calculation of the macro spin dipole-dipole energy contribution
   !> @author Jonathan Chico
   !> @todo Right now the contribution is done for every spin, maybe there is a
   !> way to be more efficient
   !-------------------------------------------------------------------------
   subroutine calc_macro_energy(i,k,field,energy,Natom,Num_macro,Mensemble,      &
      cell_index,emomM_macro,macro_nlistsize)

      implicit none

      integer, intent(in) :: i !< Atom to calculate effective field for
      integer, intent(in) :: k !< Current ensemble
      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Num_macro    !< Number of macrocells in the system
      integer, intent(in) :: Mensemble    !< Number of ensembles
      real(dblprec), dimension(3), intent(in) :: field
      integer, dimension(Natom), intent(in) :: cell_index            !< Macrocell index for each atom
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize   !< Number of atoms per macrocell
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
      real(dblprec), intent(inout) :: energy

      !.. Local variables
      integer :: icell

      icell=cell_index(i)
      energy=energy-0.5_dblprec*(field(1)*emomM_macro(1,icell,k)+field(2)*emomM_macro(2,icell,k)+&
         field(3)*emomM_macro(3,icell,k))/macro_nlistsize(icell)

   end subroutine calc_macro_energy

   !----------------------------------------------------------------------------
   ! SUBROUTINE: dipolar_field
   !> @brief Brute force calculation of the dipolar field.
   !----------------------------------------------------------------------------
   subroutine dipolar_field(i, k, field,Natom,Mensemble,emomM,Qdip)

      implicit none

      integer, intent(in) :: i !< Atom to calculate effective field for
      integer, intent(in) :: k !< Current ensemble
      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Mensemble    !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,3,Natom,Natom), intent(in) :: Qdip
      real(dblprec), dimension(3), intent(inout) :: field !< Effective field
      !
      integer :: j

#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422) && __INTEL_COMPILER < 1800
!         !$omp simd reduction(+:field)
#endif
      do j=1,Natom
         field(:) = field(:) + Qdip(1,:,j,i)*emomM(1,j,k) + Qdip(2,:,j,i)*emomM(2,j,k) + Qdip(3,:,j,i)*emomM(3,j,k)
      end do
   end subroutine dipolar_field

   !----------------------------------------------------------------------------
   ! SUBROUTINE: macro_dipolar_field
   !> @brief Calculation of the macro spin dipolar field.
   !> @details The field inside each macro cell is kept constant and it is produced
   !> by all the other macro cells a self-energy (self-demagnetizing field)
   !> is included to take into account the effect of the spins inside each macro cell.
   !> The approach is based in J. Phys.: Condens. Matter 28 (2016) 066001
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine macro_dipolar_field(i,k,field,Natom,Num_macro,Mensemble,cell_index,&
      emomM_macro,Qdip_macro)

      implicit none

      integer, intent(in) :: i !< Atom to calculate effective field for
      integer, intent(in) :: k !< Current ensemble
      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Num_macro    !< Number of macrocells in the system
      integer, intent(in) :: Mensemble    !< Number of ensembles
      integer, dimension(Natom), intent(in) :: cell_index            !< Macrocell index for each atom
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
      real(dblprec), dimension(3,3,Num_macro,Num_macro), intent(in):: Qdip_macro !< Matrix for the macro spin dipole-dipole interaction
      real(dblprec), dimension(3), intent(inout) :: field !< Effective field

      integer :: icell,jj

      ! Need to transform from atom to cell
      icell=cell_index(i)
      ! Field resulting from all the other macrocells
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422) && __INTEL_COMPILER < 1800
!         !$omp simd reduction(+:field)
#endif
      ! Contribution from the other macrocells
      do jj=1, Num_macro
         field(:)=field(:)+Qdip_macro(1,:,jj,icell)*emomM_macro(1,jj,k)&
         +Qdip_macro(2,:,jj,icell)*emomM_macro(2,jj,k)&
         +Qdip_macro(3,:,jj,icell)*emomM_macro(3,jj,k)
      enddo
      ! self-demagnetizing contribution
      field(:)=field(:)+Qdip_macro(1,:,icell,icell)*emomM_macro(1,icell,k)&
      +Qdip_macro(2,:,icell,icell)*emomM_macro(2,icell,k)&
      +Qdip_macro(3,:,icell,icell)*emomM_macro(3,icell,k)

   end subroutine macro_dipolar_field

   !----------------------------------------------------------------------------
   ! SUBROUTINE: dipole_cleanup
   !> @brief Wrapper for the deallocation of the dipole-dipole arrays
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine dipole_cleanup(do_dip,Num_macro)

      use HamiltonianData, only: allocate_dipole,allocate_macro_dipole

      implicit none

      integer, intent(in) :: Num_macro    !< Number of macrocells in the system
      integer, intent(in) :: do_dip

      if(do_dip==1) call allocate_dipole(flag=-1)
      if(do_dip==2) call allocate_macro_dipole(Num_macro,flag=-1)
      if (do_dip==3) then
#ifdef USE_FFTW
         call allocate_fft_dipole_structure(flag=-1)
#endif
#ifdef USE_MKL_FFT
         call allocate_fft_dipole_structure_MKL(flag=-1)
#endif
      endif

   end subroutine dipole_cleanup

end module DipoleManager
