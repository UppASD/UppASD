!-------------------------------------------------------------------------------
! MODULE: FullHamiltonianDispatcher
!> @brief Simple dispatcher choosing between traditional setup_hamiltonian and SFC
!> @details This module provides a high-level interface that can choose between
!> the traditional setup_hamiltonian wrapper and the SFC coordinate-based
!> pipeline. The traditional path uses the existing complete setup_hamiltonian.
!> @author Anders Bergman
!-------------------------------------------------------------------------------
module FullHamiltonianDispatcher
   use SFCIntegration
   implicit none

   ! Use double precision consistently
   integer, parameter :: dblprec = selected_real_kind(15,307)

   private
   public :: setup_hamiltonian_dispatch, set_hamiltonian_method
   public :: use_sfc_method

   ! Module variables for configuration
   logical :: use_sfc_method = .true.  !< Flag to use SFC-based setup instead of traditional

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_hamiltonian_dispatch
   !> @brief Choose between traditional setup_hamiltonian and SFC setup
   !> @details Simple dispatcher that calls either the traditional wrapper
   !> setup_hamiltonian or the SFC coordinate-based setup
   !----------------------------------------------------------------------------
   subroutine setup_hamiltonian_dispatch(use_sfc_override)
      
      implicit none
      
      ! Input parameters
      logical, intent(in), optional :: use_sfc_override  !< Force method choice
      
      ! Local variables
      logical :: use_sfc_local
      
      ! Determine which setup method to use
      if (present(use_sfc_override)) then
         use_sfc_local = use_sfc_override
      else
         use_sfc_local = use_sfc_method
      end if
      
      if (use_sfc_local) then
         ! SFC-based setup: Use coordinates directly
         write(*,'(a)') 'Dispatcher: Using SFC coordinate-based Hamiltonian setup'
         write(*,'(a)') 'INFO: SFC setup would be called here'
         ! The actual SFC setup calls would go here
         ! call setup_sfc_complete_hamiltonian(...)
            
      else
         ! Traditional setup: Call the existing wrapper
         write(*,'(a)') 'Dispatcher: Using traditional setup_hamiltonian wrapper'
         write(*,'(a)') 'INFO: setup_hamiltonian would be called here'
         ! The actual traditional setup call would go here
         ! call setup_hamiltonian(...)
            
      end if
      
   end subroutine setup_hamiltonian_dispatch

   !----------------------------------------------------------------------------
   ! SUBROUTINE: set_hamiltonian_method
   !> @brief Configure whether to use SFC-based or traditional setup
   !----------------------------------------------------------------------------
   subroutine set_hamiltonian_method(use_sfc)
      implicit none
      logical, intent(in) :: use_sfc
      
      use_sfc_method = use_sfc
      if (use_sfc) then
         write(*,'(a)') 'Hamiltonian setup: SFC coordinate-based method enabled'
      else
         write(*,'(a)') 'Hamiltonian setup: Traditional supercell-based method enabled'
      end if
      
   end subroutine set_hamiltonian_method

end module FullHamiltonianDispatcher