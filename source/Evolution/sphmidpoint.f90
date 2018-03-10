!-------------------------------------------------------------------------------
!> MODULE: SphMidpoint
!> @brief
!> The semi-implicit spherical midpoint solver for the LLG-equations
!> Documented in
!> J. Hellsvik, Progress report: Semi-implicit spherical midpoint solver
!> and constructed by combining the solvers in
!> [1] McLachlan et al, Phys. Rev. E 89, 061301(R) (2014) and
!> [2] J. H. Mentink et al, J. Phys.: Condens. Matter, 22, 176001 (2010)
!> @author Johan Hellsvik
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
!-------------------------------------------------------------------------------
module SphMidpoint
   use Parameters
   use Profiling
   use ErrorHandling

   !implicit none


contains


   !-----------------------------------------------------------------------------
   !> SUBROUTINE: sibt
   !> @brief
   !> First step of a fixed point iteration variant of the semi-implicit midpoint solver SIB [2]
   !> @author Johan Hellsvik
   !-----------------------------------------------------------------------------
   subroutine sibt(Natom, Mensemble, Landeg,bn, lambda1_array, beff, emom, emom2, &
         emomM, mmom, deltat, thermal_field)
      use Constants
      use RandomNumbers, only : ranv

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0d0)
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: beff !< Total effective field
      !from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), intent(in) :: deltat !< Time step
      !
      ! ... Local variables ...
      integer :: i, j, ij, k
      real(dblprec) :: lldamp
      !
      ! deterministic variables
      real(dblprec) :: a1x, a1y, a1z
      real(dblprec) :: b1x, b1y, b1z
      !
      ! stochastic variables
      real(dblprec) :: s1x, s1y, s1z
      real(dblprec) :: f1x, f1y, f1z
      !
      ! auxilary variables
      real(dblprec) :: Ax, Ay, Az, detAi
      real(dblprec) :: a2x, a2y, a2z
      real(dblprec) :: efpdx, efpdy, efpdz
      real(dblprec) :: efpsx, efpsy, efpsz
      !
      ! time steps
      real(dblprec) :: dt,sqrtdt
      real(dblprec) :: dtg, sqrtdtg
      !
      ! spins
      real(dblprec) :: etx, ety, etz
      real(dblprec) :: e1x, e1y, e1z

      real(dblprec) :: efpx, efpy, efpz
      real(dblprec) :: efp0x, efp0y, efp0z
      real(dblprec) :: efp1x, efp1y, efp1z
      !
      ! tolerance for fix point iteration
      real(dblprec) :: epsilon, epsilona
      real(dblprec) :: ediff, ediff2, ediff3, ediff4
      real(dblprec) :: ediffx, ediffy, ediffz
      real(dblprec) :: enorm
      !
      ! Executable statements
      call ErrorHandling_missing('Implicit solver')

   end subroutine sibt


   !-----------------------------------------------------------------------------
   !> SUBROUTINE: sibf
   !> @brief
   !> Second step of a fixed point iteration variant of the semi-implicit midpoint solver SIB [2]
   !> @author Johan Hellsvik
   !-----------------------------------------------------------------------------
   subroutine sibf(Natom, Mensemble, Landeg, bn, lambda1_array, beff, emom, emom2, deltat)
      use Constants
      use RandomNumbers, only : ranv

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0d0)
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field
      !from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), intent(in) :: deltat !< Time step
      !
      ! ... Local variables ...
      integer :: i, j, ij, k
      real(dblprec) :: lldamp
      !
      ! deterministic variables
      real(dblprec) :: a1x, a1y, a1z
      real(dblprec) :: b1x, b1y, b1z
      !
      ! stochastic variables
      real(dblprec) :: s1x, s1y, s1z
      real(dblprec) :: f1x, f1y, f1z
      !
      ! auxilary variables
      real(dblprec) :: Ax, Ay, Az, detAi
      real(dblprec) :: a2x, a2y, a2z
      real(dblprec) :: efpdx, efpdy, efpdz
      real(dblprec) :: efpsx, efpsy, efpsz
      !
      ! time steps
      real(dblprec) :: dt,sqrtdt
      real(dblprec) :: dtg, sqrtdtg
      !
      ! spins

      real(dblprec) :: e1x, e1y, e1z
      real(dblprec) :: etpx, etpy, etpz
      real(dblprec) :: efpx, efpy, efpz
      real(dblprec) :: efp0x, efp0y, efp0z
      real(dblprec) :: efp1x, efp1y, efp1z
      !
      ! tolerance for fix point iteration
      real(dblprec) :: epsilon, epsilona
      real(dblprec) :: ediff, ediff2, ediff3, ediff4
      real(dblprec) :: ediffx, ediffy, ediffz
      real(dblprec) :: enorm
      !
      !
      ! Executable statements
      call ErrorHandling_missing('Implicit solver')

   end subroutine sibf

   !-----------------------------------------------------------------------------
   !> SUBROUTINE: sispht
   !> @brief
   !> First step of a semi-implicit variant of the spherical midpoint solver [1]
   !> @author Johan Hellsvik
   !-----------------------------------------------------------------------------
   subroutine sispht(Natom, Mensemble, Landeg,bn, lambda1_array, beff, emom, emom2, &
         emomM, mmom, deltat, thermal_field)
      use Constants
      use RandomNumbers, only : ranv

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0d0)
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: beff !< Total effective field
      !from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), intent(in) :: deltat !< Time step
      !
      ! ... Local variables ...
      integer :: i, j, ij, k
      real(dblprec) :: lldamp
      !
      ! deterministic variables
      real(dblprec) :: a1x, a1y, a1z
      real(dblprec) :: b1x, b1y, b1z
      !
      ! stochastic variables
      real(dblprec) :: s1x, s1y, s1z
      real(dblprec) :: f1x, f1y, f1z
      !
      ! auxilary variables


      real(dblprec) :: efpdx, efpdy, efpdz
      real(dblprec) :: efpsx, efpsy, efpsz
      !
      ! time steps
      real(dblprec) :: dt,sqrtdt
      real(dblprec) :: dtg, sqrtdtg
      !
      ! spins

      real(dblprec) :: e1x, e1y, e1z
      real(dblprec) :: efpx, efpy, efpz
      real(dblprec) :: efp0x, efp0y, efp0z
      real(dblprec) :: efpux, efpuy, efpuz, unorm
      !
      ! tolerance for fix point iteration
      real(dblprec) :: epsilon, epsilona
      real(dblprec) :: ediff, ediff2, ediff3, ediff4
      real(dblprec) :: ediffx, ediffy, ediffz
      real(dblprec) :: enorm
      !
      !
      ! Executable statements
      call ErrorHandling_missing('Implicit solver')

   end subroutine sispht


   !-----------------------------------------------------------------------------
   !> SUBROUTINE: sisphf
   !> @brief
   !> Second step of a semi-implicit variant of the spherical midpoint solver [1]
   !> @author Johan Hellsvik
   !-----------------------------------------------------------------------------
   subroutine sisphf(Natom, Mensemble, Landeg, bn, lambda1_array, beff, emom, emom2, deltat)
      use Constants
      use RandomNumbers, only : ranv

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0d0)
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field
      !from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), intent(in) :: deltat !< Time step
      !
      ! ... Local variables ...
      integer :: i, j, ij, k
      real(dblprec) :: lldamp
      !
      ! deterministic variables
      real(dblprec) :: a1x, a1y, a1z
      real(dblprec) :: b1x, b1y, b1z
      !
      ! stochastic variables
      real(dblprec) :: s1x, s1y, s1z
      real(dblprec) :: f1x, f1y, f1z
      !
      ! auxilary variables


      real(dblprec) :: efpdx, efpdy, efpdz
      real(dblprec) :: efpsx, efpsy, efpsz
      !
      ! time steps
      real(dblprec) :: dt,sqrtdt
      real(dblprec) :: dtg, sqrtdtg
      !
      ! spins

      real(dblprec) :: e1x, e1y, e1z
      real(dblprec) :: etpx, etpy, etpz
      real(dblprec) :: efpx, efpy, efpz
      real(dblprec) :: efp0x, efp0y, efp0z
      real(dblprec) :: efpux, efpuy, efpuz, unorm
      !
      ! tolerance for fix point iteration
      real(dblprec) :: epsilon, epsilona
      real(dblprec) :: ediff, ediff2, ediff3, ediff4
      real(dblprec) :: ediffx, ediffy, ediffz
      real(dblprec) :: enorm
      !
      !
      ! Executable statements
      call ErrorHandling_missing('Implicit solver')

   end subroutine sisphf


   !> Fix-point iteration step of implicit midpoint solver
   subroutine imp_fp(Natom, Mensemble, Landeg, bn, lambda1_array, beff, &
         emom, emom2, emomM, mmom, deltat, thermal_field, SDEalgh, converged)
      use Constants
      use RandomNumbers, only : ranv

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0d0)
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: beff !< Total effective field
      !from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), intent(in) :: deltat !< Time step
      integer, intent(in) :: SDEalgh !< Solver for equations of motion
      logical, intent(inout) :: converged !< Fix-point iteration converged or not?
      !
      ! ... Local variables ...
      integer :: i, j, ij
      real(dblprec) :: lldamp
      !
      ! deterministic variables
      real(dblprec) :: a1x, a1y, a1z
      real(dblprec) :: b1x, b1y, b1z
      !
      ! stochastic variables
      real(dblprec) :: s1x, s1y, s1z
      real(dblprec) :: f1x, f1y, f1z
      !
      ! auxilary variables


      real(dblprec) :: efpdx, efpdy, efpdz
      real(dblprec) :: efpsx, efpsy, efpsz
      !
      ! time steps
      real(dblprec) :: dt,sqrtdt
      real(dblprec) :: dtg, sqrtdtg
      !
      ! spins

      real(dblprec) :: e1x, e1y, e1z
      real(dblprec) :: etpx, etpy, etpz
      real(dblprec) :: efpx, efpy, efpz
      real(dblprec) :: efpux, efpuy, efpuz, unorm
      !
      ! tolerance for fix point iteration
      real(dblprec) :: epsilon, epsilona
      real(dblprec) :: ediff, ediff2, ediff3, ediff4
      real(dblprec) :: ediffx, ediffy, ediffz
      real(dblprec) :: enorm
      !
      !
      ! Executable statements
      call ErrorHandling_missing('Implicit solver')

   end subroutine imp_fp


end module sphmidpoint
