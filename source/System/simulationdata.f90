!> Data containing simulation parameters
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module SimulationData
   use Parameters
   use Profiling
   !
   implicit none
   !
   real(dblprec) :: lambda1 !< Damping parameter
   real(dblprec) :: lambda2 !< Additional damping parameter (not used for llg=1)
   integer :: rstep !< Starting simulation step
   integer :: mstep !< Current simulation step
   real(dblprec), parameter :: bn=1.0d0   !< Scaling factor for LLG equation
   !


end module SimulationData
