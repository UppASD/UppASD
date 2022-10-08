!> Data containing simulation parameters
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
   real(dblprec), parameter :: bn=1.0_dblprec   !< Scaling factor for LLG equation
   !
   real(dblprec) :: total_energy = 0.0_dblprec


end module SimulationData
