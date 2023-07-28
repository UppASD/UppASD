!-------------------------------------------------------------------------------
! MODULE: OSO
!> @brief Othogonal Spin Optimization routines. (https://doi.org/10.1016/j.cpc.2020.107749)
!> @author Anders Bergman 
!-------------------------------------------------------------------------------
module OSO
   use Parameters
   use HamiltonianActions
   use math_functions

   implicit none

   real(dblprec), dimension(:,:,:), allocatable :: torques !! Magnetic torque on atoms
   public

contains

subroutine calculate_torques_and_gradients(Natom, Mensemble, emomM, mmom,  beff, torques, gradients)
      implicit none


      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Mensemble    !< Number of ensembles
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !! Magnitude of the magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM !! Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !! Current magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: torques !! Magnetic torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: gradients !! Magnetic gradien

      integer :: iatom, iens

      !$omp parallel do default(shared) private(iens, iatom) collapse(2)
      do iens = 1, Mensemble
            do iatom = 1, Natom
                  torques(:, iatom, iens) = f_cross_product(emomM(:, iatom, iens), beff(:, iatom, iens))
                  gradients(1, iatom, iens) =  torques(3, iatom, iens)
                  gradients(2, iatom, iens) = -torques(2, iatom, iens)
                  gradients(3, iatom, iens) =  torques(1, iatom, iens)
            end do
      end do
      !$omp end parallel do

      return

end subroutine calculate_torques_and_gradients

subroutine calculate_rloc_and_theta(Natom, Mensemble, avec_k, r_theta, r_axis)
      implicit none

      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Mensemble    !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: avec_k  !! Exponents for spin rotations
      real(dblprec), dimension(Natom*Mensemble), intent(out) :: r_theta !! Rodrigues angle (note dimension here..)
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: r_axis  !! Rodriques axis

      integer :: iatom, iens

      r_theta = f_norms(3, Natom*Mensemble, avec_k)

      do iens = 1, Mensemble
            do iatom = 1, Natom
                  r_axis(1, iatom, iens) = -avec_k(3, iatom, iens)
                  r_axis(2, iatom, iens) =  avec_k(2, iatom, iens)
                  r_axis(3, iatom, iens) = -avec_k(3, iatom, iens)
            end do
      end do

      return

end subroutine calculate_rloc_and_theta

subroutine get_cg_direction(Natom, Mensemble, gradient_c, gradient_p, pvec_c, pvec_n, is_first)
      implicit none

      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Mensemble    !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: gradient_c !! Current gradient
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: gradient_p !! Previous gradient
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: pvec_c     !! Current path-vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: pvec_n    !! Next path-vector
      logical*1 :: is_first !! Flag to see if first iteration or not

      integer :: iatom, iens
      real(dblprec) :: gnorm_c2, gnorm_p2, beta_k

      if (is_first) then
            pvec_n = -gradient_c
      else
            gnorm_c2 = sum(gradient_c*gradient_c)
            gnorm_p2 = sum(gradient_p*gradient_p)

            beta_k = sqrt(gnorm_c2/gnorm_p2)

            pvec_n = beta_k * pvec_c - gradient_c
      end if

      return

end subroutine get_cg_direction


subroutine rodrigues_evolution(Natom, Mensemble, r_theta, r_rloc, mom_in, mom_out)
      implicit none

      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Mensemble    !< Number of ensembles
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: r_theta    !! Current gradient
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: r_rloc   !! Previous gradient
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: mom_in   !! Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: mom_out !! Evolved magnetic moment vector

      integer :: iatom, iens

      do iens = 1, Mensemble
            do iatom = 1, Natom
            end do
      end do

end subroutine rodrigues_evolution

end module OSO