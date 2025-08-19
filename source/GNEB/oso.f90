!-------------------------------------------------------------------------------
! MODULE: OSO
!> @brief Othogonal Spin Optimization routines. (https://doi.org/10.1016/j.cpc.2020.107749)
!> @author Anders Bergman 
!-------------------------------------------------------------------------------
module OSO
   use Parameters
   use HamiltonianActions
   use cg_direction
   use wolfe_search
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

subroutine calculate_rloc_and_theta(Natom,Mensemble,avec_k,r_theta,r_axis)
  integer,intent(in) :: Natom,Mensemble
  real(dblprec),intent(in)  :: avec_k(3,Natom,Mensemble)      ! a12,a13,a23
  real(dblprec),intent(out) :: r_theta(Natom,Mensemble)        ! θ
  real(dblprec),intent(out) :: r_axis(3,Natom,Mensemble)       ! unit axis r
  integer :: i,j
  real(dblprec) :: a12,a13,a23,theta,invtheta

  do j=1,Mensemble
    do i=1,Natom
      a12 = avec_k(1,i,j);  a13 = avec_k(2,i,j);  a23 = avec_k(3,i,j)
      theta = sqrt(a12*a12 + a13*a13 + a23*a23)
      r_theta(i,j) = theta
      if (theta > 1.0d-14) then
        invtheta = 1.0d0/theta
        r_axis(1,i,j) = -a23*invtheta
        r_axis(2,i,j) =  a13*invtheta
        r_axis(3,i,j) = -a12*invtheta
      else
        r_axis(:,i,j) = 0.0d0
        r_axis(3,i,j) = 1.0d0  ! arbitrary unit axis when θ≈0
      end if
    end do
  end do
end subroutine

! subroutine get_cg_direction(Natom, Mensemble, gradient_c, gradient_p, pvec_c, pvec_n, is_first)
!       implicit none
! 
!       integer, intent(in) :: Natom        !< Number of atoms in system
!       integer, intent(in) :: Mensemble    !< Number of ensembles
!       real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: gradient_c !! Current gradient
!       real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: gradient_p !! Previous gradient
!       real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: pvec_c     !! Current path-vector
!       real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: pvec_n    !! Next path-vector
!       logical*1 :: is_first !! Flag to see if first iteration or not
! 
!       integer :: iatom, iens
!       real(dblprec) :: gnorm_c2, gnorm_p2, beta_k
! 
!       if (is_first) then
!             pvec_n = -gradient_c
!       else
!             gnorm_c2 = sum(gradient_c*gradient_c)
!             gnorm_p2 = sum(gradient_p*gradient_p)
! 
!             beta_k = sqrt(gnorm_c2/gnorm_p2)
! 
!             pvec_n = beta_k * pvec_c - gradient_c
!       end if
! 
!       return
! 
! end subroutine get_cg_direction


subroutine rodrigues_evolution(Natom,Mensemble,r_theta,r_axis,mom_in,mom_out)
  integer,intent(in) :: Natom,Mensemble
  real(dblprec),intent(in)  :: r_theta(Natom,Mensemble), r_axis(3,Natom,Mensemble)
  real(dblprec),intent(in)  :: mom_in(3,Natom,Mensemble)
  real(dblprec),intent(out) :: mom_out(3,Natom,Mensemble)
  integer :: i,j
  real(dblprec) :: th,ct,st,omt,rdote, r(3), e(3), rxE(3)

  do j=1,Mensemble
    do i=1,Natom
      th = r_theta(i,j);          e = mom_in(:,i,j);     r = r_axis(:,i,j)
      if (th < 1.0d-14) then
        mom_out(:,i,j) = e
      else
        ct = cos(th);  st = sin(th);  omt = 1.0d0 - ct
        ! cross and dot
        rxE(1) = r(2)*e(3) - r(3)*e(2)
        rxE(2) = r(3)*e(1) - r(1)*e(3)
        rxE(3) = r(1)*e(2) - r(2)*e(1)
        rdote  = r(1)*e(1) + r(2)*e(2) + r(3)*e(3)
        mom_out(:,i,j) = ct*e + st*rxE + omt*rdote*r
        ! (optional) re-normalize to unit length
        mom_out(:,i,j) = mom_out(:,i,j) / sqrt(sum(mom_out(:,i,j)**2))
      end if
    end do
  end do
end subroutine

end module OSO