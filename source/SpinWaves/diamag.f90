!> Data and routines for calculations of adiabatic spin wave spectra for 1q-structures based on LSWT
!> @author 
!> Anders Bergman, Manuel Pereiro
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
!!
module diamag
  use Parameters
  use Profiling
   use ErrorHandling
  !
  implicit none
  !
  !
  real(dblprec), dimension(:,:), allocatable :: magham !< M(w)
  complex(dblprec), dimension(:,:,:,:), allocatable :: magmom_qw !< M(w) in reciprocal space
  complex(dblprec), dimension(:,:), allocatable :: A_k
  complex(dblprec), dimension(:,:), allocatable :: B_k
  complex(dblprec), dimension(:,:), allocatable :: C_k
  complex(dblprec), dimension(:,:), allocatable :: h_k

  !
  character(len=1) :: do_diamag    !< Perform frequency based spin-correlation sampling (Y/N/C)
  real(dblprec)    :: diamag_mix   !< Separation between sampling steps
  real(dblprec)    :: diamag_thresh   !< Separation between sampling steps
  integer          :: diamag_niter !< Number of steps to sample

  private
  ! public subroutines
  public :: do_diamag
  public ::  setup_diamag

contains

  subroutine setup_diamag()

     implicit none

     do_diamag='Y'
     diamag_niter=1000
     diamag_mix=0.030d0
     diamag_thresh=1.0d-8

  end subroutine setup_diamag

  !> Set up the Hamiltonian for first cell
  subroutine setup_finite_hamiltonian(N1,N2,N3,NT,NA,Natom, Mensemble, conf_num, simid, emomM, &
                                            mmom, max_no_neigh, nlistsize, nlist, ncoup, kaniso )

  use Constants
    !
    implicit none
    !
    integer, intent(in) :: N1  !< Number of cell repetitions in x direction
    integer, intent(in) :: N2  !< Number of cell repetitions in y direction
    integer, intent(in) :: N3  !< Number of cell repetitions in z direction
    integer, intent(in) :: NT                                          !< Number of types of atoms
    integer, intent(in) :: NA                                          !< Number of atoms in one cell
    integer, intent(in) :: Natom     !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles
    integer, intent(in) :: conf_num  !< number of configurations for LSF
    character(len=8), intent(in) :: simid !< Name of simulation
    real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment magnitude
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector

    !
    integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
    integer, dimension(Natom), intent(in) :: nlistsize     !< Size of neighbour list for Heisenberg exchange couplings
    integer, dimension(max_no_neigh,Natom), intent(in) :: nlist       !< Neighbour list for Heisenberg exchange couplings
    real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup  !< Heisenberg exchange couplings
    real(dblprec), dimension(2,Natom), intent(in) :: kaniso  !< Anisotropy constant
    !
    integer :: i_stat, ia, ja, j, hdim, l, la
    complex(dblprec), dimension(3) :: ui, vi, uj, vj, ul, vl
    complex(dblprec) :: udot, ucdot, vdot
    complex(dblprec), dimension(NA,NA) :: dot_mat
    !
    !
            call ErrorHandling_missing('Generalized LSWT')

   end subroutine setup_finite_hamiltonian

subroutine diagonalize_quad_hamiltonian(Natom,Mensemble,NA,mmom,emomM)
   !
   use Constants
   !
   implicit none
   !
   integer, intent(in) :: Natom     !< Number of atoms in system
   integer, intent(in) :: Mensemble !< Number of ensembles
   integer, intent(in) :: NA
   real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment magnitude
   real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector
   !
   !
   complex(dblprec), dimension(:,:), allocatable :: g_mat
   complex(dblprec), dimension(:,:), allocatable :: U_mat
   complex(dblprec), dimension(:,:), allocatable :: T_mat
   complex(dblprec), dimension(:,:), allocatable :: E_mat
   complex(dblprec), dimension(:,:), allocatable :: sqE_mat
   complex(dblprec), dimension(:,:), allocatable :: L_mat
   complex(dblprec), dimension(:,:), allocatable :: K_mat
   complex(dblprec), dimension(:,:), allocatable :: iK_mat
   complex(dblprec), dimension(:,:), allocatable :: dum_mat, dum_mat2
   complex(dblprec), dimension(:,:), allocatable :: eig_vec
   complex(dblprec), dimension(:,:), allocatable :: x_mat
   complex(dblprec), dimension(:,:,:), allocatable :: bigS
   real(dblprec), dimension(:), allocatable :: eig_val
   complex(dblprec), dimension(:), allocatable :: cwork
   real(dblprec), dimension(:), allocatable :: rwork
   !
   integer :: info, lwork, hdim, i_stat, ia, ja

   complex(dblprec) :: cone, czero, fcinv, dia_eps, im

   complex(dblprec), dimension(3) :: ul, vl
   !
   !
            call ErrorHandling_missing('Generalized LSWT')

   !
end subroutine diagonalize_quad_hamiltonian
  !> Set up the Hamiltonian for first cell
  subroutine diagonalize_finite_hamiltonian(NA)

  use Constants
    !
    implicit none
    !
    integer, intent(in) :: NA                                          !< Number of atoms in one cell

    integer :: i_stat
    integer :: lwork,info, hdim


    real(dblprec), dimension(:,:), allocatable :: Hprime, eig_vec
    real(dblprec), dimension(:), allocatable :: eig_val, work
    !

            call ErrorHandling_missing('Generalized LSWT')

   end subroutine diagonalize_finite_hamiltonian

subroutine find_eigenvector(hdim,eig_vec,eig_val)

  use Constants
    !
    implicit none
    !
    integer, intent(in) :: hdim                   !< Number of atoms in one cell
    real(dblprec), dimension(hdim,hdim), intent(in) :: eig_vec
    real(dblprec), dimension(hdim), intent(in) :: eig_val
    !
    integer :: i_stat, ia, ivec, iter, i

    real(dblprec) :: enorm
    real(dblprec) :: alpha
    real(dblprec), dimension(:,:), allocatable :: Hprime, moments, err_norm
    real(dblprec), dimension(:), allocatable :: mom_vec, error
    !

            call ErrorHandling_missing('Generalized LSWT')

    !
   end subroutine find_eigenvector

subroutine moments_to_angles(moments,ndim,eig_val)
   !
   !
   implicit none
   !
   !
   integer, intent(in) :: ndim
   real(dblprec), dimension(ndim,ndim), intent(in) :: moments
   real(dblprec), dimension(ndim), intent(in) :: eig_val
   !
   real(dblprec) :: pi, mi_norm, mj_norm, mdot
   real(dblprec), dimension(ndim/3,ndim) :: angles
   real(dblprec), dimension(3) :: mom_i, mom_j

   integer :: i,j, idx
   !
            call ErrorHandling_missing('Generalized LSWT')

end subroutine moments_to_angles

subroutine normalize_moments(moments,nrow,ncol)
   !
   !
   implicit none
   !
   !
   integer, intent(in) :: nrow
   integer, intent(in) :: ncol
   real(dblprec), dimension(ncol,nrow), intent(inout) :: moments
   !
   real(dblprec) :: mnorm
   integer :: i,j
   !
            call ErrorHandling_missing('Generalized LSWT')

end subroutine normalize_moments

subroutine find_uv(u,v,mom)
   !
   !
   implicit none
   !
   complex(dblprec), dimension(3), intent(out) :: u
   complex(dblprec), dimension(3), intent(out) :: v
   real(dblprec), dimension(3), intent(in) :: mom
   !

   real(dblprec) :: mnorm
   complex(dblprec) :: im = (0.0d0,1.0d0)
   real(dblprec), dimension(3,3) :: R_prime
   real(dblprec), dimension(3) :: mom_hat
   !
   !
            u=0.0d0;v=0.0d0
            call ErrorHandling_missing('Generalized LSWT')

end subroutine find_uv

  !> Set up the Hamiltonian for first cell
  subroutine setup_tensor_hamiltonian(N1,N2,N3,NT,NA,Natom, Mensemble, conf_num, simid, emomM, &
                                             mmom, max_no_neigh, nlistsize, nlist, jtens)

  use Constants
    !
    implicit none
    !
    integer, intent(in) :: N1  !< Number of cell repetitions in x direction
    integer, intent(in) :: N2  !< Number of cell repetitions in y direction
    integer, intent(in) :: N3  !< Number of cell repetitions in z direction
    integer, intent(in) :: NT                                          !< Number of types of atoms
    integer, intent(in) :: NA                                          !< Number of atoms in one cell
    integer, intent(in) :: Natom     !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles
    integer, intent(in) :: conf_num  !< number of configurations for LSF
    character(len=8), intent(in) :: simid !< Name of simulation
    real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment magnitude
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector

    !
    integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
    integer, dimension(Natom), intent(in) :: nlistsize     !< Size of neighbour list for Heisenberg exchange couplings
    integer, dimension(max_no_neigh,Natom), intent(in) :: nlist       !< Neighbour list for Heisenberg exchange couplings
    real(dblprec), dimension(3,3,max_no_neigh,Natom), intent(in) :: jtens  !< Heisenberg exchange couplings in tensor form
    !
    integer :: i_stat, ia, ja, j, hdim, l, la
    complex(dblprec), dimension(3) :: ui, vi, uj, vj, ul, vl
    complex(dblprec) :: udot, ucdot, vdot
    complex(dblprec), dimension(NA,NA) :: dot_mat
    !
    !
            call ErrorHandling_missing('Generalized LSWT')

   end subroutine setup_tensor_hamiltonian

complex(dblprec) function sJs(u,J,v)
   !
   implicit none
   !
   complex(dblprec), dimension(3), intent(in) :: u
   real(dblprec), dimension(3,3), intent(in) :: J
   complex(dblprec), dimension(3), intent(in) :: v
   !
   complex(dblprec), dimension(3) :: dt
   !
   sJs=(0.0d0,0.0d0)
            call ErrorHandling_missing('Generalized LSWT')

end function sJs

subroutine shuffle_eig(eval,evec,ndim)
   !
   implicit none
   !
   integer, intent(in) :: ndim
   real(dblprec), dimension(ndim), intent(inout) :: eval
   complex(dblprec), dimension(ndim,ndim), intent(inout) :: evec
   !
   integer :: i, nhit, thit
   complex(dblprec), dimension(:,:), allocatable :: mtemp
   complex(dblprec), dimension(:), allocatable :: dmtemp
   real(dblprec), dimension(:), allocatable :: vtemp
   real(dblprec) :: dtemp
   !
            call ErrorHandling_missing('Generalized LSWT')

end subroutine shuffle_eig

end module diamag
