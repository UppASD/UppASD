!====================================================================!
!> @brief
!> Routines used to calculate topological properties of the magnetic system
!> Mostly related to the calculation of the skyrmion number
!
!> @author
!> Anders Bergman
!> Jonathan Chico ---> Reorganized in different modules, added site and type dependance
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
!====================================================================!
module Topology

  use Parameters
  use Profiling

  ! Parameters for the printing
  character(len=1) :: skyno      !< Perform skyrmion number measurement
  character(len=1) :: do_proj_skyno !< Perform type dependent skyrmion number measurement
  character(len=1) :: do_skyno_den  !< Perform site dependent skyrmion number measurement
  integer :: skyno_step !< Interval for sampling the skyrmion number
  integer :: skyno_buff !< Buffer size for the sampling of the skyrmion number

  public

  contains

!--------------------------------------------------------------------!
!> @brief
!> Calculates the total skyrmion number of the system
!
!> @author
!> Anders Bergman
!--------------------------------------------------------------------!
  real(dblprec) function pontryagin_no(Natom,Mensemble,emomM,grad_mom)
    use Constants

    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles
    real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
    real(dblprec), dimension(3,3,Natom, Mensemble), intent(in) :: grad_mom  !< Gradient of magnetic moment vector
 !  integer, intent(in) :: mstep !< Current time step

    integer :: iatom, k
    real(dblprec) :: thesum,cvec_x,cvec_y,cvec_z

    thesum=0.0d0

#if _OPENMP >= 201307 && __INTEL_COMPILER < 1800
    !$omp parallel do default(shared) private(iatom,k,cvec_x,cvec_y,cvec_z) reduction(+:thesum)
#endif
    do iatom=1, Natom
       do k=1, Mensemble
          cvec_x=grad_mom(2,1,iatom,k)*grad_mom(3,2,iatom,k)-grad_mom(3,1,iatom,k)*grad_mom(2,2,iatom,k)
          cvec_y=grad_mom(3,1,iatom,k)*grad_mom(1,2,iatom,k)-grad_mom(1,1,iatom,k)*grad_mom(3,2,iatom,k)
          cvec_z=grad_mom(1,1,iatom,k)*grad_mom(2,2,iatom,k)-grad_mom(2,1,iatom,k)*grad_mom(1,2,iatom,k)
          thesum=thesum+emomM(1,iatom,k)*cvec_x+emomM(2,iatom,k)*cvec_y+emomM(3,iatom,k)*cvec_z
       end do
    end do
#if _OPENMP >= 201307 && __INTEL_COMPILER < 1800
    !$omp end parallel do
#endif

    !pontryagin_no=thesum/4.0d0/pi
    pontryagin_no=thesum/2.0d0/pi/sqrt(1.0d0*Natom)/Mensemble
    !
    return
    !
  end function pontryagin_no

  !---------------------------------------------------------------------------
  !> @brief
  !> Calculation of the site dependent skyrmion number
  !
  !> @author
  !> Jonathan Chico
  !---------------------------------------------------------------------------
  real(dblprec) function pontryagin_no_density(iatom,Natom,Mensemble,emomM,grad_mom)

    use Constants

    implicit none

    !.. Input variables
    integer, intent(in) :: iatom !< Current atomic position
    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles
    real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
    real(dblprec), dimension(3,3,Natom, Mensemble), intent(in) :: grad_mom  !< Gradient of magnetic moment vector

    ! .. Local variables
    integer :: k
    real(dblprec) :: thesum,cvec_x,cvec_y,cvec_z

    thesum=0.0d0

    do k=1, Mensemble
       cvec_x=grad_mom(2,1,iatom,k)*grad_mom(3,2,iatom,k)-grad_mom(3,1,iatom,k)*grad_mom(2,2,iatom,k)
       cvec_y=grad_mom(3,1,iatom,k)*grad_mom(1,2,iatom,k)-grad_mom(1,1,iatom,k)*grad_mom(3,2,iatom,k)
       cvec_z=grad_mom(1,1,iatom,k)*grad_mom(2,2,iatom,k)-grad_mom(2,1,iatom,k)*grad_mom(1,2,iatom,k)
       thesum=thesum+emomM(1,iatom,k)*cvec_x+emomM(2,iatom,k)*cvec_y+emomM(3,iatom,k)*cvec_z
    end do

    pontryagin_no_density=thesum/4.0d0/pi

  end function pontryagin_no_density

  !---------------------------------------------------------------------------
  !> @brief
  !> Calculation of the type dependent skyrmion number
  !
  !> @author
  !> Jonathan Chico
  !---------------------------------------------------------------------------
  function proj_pontryagin_no(NT,Natom,Mensemble,atype,emomM,proj_grad_mom)
    use Constants

    implicit none

    integer, intent(in) :: NT    !< Number of types of atoms
    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles
    integer, dimension(Natom), intent(in) :: atype !< Type of atom
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
    real(dblprec), dimension(3,3,Natom,Mensemble,NT), intent(in) :: proj_grad_mom  !< Gradient of magnetic moment vector

    real(dblprec), dimension(NT) :: proj_pontryagin_no

    integer :: iatom, k,ii
    real(dblprec), dimension(NT) :: thesum,cvec_x,cvec_y,cvec_z

    thesum=0.0d0

#if _OPENMP >= 201307 && __INTEL_COMPILER < 1800
    !$omp parallel do default(shared) private(iatom,k,cvec_x,cvec_y,cvec_z,ii) reduction(+:thesum)
#endif
    do iatom=1, Natom
       do k=1, Mensemble
          ii=atype(iatom)
          cvec_x(ii)=proj_grad_mom(2,1,iatom,k,ii)*proj_grad_mom(3,2,iatom,k,ii)-proj_grad_mom(3,1,iatom,k,ii)*proj_grad_mom(2,2,iatom,k,ii)
          cvec_y(ii)=proj_grad_mom(3,1,iatom,k,ii)*proj_grad_mom(1,2,iatom,k,ii)-proj_grad_mom(1,1,iatom,k,ii)*proj_grad_mom(3,2,iatom,k,ii)
          cvec_z(ii)=proj_grad_mom(1,1,iatom,k,ii)*proj_grad_mom(2,2,iatom,k,ii)-proj_grad_mom(2,1,iatom,k,ii)*proj_grad_mom(1,2,iatom,k,ii)
          thesum(ii)=thesum(ii)+emomM(1,iatom,k)*cvec_x(ii)+emomM(2,iatom,k)*cvec_y(ii)+emomM(3,iatom,k)*cvec_z(ii)
       end do
    end do
#if _OPENMP >= 201307 && __INTEL_COMPILER < 1800
    !$omp end parallel do
#endif

    proj_pontryagin_no=thesum/4.0d0/pi
    !
    return
    !
  end function proj_pontryagin_no

end module Topology
