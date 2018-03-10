!> Spin-spiral energy minimization
!> @author
!> Anders Bergman
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module qminimizer
   use Parameters
   use Profiling
   use HamiltonianData
   use FieldData
   use HamiltonianActions
   use ErrorHandling
   !
   !
   implicit none
   !
   !
   real*8 :: fac_2d
   !
   ! Data for minimization and ss-vectors
   real*8, dimension(3,3) :: q,s
   real*8, dimension(3) :: theta, phi
   !
   integer :: nq=1
   !
   private
   !
   public :: mini_q, plot_q, qmc
   !
contains
   !
   !> Main driver
   subroutine mini_q(Natom, Mensemble, NA, coord, do_jtensor, exc_inter, do_dm, do_pd, &
         do_biqdm, do_bq, taniso, sb, do_dip, emomM, mmom, hfield, &
         OPT_flag, max_no_constellations, maxNoConstl, &
         unitCellType, constlNCoup, constellations, constellationsNeighType, mult_axis)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NA  !< Number of atoms in one cell
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      integer, intent(in) :: do_jtensor   !<  Use SKKR style exchange tensor (0=off, 1=on, 2=with biquadratic exchange)
      character(len=1),intent(in) :: exc_inter !< Interpolate Jij (Y/N)
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
      integer, intent(in) :: do_biqdm   !< Add Biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
      integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
      real(dblprec), dimension(Natom),intent(in) :: sb !< Ratio between the anisotropies
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3), intent(in) :: hfield !< Constant effective field
      character(len=1), intent(in) :: mult_axis


      !! +++ New variables due to optimization routines +++ !!
      integer, intent(in) :: max_no_constellations ! The maximum (global) length of the constellation matrix
      ! Number of entries (number of unit cells*number of atoms per unit cell) in the constellation matrix per ensemble
      integer, dimension(Mensemble), intent(in) :: maxNoConstl
      ! See OptimizationRoutines.f90 for details on classification
      integer, dimension(Natom, Mensemble), intent(in) :: unitCellType ! Array of constellation id and classification (core, boundary, or noise) per atom
      ! Matrix relating the interatomic exchanges for each atom in the constellation matrix
      real(dblprec), dimension(max_no_neigh, max_no_constellations,Mensemble), intent(in) :: constlNCoup
      ! Matrix storing all unit cells belonging to any constellation
      real(dblprec), dimension(3,max_no_constellations, Mensemble), intent(in) :: constellations
      ! Optimization flag (1 = optimization on; 0 = optimization off)
      logical, intent(in) :: OPT_flag
      ! Matrix storing the type of the neighbours within a given neighbourhood of a constellation; default is 1 outside the neighbourhood region
      ! The default is to achieve correct indexing. Note here also that constlNCoup will result in a net zero contribution to the Heissenberg exchange term
      integer, dimension(max_no_neigh,max_no_constellations,Mensemble), intent(in) :: constellationsNeighType
      ! Internal effective field arising from the optimization of the Heissenberg exchange term
      !
      !
      integer :: iq
      !
      real*8, dimension(3) :: m_i, m_j
      real*8 :: pi, qr
      integer :: i,j, k, ia, ja, lhit, nhits
      real*8 :: energy, min_energy
      !
      real*8, dimension(3,3) :: q_best, q_diff, s_diff
      real*8, dimension(3,3) :: s_save, q_save
      real*8, dimension(3) :: theta_save, theta_diff, theta_best
      real*8, dimension(3) :: phi_save, phi_diff, phi_best
      !
      integer :: niter, iter, iscale
      real*8 :: q_range, theta_range, s_range, phi_range
      !
      call ErrorHandling_missing('Q-minimization')

      !
   end subroutine mini_q

   !> Plot final configuration
   subroutine plot_q(Natom, Mensemble, coord, emom, emomM, mmom)
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      !
      real*8, dimension(3) :: srvec, m_j
      integer :: lhit, ia, nplot, k, iq
      real*8 :: pi, qr
      !
      !
      call ErrorHandling_missing('Q-minimization')

      !
   end subroutine plot_q

   !> Anisotropy
   subroutine spinspiral_ani_field(Natom,Mensemble,NA,mmom,taniso,sb,hfield,energy)
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NA  !< Number of atoms in one cell
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
      real(dblprec), dimension(Natom),intent(in) :: sb !< Ratio between the anisotropies
      real(dblprec), dimension(3), intent(in) :: hfield !< Constant effective field
      real(dblprec), intent(inout), optional :: energy !< Total energy
      !
      !
      integer :: iq,i
      real(dblprec) :: tt1, tt2, tt3, totmom, qfac
      real(dblprec), dimension(3) :: field
      !
      call ErrorHandling_missing('Q-minimization')

   end subroutine spinspiral_ani_field

   !> Rotation of vectors
   subroutine rotvec(s,m)
      !
      implicit none
      !
      real*8, dimension(3), intent(in) :: s
      real*8, dimension(3), intent(out) :: m
      !
      real*8 :: theta, dot, qnorm
      real*8, dimension(3) :: u,q_new, q
      real*8, dimension(3,3) :: I,ux,uplus, R
      !
      m=0.0d0;
      call ErrorHandling_missing('Q-minimization')

      !
   end subroutine rotvec

   !> Normalization
   subroutine normalize(v)
      !
      implicit none
      !
      real*8,dimension(3), intent(inout) :: v
      !
      real*8 :: vnorm
      !
      call ErrorHandling_missing('Q-minimization')

   end subroutine normalize

   !> Update moment
   subroutine updatrotmom_single(m_in,s_vec)
      !
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !
      !.. Formal Arguments ..
      real(dblprec), dimension(3), intent(inout) :: m_in
      real(dblprec), dimension(3), intent(in) :: s_vec
      !
      !
      !.. Local Scalars ..
      integer :: j
      real(dblprec) :: alfa, beta
      !
      !.. Local Arrays ..
      real(dblprec), dimension(3) :: v,vout, sv, mz
      real(dblprec), dimension(3,3) :: Rx, Ry, Rz
      !
      !.. Intrinsic Functions ..
      !intrinsic selected_real_kind
      !
      ! ... Executable Statements ...
      !
      !
      call ErrorHandling_missing('Q-minimization')

   end subroutine updatrotmom_single

   !> transforms cartesian (x,y,z) to spherical (Theta,Phi,R) coordinates
   subroutine car2sph(C,S)
      ! transforms cartesian (x,y,z) to spherical (Theta,Phi,R) coordinates
      !
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !
      !.. Formal Arguments ..
      real(dblprec) :: X,Y,Z, THETA, PHI, D2, R2
      real(dblprec), dimension(3), intent(in) :: C
      real(dblprec), dimension(3), intent(out) :: S
      !
      !
      ! ... Executable Statements ...
      !
      !
      S=0.0d0;
      call ErrorHandling_missing('Q-minimization')

   end subroutine car2sph

   !> Energy minimization
   subroutine qmc(Natom, Mensemble, NA, coord, do_jtensor, exc_inter, do_dm, do_pd, do_biqdm, &
         do_bq, taniso, sb, do_dip, emomM, mmom, hfield, &
         OPT_flag, max_no_constellations, maxNoConstl, &
         unitCellType, constlNCoup, constellations, constellationsNeighType, mult_axis)
      !
      use Constants, only: k_bolt, mub
      use InputData, only: Temp
      use RandomNumbers, only : rng_uniform
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NA  !< Number of atoms in one cell
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      integer, intent(in) :: do_jtensor   !<  Use SKKR style exchange tensor (0=off, 1=on, 2=with biquadratic exchange)
      character(len=1),intent(in) :: exc_inter !< Rescale Jij (Y/N)
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
      integer, intent(in) :: do_biqdm   !< Add Biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
      integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
      real(dblprec), dimension(Natom),intent(in) :: sb !< Ratio between the anisotropies
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3), intent(in) :: hfield !< Constant effective field
      character(len=1), intent(in) :: mult_axis

      !! +++ New variables due to optimization routines +++ !!
      integer, intent(in) :: max_no_constellations ! The maximum (global) length of the constellation matrix
      ! Number of entries (number of unit cells*number of atoms per unit cell) in the constellation matrix per ensemble
      integer, dimension(Mensemble), intent(in) :: maxNoConstl
      ! See OptimizationRoutines.f90 for details on classification
      integer, dimension(Natom, Mensemble), intent(in) :: unitCellType ! Array of constellation id and classification (core, boundary, or noise) per atom
      ! Matrix relating the interatomic exchanges for each atom in the constellation matrix
      real(dblprec), dimension(max_no_neigh, max_no_constellations,Mensemble), intent(in) :: constlNCoup
      ! Matrix storing all unit cells belonging to any constellation
      real(dblprec), dimension(3,max_no_constellations, Mensemble), intent(in) :: constellations
      ! Optimization flag (1 = optimization on; 0 = optimization off)
      logical, intent(in) :: OPT_flag
      ! Matrix storing the type of the neighbours within a given neighbourhood of a constellation; default is 1 outside the neighbourhood region
      ! The default is to achieve correct indexing. Note here also that constlNCoup will result in a net zero contribution to the Heissenberg exchange term
      integer, dimension(max_no_neigh,max_no_constellations,Mensemble), intent(in) :: constellationsNeighType
      ! Internal effective field arising from the optimization of the Heissenberg exchange term

      integer :: iq
      !
      real(dblprec), dimension(3) :: m_i, m_j
      real(dblprec) :: pi, qr, rn(3)
      integer :: i,j, k, ia, ja, lhit, nhits
      real(dblprec) :: energy, old_energy, flipprob(1), de, beta
      real(dblprec), dimension(3) :: avgmom
      real(dblprec) :: avgM
      !
      real(dblprec), dimension(3,3) :: q_trial, s_trial
      real(dblprec), dimension(3) :: theta_trial, phi_trial
      real(dblprec), dimension(3,Natom,Mensemble) :: emomM_tmp
      !
      integer :: niter, iter, iscale

      !
      call ErrorHandling_missing('Q-minimization')

   end subroutine qmc

end module qminimizer
