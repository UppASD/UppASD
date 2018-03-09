!====================================================================!
!> @brief
!> Intended to use for the treatment of induced magnetic moments
!
!> @author
!> Jonathan Chico
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
!====================================================================!

! Routines needed to calculate the effect of induced moments via linear response
! treatment as outlined by S. Polesya et al. in Phys. Rev. B 82, 214409.
! In this approach the magnitude and direction of the induced moments is determined
! by the nearest neighbours fixed moments
module InducedMoments

  use Parameters
  use Profiling
  use ErrorHandling

  implicit none

  ! Interestingly enough, one cannot define this arrays with different names, if one does
  ! the routine setup_nm complains
  integer, dimension(:,:,:), allocatable :: nm_ind !< Neighbour map for nearest neighbours
  integer, dimension(:,:), allocatable :: nmdim_ind !< Dimension of neighbour map for the nearest neighbours
  real(dblprec), dimension(:,:), allocatable :: redcoord_mod
  real(dblprec), dimension(:,:,:), allocatable :: ind_redcoord

  private
  public :: induced_mapping, mc_update_ind_mom, calculate_ind_mom_linear

contains


  !---------------------------------------------------------------------------
  !> @brief
  !> Wrapper routine for the creation of the needed lists and arrays for the
  !> treatment of induced moments
  !>
  !> @author
  !> Jonathan Chico
  !---------------------------------------------------------------------------
  subroutine induced_mapping(Natom,NT,NA,N1,N2,N3,sym,max_no_shells,nn,atype,&
             ind_nlistsize,ind_nlist,do_sortcoup,Nchmax,do_ralloy,Natom_full,&
             atype_ch,acellnumb,C1,C2,C3,Bas,BC1,BC2,BC3,&
             ind_tol,redcoord,ind_mom)

    use NeighbourMap, only : setup_nm

    implicit none
    ! .. Input variables
    integer, intent(in) :: NT  !< Number of types of atoms
    integer, intent(in) :: NA  !< Number of atoms in one cell
    integer, intent(in) :: N1  !< Number of cell repetitions in x direction
    integer, intent(in) :: N2  !< Number of cell repetitions in y direction
    integer, intent(in) :: N3  !< Number of cell repetitions in z direction
    integer, intent(in) :: sym !< Symmetry of system (0-3)
    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
    integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
    integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
    integer, intent(in) :: max_no_shells !< Calculated maximum of shells for exchange
    integer, dimension(:), intent(in):: nn !< Number of neighbour shells
    integer, dimension(Natom), intent(in) :: atype !< Type of atom
    integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual site of atom for dilute system
    integer, dimension(Natom_full), optional ,intent(in) :: acellnumb !< List for translating atom no. in full cell to actual cell
    integer, dimension(NA,Nchmax), intent(in) :: ind_mom !< Indication of whether a given moment is induced (1) or fixed (0)
    character :: do_sortcoup !< Sort the entries of ncoup arrays (Y/N)
    character(len=1), intent(in) :: BC1  !< Boundary conditions in x-direction
    character(len=1), intent(in) :: BC2  !< Boundary conditions in y-direction
    character(len=1), intent(in) :: BC3  !< Boundary conditions in z-direction
    real(dblprec), intent(in) :: ind_tol !< Value for the tolerance between neighbouring shells
    real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
    real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
    real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
    real(dblprec), dimension(NT,max_no_shells), intent(in) :: redcoord
    ! .. Inout variables
    integer, dimension(:), allocatable, intent(inout) :: ind_nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, dimension(:,:), allocatable, intent(inout) :: ind_nlist !< Neighbour list for Heisenberg exchange couplings
    real(dblprec), dimension(3,NA), intent(inout) :: Bas !< Coordinates for basis atoms
    ! .. Local variables
    integer :: i_stat, i_all
    integer :: max_no_neigh !< Calculated maximum of neighbours for exchange
    integer :: max_no_equiv !< Calculated maximum of neighbours in one shell for exchange

    call ErrorHandling_missing('Induced moments')

  end subroutine induced_mapping

  !---------------------------------------------------------------------------
  !> @brief
  !> Creation of the bonding vector between neighbours for ony the nearest neighbours
  !>
  !> @author
  !> Jonathan Chico
  !---------------------------------------------------------------------------
  subroutine setup_induced_nearest(NT,max_no_shells,ind_tol,redcoord)

    implicit none

    ! .. Input variables
    integer, intent(in) :: NT !< Number of types of atoms
    integer, intent(in) :: max_no_shells !< Calculated maximum of shells for exchange
    real(dblprec), intent(in) :: ind_tol !< Value for the tolerance between neighbouring shells
    real(dblprec), dimension(NT,max_no_shells,3), intent(in) :: redcoord !< Coordinates for Heisenberg exchange couplings

    ! .. Local variables
    integer :: IT, shells, i_stat,i_all
    integer :: temp_counter, temp_max_counter
    real(dblprec) :: temp_min,diff
    real(dblprec), dimension(NT) :: redcoord_dist
    real(dblprec), dimension(NT,3) :: redcoord_tmp

    call ErrorHandling_missing('Induced moments')

  end subroutine  setup_induced_nearest

  !---------------------------------------------------------------------------
  !> @brief
  !> Setup auxiliarly neighbour list to treat induced moments
  !> it creates and the list for the induced atom neighbours
  !> @author
  !> Jonathan Chico
  !---------------------------------------------------------------------------
  subroutine setup_induced_list(Natom,NT,NA,atype,max_no_neigh,max_no_equiv,max_no_shells,&
             do_sortcoup,Natom_full,atype_ch,Nchmax,do_ralloy,nm_ind,nn,nmdim_ind,ind_nlistsize,&
             ind_nlist,ind_mom)

    use Sorting, only : MergeSortIR
    !
    implicit none
    ! .. Input variables
    integer, intent(in) :: NT  !< Number of types of atoms
    integer, intent(in) :: NA  !< Number of atoms in one cell
    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
    integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
    integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
    integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
    integer, intent(in) :: max_no_equiv !< Calculated maximum of neighbours in one shell for exchange
    integer, intent(in) :: max_no_shells !< Calculated maximum of shells for exchange
    integer, dimension(NT), intent(in):: nn !< Number of neighbour shells
    integer, dimension(Natom), intent(in) :: atype !< Type of atom
    integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual site of atom for dilute system
    integer, dimension(NA,Nchmax), intent(in) :: ind_mom !< Indication of whether a given moment is induced (1) or fixed (0)
    integer, dimension(max_no_shells,Natom), intent(in) :: nmdim_ind !< Dimension of neighbour map
    integer, dimension(Natom,max_no_shells,max_no_equiv), intent(in) :: nm_ind !< Neighbour map
    character :: do_sortcoup !< Sort the entries of ncoup arrays (Y/N)
    ! .. Output variables
    integer, dimension(Natom), intent(out):: ind_nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, dimension(max_no_neigh,Natom), intent(out) :: ind_nlist !< Neighbour list for Heisenberg exchange couplings
    ! .. Local variables
    integer :: i, j, k, l
    integer :: jchem, ichem
    integer :: tempn
    integer :: ncount,ncounter
    logical :: exis

    call ErrorHandling_missing('Induced moments')

  end subroutine setup_induced_list

  !---------------------------------------------------------------------------
  !> @brief
  !> Calculation of the weighting/renormalization factor for the moments from
  !> linear response, this array will be then used to update the magnitude of the
  !> induced moments.
  !> @author
  !> Jonathan Chico
  !---------------------------------------------------------------------------
  subroutine calculate_ind_mom_linear(NA,Natom,Nchmax,Mensemble,do_ralloy,Natom_full,&
             max_no_neigh,atype,ind_nlistsize,atype_ch,ind_nlist,ind_mom,mmom,&
             emomM,sus_ind)

    implicit none

    ! .. Input variables
    integer, intent(in) :: NA  !< Number of atoms in one cell
    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
    integer, intent(in) :: Mensemble !< Number of ensembles
    integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
    integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
    integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
    integer, dimension(Natom), intent(in) :: atype !< Type of atom
    integer, dimension(Natom), intent(in):: ind_nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual site of atom for dilute system
    integer, dimension(NA,Nchmax), intent(in) :: ind_mom !< Indication of whether a given moment is induced (1) or fixed (0)
    integer, dimension(max_no_neigh,Natom), intent(in) :: ind_nlist !< Neighbour list for Heisenberg exchange couplings
    real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
    ! .. Output variables
    real(dblprec), dimension(Natom), intent(out) :: sus_ind !< Scaling factor for the magneitc moment of the induced moments
    ! .. Local variables
    integer :: ichem,iatom, ens, ii
    real(dblprec) :: mod_sum

    call ErrorHandling_missing('Induced moments')

  end subroutine calculate_ind_mom_linear

  !---------------------------------------------------------------------------
  !> @brief
  !> Metropolis algorithm MonteCarlo update for a system with both fixed and induced
  !> magnetic moments. The approach is based in the work presented by Polesya et al.
  !> in Phys. Rev. B 82, 214409.
  !> @author
  !> Jonathan Chico
  !---------------------------------------------------------------------------
  subroutine mc_update_ind_mom(NA,Natom,Nchmax,Mensemble,do_ralloy,Natom_full,&
             iflip_a,atype,atype_ch,temperature,temprescale,mode,max_no_neigh,&
             nlistsize,nlist,ncoup,ncoupD,conf_num,exc_inter,do_dm,max_no_dmneigh,&
             dmlistsize,dmlist,dm_vect,do_pd,nn_pd_tot,pdlistsize,pdlist,pd_vect,&
             do_biqdm,nn_biqdm_tot,biqdmlistsize,biqdmlist,biqdm_vect,do_bq,nn_bq_tot,&
             bqlistsize,bqlist,j_bq,taniso,taniso_diff,eaniso,eaniso_diff,kaniso,&
             kaniso_diff,sb,sb_diff,mult_axis,mmom,emomM,emom,extfield,do_dip,Qdip,&
             ind_nlistsize,ind_nlist,ind_mom,sus_ind,do_lsf,lsf_metric,ind_mom_flag)
    !
    use RandomNumbers, only: rng_uniform,rng_gaussian
    use montecarlo_common

    !.. Input variables
    implicit none

    integer, intent(in) :: NA
    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Nchmax
    integer, intent(in) :: Mensemble !< Number of ensembles
    integer, intent(in) :: do_ralloy
    integer, intent(in) :: Natom_full
    integer, dimension(Natom),intent(in) :: iflip_a !< Flipping pattern
    integer, dimension(Natom), intent(in) :: atype !< Type of atom
    integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual site of atom for dilute system
    real(dblprec), intent(in) :: temperature !< Temperature
    real(dblprec), intent(in) :: temprescale !< Temperature rescaling according to QHB
    character(len=1) :: mode !< Simulation mode (M=MC, H=MC Heat Bath)
    ! LSF variables
    integer, intent(in) :: conf_num   !< Number of configurations for LSF
    integer,intent(in) :: lsf_metric !< LSF metric or phase space measure
    character(len=1), intent(in) :: do_lsf     !< Including LSF energy
    character(len=1), intent(in) :: exc_inter !< Interpolation of Jij between FM/DLM (Y/N)
    ! Heisenberg exchange variables
    integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
    integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
    real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
    real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoupD !< Heisenberg exchange couplings (DLM)
    ! DMI  variables
    integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
    integer, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
    integer, dimension(Natom),intent(in) :: dmlistsize !< Size of neighbour list for DM
    integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist   !< List of neighbours for DM
    real(dblprec), dimension(3,max_no_dmneigh,Natom), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
    ! PD interactions variables
    integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
    integer, intent(in) :: nn_pd_tot !< Calculated number of neighbours with PD interactions
    integer, dimension(Natom),intent(in) :: pdlistsize !< Size of neighbour list for PD
    integer, dimension(nn_pd_tot,Natom), intent(in) :: pdlist   !< List of neighbours for PD
    real(dblprec), dimension(6,nn_pd_tot,Natom), intent(in) :: pd_vect !< Pseudo-Dipolar exchange vector
    ! BIQDM variables
    integer, intent(in) :: do_biqdm   !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
    integer, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
    integer, dimension(Natom),intent(in) :: biqdmlistsize !< Size of neighbour list for BIQDM
    integer, dimension(nn_biqdm_tot,Natom), intent(in) :: biqdmlist   !< List of neighbours for BIQDM
    real(dblprec), dimension(1,nn_biqdm_tot,Natom), intent(in) :: biqdm_vect !< BIQDM exchange coupling
    ! BQ variables
    integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
    integer, intent(in) :: nn_bq_tot !< Calculated number of neighbours with BQ interactions
    integer, dimension(Natom),intent(in) :: bqlistsize !< Size of neighbour list for BQ
    integer, dimension(nn_bq_tot,Natom), intent(in) :: bqlist   !< List of neighbours for BQ
    real(dblprec), dimension(nn_bq_tot,Natom), intent(in) :: j_bq !< Biquadratic exchange couplings
    ! Anisotropy variables
    integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
    integer, dimension(Natom),intent(in) :: taniso_diff !< Type of anisotropy (0-2)
    real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
    real(dblprec), dimension(3,Natom), intent(in) :: eaniso_diff !< Unit anisotropy vector
    real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
    real(dblprec), dimension(2,Natom), intent(in) :: kaniso_diff !< Anisotropy constant
    real(dblprec), dimension(Natom), intent(in) :: sb!< Ratio between the Anisotropy constants
    real(dblprec), dimension(Natom), intent(in) :: sb_diff!< Ratio between the Anisotropy constants
    character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
    ! Moments variables
    real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
    ! External magnetic fields
    real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
    ! Dipolar interaction variables
    integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
    real(dblprec), dimension(3,3,Natom,Natom), intent(in), optional :: Qdip !< Matrix for dipole-dipole interaction
    ! Induced moment variables
    integer, dimension(Natom),intent(in) :: ind_nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, dimension(NA,Nchmax), intent(in) :: ind_mom !< Indication of whether a given moment is induced (1) or fixed (0)
    integer, dimension(max_no_neigh,Natom), intent(in) :: ind_nlist !< Neighbour list for Heisenberg exchange couplings
    real(dblprec), dimension(Natom), intent(in) :: sus_ind
    character(len=1), intent(in) :: ind_mom_flag

    !.. Local scalars
    !
    integer :: iatom,k,ichem
    real(dblprec) :: de !< New trial magnitude of moment
    real(dblprec),dimension(1):: fran
    !
    !.. Local arrays
    !
    real(dblprec), dimension(3) :: newmom
    real(dblprec), dimension(natom) :: flipprob_a ,newmmom_a
    real(dblprec), dimension(3,natom) :: newmom_a

    call ErrorHandling_missing('Induced moments')

  end subroutine mc_update_ind_mom

  !---------------------------------------------------------------------------
  !> @brief
  !> Calculate the total energy of a single fixed magnetic moment plus including the
  !> indirect energy contributions of the induced magnetic moments
  !>
  !> @author
  !> Jonathan Chico
  !---------------------------------------------------------------------------
  subroutine calculate_energy_wIND(k,NA,Natom,Nchmax,do_ralloy,Mensemble,Natom_full,atype,&
             atype_ch,max_no_neigh,nlistsize,nlist,ncoup,ncoupD,conf_num,exc_inter,&
             do_dm,max_no_dmneigh,dmlistsize,dmlist,dm_vect,do_pd,nn_pd_tot,&
             pdlistsize,pdlist,pd_vect,do_biqdm,nn_biqdm_tot,biqdmlistsize,biqdmlist,biqdm_vect,&
             do_bq,nn_bq_tot,bqlistsize,bqlist,j_bq,taniso,taniso_diff,eaniso,eaniso_diff,kaniso,kaniso_diff,&
             sb,sb_diff,mult_axis,iflip,newmom,mmom,emomM,emom,extfield,do_dip,Qdip,ind_nlistsize,&
             ind_nlist,ind_mom,sus_ind,de)

     use Constants, only : mub

     !.. Implicit declarations
     implicit none

     ! System variables
     integer, intent(in) :: k !< Current ensemble
     integer, intent(in) :: NA
     integer, intent(in) :: Natom !< Number of atoms in system
     integer, intent(in) :: Nchmax
     integer, intent(in) :: do_ralloy
     integer, intent(in) :: Mensemble !< Number of ensembles
     integer, intent(in) :: Natom_full
     integer, dimension(Natom), intent(in) :: atype !< Type of atom
     integer, dimension(Natom_full), intent(in) :: atype_ch
     ! LSF variables
     integer, intent(in) :: conf_num   !< Number of configurations for LSF
     character(len=1), intent(in) :: exc_inter !> Interpolation of Jij between FM/DLM (Y/N)
     ! Heisenberg exchange variables
     integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
     integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
     integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
     real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
     real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoupD !< Heisenberg exchange couplings (DLM)
     ! DMI  variables
     integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
     integer, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
     integer, dimension(Natom),intent(in) :: dmlistsize !< Size of neighbour list for DM
     integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist   !< List of neighbours for DM
     real(dblprec), dimension(3,max_no_dmneigh,Natom), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
     ! PD interactions variables
     integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
     integer, intent(in) :: nn_pd_tot !< Calculated number of neighbours with PD interactions
     integer, dimension(Natom),intent(in) :: pdlistsize !< Size of neighbour list for PD
     integer, dimension(nn_pd_tot,Natom), intent(in) :: pdlist   !< List of neighbours for PD
     real(dblprec), dimension(6,nn_pd_tot,Natom), intent(in) :: pd_vect !< Pseudo-Dipolar exchange vector
     ! BIQDM variables
     integer, intent(in) :: do_biqdm   !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
     integer, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
     integer, dimension(Natom),intent(in) :: biqdmlistsize !< Size of neighbour list for BIQDM
     integer, dimension(nn_biqdm_tot,Natom), intent(in) :: biqdmlist   !< List of neighbours for BIQDM
     real(dblprec), dimension(1,nn_biqdm_tot,Natom), intent(in) :: biqdm_vect !< BIQDM exchange coupling
     ! BQ variables
     integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
     integer, intent(in) :: nn_bq_tot !< Calculated number of neighbours with BQ interactions
     integer, dimension(Natom),intent(in) :: bqlistsize !< Size of neighbour list for BQ
     integer, dimension(nn_bq_tot,Natom), intent(in) :: bqlist   !< List of neighbours for BQ
     real(dblprec), dimension(nn_bq_tot,Natom), intent(in) :: j_bq !< Biquadratic exchange couplings
     ! Anisotropy variables
     integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
     integer, dimension(Natom),intent(in) :: taniso_diff !< Type of anisotropy (0-2)
     real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
     real(dblprec), dimension(3,Natom), intent(in) :: eaniso_diff !< Unit anisotropy vector
     real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
     real(dblprec), dimension(2,Natom), intent(in) :: kaniso_diff !< Anisotropy constant
     real(dblprec), dimension(Natom), intent(in) :: sb!< Ratio between the Anisotropy constants
     real(dblprec), dimension(Natom), intent(in) :: sb_diff!< Ratio between the Anisotropy constants
     character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
     ! Moments variables
     integer, intent(in) :: iflip !< Atom to flip spin for
     real(dblprec), dimension(3), intent(in) :: newmom !< New trial moment
     real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
     real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
     real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
     ! External magnetic fields
     real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
     ! Dipolar interaction variables
     integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
     real(dblprec), dimension(3,3,Natom,Natom), intent(in), optional :: Qdip !< Matrix for dipole-dipole interaction
     ! Induced moment variables
     integer, dimension(Natom),intent(in) :: ind_nlistsize !< Size of neighbour list for Heisenberg exchange couplings
     integer, dimension(NA,Nchmax), intent(in) :: ind_mom !< Indication of whether a given moment is induced (1) or fixed (0)
     integer, dimension(max_no_neigh,Natom), intent(in) :: ind_nlist !< Neighbour list for Heisenberg exchange couplings
     real(dblprec), dimension(Natom), intent(in) :: sus_ind

     ! Output variables
     real(dblprec), intent(out):: de  !< Energy difference

     !.. Local scalars
     integer :: j,fix,curr_fix,neighchem,neigh_test
     real(dblprec) :: tt, tta, ttb, e_c, e_t
     real(dblprec) :: bqmdot, excscale
     real(dblprec) :: aw1,aw2

     !.. Local arrays
     real(dblprec), dimension(3) :: beff_t, trialmom,ave_mom

    call ErrorHandling_missing('Induced moments')
     return
   end subroutine calculate_energy_wIND

end module InducedMoments
