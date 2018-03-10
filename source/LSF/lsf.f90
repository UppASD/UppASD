!-------------------------------------------------------------------------------
! MODULE: LSF
!> @brief
!> Collection of LSF routines including Interpolation to obtain the moment size dependent LSF energy and exchange parameters
!> @author
!> Fan Pan and Lars Bergqvist
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License.
!! See http://www.gnu.org/copyleft/gpl.txt
!-------------------------------------------------------------------------------
module LSF
   use Parameters
   use Profiling
   use inputdata, only : ammom_inp
   use ChemicalData, only : achtype
   use ErrorHandling

   implicit none

   real(dblprec), allocatable, dimension(:,:) :: ammom_int,mom_int,LSF_energy            !< trimed variable from ammom_inp
   real(dblprec), allocatable, dimension(:) :: ammom_llim, ammom_hlim   !< the the max/min value in ammom
   integer, allocatable, dimension(:) :: ngrids            !< number of grids for each components (from modified momfile)

   private :: ngrids, ammom_llim, ammom_hlim,ammom_int,LSF_energy,mom_int

   public :: allocate_lsfdata, LSF_datareshape,mc_update_LSF, totalenergy_LSF, read_LSF

contains
   !---------------------------------------------------------------------------
   !> @brief
   !! Driver routine for Monte Carlo LSF
   !> @author
   !! Lars Bergqvist
   !--------------------------------------------------------------------------
   subroutine mc_update_LSF(Natom,Nchmax,Mensemble,max_no_neigh, conf_num, ncoup, ncoupD,do_lsf, nlist, nlistsize , &
         taniso, eaniso, kaniso,sb,emomM, emom, mmom, temperature, temprescale, extfield, &
         mult_axis, taniso_diff, eaniso_diff, kaniso_diff,sb_diff,mode,fs_nlist, fs_nlistsize, &
         nind,lsf_interpolate,lsf_field,lsf_window,lsf_metric,exc_inter,iflip_a,ind_nlistsize,&
         ind_nlist,sus_ind,ind_mom_flag)
      !
      use RandomNumbers, only: rng_uniform,rng_uniformP, rng_gaussian, use_vsl
      use montecarlo_common

      !.. Implicit declarations
      implicit none
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Nchmax  !< Number of chemical type

      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, intent(in) :: conf_num  !< number of configurations for LSF
      real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
      real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoupD !< Heisenberg exchange couplings (DLM)
      character(len=1), intent(in)  ::  do_lsf     !< Including LSF energy
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
      integer, dimension(Natom),intent(in) :: taniso_diff !< Type of anisotropy (0-2)
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso_diff !< Unit anisotropy vector
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso_diff !< Anisotropy constant
      real(dblprec), dimension(Natom), intent(in) :: sb !< Ratio between anisotropy constants
      real(dblprec), dimension(Natom), intent(in) :: sb_diff !< Ratio between anisotropy constants
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
      real(dblprec), intent(in) :: temperature !< Temperature
      real(dblprec), intent(in) :: temprescale !< Temperature rescaling according to QHB
      real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
      character(len=1), intent(in) :: mult_axis !< Multiple uniaxial anisotropies
      character(len=1) :: mode !< Simulation mode (M=MC, H=MC Heat Bath)
      integer, dimension(:,:), intent(in) :: fs_nlist !< First shell Neighbouring list for centered atom
      integer, dimension(:), intent(in) :: fs_nlistsize  !< Size of first shell neighbouring list for centered atom
      integer, dimension(:,:), intent(in) :: nind !< index of firstshell-neighbour-list corresponds to neighbour-list
      character(len=1), intent(in)  ::  lsf_interpolate     !< Interpolate LSF or not
      character(len=1), intent(in)  ::  lsf_field           !< LSF field contribution (Local/Total)
      real(dblprec),intent(in) :: lsf_window                !< Range of moment variation in LSF
      integer,intent(in) :: lsf_metric                !< LSF metric in phase space integration (1=Murata-Doniach,2=Jacobian)
      character(len=1), intent(in)  ::  exc_inter           !< Exchange interpolation between FM/DLM (Y/N)
      integer,dimension(natom),intent(in) :: iflip_a !< Flipping pattern

      integer, dimension(Natom),intent(in) :: ind_nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(in) :: ind_nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(Natom), intent(in) :: sus_ind
      character(len=1), intent(in) :: ind_mom_flag
      !.. Local scalars
      !
      integer :: i,ip,k
      real(dblprec) :: de, newmmom !< New trial magnitude of moment
      !
      !.. Local arrays
      !
      integer, dimension(Natom,mensemble) :: iin
      real(dblprec), dimension(3) :: newmom !< New trial moment
      real(dblprec), dimension(3) :: totfield  !<Total effective field acting on each moment

      real(dblprec), dimension(Natom,mensemble) :: rn, newmmom_a,flipprob_a
      real(dblprec),dimension(3,natom,mensemble) :: newmom_a

      call ErrorHandling_missing('Local spin fluctuations')

   end subroutine mc_update_LSF
   !---------------------------------------------------------------------------
   !> @brief
   !! Metropolis Monte Carlo inclsuing LSF
   !> @author
   !! Lars Bergqvist and Fan Pan
   !--------------------------------------------------------------------------
   subroutine calculate_energy_wLSF(Natom,Nchmax, Mensemble, max_no_neigh, conf_num,ncoup,ncoupD,nlist, nlistsize, fs_nlist, &
         fs_nlistsize,nind, taniso, eaniso, kaniso,sb,emomM, emom, mmom, iflip, newmom,newmmom, extfield, de, &
         k,lsf_interpolate,lsf_field,exc_inter)

      use Constants, only : mub,mry
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Nchmax  !< Number of chemical type
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, intent(in) :: conf_num   !< Number of configurations for LSF
      real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
      real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoupD !< Heisenberg exchange couplings (DLM)
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist    !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom),intent(in) :: nlistsize              !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(:,:), intent(in) :: fs_nlist !< First shell Neighbouring list for centered atom
      integer, dimension(:), intent(in) :: fs_nlistsize          !< Size of first shell neighbouring list for centered atom
      integer, dimension(:,:), intent(in) :: nind !< index of firstshell-neighbour-list corresponds to neighbour-list
      integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
      real(dblprec), dimension(Natom), intent(in) :: sb!< Ratio between the Anisotropy constants
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      integer, intent(in) :: iflip !< Atom to flip spin for
      real(dblprec), dimension(3), intent(in) :: newmom !< New trial unit moment
      real(dblprec), intent(in) :: newmmom !< New trial magnitude of moment
      real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
      real(dblprec), intent(out):: de  !< Energy difference
      integer, intent(in) :: k !< Current ensemble
      character(len=1),intent(in) :: lsf_interpolate !< Interpolate LSF or not
      character(len=1),intent(in) :: lsf_field       !< LSF field contribution (Local/Total)
      character(len=1),intent(in) :: exc_inter    !< Rescaling of exchange interactions (Y/N)

      real(dblprec) :: aw1,aw2 !,tmp1,tmp2,tmp3,tmp4
      real(dblprec),dimension(2) :: lsf_c, lsf_t   !< lsf energy for a single site
      real(dblprec),dimension(Nchmax) :: nbsum     !< sum of the moments in all neighbours
      real(dblprec),dimension(3) :: nbsumfield     !< sum of the moments in all neighbours
      real(dblprec), dimension(Nchmax) :: ammom_inter  !< moments on grids of  each of the ch_type
      real(dblprec), dimension(max_no_neigh) :: ncoup_c,ncoup_t  !<interpolated ncoup only at flipping site
      real(dblprec), dimension(max_no_neigh) :: fs_ncoup_c, fs_ncoup_t  !< first-shell  ncoup  on  flipping site
      !.. Local scalars
      integer :: inttype              !< (0/1) pick the nearest grids or do interpolation
      integer :: i, j
      real(dblprec) :: tt, tta, ttb, e_c, e_t, fc, excscale

      !.. Local arrays
      integer(dblprec),dimension(nchmax) :: counter     ! counting for every Nchtype

      real(dblprec), dimension(3) :: beff_c,beff_t, trialmom
      real(dblprec), dimension(3) :: fs_beff_c, fs_beff_t  !< beff generated from first shell

      de=0.0d0
      call ErrorHandling_missing('Local spin fluctuations')


   end subroutine calculate_energy_wLSF

   !---------------------------------------------------------------------------
   !> @brief
   !> Vary the magnitude of the moment
   !> @author
   !! Lars Bergqvist and Fan Pan
   !---------------------------------------------------------------------------
   subroutine vary_moment_magnitude(Nchmax,ip,oldmmom,newmmom,ammom_hlim,ammom_llim,wind,temperature)

      use RandomNumbers, only : rng_gaussian,rng_uniform
      use Constants

      implicit none

      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, intent(in) :: ip  !< site number
      real(dblprec), intent(in) :: oldmmom !< chosen magnitude of moment to be varied
      real(dblprec), intent(in), dimension(Nchmax) :: ammom_hlim,ammom_llim
      real(dblprec), intent(out) :: newmmom !< New trial magnitude of moment
      real(dblprec), intent(in) :: wind !< window for changing moment
      real(dblprec), intent(in):: temperature !< Temperature

      real(dblprec) :: locwind,delta
      real(dblprec), dimension(1) :: gasrun,tempmom !< Gaussian RN
      integer :: check

      newmmom=0.0d0
      call ErrorHandling_missing('Local spin fluctuations')

   end subroutine vary_moment_magnitude


   !---------------------------------------------------------------------------
   !> @brief
   !> Extract total energy in LSF
   !> @author
   !! Lars Bergqvist
   !> @todo
   !> Reinstate site resolved energies
   !---------------------------------------------------------------------------
   subroutine totalenergy_LSF(Natom, Nchmax, Mensemble, conf_num, emom, emomM, mmom,simid, &
         plotenergy, mstep, extfield,eenergy, aenergy, fenergy, lsfenergy, &
         max_no_neigh, nlistsize, nlist, ncoup, ncoupD,exc_inter, &
         taniso, eaniso, kaniso,sb,do_lsf,fs_nlist,fs_nlistsize,nind,inttype,lsf_field,Temp)
      use Constants

      implicit none


      integer, intent(in) :: mstep   !< Current simulation step
      integer, intent(in) :: Natom   !< Number of atoms in system
      integer, intent(in) :: Nchmax  !< Number of chemical type
      integer, intent(in) :: Mensemble    !< Number of ensembles
      integer, intent(in) :: conf_num  !< number of configurations for LSF
      integer, intent(in) :: plotenergy   !< Calculate and plot energy (0/1)
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(Natom),intent(in) :: taniso         !< Type of anisotropy (0-2)
      integer, dimension(Natom), intent(in) :: nlistsize     !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist       !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(Natom),intent(in) :: sb !< Ratio between Cubic and Uniaxial anisotropy
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso      !< Anisotropy constant
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom     !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom    !< Current magnetic moment
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: extfield !< External magnetic field
      real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup  !< Heisenberg exchange couplings
      real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoupD !< Heisenberg exchange couplings (DLM)
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in)  ::  do_lsf     !< Including LSF energy
      character(len=1),intent(in) :: exc_inter !< Interpolation of Jij between FM/DLM
      integer, dimension(:,:), intent(in) :: fs_nlist !< First shell Neighbouring list for centered atom
      integer, dimension(:), intent(in) :: fs_nlistsize  !< Size of first shell neighbouring list for centered atom
      integer, dimension(:,:), intent(in) :: nind !< index of firstshell-neighbour-list corresponds to neighbour-list
      integer, intent(in) :: inttype !< Interpolatation type in LSF
      character(len=1), intent(in) :: lsf_field       !< LSF field contribution (Local/Total)
      real(dblprec), intent(in) :: Temp               !< Temperature

      !.. Subroutine output
      real(dblprec), dimension(Mensemble), intent(out) :: aenergy !< Anisotropy energy
      real(dblprec), dimension(Mensemble), intent(out) :: eenergy !< Total exchange (transversal) energy
      real(dblprec), dimension(Mensemble), intent(out) :: fenergy !< Total external energy
      real(dblprec), dimension(Mensemble), intent(out) :: lsfenergy !< Total LSF (longitudinal) energy

      !...Local variables
      integer :: i, j, k
      real(dblprec) :: fcinv,fc
      real(dblprec) :: excscale !< Interpolation parameter FM/DLM
      real(dblprec) :: ieenergy, exptemp
      real(dblprec) :: xu1,xu2


      !...Local arrays
      real(dblprec), dimension(Mensemble) :: tt
      real(dblprec), dimension(3,Mensemble) :: ttv
      real(dblprec), dimension(Mensemble) :: tempk1, tempk2
      real(dblprec), dimension(Mensemble) :: aeatom
      real(dblprec),dimension(2) :: lsfE
      real(dblprec),dimension(Nchmax) :: nbsum     !< sum of the moments in all neighbours
      real(dblprec),dimension(3) :: nbsumfield
      real(dblprec), dimension(Nchmax) :: ammom_inter  !< moments on grids of  each of the ch_type
      real(dblprec), dimension(max_no_neigh) :: ncoup_c  !<interpolated ncoup only at flipping site
      real(dblprec), dimension(max_no_neigh) :: fs_ncoup  !< first-shell  ncoup  on  flipping site
      real(dblprec), dimension(3) :: fs_beff
      integer(dblprec),dimension(nchmax) :: counter     ! counting for every Nchtype
      !    real(dblprec),dimension(Natom,Nchmax,Mensemble) :: xiS,xiJ     !< Statistical weights, standard and Jacobian
      !    real(dblprec),dimension(Nchmax,Mensemble) :: zS,zJ     !< Partition functions with and without Jacobian
      !    real(dblprec),dimension(Nchmax,Mensemble,5) :: mloc      !< Local magnetic moment
      call ErrorHandling_missing('Local spin fluctuations')
      eenergy=0.0d0;fenergy=0.0d0;lsfenergy=0.0d0;aenergy=0.0d0

   end subroutine totalenergy_LSF



   !---------------------------------------------------------------------------
   !> @brief
   !> Interpolation to given values from precalculated grid of lsf energies
   !! and exchange parameters
   !---------------------------------------------------------------------------
   subroutine allocate_lsfdata(NA, Nchmax, conf_num, flag)
      !
      implicit none
      !
      integer, intent(in) :: NA       !< Number of atoms in one cell
      integer, intent(in) :: Nchmax   !< Max number of chemical components on each site in cell
      integer, intent(in) :: conf_num !< Number of configurations for LSF
      integer, intent(in) :: flag     !< Allocate or deallocate (1/-1)
      !
      integer :: i_all, i_stat        !< for memory allocation


      if(flag >0) then
         allocate(LSF_energy(conf_num,2),stat=i_stat)
         call memocc(i_stat,product(shape(LSF_energy))*kind(LSF_energy),'LSF_energy','allocate_lsfdata')
         !
         allocate(ammom_int(Nchmax,conf_num),stat=i_stat)
         call memocc(i_stat,product(shape(ammom_int))*kind(ammom_int),'ammom_int','allocate_lsfdata')
         !
         allocate(ngrids(Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(Nchmax))*kind(Nchmax),'ngrids','allocate_lsfdata')
         !
         allocate(ammom_llim(Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(ammom_llim))*kind(ammom_llim),'ammom_llim','allocate_lsfdata')
         !
         allocate(ammom_hlim(Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(ammom_hlim))*kind(ammom_hlim),'ammom_hlim','allocate_lsfdata')
      else
         i_all=-product(shape(ammom_int))*kind(ammom_int)
         deallocate(ammom_int,stat=i_stat)
         call memocc(i_stat,i_all,'ammom_int','allocate_lsfdata')
         !
         i_all=-product(shape(ngrids))*kind(ngrids)
         deallocate(ngrids,stat=i_stat)
         call memocc(i_stat,i_all,'ngrids','allocate_lsfdata')
         !
         i_all=-product(shape(ammom_llim))*kind(ammom_llim)
         deallocate(ammom_llim,stat=i_stat)
         call memocc(i_stat,i_all,'ammom_llim','allocate_lsfdata')
         !
         i_all=-product(shape(ammom_hlim))*kind(ammom_hlim)
         deallocate(ammom_hlim,stat=i_stat)
         call memocc(i_stat,i_all,'ammom_hlim','allocate_lsfdata')
         !
         i_all=-product(shape(mom_int))*kind(mom_int)
         deallocate(mom_int,stat=i_stat)
         call memocc(i_stat,i_all,'mom_int','allocate_lsfdata')
         !
         !> Deallocate LSF energy
         i_all=-product(shape(LSF_energy))*kind(LSF_energy)
         deallocate(LSF_energy,stat=i_stat)
         call memocc(i_stat,i_all,'LSF_energy','allocate_lsfdata')
      endif
   end subroutine allocate_lsfdata

   !> Reshape of moment input to more suitable format
   !! Trim ammom_inp into ammom_int with dimension (Nchmax,conf_num)
   subroutine LSF_datareshape(NA, Nchmax, conf_num)
      !
      implicit none
      !
      integer, intent(in) :: NA        !< Number of atoms in one cell
      integer, intent(in) :: Nchmax    !< Max number of chemical components on each site in cell
      integer, intent(in) :: conf_num  !< Number of configurations for LSF
      !
      integer :: i_stat  !<, for, memory, allocation
      integer :: i, j, k, tmp, num !< loop index
      logical,dimension(conf_num) :: mask
      real(dblprec),dimension(conf_num,nchmax) :: mint_tmp

      call ErrorHandling_missing('Local spin fluctuations')

   end subroutine LSF_datareshape

   !> Interpolation of LSF energy and exchange interactions to given moment size from existing grid
   subroutine do_interpolation_ncoup_and_lsf(Natom,Nchmax,max_no_neigh,nlist,nlistsize,conf_num,ammom_inter, &
         ncoup,ncoupD,iflip,itp_noup,obj,k,inttype,excscale,exc_inter)

      !
      implicit none

      real(dblprec), intent(in), dimension(:) :: ammom_inter !< moments at grids
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Nchmax  !< Number of chemical type
      integer, intent(in) :: k !< Current ensemble
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, intent(in) :: conf_num   !< Number of configurations for LSF
      real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
      real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoupD !< Heisenberg exchange couplings (DLM)
      integer, intent(in) :: iflip !< Atom to flip spin for
      real(dblprec), dimension(max_no_neigh),intent(out) :: itp_noup  !< Interpolated ncoup
      real(dblprec),dimension(2), intent(out) :: obj !< the thing to be interpolated (possibly LSF_energy)
      integer, intent(in) :: inttype !< (0/1) pick the nearest grids or do interpolation
      real(dblprec),intent(in) :: excscale !< Interpolaton parameter FM/DLM
      character(len=1),intent(in) :: exc_inter !< Exchange rescaling
      !local variables
      integer, dimension(nchmax) :: ind  !< the lower-limit to the moment on grids
      integer, dimension(2**nchmax) :: ac  !< array of conf-num for the interpolating mash
      integer :: i, j, iconf, ii, jj
      real(dblprec),dimension(4) :: temp
      real(dblprec) :: invtemp, invtemp2
      real(dblprec), dimension(max_no_neigh) :: tmp_noup  !< Interpolated ncoup

      itp_noup=0.0d0;obj=0.0d0
      call ErrorHandling_missing('Local spin fluctuations')

   end subroutine do_interpolation_ncoup_and_lsf

   !---------------------------------------------------------------------------
   !> @brief
   !> Reading the LSF energy
   !>
   !> @author
   !> Fan Pan and someguy
   !---------------------------------------------------------------------------
   !
   subroutine read_LSF(lsffile,nchmax,conf_num,mconf)
      !
      implicit none
      !
      character(len=35),intent(in) :: lsffile
      integer,intent(in) :: nchmax
      integer,intent(in) :: conf_num
      integer,intent(inout) :: mconf
      integer :: m, iconf
      integer,dimension(nchmax) :: mconftmp

      call ErrorHandling_missing('Local spin fluctuations')

   end subroutine read_LSF


end module LSF
