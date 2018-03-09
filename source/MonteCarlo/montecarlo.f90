!> Data and routines for MonteCarlo simulations of either Ising or Heisenberg models
!> @author
!> A. Bergman, L. Bergqvist, Fan Pan, Jonathan Chico, J. Hellsvik + ...
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
!> -----------------------------------------------------------------------------------
module MonteCarlo
  use Parameters
  use Profiling
  implicit none


  integer :: mcmstep !< Current Monte Carlo step
  integer, dimension(:), allocatable :: indxb_mcavrg !< Buffer index for averages
  real(dblprec), dimension(:,:), allocatable :: mcmavg_buff !< Buffer for magnetic moments
  !
  integer, dimension(:), allocatable :: iflip_a !< Array of atoms to flip spins for
  !

  private

  public :: indxb_mcavrg, mcmavg_buff ,iflip_a

  ! public subroutines
  public :: mc_evolve, prn_mcmoments, choose_random_atom_x , allocate_mcdata

contains

  !--------------------------------------------------------------------------
  ! DESCRIPTION
  !> @brief
  !> Main driver for Monte Carlo simulations using various implementations.
  !! For Heisenberg spin systems the options are either Metropolis,Heat bath or Wolff cluster.
  !! Based on the Hamiltonian \f$H=-\frac{1}{2} \sum_{ij} J_{ij} \mathbf{m}_i \cdot \mathbf{m}_j \f$.
  !! Ising model and spin-ice models implemented in Metropolis or Loop Algorithm.
  !---------------------------------------------------------------------------------
subroutine mc_evolve(NA,Natom,Nchmax,Mensemble,do_ralloy,Natom_full,atype,atype_ch,temperature,temprescale,&
           mode,conf_num,lsf_metric,fs_nlistsize,fs_nlist,nind,lsf_window,do_lsf,lsf_field,exc_inter,&
           lsf_interpolate,do_jtensor,max_no_neigh,nlistsize,nlist,ncoup,ncoupD,j_tens,do_dm,max_no_dmneigh,&
           dmlistsize,dmlist,dm_vect,do_pd,nn_pd_tot,pdlistsize,pdlist,pd_vect,do_biqdm,nn_biqdm_tot,&
           biqdmlistsize,biqdmlist,biqdm_vect,do_bq,nn_bq_tot,bqlistsize,bqlist,j_bq,taniso,taniso_diff,&
           eaniso,eaniso_diff,kaniso,kaniso_diff,sb,sb_diff,mult_axis,iflip_a,emomM,emom,mmom,ind_nlistsize,&
           ind_nlist,ind_mom,sus_ind,ind_mom_flag,extfield,do_dip,Qdip)
    !
    use RandomNumbers, only: rng_uniform,rng_uniformP,rng_gaussian, use_vsl
    use LSF, only : mc_update_LSF
    use SpinIce , only: mc_update_spinice
    use montecarlo_common
    use InducedMoments, only : mc_update_ind_mom
    !
    implicit none

    !.. Input variables
    ! System variables
    integer, intent(in) :: NA
    integer, intent(in) :: Natom     !< Number of atoms in system
    integer, intent(in) :: Nchmax    !< Number of chemical type
    integer, intent(in) :: Mensemble !< Number of ensembles
    integer, intent(in) :: do_ralloy
    integer, intent(in) :: Natom_full
    integer, dimension(Natom), intent(in) :: atype !< Type of atom
    integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual site of atom for dilute system
    real(dblprec), intent(in) :: temperature !< Temperature
    real(dblprec), intent(in) :: temprescale !< Temperature rescaling according to QHB
    character(len=1) :: mode !< Simulation mode (M=MC, H=MC Heat Bath,D=MC Glauber)
    ! LSF variables
    integer, intent(in) :: conf_num  !< number of configurations for LSF
    integer, intent(in) :: lsf_metric !< LSF metric in phase space integration (1=Murata-Doniach,2=Jacobian)
    integer, dimension(:), intent(in) :: fs_nlistsize  !< Size of first shell neighbouring list for centered atom
    integer, dimension(:,:), intent(in) :: nind        !< index of firstshell-neighbour-list corresponds to neighbour-list
    integer, dimension(:,:), intent(in) :: fs_nlist    !< First shell Neighbouring list for centered atom
    real(dblprec),intent(in) :: lsf_window             !< Range of moment variation in LSF
    character(len=1), intent(in)  ::  do_lsf           !< Including LSF energy
    character(len=1), intent(in)  ::  lsf_field        !< LSF field contribution (Local/Total)
    character(len=1), intent(in)  ::  exc_inter        !< Exchange interpolation between FM/DLM (Y/N)
    character(len=1), intent(in)  ::  lsf_interpolate  !< Interpolate LSF or not
    ! Heisenberg exchange variables
    integer, intent(in) ::  do_jtensor  !<  Use SKKR style exchange tensor (0=off, 1=on, 2=with biquadratic exchange)
    integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
    integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
    real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
    real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoupD !< Heisenberg exchange couplings (DLM)
    real(dblprec), dimension(3,3,max_no_neigh,Natom), intent(in) :: j_tens !< Tensorial exchange couplings (SKKR)
    ! DMI variables
    integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
    integer, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
    integer, dimension(Natom),intent(in) :: dmlistsize !< Size of neighbour list for DM
    integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist   !< List of neighbours for DM
    real(dblprec), dimension(3,max_no_dmneigh,Natom), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
    ! PD variables
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
    real(dblprec), dimension(Natom), intent(in) :: sb !< Ratio between anisotropy constants
    real(dblprec), dimension(Natom), intent(in) :: sb_diff !< Ratio between anisotropy constants
    character(len=1), intent(in) :: mult_axis
    ! Moments variables
    integer, dimension(Natom),intent(in) :: iflip_a !< Flipping pattern
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
    real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
    ! Induced moments variables
    integer, dimension(Natom),intent(in) :: ind_nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, dimension(NA,Nchmax), intent(in) :: ind_mom !< Indication of whether a given moment is induced (1) or fixed (0)
    integer, dimension(max_no_neigh,Natom), intent(in) :: ind_nlist !< Neighbour list for Heisenberg exchange couplings
    real(dblprec), dimension(Natom), intent(in) :: sus_ind
    character(len=1), intent(in) :: ind_mom_flag
    ! External fields variables
    real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
    ! Dipolar interactions variables
    integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
    real(dblprec), dimension(3,3,Natom,Natom), intent(in), optional :: Qdip !< Matrix for dipole-dipole interaction

    !.. Local variables
    integer :: i, k, ip

    real(dblprec) :: de !< Energy difference
    real(dblprec) :: cluster_size
    !.. Local arrays
    integer, dimension(Natom) :: visited_atoms
    real(dblprec), dimension(3) :: newmom !< New trial moment
    real(dblprec), dimension(3) :: totfield  !<Total effective field acting on each moment
    real(dblprec),dimension(natom,mensemble) :: flipprob_a ,newmmom_a
    real(dblprec),dimension(3,natom,mensemble) :: newmom_a

    integer, external :: omp_get_thread_num

    cluster_size=0.0D0
    visited_atoms=0

    if (do_lsf=='Y') then
       call mc_update_LSF(Natom,Nchmax,Mensemble,max_no_neigh, conf_num, ncoup, ncoupD,do_lsf, nlist, nlistsize , &
           taniso, eaniso, kaniso,sb,emomM, emom, mmom, temperature, temprescale, extfield, &
           mult_axis, taniso_diff, eaniso_diff, kaniso_diff,sb_diff,mode,fs_nlist, fs_nlistsize, &
           nind,lsf_interpolate,lsf_field,lsf_window,lsf_metric,exc_inter,iflip_a,ind_nlistsize,&
           ind_nlist,sus_ind,ind_mom_flag)

    else if (ind_mom_flag=='Y') then
       call mc_update_ind_mom(NA,Natom,Nchmax,Mensemble,do_ralloy,Natom_full,iflip_a,atype,atype_ch,&
            temperature,temprescale,mode,max_no_neigh,nlistsize,nlist,ncoup,ncoupD,conf_num,exc_inter,&
            do_dm,max_no_dmneigh,dmlistsize,dmlist,dm_vect,do_pd,nn_pd_tot,pdlistsize,&
            pdlist,pd_vect,do_biqdm,nn_biqdm_tot,biqdmlistsize,biqdmlist,biqdm_vect,&
            do_bq,nn_bq_tot,bqlistsize,bqlist,j_bq,taniso,taniso_diff,eaniso,eaniso_diff,kaniso,&
            kaniso_diff,sb,sb_diff,mult_axis,mmom,emomM,emom,extfield,do_dip,Qdip,&
            ind_nlistsize,ind_nlist,ind_mom,sus_ind,do_lsf,lsf_metric,ind_mom_flag)

    else
    ! Set up trial directions of magnetic moments
        if (mode=='L') then
           call mc_update_spinice(Natom, Mensemble, max_no_neigh, conf_num, ncoup, ncoupD, nlist, nlistsize, &
                do_dm, max_no_dmneigh, dm_vect, dmlist, dmlistsize, do_pd, nn_pd_tot, pd_vect, pdlist, pdlistsize, &
                do_biqdm, nn_biqdm_tot, biqdm_vect, biqdmlist, biqdmlistsize, do_bq, nn_bq_tot, j_bq, bqlist, bqlistsize, &
                taniso, eaniso, kaniso,sb,emomM, emom, mmom, iflip_a, extfield, &
                mult_axis,taniso_diff, eaniso_diff, kaniso_diff,sb_diff, &
                do_dip, Qdip,exc_inter,temperature,temprescale,ind_nlistsize,ind_nlist,sus_ind,ind_mom_flag)
        else
            if (mode=='M'.or.mode=='H'.or.mode=='D') then
                 call rng_uniformP(flipprob_a(:,:),natom*mensemble)

                 if(use_vsl) then
#ifdef VSL
                    !$omp parallel do default(shared),private(i,k,newmom),schedule(auto),collapse(2)
#endif
                    do i=1,Natom
                       do k=1,mensemble
                          call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,temperature,flipprob_a(i,k))
                          newmom_a(1:3,i,k)=newmom(1:3)
                        enddo
                    enddo
#ifdef VSL
                   !$omp end parallel do
#endif
                 else
                    do i=1,Natom
                       do k=1,mensemble
                         call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,temperature,flipprob_a(i,k))
                         newmom_a(1:3,i,k)=newmom(1:3)
                       enddo
                    enddo
                 end if
            else if (mode=='I') then
               ! First Ising random atom
              !$omp parallel do default(shared),private(i,k,newmom),schedule(auto),collapse(2)
                do i=1,Natom
                    do k=1,mensemble
                       call Ising_random_flip(emom,newmom,iflip_a(i),k,Natom,Mensemble)
                       newmom_a(:,i,k)=newmom(:)
                    enddo
                enddo
              !$omp end parallel do
            endif

            call rng_uniformP(flipprob_a(:,:),natom*mensemble)

            ! Calculate energy and flip spin if preferred
            !$omp parallel do default(shared), private(i,k,de,totfield), schedule(auto),collapse(2)
            do i=1, Natom
               do k=1,mensemble
               if(do_jtensor==1) then
                  call calculate_efield_tensor(Natom, Mensemble, max_no_neigh, j_tens, nlist, nlistsize, &
                      do_biqdm, nn_biqdm_tot, biqdm_vect, biqdmlist, biqdmlistsize, do_bq, nn_bq_tot, j_bq, bqlist, bqlistsize, &
                      taniso, eaniso, kaniso,sb,emomM, iflip_a(i),extfield,k, totfield, do_dip, Qdip)

                  call flip_h(Natom, Mensemble, emom, emomM, mmom(iflip_a(i),k),mmom(iflip_a(i),k),iflip_a(i), &
                      temperature,temprescale,k,flipprob_a(i,k),totfield)
               else
                  if (mode=='H') then
                     ! Heat bath algorithm
                     call calculate_efield(Natom, Mensemble, max_no_neigh, conf_num, ncoup, ncoupD, nlist, nlistsize, &
                          do_dm, max_no_dmneigh, dm_vect, dmlist, dmlistsize, &
                          do_pd, nn_pd_tot, pd_vect, pdlist, pdlistsize, &
                          do_biqdm, nn_biqdm_tot, biqdm_vect, biqdmlist, biqdmlistsize, &
                          do_bq, nn_bq_tot, j_bq, bqlist, bqlistsize, &
                          taniso, eaniso, kaniso,sb,emomM, emom, &
                          mult_axis,taniso_diff, eaniso_diff, kaniso_diff,sb_diff,&
                          iflip_a(i),extfield, do_lsf, k, totfield,exc_inter)
!
                     call flip_h(Natom, Mensemble, emom, emomM, mmom(iflip_a(i),k), mmom(iflip_a(i),k), &
                         iflip_a(i),temperature,temprescale, k,flipprob_a(i,k),totfield)
                  else
! Metropolis algorithm, either in Ising or Loop Algorithm form
                     call calculate_energy(Natom, Mensemble, max_no_neigh, conf_num, ncoup, ncoupD, nlist, nlistsize, &
                          do_dm, max_no_dmneigh, dm_vect, dmlist, dmlistsize, do_pd, nn_pd_tot, pd_vect, pdlist, pdlistsize, &
                          do_biqdm, nn_biqdm_tot, biqdm_vect, biqdmlist, biqdmlistsize, do_bq, nn_bq_tot, j_bq, bqlist, bqlistsize, &
                          taniso, eaniso, kaniso,sb,emomM, emom, mmom, iflip_a(i), newmom_a(1:3,iflip_a(i),k), extfield, de, k, &
                          mult_axis,taniso_diff, eaniso_diff, kaniso_diff,sb_diff, &
                          do_dip, Qdip,exc_inter)
                     if(mode=='D') then
                        call flip_g(Natom, Mensemble, emom, emomM, mmom, iflip_a(i),newmom_a(1:3,iflip_a(i),k),newmmom_a(iflip_a(i),k), &
                             de,temperature,temprescale,do_lsf,k,flipprob_a(i,k),lsf_metric,ind_nlistsize,&
                             ind_nlist,ind_mom_flag,max_no_neigh,sus_ind)
                     else
                        call flip_a(Natom, Mensemble, emom, emomM, mmom, iflip_a(i),newmom_a(1:3,iflip_a(i),k),newmmom_a(iflip_a(i),k), &
                          de,temperature,temprescale,do_lsf,k,flipprob_a(i,k),lsf_metric,ind_nlistsize,&
                          ind_nlist,ind_mom_flag,max_no_neigh,sus_ind)
                     endif
                  endif
               endif   !jtensor
               enddo    !ensemble
            enddo     !atom
           !$omp end parallel do
!           enddo    !ensemble
        endif       !mode
    endif         !lsf
   return

  end subroutine mc_evolve



  !> Sets up an array with moments for trial spin flips
  subroutine choose_random_atom_x(Natom,iflip_a)

    use RandomNumbers, only : rng_uniform

    !.. Implicit declarations
    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer,dimension(natom),intent(out) :: iflip_a !< FLip pattern

    real(dblprec) :: dshift(natom)
    integer :: i,ishift,itmp

     do i=1,natom
       iflip_a(i)=i
     end do
    call rng_uniform(dshift,natom)
    do i=1,natom
       ishift=int(dshift(i)*Natom)
       itmp=iflip_a(i)
       iflip_a(i)=iflip_a(mod(ishift,Natom)+1)
       iflip_a(mod(ishift,Natom)+1)=itmp
    end do


  end subroutine choose_random_atom_x



 !> Calculate total field for a single spin
  subroutine calculate_efield(Natom, Mensemble, max_no_neigh, conf_num, ncoup, ncoupD, nlist, nlistsize, &
             do_dm, max_no_dmneigh, dm_vect, dmlist, dmlistsize, &
             do_pd, nn_pd_tot, pd_vect, pdlist, pdlistsize, &
             do_biqdm, nn_biqdm_tot, biqdm_vect, biqdmlist, biqdmlistsize, &
             do_bq, nn_bq_tot, j_bq, bqlist, bqlistsize, &
             taniso, eaniso, kaniso,sb,emomM, emom, &
             mult_axis,taniso_diff, eaniso_diff, kaniso_diff,sb_diff,&
             iflip,extfield, do_lsf, k, totfield,exc_inter)

    !.. Implicit declarations
    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles
    integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
    integer, intent(in) :: conf_num    !< Number of configurations for LSF
    real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
    real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoupD !< Heisenberg exchange couplings (DLM)
    integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
    integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
    integer, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
    real(dblprec), dimension(3,max_no_dmneigh,Natom), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
    integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist   !< List of neighbours for DM
    integer, dimension(Natom),intent(in) :: dmlistsize !< Size of neighbour list for DM
    integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
    integer, intent(in) :: nn_pd_tot !< Calculated number of neighbours with PD interactions
    real(dblprec), dimension(6,nn_pd_tot,Natom), intent(in) :: pd_vect !< Pseudo-Dipolar exchange vector
    integer, dimension(nn_pd_tot,Natom), intent(in) :: pdlist   !< List of neighbours for PD
    integer, dimension(Natom),intent(in) :: pdlistsize !< Size of neighbour list for PD
    integer, intent(in) :: do_biqdm   !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
    integer, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
    real(dblprec), dimension(1,nn_biqdm_tot,Natom), intent(in) :: biqdm_vect !< BIQDM exchange coupling
    integer, dimension(nn_biqdm_tot,Natom), intent(in) :: biqdmlist   !< List of neighbours for BIQDM
    integer, dimension(Natom),intent(in) :: biqdmlistsize !< Size of neighbour list for BIQDM
    integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
    integer, intent(in) :: nn_bq_tot !< Calculated number of neighbours with BQ interactions
    real(dblprec), dimension(nn_bq_tot,Natom), intent(in) :: j_bq !< Biquadratic exchange couplings
    integer, dimension(nn_bq_tot,Natom), intent(in) :: bqlist   !< List of neighbours for BQ
    integer, dimension(Natom),intent(in) :: bqlistsize !< Size of neighbour list for BQ
    integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
    integer, dimension(Natom),intent(in) :: taniso_diff !< Type of anisotropy (0-2)
    real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
    real(dblprec), dimension(3,Natom), intent(in) :: eaniso_diff !< Unit anisotropy vector
    real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
    real(dblprec), dimension(2,Natom), intent(in) :: kaniso_diff !< Anisotropy constant
    real(dblprec), dimension(Natom), intent(in) :: sb!< Ratio between the Anisotropy constants
    real(dblprec), dimension(Natom), intent(in) :: sb_diff!< Ratio between the Anisotropy constants
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
    integer, intent(in) :: iflip !< Atom to flip spin for
    real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
    character(len=1), intent(in)  ::  do_lsf     !< Including LSF energy
    integer, intent(in) :: k !< Current ensemble
    real(dblprec), dimension(3),intent(out) :: totfield  !<Total effective field
    character(len=1), intent(in) :: mult_axis
    character(len=1), intent(in) :: exc_inter !> Interpolation of Jij between FM/DLM (Y/N)

    !.. Local scalars
    integer :: j
    real(dblprec) :: tta, aw1,aw2
    real(dblprec) :: bqmdot,excscale
    real(dblprec) :: sxy, syz, szx

    !.. Executable statements

    ! First calculate effective field
    totfield(:) = 0.d0

    if (exc_inter=='N') then
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
       !$omp simd reduction(+:totfield)
#endif
       do j=1,nlistsize(iflip)
          totfield(:) = totfield(:)+ ncoup(j,iflip,1)*emomM(:,nlist(j,iflip),k)
       end do
    else
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
       !$omp simd reduction(+:totfield)
#endif
       do j=1,nlistsize(iflip)
          totfield=totfield+emom(:,nlist(j,iflip),k)
       enddo
       excscale=sqrt(totfield(1)**2+totfield(2)**2+totfield(3)**2)/nlistsize(iflip)
       totfield=0.d0
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
       !$omp simd reduction(+:totfield)
#endif
       do j=1,nlistsize(iflip)
          totfield(:)=totfield(:)+((excscale*ncoup(j,iflip,1)+(1.d0-excscale)*ncoupD(j,iflip,1)))*emomM(:,nlist(j,iflip),k)
       end do
    endif

    ! Anisotropy
   ! Uniaxial anisotropy  scaled down to match heatbath
    if (taniso(iflip)==1) then
       tta=sum(emomM(:,iflip,k)*eaniso(:,iflip))

       ! K1*(sin theta)^2
       totfield(1:3) = totfield(1:3)  &
            - 1.0d0*kaniso(1,iflip)*tta*eaniso(1:3,iflip) &
                                ! K2*(sin theta)^4
            - 1.0d0*kaniso(2,iflip)*(tta**2)*tta*eaniso(1:3,iflip)
       ! Cubic anisotropy
    elseif (taniso(iflip)==2) then
       ! K1*(sin theta)^2
       totfield(1) = totfield(1) &
            + 2*kaniso(1,iflip)*emomM(1,iflip,k)*(emomM(2,iflip,k)**2+emomM(3,iflip,k)**2) &
            + 2*kaniso(2,iflip)+emomM(1,iflip,k)*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2

       totfield(2) = totfield(2) &
            + 2*kaniso(1,iflip)*emomM(2,iflip,k)*(emomM(3,iflip,k)**2+emomM(1,iflip,k)**2) &
            + 2*kaniso(2,iflip)+emomM(2,iflip,k)*emomM(3,iflip,k)**2*emomM(1,iflip,k)**2

       totfield(3) = totfield(3) &
            + 2*kaniso(1,iflip)*emomM(3,iflip,k)*(emomM(1,iflip,k)**2+emomM(2,iflip,k)**2) &
            + 2*kaniso(2,iflip)+emomM(3,iflip,k)*emomM(1,iflip,k)**2*emomM(2,iflip,k)**2

    endif
    ! When both Cubic and Uniaxial are switched on
    if (taniso(iflip)==7) then

       ! Uniaxial anisotropy
       tta=sum(emomM(:,iflip,k)*eaniso(:,iflip))

       ! K1*(sin theta)^2
       totfield(1:3) = totfield(1:3)  &
            - 1.0d0*kaniso(1,iflip)*tta*eaniso(1:3,iflip) &
                                ! K2*(sin theta)^4
            - 1.0d0*kaniso(2,iflip)*(tta**2)*tta*eaniso(1:3,iflip)
       ! Cubic anisotropy
       ! K1*(sin theta)^2
       aw1=kaniso(1,iflip)*sb(iflip)
       aw2=kaniso(2,iflip)*sb(iflip)
       totfield(1) = totfield(1) &
            + 2*aw1*emomM(1,iflip,k)*(emomM(2,iflip,k)**2+emomM(3,iflip,k)**2) &
            + 2*aw2+emomM(1,iflip,k)*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2

       totfield(2) = totfield(2) &
            + 2*aw1*emomM(2,iflip,k)*(emomM(3,iflip,k)**2+emomM(1,iflip,k)**2) &
            + 2*aw2+emomM(2,iflip,k)*emomM(3,iflip,k)**2*emomM(1,iflip,k)**2

       totfield(3) = totfield(3) &
            + 2*aw1*emomM(3,iflip,k)*(emomM(1,iflip,k)**2+emomM(2,iflip,k)**2) &
            + 2*aw2+emomM(3,iflip,k)*emomM(1,iflip,k)**2*emomM(2,iflip,k)**2
    endif


   if (mult_axis=='Y') then
      ! Uniaxial anisotropy
      if (taniso_diff(iflip)==1) then
         tta=sum(emomM(:,iflip,k)*eaniso_diff(:,iflip))

         ! K1*(sin theta)^2
         totfield(1:3) = totfield(1:3)  &
              - 2.0d0*kaniso_diff(1,iflip)*tta*eaniso_diff(1:3,iflip) &
                                  ! K2*(sin theta)^4
              - 4.0d0*kaniso_diff(2,iflip)*(tta**2)*tta*eaniso_diff(1:3,iflip)
         ! Cubic anisotropy
      elseif (taniso_diff(iflip)==2) then
         ! K1*(sin theta)^2
         totfield(1) = totfield(1) &
              + 2*kaniso_diff(1,iflip)*emomM(1,iflip,k)*(emomM(2,iflip,k)**2+emomM(3,iflip,k)**2) &
              + 2*kaniso_diff(2,iflip)+emomM(1,iflip,k)*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2

         totfield(2) = totfield(2) &
              + 2*kaniso_diff(1,iflip)*emomM(2,iflip,k)*(emomM(3,iflip,k)**2+emomM(1,iflip,k)**2) &
              + 2*kaniso_diff(2,iflip)+emomM(2,iflip,k)*emomM(3,iflip,k)**2*emomM(1,iflip,k)**2

         totfield(3) = totfield(3) &
              + 2*kaniso_diff(1,iflip)*emomM(3,iflip,k)*(emomM(1,iflip,k)**2+emomM(2,iflip,k)**2) &
              + 2*kaniso_diff(2,iflip)+emomM(3,iflip,k)*emomM(1,iflip,k)**2*emomM(2,iflip,k)**2

      endif
      ! When both Cubic and Uniaxial are switched on
      if (taniso_diff(iflip)==7) then

         ! Uniaxial anisotropy
         tta=sum(emomM(:,iflip,k)*eaniso_diff(:,iflip))

         ! K1*(sin theta)^2
         totfield(1:3) = totfield(1:3)  &
              - 2.0d0*kaniso_diff(1,iflip)*tta*eaniso_diff(1:3,iflip) &
                                  ! K2*(sin theta)^4
              - 4.0d0*kaniso_diff(2,iflip)*(tta**2)*tta*eaniso_diff(1:3,iflip)
         ! Cubic anisotropy
         ! K1*(sin theta)^2
         aw1=kaniso_diff(1,iflip)*sb_diff(iflip)
         aw2=kaniso_diff(2,iflip)*sb_diff(iflip)
         totfield(1) = totfield(1) &
              + 2*aw1*emomM(1,iflip,k)*(emomM(2,iflip,k)**2+emomM(3,iflip,k)**2) &
              + 2*aw2+emomM(1,iflip,k)*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2

         totfield(2) = totfield(2) &
              + 2*aw1*emomM(2,iflip,k)*(emomM(3,iflip,k)**2+emomM(1,iflip,k)**2) &
              + 2*aw2+emomM(2,iflip,k)*emomM(3,iflip,k)**2*emomM(1,iflip,k)**2

         totfield(3) = totfield(3) &
              + 2*aw1*emomM(3,iflip,k)*(emomM(1,iflip,k)**2+emomM(2,iflip,k)**2) &
              + 2*aw2+emomM(3,iflip,k)*emomM(1,iflip,k)**2*emomM(2,iflip,k)**2
      endif

   endif

    ! DM interaction
    if(do_dm==1) then
       do j=1,dmlistsize(iflip)
          totfield(1) = totfield(1) - dm_vect(3,j,iflip)*emomM(2,dmlist(j,iflip),k) +&
               dm_vect(2,j,iflip)*emomM(3,dmlist(j,iflip),k)
          totfield(2) = totfield(2) - dm_vect(1,j,iflip)*emomM(3,dmlist(j,iflip),k) +&
               dm_vect(3,j,iflip)*emomM(1,dmlist(j,iflip),k)
          totfield(3) = totfield(3) - dm_vect(2,j,iflip)*emomM(1,dmlist(j,iflip),k) +&
               dm_vect(1,j,iflip)*emomM(2,dmlist(j,iflip),k)
       end do
    end if

    ! PD interaction
    if(do_pd==1) then
       do j=1,pdlistsize(iflip)
          totfield(1) = totfield(1) + pd_vect(1,j,iflip)*emomM(1,pdlist(j,iflip),k) +&
               pd_vect(4,j,iflip)*emomM(2,pdlist(j,iflip),k) +&
               pd_vect(5,j,iflip)*emomM(3,pdlist(j,iflip),k)
          totfield(2) = totfield(2) + pd_vect(4,j,iflip)*emomM(1,pdlist(j,iflip),k) +&
               pd_vect(2,j,iflip)*emomM(2,pdlist(j,iflip),k) +&
               pd_vect(6,j,iflip)*emomM(3,pdlist(j,iflip),k)
          totfield(3) = totfield(3) + pd_vect(5,j,iflip)*emomM(1,pdlist(j,iflip),k) +&
               pd_vect(6,j,iflip)*emomM(2,pdlist(j,iflip),k) +&
               pd_vect(3,j,iflip)*emomM(3,pdlist(j,iflip),k)
       end do
    end if

    ! BIQDM interaction
    if(do_biqdm==1) then
       do j=1,biqdmlistsize(iflip)
          sxy = emomM(1,iflip,k)*emomM(2,biqdmlist(j,iflip),k)-&
               emomM(2,iflip,k)*emomM(1,biqdmlist(j,iflip),k)
          syz = emomM(2,iflip,k)*emomM(3,biqdmlist(j,iflip),k)-&
               emomM(3,iflip,k)*emomM(2,biqdmlist(j,iflip),k)
          szx = emomM(3,iflip,k)*emomM(1,biqdmlist(j,iflip),k)-&
               emomM(1,iflip,k)*emomM(3,biqdmlist(j,iflip),k)
          totfield(1) = totfield(1) + 2.0d0*biqdm_vect(1,j,iflip)*(&
               szx*emomM(3,biqdmlist(j,iflip),k)-&
               sxy*emomM(2,biqdmlist(j,iflip),k))
          totfield(2) = totfield(2) + 2.0d0*biqdm_vect(1,j,iflip)*(&
               sxy*emomM(1,biqdmlist(j,iflip),k)-&
               syz*emomM(3,biqdmlist(j,iflip),k))
          totfield(3) = totfield(3) + 2.0d0*biqdm_vect(1,j,iflip)*(&
               syz*emomM(2,biqdmlist(j,iflip),k)-&
               szx*emomM(1,biqdmlist(j,iflip),k))
       end do
    end if

    ! BQ interaction scaled down to match Heatbath
    if(do_bq==1) then
       do j=1,bqlistsize(iflip)
          ! current spin
          !!!!!++ Lars changed the last emom to emomM
          bqmdot=emomM(1,bqlist(j,iflip),k)*emomM(1,iflip,k)+&
               emomM(2,bqlist(j,iflip),k)*emomM(2,iflip,k)+&
               emomM(3,bqlist(j,iflip),k)*emomM(3,iflip,k)
          totfield(:) = totfield(:) + 1.0d0*j_bq(j,iflip)*bqmdot*emomM(:,bqlist(j,iflip),k)
       end do
    end if

    ! Add static field
    totfield(1:3) = totfield(1:3)+extfield(1:3)
    !
  end subroutine calculate_efield


  !> Calculate total field of a single spin of tensorial exchange
  subroutine calculate_efield_tensor(Natom, Mensemble, max_no_neigh, j_tens, nlist, nlistsize, &
             do_biqdm, nn_biqdm_tot, biqdm_vect, biqdmlist, biqdmlistsize, do_bq, nn_bq_tot, j_bq, bqlist, bqlistsize, &
             taniso, eaniso, kaniso,sb,emomM, iflip,extfield,k, totfield, do_dip, Qdip)

    !.. Implicit declarations
    implicit none

    integer, intent(in) :: Natom      !< Number of atoms in system
    integer, intent(in) :: Mensemble  !< Number of ensembles
    integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
    real(dblprec), dimension(3,3,max_no_neigh,Natom), intent(in) :: j_tens !< Tensor exchange couplings (SKKR)
    integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
    integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, intent(in) :: do_biqdm     !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
    integer, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
    real(dblprec), dimension(1,nn_biqdm_tot,Natom), intent(in) :: biqdm_vect !< BIQDM exchange coupling
    integer, dimension(nn_biqdm_tot,Natom), intent(in) :: biqdmlist   !< List of neighbours for BIQDM
    integer, dimension(Natom),intent(in) :: biqdmlistsize !< Size of neighbour list for BIQDM
    integer, intent(in) :: do_bq     !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
    integer, intent(in) :: nn_bq_tot !< Calculated number of neighbours with BQ interactions
    real(dblprec), dimension(nn_bq_tot,Natom), intent(in) :: j_bq !< Biquadratic exchange couplings
    integer, dimension(nn_bq_tot,Natom), intent(in) :: bqlist   !< List of neighbours for BQ
    integer, dimension(Natom),intent(in) :: bqlistsize !< Size of neighbour list for BQ
    integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
    real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
    real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
    real(dblprec), dimension(Natom), intent(in) :: sb !< Ratio between anisotropy constants
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
    integer, intent(in) :: iflip !< Atom to flip spin for
    real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
    integer, intent(in) :: k !< Current ensemble
    real(dblprec), dimension(3), intent(out) :: totfield !< Total effective field
    real(dblprec) :: aw1,aw2
    integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
    real(dblprec), dimension(3,3,Natom,Natom), intent(in), optional :: Qdip !< Matrix for dipole-dipole interaction

    !.. Local scalars
    integer :: j
    real(dblprec) :: tta
    real(dblprec) :: bqmdot
    real(dblprec) :: sxy, syz, szx

    !.. Local arrays
    real(dblprec), dimension(3) :: befftemp

    !.. Executable statements


    ! First calculate effective field
    befftemp(1) = 0d0
    befftemp(2) = 0d0
    befftemp(3) = 0d0

    ! Exchange
    do j=1,nlistsize(iflip)
       befftemp(1) = befftemp(1)+ j_tens(1,1,j,iflip)*emomM(1,nlist(j,iflip),k) + &
                     j_tens(1,2,j,iflip)*emomM(2,nlist(j,iflip),k) + j_tens(1,3,j,iflip)*emomM(3,nlist(j,iflip),k)
       befftemp(2) = befftemp(2)+ j_tens(2,1,j,iflip)*emomM(1,nlist(j,iflip),k) + &
                     j_tens(2,2,j,iflip)*emomM(2,nlist(j,iflip),k) + j_tens(2,3,j,iflip)*emomM(3,nlist(j,iflip),k)
       befftemp(3) = befftemp(3)+ j_tens(3,1,j,iflip)*emomM(1,nlist(j,iflip),k) + &
                     j_tens(3,2,j,iflip)*emomM(2,nlist(j,iflip),k) + j_tens(3,3,j,iflip)*emomM(3,nlist(j,iflip),k)
    end do

    ! BIQDM interaction
    if(do_biqdm==1) then
       do j=1,biqdmlistsize(iflip)
          sxy = emomM(1,iflip,k)*emomM(2,biqdmlist(j,iflip),k)-&
                emomM(2,iflip,k)*emomM(1,biqdmlist(j,iflip),k)
          syz = emomM(2,iflip,k)*emomM(3,biqdmlist(j,iflip),k)-&
                emomM(3,iflip,k)*emomM(2,biqdmlist(j,iflip),k)
          szx = emomM(3,iflip,k)*emomM(1,biqdmlist(j,iflip),k)-&
                emomM(1,iflip,k)*emomM(3,biqdmlist(j,iflip),k)
          befftemp(1) = befftemp(1) + 2.0d0*biqdm_vect(1,j,iflip)*(&
                        szx*emomM(3,biqdmlist(j,iflip),k)-&
                        sxy*emomM(2,biqdmlist(j,iflip),k))
          befftemp(2) = befftemp(2) + 2.0d0*biqdm_vect(1,j,iflip)*(&
                        sxy*emomM(1,biqdmlist(j,iflip),k)-&
                        syz*emomM(3,biqdmlist(j,iflip),k))
          befftemp(3) = befftemp(3) + 2.0d0*biqdm_vect(1,j,iflip)*(&
                        syz*emomM(2,biqdmlist(j,iflip),k)-&
                        szx*emomM(1,biqdmlist(j,iflip),k))
       end do
    end if

    ! BQ interaction
    if(do_bq==1) then
       do j=1,bqlistsize(iflip)
          bqmdot=emomM(1,bqlist(j,iflip),k)*emomM(1,iflip,k)+&
                 emomM(2,bqlist(j,iflip),k)*emomM(2,iflip,k)+&
                 emomM(3,bqlist(j,iflip),k)*emomM(3,iflip,k)

          !!!!!++ Jonathan used + instead of original -
          !!!!!++ Lars rescaled the 2.0 to 1.0
          befftemp(1) = befftemp(1)- 1.0d0*j_bq(j,iflip)*bqmdot*emomM(1,bqlist(j,iflip),k)
          befftemp(2) = befftemp(2)- 1.0d0*j_bq(j,iflip)*bqmdot*emomM(2,bqlist(j,iflip),k)
          befftemp(3) = befftemp(3)- 1.0d0*j_bq(j,iflip)*bqmdot*emomM(3,bqlist(j,iflip),k)
       end do
    end if

    ! Anisotropy
    ! Uniaxial anisotropy
    if (taniso(iflip)==1) then
       tta=(emomM(1,iflip,k)*eaniso(1,iflip)+emomM(2,iflip,k)*eaniso(2,iflip)+emomM(3,iflip,k)*eaniso(3,iflip))

       !!!!!++ Lars rescaled 2.0 and 4.0 to 1.0
       ! K1*(sin theta)^2
       befftemp(1:3) = befftemp(1:3)  &
            - 1.0d0*kaniso(1,iflip)*&
            tta&
            *eaniso(1:3,iflip) &

                                ! K2*(sin theta)^4
            - 1.0d0*kaniso(2,iflip)*&
            (tta**2)*tta*eaniso(1:3,iflip)

       ! Cubic anisotropy
    elseif (taniso(iflip)==2) then
       ! K1*(sin theta)^2
       befftemp(1) = befftemp(1) &
            + 2*kaniso(1,iflip)*emomM(1,iflip,k)*(emomM(2,iflip,k)**2+emomM(3,iflip,k)**2) &
            + 2*kaniso(2,iflip)+emomM(1,iflip,k)*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2

       befftemp(2) = befftemp(2) &
            + 2*kaniso(1,iflip)*emomM(2,iflip,k)*(emomM(3,iflip,k)**2+emomM(1,iflip,k)**2) &
            + 2*kaniso(2,iflip)+emomM(2,iflip,k)*emomM(3,iflip,k)**2*emomM(1,iflip,k)**2

       befftemp(3) = befftemp(3) &
            + 2*kaniso(1,iflip)*emomM(3,iflip,k)*(emomM(1,iflip,k)**2+emomM(2,iflip,k)**2) &
            + 2*kaniso(2,iflip)+emomM(3,iflip,k)*emomM(1,iflip,k)**2*emomM(2,iflip,k)**2

    endif

    ! When both Cubic and Uniaxial are switched on
    if (taniso(iflip)==7) then

       ! Uniaxial anisotropy
       tta=(emomM(1,iflip,k)*eaniso(1,iflip)+emomM(2,iflip,k)*eaniso(2,iflip)+emomM(3,iflip,k)*eaniso(3,iflip))
       ! K1*(sin theta)^2
       befftemp(1:3) = befftemp(1:3)  &
            - 2.0d0*kaniso(1,iflip)*&
            tta&
            *eaniso(1:3,iflip) &
                                ! K2*(sin theta)^4
            -4.0d0*kaniso(2,iflip)*&
            (tta**2)*tta*eaniso(1:3,iflip)

       ! Cubic anisotropy
       ! K1*(sin theta)^2
       aw1=kaniso(1,iflip)*sb(iflip)
       aw2=kaniso(2,iflip)*sb(iflip)
       befftemp(1) = befftemp(1) &
            + 2*aw1*emomM(1,iflip,k)*(emomM(2,iflip,k)**2+emomM(3,iflip,k)**2) &
            + 2*aw2+emomM(1,iflip,k)*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2

       befftemp(2) = befftemp(2) &
            + 2*aw1*emomM(2,iflip,k)*(emomM(3,iflip,k)**2+emomM(1,iflip,k)**2) &
            + 2*aw2+emomM(2,iflip,k)*emomM(3,iflip,k)**2*emomM(1,iflip,k)**2

       befftemp(3) = befftemp(3) &
            + 2*aw1*emomM(3,iflip,k)*(emomM(1,iflip,k)**2+emomM(2,iflip,k)**2) &
            + 2*aw2+emomM(3,iflip,k)*emomM(1,iflip,k)**2*emomM(2,iflip,k)**2

    endif

  ! Dipolar Interaction Jonathan 19-07-2012
    if(present(Qdip)) then
      if(do_dip==1) then
          do j=1,Natom
             !          do k=1,Mensemble
             befftemp(1) = befftemp(1) + Qdip(1,1,j,iflip)*emomM(1,j,k) + Qdip(2,1,j,iflip)*emomM(2,j,k) + Qdip(3,1,j,iflip)*emomM(3,j,k)
             befftemp(2) = befftemp(2) + Qdip(1,2,j,iflip)*emomM(1,j,k) + Qdip(2,2,j,iflip)*emomM(2,j,k) + Qdip(3,2,j,iflip)*emomM(3,j,k)
             befftemp(3) = befftemp(3) + Qdip(1,3,j,iflip)*emomM(1,j,k) + Qdip(2,3,j,iflip)*emomM(2,j,k) + Qdip(3,3,j,iflip)*emomM(3,j,k)
             !          end do
          end do
       end if
    end if

!-----------------------------------------------------------------------------------------------------------------------

    ! Add static field
    totfield(1:3) = befftemp(1:3)+extfield(1:3)

  end subroutine calculate_efield_tensor



  !> Print out spin configuration of MC run
  subroutine prn_mcmoments(Natom,Mensemble,simid,emom)


    !.. Implicit declarations
    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles
    character(len=8),intent(in) :: simid !< Name of simulation
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
    !
    integer :: i, j
    character(len=30) :: filn

    !.. Executable statements
    write (filn,'(''moment.'',a8,''.out'')') simid
    open(ofileno, file=filn, status='unknown')
    do i=1, Mensemble
       do j=1, Natom
          write (ofileno,10003) i, j, emom(1,j,i), emom(2,j,i), emom(3,j,i), &
               emom(1,j,i)**2+emom(2,j,i)**2+emom(3,j,i)**2
       end do
    end do
    close(ofileno)

10003 format (i8,i8,2x,es16.8,es16.8,es16.8,es16.8)

  end subroutine prn_mcmoments

!> Allocates MC arrays
  subroutine allocate_mcdata(Natom,flag)

    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

    integer :: i_all, i_stat

    if(flag>0) then
      allocate(iflip_a(natom),stat=i_stat)
      call memocc(i_stat,product(shape(iflip_a))*kind(iflip_a),'iflip_a','allocate_mcdata')
    else
      i_all=-product(shape(iflip_a))*kind(iflip_a)
      deallocate(iflip_a,stat=i_stat)
      call memocc(i_stat,i_all,'iflip_a','allocate_mcdata')
    endif
  end subroutine allocate_mcdata

end module MonteCarlo
