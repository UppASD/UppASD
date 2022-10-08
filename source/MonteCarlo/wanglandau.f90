!-------------------------------------------------------------------------------
! MODULE: WangLandau
!> @brief
!> Data and routines for Wang-Landau simulations of Heisenberg models
!> @author
!> A. Bergman
!> @copyright
!> GNU Public License.
!> -----------------------------------------------------------------------------
module WangLandau
   use Parameters
   use Profiling
   use HamiltonianData, only : ham
   implicit none

   !
   real(dblprec), dimension(:), allocatable :: wl_dos !< DOS for Wang-Landau
   real(dblprec), dimension(:), allocatable :: wl_ene !< energy for Wang-Landau
   real(dblprec), dimension(:), allocatable :: wl_hist !< Histogram for Wang-Landau
   real(dblprec), dimension(:), allocatable :: wl_mavg !< Average magnetization
   real(dblprec), dimension(:,:), allocatable :: wl_mhist !< Histogram for magnetization
   real(dblprec), dimension(:,:,:), allocatable :: wl_ran !< Random work array
   !$omp threadprivate(wl_ran)
   real(dblprec) :: wl_totenergy
   real(dblprec) :: wl_emin
   real(dblprec) :: wl_emax
   real(dblprec) :: wl_estep
   real(dblprec) :: wl_stepsize
   real(dblprec) :: wl_gfac
   real(dblprec) :: wl_lcut
   real(dblprec) :: wl_hcut
   real(dblprec) :: wl_maxT
   real(dblprec) :: wl_minT
   integer :: wl_nhist
   integer :: wl_nloop
   integer :: wl_ntemp
   real(dblprec) :: wl_sigma
   real(dblprec) :: wl_hit
   real(dblprec) :: wl_miss
   real(dblprec) :: wl_fac
   real(dblprec) :: wl_direction
   !  deltaT=(wl_maxT-wl_minT)/wl_ntemp
   !character*1 :: do_wanglandau = 'Y'
   !
   public

   !public :: wl_dos,wl_hist, wl_emin,wl_nhist,wl_nloop, wl_estep, wl_fac, wl_direction, wl_totenergy, wl_ene
   !public :: wl_dos,wl_hist, wl_emin,wl_nhist,wl_nloop, wl_estep, wl_fac, wl_direction, wl_totenergy, wl_ene
   !public ::  wl_evolve, allocate_wldata


contains

   !--------------------------------------------------------------------------
   ! SUBROUTINE: wl_evolve
   ! DESCRIPTION
   !> @brief
   !> Main driver for Wang-Landau simulations
   !> For Heisenberg spin systems based on the Hamiltonian
   !>  \f$H=-\frac{1}{2} \sum_{ij} J_{ij} \mathbf{m}_i \cdot \mathbf{m}_j \f$.
   !> Ising model and spin-ice models to be implemented
   !> @author
   !> Anders Bergman (cloned from mc_evolve by Lars Bergqvist)
   !---------------------------------------------------------------------------------

   subroutine wl_evolve(Natom,Nchmax,Mensemble,nHam,&
         mode,conf_num,lsf_metric,lsf_window,do_lsf,lsf_field,exc_inter,&
         lsf_interpolate,do_jtensor,do_dm, do_pd, do_biqdm,do_bq,do_ring,do_chir,do_sa,&
         mult_axis,iflip_a,emomM,emom,mmom,ind_mom_flag,&
         extfield,do_dip,Num_macro,max_num_atom_macro_cell,&
         cell_index,macro_nlistsize,macro_atom_nlist,emomM_macro,emom_macro,mmom_macro,&
         hist,dos,wl_totenergy,m_hist,m_avg,wl_lhist_min,wl_lhist_max,do_anisotropy,stepsize)
         !hist,dos,wl_totenergy,ran_w,m_hist,m_avg)
      !
      use RandomNumbers, only: rng_uniform,rng_uniformP,rng_gaussian
      use LSF, only : mc_update_LSF
      use montecarlo_common
      use optimizationRoutines
      use omp_lib

      !
      implicit none

      !.. Input variables
      ! System variables
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: nHam      !< Number of atoms in Hamiltonian
      character(len=1) :: mode !< Simulation mode (M=MC, H=MC Heat Bath,D=MC Glauber)
      ! LSF variables
      integer, intent(in) :: conf_num  !< number of configurations for LSF
      integer, intent(in) :: lsf_metric !< LSF metric in phase space integration (1=Murata-Doniach,2=Jacobian)
      real(dblprec),intent(in) :: lsf_window             !< Range of moment variation in LSF
      character(len=1), intent(in)  ::  do_lsf           !< Including LSF energy
      character(len=1), intent(in)  ::  lsf_field        !< LSF field contribution (Local/Total)
      character(len=1), intent(in)  ::  exc_inter        !< Exchange interpolation between FM/DLM (Y/N)
      character(len=1), intent(in)  ::  lsf_interpolate  !< Interpolate LSF or not
      ! Heisenberg exchange variables
      integer, intent(in) ::  do_jtensor  !<  Use SKKR style exchange tensor (0=off, 1=on, 2=with biquadratic exchange)
      ! DMI variables
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      ! PD variables
      integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
      !integer, intent(in) :: nn_pd_tot !< Calculated number of neighbours with PD interactions
      ! BIQDM variables
      integer, intent(in) :: do_biqdm   !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      ! BQ variables
      integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      ! Four-spin ring variables
      integer, intent(in) :: do_ring   !< Add four-spin ring (4SR) term to Hamiltonian (0/1)
      ! CHIR variables
      integer, intent(in) :: do_chir !< Add scalar chirality (CHIR) term to Hamiltonian (0/1)
      ! SA  variables
      integer, intent(in) :: do_sa   !< Add Symmetric anisotropic (SA) term to Hamiltonian (0/1)
      ! Anisotropy variables
      character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
      integer, intent(in) :: do_anisotropy
      ! Moments variables
      integer, dimension(Natom),intent(in) :: iflip_a !< Flipping pattern
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      ! Induced moments variables
      character(len=1), intent(in) :: ind_mom_flag
      ! External fields variables
      real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
      ! Dipolar interactions variables
      integer, intent(in) :: do_dip  !< Calculate dipole-dipole contribution (0=Off, 1=Brute Force, 2=macrocell)
      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, dimension(Natom), intent(in) :: cell_index !< Macrocell index for each atom
      integer, intent(in) :: max_num_atom_macro_cell !< Maximum number of atoms in  a macrocell
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
      integer, dimension(Num_macro,max_num_atom_macro_cell), intent(in) :: macro_atom_nlist !< List containing the information of which atoms are in a given macrocell
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emomM_macro !< The full vector of the macrocell magnetic moment
      real(dblprec), dimension(Num_macro,Mensemble), intent(inout) :: mmom_macro !< Magnitude of the macrocell magnetic moments
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emom_macro !< Unit vector of the macrocell magnetic moment
      real(dblprec), dimension(wl_nhist), intent(inout) :: hist
      real(dblprec), dimension(wl_nhist), intent(inout) :: dos
      real(dblprec), intent(inout) :: wl_totenergy
!     real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: ran_w
      real(dblprec), dimension(3,wl_nhist), intent(inout) :: m_hist
      real(dblprec), dimension(3), intent(inout) :: m_avg
      integer ,intent(in) :: wl_lhist_min
      integer ,intent(in) :: wl_lhist_max
      real(dblprec),intent(in) :: stepsize
      !.. Local variables
      integer :: i, k, icell

      real(dblprec) :: de !< Energy difference
      real(dblprec) :: macro_mag_trial
      !.. Local arrays
      real(dblprec), dimension(3) :: newmom !< New trial moment
      real(dblprec), dimension(3) :: macro_trial,halfarray
      real(dblprec),dimension(natom,mensemble) :: flipprob_a ,newmmom_a
      real(dblprec),dimension(3,natom,mensemble) :: newmom_a
!     real(dblprec),dimension(3,natom,mensemble) :: grot
!      real(dblprec),dimension(3) :: grot

!     integer, external :: omp_get_thread_num
      real(dblprec) :: temp_ene

      halfarray=0.5_dblprec


      !!!       if (do_lsf=='Y') then
      !!!          call mc_update_LSF(Natom,Nchmax,Mensemble,nHam, conf_num,do_lsf,emomM, emom, mmom, temperature, temprescale,  &
      !!!             extfield,mult_axis,mode,lsf_interpolate,lsf_field,lsf_window,lsf_metric,exc_inter,iflip_a,&
      !!!             ind_mom_flag,do_dip,Num_macro,mmom_macro,emom_macro,emomM_macro)
      !!!
      !!!       ! Evolution of induced moments as described by Polesya et al.
      !!!       else if (ind_mom_flag=='Y') then
      !!!          call mc_update_ind_mom(Natom,Mensemble,iflip_a,&
      !!!             temperature,temprescale,mode,ham%max_no_neigh,ham%nlistsize,ham%nlist,ham%ncoup,ham%ncoupD,conf_num,exc_inter,&
      !!!             do_dm,ham%max_no_dmneigh,ham%dmlistsize,ham%dmlist,ham%dm_vect,do_pd,ham%nn_pd_tot,ham%pdlistsize,&
      !!!             ham%pdlist,ham%pd_vect,do_biqdm,ham%nn_biqdm_tot,ham%biqdmlistsize,ham%biqdmlist,ham%biqdm_vect,&
      !!!             do_bq,ham%nn_bq_tot,ham%bqlistsize,ham%bqlist,ham%j_bq,ham%taniso,ham%taniso_diff,ham%eaniso,ham%eaniso_diff,ham%kaniso,&
      !!!             ham%kaniso_diff,ham%sb,ham%sb_diff,mult_axis,mmom,emomM,emom,extfield,do_dip,ham%Qdip,&
      !!!             Num_macro,max_num_atom_macro_cell,&
      !!!             cell_index,macro_nlistsize,macro_atom_nlist,emomM_macro,ham%Qdip_macro,emom_macro,mmom_macro,&
      !!!             ham%ind_nlistsize,ham%ind_nlist,ham%ind_list_full,ham%sus_ind,do_lsf,lsf_metric,ind_mom_flag,ham%max_no_neigh_ind)
      !!!
      ! Here mode==E
      ! Random numbers for new spins
!      call rng_uniformP(flipprob_a(:,:),natom*mensemble)
!       print *,'allocated wl_ran?',allocated(wl_ran), omp_get_thread_num()
!       print *,'shape wl_ran',shape(wl_ran), omp_get_thread_num()
       call rng_gaussian(wl_ran(:,:,:),3*natom*mensemble,1.0_dblprec)

         do i=1,Natom
            do k=1,mensemble
!              call rng_gaussian(grot(:),3,1.0_dblprec)
               !call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,1e3_dblprec,0.5_dblprec)
               !call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,temperature,0.5_dblprec)
!              print '(3f12.6)',wl_ran(:,i,k)
                call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,stepsize,halfarray,wl_ran(:,i,k))
!              call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,7e3_dblprec,0.5_dblprec,grot(:))
               !print '(6f10.5)',emom(1:3,i,k),newmom(1:3)
               newmom_a(1:3,i,k)=newmom(1:3)
            enddo
         enddo

      ! Random numbers for probability
      ! print *,'aft choose_random_flip',omp_get_thread_num()
      call rng_uniform(flipprob_a(:,:),natom*mensemble)
      ! print *,'aft choose_random_flip',omp_get_thread_num()

      ! Calculate energy and flip spin if preferred
      temp_ene=0.0_dblprec
      !!$omp parallel do default(shared), private(i,k,de,totfield), schedule(auto),collapse(2) &
      !!$omp reduction(+:temp_ene)
      do i=1, Natom
         do k=1,mensemble

            ! Calculate the energy
            call calculate_energy(Natom, Mensemble, nHam, conf_num, do_dm , do_pd, do_biqdm, do_bq, do_ring, do_chir,do_sa,&
               emomM, emom, mmom, iflip_a(i), newmom_a(1:3,iflip_a(i),k), extfield, de, k, &
               mult_axis, do_dip,Num_macro,max_num_atom_macro_cell,cell_index,macro_nlistsize,&
               macro_atom_nlist,emomM_macro,icell,macro_mag_trial,macro_trial,exc_inter,do_anisotropy, do_jtensor)
      !     call effective_field_extralite(Natom,Mensemble,iflip_a(i),iflip_a(i),emomM,mmom,temp_ene,beff)
      !      call effective_field_lite(Natom,Mensemble,i,i,do_jtensor,      &
      !         do_anisotropy,exc_inter,do_dm,do_pd,do_biqdm,do_bq,do_chir,do_dip,emomM,mmom, &
      !         external_field,time_external_field,beff,beff1,beff2,mult_axis,totenergy,NA,N1,N2,N3)

            call flip_wl(Natom, Mensemble, emom, emomM, mmom, iflip_a(i),newmom_a(1:3,iflip_a(i),k),newmmom_a(iflip_a(i),k), &
               de,do_lsf,k,flipprob_a(i,k),lsf_metric,ham%ind_nlistsize,&
               ham%ind_nlist,ind_mom_flag,ham%max_no_neigh_ind,ham%sus_ind,ham%ind_list_full,&
               do_dip,Num_macro,icell,macro_mag_trial,macro_trial,mmom_macro,emom_macro,emomM_macro,&
               wl_nhist,wl_lhist_min,wl_lhist_max,dos,hist,m_hist,m_avg,wl_estep,wl_emin,wl_fac,wl_totenergy)

            !temp_ene=temp_ene+de
            !!$omp atomic
            wl_totenergy=wl_totenergy+de!/mry
            !print *,'E= ',wl_totenergy+temp_ene
            !!!                     endif
            !!!                  endif   !jtensor
         enddo    !ensemble
      enddo     !atom
      !!$omp end parallel do
      !wl_totenergy=wl_totenergy+temp_ene!/mry
      !print *,wl_totenergy
      !!         endif       !mode
      !!      endif         !lsf
      return

   end subroutine wl_evolve


   !--------------------------------------------------------------------------
   ! SUBROUTINE: wl_warmup
   ! DESCRIPTION
   !> @brief
   !> Main driver for Wang-Landau initialization
   !> For Heisenberg spin systems based on the Hamiltonian
   !>  \f$H=-\frac{1}{2} \sum_{ij} J_{ij} \mathbf{m}_i \cdot \mathbf{m}_j \f$.
   !> Ising model and spin-ice models to be implemented
   !> @author
   !> Anders Bergman (cloned from mc_evolve by Lars Bergqvist)
   !---------------------------------------------------------------------------------
   !call wl_warmup(Natom,Nchmax,Mensemble,nHam,Temp,temprescale,&
   !   mode,conf_num,lsf_metric,lsf_window,do_lsf,lsf_field,exc_inter,&
   !   lsf_interpolate,do_jtensor,do_dm, do_pd, do_biqdm,do_bq,&
   !   mult_axis,iflip_a,emomM,emom,mmom,ind_mom_flag,hfield,&
   !   do_dip,Num_macro,max_num_atom_macro_cell,&
   !   cell_index,macro_nlistsize,macro_atom_nlist,emomM_macro,emom_macro,mmom_macro)

   subroutine wl_warmup(Natom,Nchmax,Mensemble,nHam,&
         mode,conf_num,lsf_metric,lsf_window,do_lsf,lsf_field,exc_inter,&
         lsf_interpolate,do_jtensor,do_dm, do_pd, do_biqdm,do_bq,do_ring,do_chir,do_sa,&
         mult_axis,iflip_a,emomM,emom,mmom,ind_mom_flag,&
         extfield,do_dip,Num_macro,max_num_atom_macro_cell,&
         cell_index,macro_nlistsize,macro_atom_nlist,emomM_macro,emom_macro,mmom_macro,&
         do_anisotropy)
      !
      use RandomNumbers, only: rng_uniform,rng_uniformP,rng_gaussianP,use_vsl
      use LSF, only : mc_update_LSF
      use montecarlo_common
      !
      implicit none

      !.. Input variables
      ! System variables
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: nHam      !< Number of atoms in Hamiltonian
      character(len=1) :: mode !< Simulation mode (M=MC, H=MC Heat Bath,D=MC Glauber)
      ! LSF variables
      integer, intent(in) :: conf_num  !< number of configurations for LSF
      integer, intent(in) :: lsf_metric !< LSF metric in phase space integration (1=Murata-Doniach,2=Jacobian)
      real(dblprec),intent(in) :: lsf_window             !< Range of moment variation in LSF
      character(len=1), intent(in)  ::  do_lsf           !< Including LSF energy
      character(len=1), intent(in)  ::  lsf_field        !< LSF field contribution (Local/Total)
      character(len=1), intent(in)  ::  exc_inter        !< Exchange interpolation between FM/DLM (Y/N)
      character(len=1), intent(in)  ::  lsf_interpolate  !< Interpolate LSF or not
      ! Heisenberg exchange variables
      integer, intent(in) ::  do_jtensor  !<  Use SKKR style exchange tensor (0=off, 1=on, 2=with biquadratic exchange)
      ! DMI variables
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      ! PD variables
      integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
      !integer, intent(in) :: nn_pd_tot !< Calculated number of neighbours with PD interactions
      ! BIQDM variables
      integer, intent(in) :: do_biqdm   !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      ! BQ variables
      integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      ! Four-spin ring (4SR) variables
      integer, intent(in) :: do_ring   !< Add four-spin ring (4SR) term to Hamiltonian (0/1)
      ! CHIR variables
      integer, intent(in) :: do_chir !< Add scalar chirality (CHIR) term to Hamiltonian (0/1)
      ! SA  variables
      integer, intent(in) :: do_sa   !< Add Symmetric anisotropic (SA) term to Hamiltonian (0/1)
      ! Anisotropy variables
      character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
      integer, intent(in) :: do_anisotropy
      ! Moments variables
      integer, dimension(Natom),intent(in) :: iflip_a !< Flipping pattern
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      ! Induced moments variables
      character(len=1), intent(in) :: ind_mom_flag
      ! External fields variables
      real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
      ! Dipolar interactions variables
      integer, intent(in) :: do_dip  !< Calculate dipole-dipole contribution (0=Off, 1=Brute Force, 2=macrocell)
      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, dimension(Natom), intent(in) :: cell_index !< Macrocell index for each atom
      integer, intent(in) :: max_num_atom_macro_cell !< Maximum number of atoms in  a macrocell
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
      integer, dimension(Num_macro,max_num_atom_macro_cell), intent(in) :: macro_atom_nlist !< List containing the information of which atoms are in a given macrocell
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emomM_macro !< The full vector of the macrocell magnetic moment
      real(dblprec), dimension(Num_macro,Mensemble), intent(inout) :: mmom_macro !< Magnitude of the macrocell magnetic moments
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emom_macro !< Unit vector of the macrocell magnetic moment
      !.. Local variables
      integer :: i, k,icell

      real(dblprec) :: de !< Energy difference
      real(dblprec) :: macro_mag_trial
      !.. Local arrays
      real(dblprec), dimension(3) :: newmom !< New trial moment
      real(dblprec), dimension(3) :: macro_trial
      real(dblprec), dimension(3) :: totfield  !<Total effective field acting on each moment
      real(dblprec),dimension(natom,mensemble) :: flipprob_a ,newmmom_a
      real(dblprec),dimension(3,natom,mensemble) :: newmom_a
!      real(dblprec),dimension(3,natom,mensemble) :: grot
      real(dblprec),dimension(3) :: grot, guni

      integer, external :: omp_get_thread_num
      real(dblprec) :: temp_ene

      grot=0.0_dblprec ; temp_ene=0.0_dblprec ; guni=0.0_dblprec

      ! Here mode==E
      ! Random numbers for new spins
      !call rng_uniformP(flipprob_a(:,:),natom*mensemble)
!      call rng_gaussianP(grot,3*natom*mensemble,1.0_dblprec)

      if(use_vsl) then
#ifdef VSL
        !$omp parallel do default(shared),private(i,k,newmom,guni,grot),schedule(auto),collapse(2)
#endif
         do i=1,Natom
            do k=1,mensemble
               call rng_uniform(guni(2:3),2)
               !call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,temperature,flipprob_a(i,k))
               call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,0.0_dblprec,guni(:),grot(:))
               newmom_a(1:3,i,k)=newmom(1:3)
            enddo
         enddo
#ifdef VSL
        !$omp end parallel do
#endif
      else
         do i=1,Natom
            do k=1,mensemble
               call rng_uniform(guni(2:3),2)
               !call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,temperature,flipprob_a(i,k))
               call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,0.0_dblprec,guni(:),grot(:))
               newmom_a(1:3,i,k)=newmom(1:3)
            enddo
         enddo
      end if

      ! Random numbers for probability
      call rng_uniformP(flipprob_a(:,:),natom*mensemble)

      ! Calculate energy and flip spin if preferred
      !$omp parallel do default(shared), private(i,k,de,totfield), schedule(auto),collapse(2)
      !!$omp reduction(+:temp_ene)
      do i=1, Natom
         do k=1,mensemble
            ! Calculate the energy
            call calculate_energy(Natom, Mensemble, nHam, conf_num, do_dm , do_pd, do_biqdm, do_bq, do_ring, do_chir,do_sa,&
               emomM, emom, mmom, iflip_a(i), newmom_a(1:3,iflip_a(i),k), extfield, de, k, &
               mult_axis, do_dip,Num_macro,max_num_atom_macro_cell,cell_index,macro_nlistsize,&
               macro_atom_nlist,emomM_macro,icell,macro_mag_trial,macro_trial,exc_inter,do_anisotropy, do_jtensor)

            call minimaxi(Natom, Mensemble, emom, emomM, mmom, iflip_a(i),newmom_a(1:3,iflip_a(i),k),newmmom_a(iflip_a(i),k), &
               de,do_lsf,k,flipprob_a(i,k),lsf_metric,ham%ind_nlistsize,&
               ham%ind_nlist,ind_mom_flag,ham%max_no_neigh_ind,ham%sus_ind,ham%ind_list_full,&
               do_dip,Num_macro,icell,macro_mag_trial,macro_trial,mmom_macro,emom_macro,emomM_macro,&
               wl_nhist,wl_dos,wl_hist,wl_estep,wl_emin,wl_fac,wl_totenergy,wl_direction)
               !wl_nhist,wl_dos,wl_hist,wl_estep,wl_emin,wl_fac,wl_totenergy+temp_ene,wl_direction)

!           temp_ene=temp_ene+de
           !!$omp atomic
            wl_totenergy=wl_totenergy+de
            !print *,de

         enddo    !ensemble
      enddo     !atom
     !$omp end parallel do
!    wl_totenergy=wl_totenergy+temp_ene!/mry
!      print *,wl_totenergy
     !write(777,*) wl_totenergy
      return

   end subroutine wl_warmup


   subroutine allocate_wldata(Natom,Mensemble,flag)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles in system
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(wl_hist(wl_nhist),stat=i_stat)
         call memocc(i_stat,product(shape(wl_hist))*kind(wl_hist),'wl_hist','allocate_wldata')
         allocate(wl_dos(wl_nhist),stat=i_stat)
         call memocc(i_stat,product(shape(wl_dos))*kind(wl_dos),'wl_dos','allocate_wldata')
         allocate(wl_ene(wl_nhist),stat=i_stat)
         call memocc(i_stat,product(shape(wl_ene))*kind(wl_ene),'wl_ene','allocate_wldata')
         !$omp parallel
         allocate(wl_ran(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(wl_ran))*kind(wl_ran),'wl_ran','allocate_wldata')
         !$omp end parallel
         allocate(wl_mhist(3,wl_nhist),stat=i_stat)
         call memocc(i_stat,product(shape(wl_mhist))*kind(wl_mhist),'wl_mhist','allocate_wldata')
         allocate(wl_mavg(3),stat=i_stat)
         call memocc(i_stat,product(shape(wl_mavg))*kind(wl_mavg),'wl_mavg','allocate_wldata')
      else
         i_all=-product(shape(wl_hist))*kind(wl_hist)
         deallocate(wl_hist,stat=i_stat)
         call memocc(i_stat,i_all,'wl_hist','allocate_wldata')
         i_all=-product(shape(wl_dos))*kind(wl_dos)
         deallocate(wl_dos,stat=i_stat)
         call memocc(i_stat,i_all,'wl_dos','allocate_wldata')
         i_all=-product(shape(wl_ene))*kind(wl_ene)
         deallocate(wl_ene,stat=i_stat)
         call memocc(i_stat,i_all,'wl_ene','allocate_wldata')
         !$omp parallel
         i_all=-product(shape(wl_ran))*kind(wl_ran)
         deallocate(wl_ran,stat=i_stat)
         call memocc(i_stat,i_all,'wl_ran','allocate_wldata')
         i_all=-product(shape(wl_mhist))*kind(wl_mhist)
         !$omp end parallel
         deallocate(wl_mhist,stat=i_stat)
         call memocc(i_stat,i_all,'wl_mhist','allocate_wldata')
         i_all=-product(shape(wl_mavg))*kind(wl_mavg)
         deallocate(wl_mavg,stat=i_stat)
         call memocc(i_stat,i_all,'wl_mavg','allocate_wldata')
      endif
   end subroutine allocate_wldata

   !---------------------------------------------------------------------------
   !> @brief
   !> Flip selected spins according to Wang-Landau DOS
   !>
   !> @author
   !> Anders Bergman (cloned from flip_a)
   !>
   !---------------------------------------------------------------------------
   subroutine flip_wl(Natom, Mensemble, emom, emomM, mmom, iflip,newmom,newmmom,de,&
         do_lsf,k,flipprob,lsf_metric,ind_nlistsize,&
         ind_nlist,ind_mom_flag,max_no_neigh_ind,sus_ind,ind_list_full,&
         do_dip,Num_macro,icell,macro_mag_trial,macro_trial,mmom_macro,emom_macro,emomM_macro,&
         nhist,lhist_min,lhist_max,dos,hist,m_hist,m_avg,estep,emin,fac,wl_totenergy)

      use Constants
      !use MonteCarlo

      !.. Implicit declarations
      implicit none

      !.. Input variables
      ! System/geometry variables
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      ! Moment variables
      integer, intent(in) :: iflip !< Atom to flip spin for
      real(dblprec), intent(in) :: newmmom !< New trial magnitude of moment
      real(dblprec), dimension(3), intent(in) :: newmom !< New trial moment
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      ! Probabilities variables
      real(dblprec), intent(inout):: de  !< Energy difference
      real(dblprec), intent(in) :: flipprob !< Probability of flipping spin
      ! LSF variables
      integer, intent(in) :: lsf_metric !< LSF metric or phase space measure
      character(len=1), intent(in)  ::  do_lsf     !< Including LSF energy
      ! Induced moment variables
      integer, intent(in) :: max_no_neigh_ind !< Calculated maximum of neighbours for induced moments
      integer, dimension(Natom), intent(in) :: ind_list_full
      integer, dimension(Natom), intent(in) :: ind_nlistsize !< Size of neighbour list for induced moments
      integer, dimension(max_no_neigh_ind,Natom), intent(in) :: ind_nlist !< Neighbour list for induced moments
      real(dblprec), dimension(Natom), intent(in) :: sus_ind
      character(len=1), intent(in) :: ind_mom_flag
      integer, intent(in) :: do_dip
      integer, intent(in) :: Num_macro
      integer, intent(in) :: icell
      real(dblprec), intent(in) :: macro_mag_trial
      real(dblprec), dimension(3), intent(in) :: macro_trial
      real(dblprec), dimension(Num_macro,Mensemble), intent(inout) :: mmom_macro
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emom_macro
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emomM_macro

      integer, intent(in) :: nhist !< Number of histogram bins
      integer, intent(in) :: lhist_min !< Lower local histogram range
      integer, intent(in) :: lhist_max !< Upper local histogram range
      real(dblprec), dimension(nhist), intent(inout) :: dos !< WL dos
      real(dblprec), dimension(nhist), intent(inout) :: hist !< WL E-histogram
      real(dblprec), dimension(3,nhist), intent(inout) :: m_hist !< WL M-histogram
      real(dblprec), dimension(3), intent(inout) :: m_avg  !< avergage magnetization
      real(dblprec), intent(in) :: estep  !< WL bin width
      real(dblprec), intent(in) :: emin   !< WL minimum energy
      real(dblprec), intent(in) :: fac   !< WL histogram factor
      real(dblprec), intent(in) :: wl_totenergy !< Total energy

      integer :: k !< Current ensemble
      integer :: idx_0, idx_1, idx_s , delta
      real(dblprec) :: edist, gam_fac
      real(dblprec) :: g1,g0, prob, conv_de

      real(dblprec) :: e_trial

      !e_min=-1589.000_dblprec
      !e_max=664.000_dblprec
      gam_fac=1.0_dblprec
      delta=floor(nhist*wl_sigma)

      conv_de=de/mry
      !if(wl_totenergy+conv_de>e_max.or.wl_totenergy+conv_de<e_min)  direction=-1.0_dblprec*direction
      !print '(3f10.4,2i6)',wl_totenergy,de ,conv_de,int((wl_totenergy-emin)/estep),int((wl_totenergy+conv_de-emin)/estep)
      !if(wl_totenergy>emin.and.wl_totenergy+conv_de<emax.and.wl_totenergy+conv_de>emin) then
      idx_0=int((wl_totenergy-emin)/(estep))+1
      idx_1=int((wl_totenergy+conv_de-emin)/(estep))+1
      e_trial=wl_totenergy+conv_de
      !print *,'evolve: 1  ',wl_totenergy-emin,wl_totenergy+de-emin
      !print *,'evolve: 2  ',wl_totenergy,wl_totenergy+de
      !print *,'evolve: 3  ',dos(idx_0)-dos(idx_1),log(flipprob)
      !print *,'evolve: 4  ',idx_0,idx_1,nhist
      !print *,wl_totenergy,direction,wl_totenergy+conv_de<e_min,wl_totenergy+conv_de>e_max
      !if(idx_0>=1.and.idx_0<=nhist.and.idx_1>=1.and.idx_1<=nhist) then
      if(idx_0>=lhist_min.and.idx_0<=lhist_max.and.idx_1>=lhist_min.and.idx_1<=lhist_max) then
         !!print *,idx_0,idx_1
         g0=dos(idx_0)
         g1=dos(idx_1)
         prob=min(0.0_dblprec,wl_gfac*(g0-g1))
         !beta=1.0_dblprec/k_bolt/(temprescale*Temperature+1.0d-15)
         !!!if (do_lsf=='N'.and.ind_mom_flag=='N') then
         !if(direction*de>0.0_dblprec) then
         !if(e_trial>-120.0_dblprec.and.e_trial<120.0_dblprec.and.log(flipprob)<prob) then
         !print *,idx_0,idx_1,log(flipprob)<=prob
         if(log(flipprob)<=prob) then
            ! Updating magnetic average
            m_avg=m_avg-emomM(:,iflip,k)+mmom(iflip,k)*newmom(:)
!           if(idx_1==1) print '(a,5f20.6)','MOM',m_avg,sqrt(sum(m_avg*m_avg)),e_trial
            !m_avg=m_avg-emomM(:,iflip,k)+mmom(iflip,k)*newmom(:)
            ! Making the flip
            emom(1:3,iflip,k)=newmom(1:3)
            emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
            if (do_dip==2) then
               emomM_macro(:,icell,k)=macro_trial(:)
               emom_macro(:,icell,k)=macro_trial(:)/macro_mag_trial
               mmom_macro(icell,k)=macro_mag_trial
            endif
            wl_hit=wl_hit+1.0_dblprec
            ! Possible smearing below
!!$omp simd private(idx_s,edist)
            do idx_s=max(idx_1-delta,1),min(idx_1+delta,nhist)
!                edist=exp(-(( (idx_1-idx_s)**2*estep)/(4.0_dblprec*delta)))   !Gaussian
!                edist=1.0_dblprec-(((idx_1-idx_s)*1.0_dblprec)/delta)**2   !Epanchnikov - gives oscillations at E=0
                edist=1.0_dblprec-abs((1.0_dblprec*(idx_1-idx_s)/delta)) !Triangular

!                write(*,'(i5,i5,1f10.7)') idx_1,idx_s,edist
                dos(idx_s)=dos(idx_s)+gam_fac*edist*fac
                hist(idx_s)=hist(idx_s)+edist
            end do
            m_hist(:,idx_1)=m_hist(:,idx_1)+m_avg
            de=conv_de
            !!$omp atomic
!            dos(idx_1)=dos(idx_1)+fac
            !!$omp atomic
!            hist(idx_1)=hist(idx_1)+1.0_dblprec !prefac
            !!$omp atomic
!           m_hist(:,idx_1)=m_hist(:,idx_1)+m_avg
         else
            de=0.0_dblprec
            hist(idx_0)=hist(idx_0)+1.0_dblprec
            m_hist(:,idx_0)=m_hist(:,idx_0)+m_avg
            dos(idx_0)=dos(idx_0)+fac
            wl_miss=wl_miss+1.0_dblprec
         endif
         !!!else if (do_lsf=='N'.and.ind_mom_flag=='Y') then
         !!!   if(de<=0.0_dblprec .or. flipprob<exp(-beta*de)) then
         !!!      emom(:,iflip,k)=newmom(:)
         !!!      emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
         !!!      ! If the flip is accepted then the magnitude of the induced moments must change
         !!!      do ind_neigh=1,ind_nlistsize(iflip)
         !!!         ! Select the current induced moment
         !!!         curr_ind=ind_nlist(ind_neigh,iflip)
         !!!         if (ind_list_full(curr_ind).eq.1) then

         !!!            ! Sum over the nearest neighbours that are fixed
         !!!            ave_mom=0.0_dblprec
         !!!            do fix_neigh=1, ind_nlistsize(curr_ind)
         !!!               ave_mom(1)=ave_mom(1)+emomM(1,ind_nlist(fix_neigh,curr_ind),k)
         !!!               ave_mom(2)=ave_mom(2)+emomM(2,ind_nlist(fix_neigh,curr_ind),k)
         !!!               ave_mom(3)=ave_mom(3)+emomM(3,ind_nlist(fix_neigh,curr_ind),k)
         !!!            enddo
         !!!            ! Vary the magnitude of the induced moments
         !!!            emomM(:,curr_ind,k)=ave_mom(:)*sus_ind(curr_ind)
         !!!            mmom(curr_ind,k)=sqrt(emomM(1,curr_ind,k)**2+emomM(2,curr_ind,k)**2+emomM(3,curr_ind,k)**2)
         !!!            emom(:,curr_ind,k)=emomM(:,curr_ind,k)/mmom(curr_ind,k)
         !!!         endif
         !!!      enddo
         !!!   endif
         !!!else
         !!!   if(lsf_metric==1) then
         !!!      lmetric=1.0_dblprec
         !!!   elseif(lsf_metric==2) then
         !!!      lmetric=(mmom(iflip,k)/newmmom)**2
         !!!   else
         !!!      lmetric=(mmom(iflip,k)/newmmom)
         !!!   endif
         !!!   des=-log(lmetric)/beta
         !!!   if(de<=des .or. flipprob<exp(-beta*de)/lmetric) then
         !!!      mmom(iflip,k) = newmmom
         !!!      emom(:,iflip,k)=newmom(:)
         !!!      emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
         !!!   endif
         !!!endif
      else
         de=0.0_dblprec
         dos(idx_0)=dos(idx_0)+fac
         hist(idx_0)=hist(idx_0)+1.0_dblprec !prefac
         wl_miss=wl_miss+1.0_dblprec
         !print *,idx_0,idx_1,nhist
      end if
   end subroutine flip_wl

   !---------------------------------------------------------------------------
   !> @brief
   !> Flip selected spins according to Wang-Landau DOS
   !>
   !> @author
   !> Anders Bergman (cloned from flip_a)
   !>
   !---------------------------------------------------------------------------
   !call minimaxi(Natom, Mensemble, emom, emomM, mmom, iflip_a(i),newmom_a(1:3,iflip_a(i),k),newmmom_a(iflip_a(i),k), &
   !   de,do_lsf,k,flipprob_a(i,k),lsf_metric,ham%ind_nlistsize,&
   !   ham%ind_nlist,ind_mom_flag,ham%max_no_neigh_ind,ham%sus_ind,ham%ind_list_full,&
   !   do_dip,Num_macro,icell,macro_mag_trial,macro_trial,mmom_macro,emom_macro,emomM_macro,&
   !   wl_nhist,wl_dos,wl_hist,wl_estep,wl_emin,wl_fac,wl_totenergy+temp_ene,wl_direction)

   subroutine minimaxi(Natom, Mensemble, emom, emomM, mmom, iflip,newmom,newmmom,de,&
         do_lsf,k,flipprob,lsf_metric,ind_nlistsize,&
         ind_nlist,ind_mom_flag,max_no_neigh_ind,sus_ind,ind_list_full,&
         do_dip,Num_macro,icell,macro_mag_trial,macro_trial,mmom_macro,emom_macro,emomM_macro,&
         wl_nhist,wl_dos,wl_hist,wl_estep,wl_emin,wl_fac,wl_totenergy,direction)

      use Constants
      !use MonteCarlo

      !.. Implicit declarations
      implicit none

      !.. Input variables
      ! System/geometry variables
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      ! Moment variables
      integer, intent(in) :: iflip !< Atom to flip spin for
      real(dblprec), intent(in) :: newmmom !< New trial magnitude of moment
      real(dblprec), dimension(3), intent(in) :: newmom !< New trial moment
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      ! Probabilities variables
      real(dblprec), intent(inout):: de  !< Energy difference
      real(dblprec), intent(in) :: flipprob !< Probability of flipping spin
      ! LSF variables
      integer, intent(in) :: lsf_metric !< LSF metric or phase space measure
      character(len=1), intent(in)  ::  do_lsf     !< Including LSF energy
      ! Induced moment variables
      integer, intent(in) :: max_no_neigh_ind !< Calculated maximum of neighbours for induced moments
      integer, dimension(Natom), intent(in) :: ind_list_full
      integer, dimension(Natom), intent(in) :: ind_nlistsize !< Size of neighbour list for induced moments
      integer, dimension(max_no_neigh_ind,Natom), intent(in) :: ind_nlist !< Neighbour list for induced moments
      real(dblprec), dimension(Natom), intent(in) :: sus_ind
      character(len=1), intent(in) :: ind_mom_flag
      integer, intent(in) :: do_dip
      integer, intent(in) :: Num_macro
      integer, intent(in) :: icell
      real(dblprec), intent(in) :: macro_mag_trial
      real(dblprec), dimension(3), intent(in) :: macro_trial
      real(dblprec), dimension(Num_macro,Mensemble), intent(inout) :: mmom_macro
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emom_macro
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emomM_macro

      integer, intent(in) :: wl_nhist !< Number of histogram bins
      real(dblprec), dimension(wl_nhist), intent(inout) :: wl_dos !< WL dos
      real(dblprec), dimension(wl_nhist), intent(inout) :: wl_hist !< WL histogram
      real(dblprec), intent(in) :: wl_estep  !< WL bin width
      real(dblprec), intent(in) :: wl_emin   !< WL minimum energy
      real(dblprec), intent(in) :: wl_fac   !< WL histogram factor
      real(dblprec), intent(in) :: wl_totenergy !< Total energy
      real(dblprec), intent(in) :: direction  !< Direction of energy sweep

      integer :: k !< Current ensemble


      !e_min=-1589.000_dblprec
      !e_max=664.000_dblprec

      !conv_de=de/mry
      !!! !if(wl_totenergy+conv_de>e_max.or.wl_totenergy+conv_de<e_min)  wl_direction=-1.0_dblprec*wl_direction
      !!! !print *,wl_totenergy,de ,conv_de!int((wl_totenergy-wl_emin)/wl_estep),wl_nhist
      !!! idx_0=int((wl_totenergy-wl_emin)/wl_estep)+1
      !!! idx_1=int((wl_totenergy+conv_de-wl_emin)/wl_estep)+1
      !!! e_trial=wl_totenergy+conv_de
      !!! !!print *,wl_totenergy-wl_emin,wl_totenergy+de-wl_emin
      !!! !print *,wl_totenergy,wl_direction,wl_totenergy+conv_de<e_min,wl_totenergy+conv_de>e_max
      !!! !!print *,idx_0,idx_1
      !!! g0=wl_dos(idx_0)
      !!! g1=wl_dos(idx_1)
      !!! prob=min(0.0_dblprec,g0-g1)
      !!! !beta=1.0_dblprec/k_bolt/(temprescale*Temperature+1.0d-15)
      !!! !!!if (do_lsf=='N'.and.ind_mom_flag=='N') then
      !!! !if(wl_direction*de>0.0_dblprec) then
      !!! !if(e_trial>-120.0_dblprec.and.e_trial<120.0_dblprec.and.log(flipprob)<prob) then
      !print '(3f10.5,a4,3f10.5)',emom(:,iflip,k),newmom(:)
      if(direction*de>0.0_dblprec) then
         emom(:,iflip,k)=newmom(:)
         emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
         if (do_dip==2) then
            emomM_macro(:,icell,k)=macro_trial(:)
            emom_macro(:,icell,k)=macro_trial(:)/macro_mag_trial
            mmom_macro(icell,k)=macro_mag_trial
         endif
         !print *,'wl_hit!',idx_0,idx_1
         ! do idx_s=max(idx_1-1,1),min(idx_1+1,wl_nhist)
         !    edist=wl_totenergy+conv_de-((idx_s-1.0_dblprec)*wl_estep+wl_emin)
         !    edist=edist/(wl_estep)
         !    prefac=3.0_dblprec/4.0_dblprec*(1.0_dblprec-edist**2)
         !    !print '(2x,3f12.5,2i6,f12.5)',wl_totenergy+conv_de,edist,prefac,idx_1,idx_s,wl_estep
         !    !print *,wl_totenergy,edist,prefac
         !    wl_dos(idx_s)=wl_dos(idx_s)+wl_fac*prefac
         !    wl_hist(idx_s)=wl_hist(idx_s)+prefac
         ! end do
         !wl_hist(idx_1)=wl_hist(idx_1)+1
         !!! !wl_dos(idx_1)=wl_dos(idx_1)+wl_fac*prefac
         !!! edist=wl_totenergy+conv_de-((idx_1-1.0_dblprec)*wl_estep+wl_emin)
         !!! edist=edist/(wl_estep)
         !!! prefac=(1.0_dblprec-edist**2)
         !!! !print *,edist,prefac
         !!! wl_dos(idx_1)=wl_dos(idx_1)+wl_fac*prefac
         !!! wl_hist(idx_1)=wl_hist(idx_1)+1.0_dblprec !prefac
         de=de/mry
      else
         de=0.0_dblprec
         !wl_hist(idx_0)=wl_hist(idx_0)+1
         !wl_dos(idx_0)=wl_dos(idx_0)+wl_fac
      endif
      !!!else if (do_lsf=='N'.and.ind_mom_flag=='Y') then
      !!!   if(de<=0.0_dblprec .or. flipprob<exp(-beta*de)) then
      !!!      emom(:,iflip,k)=newmom(:)
      !!!      emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
      !!!      ! If the flip is accepted then the magnitude of the induced moments must change
      !!!      do ind_neigh=1,ind_nlistsize(iflip)
      !!!         ! Select the current induced moment
      !!!         curr_ind=ind_nlist(ind_neigh,iflip)
      !!!         if (ind_list_full(curr_ind).eq.1) then

      !!!            ! Sum over the nearest neighbours that are fixed
      !!!            ave_mom=0.0_dblprec
      !!!            do fix_neigh=1, ind_nlistsize(curr_ind)
      !!!               ave_mom(1)=ave_mom(1)+emomM(1,ind_nlist(fix_neigh,curr_ind),k)
      !!!               ave_mom(2)=ave_mom(2)+emomM(2,ind_nlist(fix_neigh,curr_ind),k)
      !!!               ave_mom(3)=ave_mom(3)+emomM(3,ind_nlist(fix_neigh,curr_ind),k)
      !!!            enddo
      !!!            ! Vary the magnitude of the induced moments
      !!!            emomM(:,curr_ind,k)=ave_mom(:)*sus_ind(curr_ind)
      !!!            mmom(curr_ind,k)=sqrt(emomM(1,curr_ind,k)**2+emomM(2,curr_ind,k)**2+emomM(3,curr_ind,k)**2)
      !!!            emom(:,curr_ind,k)=emomM(:,curr_ind,k)/mmom(curr_ind,k)
      !!!         endif
      !!!      enddo
      !!!   endif
      !!!else
      !!!   if(lsf_metric==1) then
      !!!      lmetric=1.0_dblprec
      !!!   elseif(lsf_metric==2) then
      !!!      lmetric=(mmom(iflip,k)/newmmom)**2
      !!!   else
      !!!      lmetric=(mmom(iflip,k)/newmmom)
      !!!   endif
      !!!   des=-log(lmetric)/beta
      !!!   if(de<=des .or. flipprob<exp(-beta*de)/lmetric) then
      !!!      mmom(iflip,k) = newmmom
      !!!      emom(:,iflip,k)=newmom(:)
      !!!      emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
      !!!   endif
      !!!endif
   end subroutine minimaxi

   subroutine wl_integration(nhist,ene,dos,hist,simid)
      !
      use Constants
      use QHB
      use InputData, only : Natom
      !
      implicit none
      !
      integer, intent(in) :: nhist
      real(dblprec), dimension(nhist), intent(in) :: ene
      real(dblprec), dimension(nhist), intent(in) :: dos
      real(dblprec), dimension(nhist), intent(in) :: hist
      character(len=8), intent(in) :: simid  !< Name of simulation


      integer :: ihist, itemp
      character(len=30) :: filn
      !real(dblprec), dimension(nhist) :: dos, hist, ene, dosS, histS
      !real(dblprec), dimension(-3:3) :: dsmear
      real(dblprec) :: lambda, U, E2, Z, F, iS, C, deltaT, temp, tempC, kb
      real(dblprec) :: lambdaQ, UQ, E2Q, ZQ, FQ, iSQ, CQ
      real(dblprec) :: mavg, temprescale,temprescalegrad
      !
      kb=k_bolt/mry
      !
      !
      write(filn,'(''wloutput.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position="append")
      write(ofileno,'(a)') "#   T (K)    x      U(classical)      U(quantum)     Cv (classical)    Cv(quantum)  F(classical)     F(quantum)     S(classical)       S(quantum)"
      iS=0.0_dblprec ; iSQ=0.0_dblprec
      deltaT=(wl_maxT-wl_minT)/wl_ntemp
      do itemp=0,wl_ntemp
         tempC=1.0d-6+wl_minT+deltat*itemp
         ! Rescaling of temperature according to Quantum Heat bath
         temprescale=1.0_dblprec
         temprescalegrad=0.0_dblprec
         if (do_qhb=='Q' .or. do_qhb=='R' .or. do_qhb=='P' .or. do_qhb=='T') then
            if(qhb_mode=='MT') then
               !!call calc_mavrg(Natom,Mensemble,emomM,mavg)
               !mavg=
               !call qhb_rescale(temp,temprescale,do_qhb,qhb_mode,mavg)
            else
               call qhb_rescale(tempC,temprescale,temprescalegrad,do_qhb,qhb_mode,mavg)
            endif
         endif
         temp=tempC*temprescale

         lambda=maxval(dos-ene/(tempC*kb))
         lambdaQ=maxval(dos-ene/(temp*kb))
         U=0.0_dblprec
         UQ=0.0_dblprec
         Z=0.0_dblprec
         ZQ=0.0_dblprec
         E2=0.0_dblprec
         E2Q=0.0_dblprec
         F=0.0_dblprec
         FQ=0.0_dblprec
         do ihist=1,nhist
            Z=Z+exp(dos(ihist)-ene(ihist)/(tempC*kb)-lambda)
            U=U+ene(ihist)*exp(dos(ihist)-ene(ihist)/(tempC*kb)-lambda)
            E2=E2+ene(ihist)**2*exp(dos(ihist)-ene(ihist)/(tempC*kb)-lambda)
!            F=F+exp(dos(ihist)-ene(ihist)/(tempC*kb)-lambda)
            ZQ=ZQ+exp(dos(ihist)-ene(ihist)/(temp*kb)-lambdaQ)
            UQ=UQ+ene(ihist)*exp(dos(ihist)-ene(ihist)/(temp*kb)-lambdaQ)
            E2Q=E2Q+ene(ihist)**2*exp(dos(ihist)-ene(ihist)/(temp*kb)-lambdaQ)
!            FQ=FQ+exp(dos(ihist)-ene(ihist)/(tempC*kb)-lambda)
         end do
         U=U/Z
         E2=E2/Z
         C=((E2-U**2)/(kb*tempC**2*natom))/kb
         iS=iS+C/tempC*deltaT
         !F=-kb*temp*log(F)/lambda
         F=U-tempC*iS*natom*kb

         UQ=UQ/ZQ
         E2Q=E2Q/ZQ
         CQ=(temprescalegrad*tempC+temprescale)*((E2Q-UQ**2)/(kb*temp**2*natom))/kb
!         if (itemp==0) umin=U
!         UQ=(U-umin)*temprescale+umin
!         CQ=C*temprescale
         iSQ=iSQ+CQ/tempC*deltaT
         FQ=UQ-tempC*iSQ*natom*kb
         write(ofileno,'(1x,f8.2,2x,f6.4,1x,8g16.6)') tempC,temprescale,U,UQ,C,CQ,F,FQ,iS,iSQ
         !write(ofileno,'(1x,f12.4,10g20.6)') temp,U,U**2,E2,(E2-U**2)/temp**2,F,iS,lambda
         !print *,temp,lambda,maxval(dos)
      end do
      write(ofileno,'(1x,a)') " "
      !
      close(ofileno)
      !
   end subroutine wl_integration


   !---------------------------------------------------------------------------
   !> @brief
   !> Read input parameters.
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_parameters_wanglandau(ifile)

      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword
      integer :: rd_len, i_err, i_errb
      logical :: comment

      do
         10     continue
         ! Read file character for character until first whitespace
         keyword=""
         call bytereader(keyword,rd_len,ifile,i_errb)

         ! converting Capital letters
         call caps2small(keyword)

         ! check for comment markers (currently % and #)
         comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&
            (scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1.or.&
            (scan(trim(keyword),'!')==1))

         if (comment) then
            read(ifile,*)
         else
            ! Parse keyword
            keyword=trim(keyword)
            select case(keyword)
            ! This is the flags for the S(q,w)
            !> - wl_emin
            !! Lowest energy used for Wang-Landau sampling.
            !! Can be set in input or found automatically
         case('wl_emin')
            read(ifile,*,iostat=i_err) wl_emin
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - wl_emax
            !! Highest energy used for Wang-Landau sampling.
            !! Can be set in input or found automatically
         case('wl_emax')
            read(ifile,*,iostat=i_err) wl_emax
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - wl_gfac
            !! Prefactor for the acceptance rate (default=1.0)
         case('wl_gfac')
            read(ifile,*,iostat=i_err) wl_gfac
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - wl_sigma
            !! Broadening of kernel function (number of bin-widths)
         case('wl_sigma')
            read(ifile,*,iostat=i_err) wl_sigma
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - wl_lcut
            !! Lower cutoff of energy window
         case('wl_lcut')
            read(ifile,*,iostat=i_err) wl_lcut
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - wl_hcut
            !! Higher cutoff of energy window
         case('wl_hcut')
            read(ifile,*,iostat=i_err) wl_hcut
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !! Number of bins in Wang-Landau histogram
            !! If unset, the number is chosen so the bin width is 1 mRy (1 J if aunits Y)
         case('wl_nhist')
            read(ifile,*,iostat=i_err) wl_nhist
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - wl_maxT
            !! Maximum temperature for thermodynamic integration
         case('wl_maxt')
            read(ifile,*,iostat=i_err) wl_maxT
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - wl_minT
            !! Minimum temperature for thermodynamic integration
         case('wl_mint')
            read(ifile,*,iostat=i_err) wl_minT
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - wl_nloop
            !! Number of Wang-Landau histogram loops
            !! For every loop, scaling factor f is decreased to the square root of the
            !! previous value. f_0=e so f_10=1e-3, f_20=1e-6 and so on.
            !! Idealy, the error in the DOS scales as log(f)
         case('wl_nloop')
            read(ifile,*,iostat=i_err) wl_nloop
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         end select
      end if

      ! End of file
      if (i_errb==20) goto 20
      ! End of row
      if (i_errb==10) goto 10
   end do

   20  continue

   rewind(ifile)
   return
end subroutine read_parameters_wanglandau

subroutine wanglandau_init()
   !
   implicit none
   !
   !
   wl_gfac=1.0_dblprec
   wl_nloop=20
   wl_sigma=0.005
   wl_ntemp=100
   wl_lcut=0.975_dblprec
   wl_hcut=0.950_dblprec
   wl_minT=1.0_dblprec
   wl_maxT=1500.0_dblprec
   wl_emin=0.0_dblprec
   wl_emax=0.0_dblprec
   wl_nhist=0
   !
end subroutine wanglandau_init

end module WangLandau
