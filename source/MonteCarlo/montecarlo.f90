!-------------------------------------------------------------------------------
! MODULE: MonteCarlo
!> @brief
!> Data and routines for MonteCarlo simulations of either Ising or Heisenberg models
!> @author
!> A. Bergman, L. Bergqvist, Fan Pan, Jonathan Chico, J. Hellsvik + ...
!> @copyright
!> GNU Public License.
!> -----------------------------------------------------------------------------
module MonteCarlo
   use Parameters
   use Profiling
   use HamiltonianData, only : ham
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
   ! SUBROUTINE: mc_evolve
   ! DESCRIPTION
   !> @brief
   !> Main driver for Monte Carlo simulations using various implementations.
   !> For Heisenberg spin systems the options are either Metropolis,Heat bath or Wolff cluster.
   !> Based on the Hamiltonian \f$H=-\frac{1}{2} \sum_{ij} J_{ij} \mathbf{m}_i \cdot \mathbf{m}_j \f$.
   !> Ising model and spin-ice models implemented in Metropolis or Loop Algorithm.
   !> @author
   !> Lars Bergqvist
   !---------------------------------------------------------------------------------

   subroutine mc_evolve(Natom,Nchmax,Mensemble,nHam,temperature,temprescale,&
         mode,conf_num,lsf_metric,lsf_window,do_lsf,lsf_field,exc_inter,&
         lsf_interpolate,do_jtensor,do_dm, do_pd, do_biqdm,do_bq,do_ring,do_chir,do_sa,&
         mult_axis,iflip_a,emomM,emom,mmom,ind_mom_flag,&
         extfield,do_dip,Num_macro,max_num_atom_macro_cell,&
         cell_index,macro_nlistsize,macro_atom_nlist,emomM_macro,emom_macro,mmom_macro,do_anisotropy)
      !
      use RandomNumbers, only: rng_uniform,rng_uniformP,rng_gaussian, rng_gaussianP, use_vsl
      use LSF, only : mc_update_LSF
      use SpinIce , only: mc_update_spinice
      use montecarlo_common
      use InducedMoments, only : mc_update_ind_mom
      use Constants, only : mub,k_bolt
      use FieldData !,             only : allocation_field_time
      use OptimizationRoutines
      use AdaptiveTimeStepping
      use HamiltonianActions, only : effective_field
      use InputData, only : NA, N1, N2, N3, demag, ham_inp
      use DemagField
      use DipoleManager, only : dipole_field_calculation

            !
      implicit none

      !.. Input variables
      ! System variables
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: nHam      !< Number of atoms in Hamiltonian
      real(dblprec), intent(in) :: temperature !< Temperature
      real(dblprec), intent(in) :: temprescale !< Temperature rescaling according to QHB
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
      integer, intent(in) :: do_ring !< Add four-spin ring (4SR) term to Hamiltonian (0/1)
      ! CHIR variables
      integer, intent(in) :: do_chir   !< Add scalar chirality term (CHIR) to Hamiltonian (0/1)
      ! SA variables
      integer, intent(in) :: do_sa   !< Add Symmetric anisotropic (SA) term to Hamiltonian (0/1)
      ! Anisotropy variables
      integer, intent(in) :: do_anisotropy
      character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
      ! Moments variables
      integer, dimension(Natom),intent(in) :: iflip_a !< Flipping pattern
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      ! Induced moments variables
      character(len=1), intent(in) :: ind_mom_flag
      ! External fields variables
      real(dblprec), dimension(3), intent(inout) :: extfield !< External magnetic field
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
      integer :: i, k, icell

      real(dblprec) :: de !< Energy difference
      real(dblprec) :: cluster_size,macro_mag_trial,delta
      !.. Local arrays
      integer, dimension(Natom) :: visited_atoms
      real(dblprec), dimension(3) :: newmom !< New trial moment
      real(dblprec), dimension(3) :: macro_trial
      real(dblprec), dimension(3) :: totfield  !<Total effective field acting on each moment
      real(dblprec),dimension(natom,mensemble) :: flipprob_a ,newmmom_a,mflip
      real(dblprec),dimension(3,natom,mensemble) :: newmom_a,flipprob_g,flipprob_m

      integer, external :: omp_get_thread_num
      real(dblprec) :: henergy
      real(dblprec), dimension(3) :: loc_demag_fld, loc_mag_fld

      cluster_size=0.0_dblprec
      visited_atoms=0
      delta=(2.0/25.0)*(k_bolt*temperature/mub)**(0.20_dblprec)

      if (do_lsf=='Y') then
         call mc_update_LSF(Natom,Nchmax,Mensemble, do_lsf,emomM, emom, mmom, temperature, temprescale,  &
            extfield,mode,lsf_interpolate,lsf_field,lsf_window,lsf_metric,exc_inter,iflip_a,&
            ind_mom_flag,do_dip,Num_macro,mmom_macro,emom_macro,emomM_macro,do_anisotropy)

      ! Evolution of induced moments as described by Polesya et al.
      else if (ind_mom_flag=='Y') then
         call mc_update_ind_mom(Natom,Mensemble,iflip_a,&
            temperature,temprescale,mode,ham%max_no_neigh,ham%nlistsize,ham%nlist,ham%ncoup,conf_num,&
            mmom,emomM,emom,extfield,do_dip,Num_macro, emomM_macro,emom_macro,mmom_macro,&
            ham%ind_nlistsize,ham%ind_nlist,ham%ind_list_full,ham%sus_ind,&
            do_lsf,lsf_metric,ind_mom_flag,ham%max_no_neigh_ind)

      else
         ! Set up trial directions of magnetic moments
         if (mode=='L') then
            call mc_update_spinice(Natom, Mensemble, nHam, ham%max_no_neigh, conf_num, ham%ncoup, ham%ncoupD, ham%nlist, ham%nlistsize, ham%aham, &
               do_dm, ham%max_no_dmneigh, ham%dm_vect, ham%dmlist, ham%dmlistsize, do_pd, ham%nn_pd_tot, ham%pd_vect, ham%pdlist, ham%pdlistsize, &
               do_biqdm, ham%nn_biqdm_tot, ham%biqdm_vect, ham%biqdmlist, ham%biqdmlistsize, do_bq, ham%nn_bq_tot, ham%j_bq, ham%bqlist, ham%bqlistsize, &
               do_chir, ham%nn_chir_tot, ham%chir_coup, ham%chirlist, ham%chirlistsize, &
               do_sa, ham%max_no_saneigh, ham%sa_vect, ham%salist, ham%salistsize, &
               ham%taniso, ham%eaniso, ham%kaniso,ham%sb,emomM, emom, mmom, iflip_a, extfield, &
               mult_axis,ham%taniso_diff, ham%eaniso_diff, ham%kaniso_diff,ham%sb_diff, &
               do_dip, ham%Qdip,exc_inter,temperature,temprescale,ham%ind_nlistsize,ham%ind_nlist,ham%sus_ind,ind_mom_flag,ham%ind_list_full,ham%max_no_neigh_ind,&
               Num_macro,max_num_atom_macro_cell,cell_index,macro_nlistsize,macro_atom_nlist,&
               mmom_macro,emom_macro,emomM_macro,ham%Qdip_macro,do_anisotropy)
         else
            if (mode=='M'.or.mode=='H'.or.mode=='D') then
               call rng_uniformP(flipprob_m(:,:,:),3*natom*mensemble)
               call rng_gaussianP(flipprob_g(:,:,:),3*natom*mensemble,1.0_dblprec)

               if(use_vsl) then
#ifdef VSL
                  !$omp parallel do default(shared),private(i,k,newmom),schedule(auto),collapse(2)
#endif
                  do i=1,Natom
                     do k=1,mensemble
                        call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,delta,flipprob_m(:,i,k),flipprob_g(:,i,k))
                        newmom_a(1:3,i,k)=newmom(1:3)
                     enddo
                  enddo
#ifdef VSL
                  !$omp end parallel do
#endif
               else
                  do i=1,Natom
                     do k=1,mensemble
                        call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,delta,flipprob_m(:,i,k),flipprob_g(:,i,k))
                        newmom_a(1:3,i,k)=newmom(1:3)
                     enddo
                  enddo
               end if
               if(mode=='H') call rng_uniformP(mflip(:,:),natom*mensemble)
            else if (mode=='I') then
               ! First Ising random atom
               !$omp parallel do default(shared),private(i,k,newmom),schedule(auto),collapse(2)
               do i=1,Natom
                  do k=1,mensemble
                     call Ising_random_flip(emom,newmom,iflip_a(i),k,Natom,Mensemble)
                     newmom_a(:,iflip_a(i),k)=newmom(:)
                  enddo
               enddo
               !$omp end parallel do
            endif

            call rng_uniformP(flipprob_a(:,:),natom*mensemble)

            loc_demag_fld = 0.0_dblprec
            ! Wrapper for the calculation of the dipole-dipole interaction field
            ! The field is stored in the bfield array which then is passed to the main loop
            ! This is inefficient for the brute-force methods, but it the best way to ensure
            ! that the FFT approaches can be used in an appropriate way
            if (ham_inp%do_dip>0) then
               call timing(0,'Hamiltonian   ','OF')
               call timing(0,'Dipolar Int.  ','ON')
               beff = 0.0_dblprec
               call dipole_field_calculation(NA,N1,N2,N3,Natom,ham_inp%do_dip,Num_macro,          &
                  Mensemble,Natom,1,cell_index,macro_nlistsize,emomM,        &
                  emomM_macro,ham%Qdip,ham%Qdip_macro,henergy,beff)
               call timing(0,'Dipolar Int.  ','OF')
               call timing(0,'Hamiltonian   ','ON')
            endif
            ! Calculate energy and flip spin if preferred
            !$omp parallel do default(shared), private(i,k,de,totfield,loc_mag_fld), firstprivate(loc_demag_fld) schedule(auto),collapse(2)
            do i=1, Natom
               do k=1,mensemble
                  if (demag=='Y') loc_demag_fld = loc_demag_fld + demag_update(emomM(:,iflip_a(i),k))
                  loc_mag_fld = extfield + beff(:,iflip_a(i),k) + loc_demag_fld
                  if (mode=='H') then
                     ! Heat bath algorithm
                     call effective_field(iflip_a(i), k, totfield)
                     
                     call flip_h(Natom, Mensemble, emom, emomM, mmom(iflip_a(i),k), mmom(iflip_a(i),k), &
                        iflip_a(i),temperature,temprescale, k,flipprob_a(i,k),totfield, mflip(i,k))

                     if (demag=='Y') loc_demag_fld = loc_demag_fld - demag_update(emomM(:,iflip_a(i),k))
                  else
                     ! Metropolis algorithm, either in Ising or Loop Algorithm form
                     call calculate_energy(Natom, Mensemble, nHam, conf_num, do_dm , do_pd, do_biqdm, do_bq, do_ring, do_chir, do_sa,&
                         emomM, emom, mmom, iflip_a(i), newmom_a(1:3,iflip_a(i),k), loc_mag_fld, de, k, &
                         mult_axis, do_dip,Num_macro,max_num_atom_macro_cell,cell_index,macro_nlistsize,&
                         macro_atom_nlist,emomM_macro,icell,macro_mag_trial,macro_trial,exc_inter,do_anisotropy,do_jtensor)
                      !!!  call effective_field(Natom,Mensemble,iflip_a(i),iflip_a(i),emomM,   &
                      !!!     mmom,            &
                      !!!     external_field,time_external_field,beff,beff1,beff2,OPT_flag,     &
                      !!!     max_no_constellations,maxNoConstl,unitCellType,constlNCoup,       &
                      !!!     constellations,constellationsNeighType,de,    &
                      !!!     Num_macro,cell_index,emomM_macro,macro_nlistsize,  &
                      !!!     NA,N1,N2,N3)

                     if(mode=='D') then
                        call flip_g(Natom, Mensemble, emom, emomM, mmom, iflip_a(i),newmom_a(1:3,iflip_a(i),k),newmmom_a(iflip_a(i),k), &
                           de,temperature,temprescale,do_lsf,k,flipprob_a(i,k),lsf_metric,ham%ind_nlistsize,&
                           ham%ind_nlist,ind_mom_flag,ham%max_no_neigh,ham%sus_ind,&
                           do_dip,Num_macro,icell,macro_mag_trial,macro_trial,mmom_macro,emom_macro,emomM_macro)
                     else
                        call flip_a(Natom, Mensemble, emom, emomM, mmom, iflip_a(i),newmom_a(1:3,iflip_a(i),k),newmmom_a(iflip_a(i),k), &
                           de,temperature,temprescale,do_lsf,k,flipprob_a(i,k),lsf_metric,ham%ind_nlistsize,&
                           ham%ind_nlist,ind_mom_flag,ham%max_no_neigh_ind,ham%sus_ind,ham%ind_list_full,&
                           do_dip,Num_macro,icell,macro_mag_trial,macro_trial,mmom_macro,emom_macro,emomM_macro)
                     endif
                  endif
               enddo    !ensemble
            enddo     !atom
            !$omp end parallel do
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

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calculate_efield
   !> @brief
   !> Calculate total field for a single spin, used in the Heatbath algorithm
   !> @authors
   !> Lars Bergqvist
   !---------------------------------------------------------------------------------
   subroutine calculate_efield(Natom, Mensemble, conf_num,  do_dm,  do_pd, do_biqdm, do_bq, do_ring, do_chir, do_sa,&
         emomM, emom, mult_axis, iflip, extfield, do_lsf, k, totfield, exc_inter,do_anisotropy)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: k !< Current ensemble
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: conf_num    !< Number of configurations for LSF
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
      integer, intent(in) :: do_biqdm   !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      integer, intent(in) :: do_ring   !< Add four-spin ring exchange (4SR) term to Hamiltonian (0/1)      
      integer, intent(in) :: do_chir !< Add scalar chirality exchange to Hamiltonian (0/1)
      integer, intent(in) :: do_sa   !< Add symmetric anisotropic (SA) term to Hamiltonian (0/1)
      integer, intent(in) :: do_anisotropy
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      integer, intent(in) :: iflip !< Atom to flip spin for
      real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
      character(len=1), intent(in)  ::  do_lsf     !< Including LSF energy
      real(dblprec), dimension(3),intent(inout) :: totfield  !<Total effective field
      character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
      character(len=1), intent(in) :: exc_inter !< Interpolation of Jij between FM/DLM (Y/N)

      !.. Local scalars
      integer :: j, iflip_ham
      real(dblprec) :: tta, aw1,aw2
      real(dblprec) :: bqmdot,ringmdotkl,ringmdotkj,ringmdotjl,excscale
      real(dblprec) :: sxy, syz, szx
      integer :: im1,im2,ip1,ip2

      real(dblprec), dimension(3) :: fbef,faft
      !.. Executable statements

      ! First calculate effective field
      totfield(:) = 0.0_dblprec
      iflip_ham=ham%aham(iflip)

      if (exc_inter=='N') then
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422) && __INTEL_COMPILER < 1800
         !$omp simd reduction(+:totfield)
#endif
         do j=1,ham%nlistsize(iflip_ham)
            totfield(:) = totfield(:)+ ham%ncoup(j,iflip_ham,1)*emomM(:,ham%nlist(j,iflip),k)
         end do
      else
         totfield=0.0_dblprec
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422) && __INTEL_COMPILER < 1800
         !$omp simd private(excscale) reduction(+:totfield)
#endif
         do j=1,ham%nlistsize(iflip_ham)
            excscale=abs(sum(emom(:,ham%nlist(j,iflip),k)*emom(:,iflip,k)))
            totfield(:)=totfield(:)+((excscale*ham%ncoup(j,iflip_ham,1)+(1.0_dblprec-excscale)*ham%ncoupD(j,iflip_ham,1)))*emomM(:,ham%nlist(j,iflip),k)
         end do
      endif

      if (do_anisotropy==1) then
         ! Anisotropy
         ! Uniaxial anisotropy  scaled down to match heatbath
         ! Corrected Anisotropy 
         if (ham%taniso(iflip)==1) then
            tta=sum(emomM(:,iflip,k)*ham%eaniso(:,iflip))
            ! K1*(sin theta)^2
            totfield(1:3) = totfield(1:3)  &
               - 2.0_dblprec*ham%kaniso(1,iflip)*tta*ham%eaniso(1:3,iflip) &
               ! K2*(sin theta)^4
               - 4.0_dblprec*ham%kaniso(2,iflip)*(tta**2)*tta*ham%eaniso(1:3,iflip)
               ! Cubic anisotropy
         elseif (ham%taniso(iflip)==2) then
            ! K1*(sin theta)^2
            totfield(1) = totfield(1) &
            + 2*ham%kaniso(1,iflip)*emomM(1,iflip,k)*(emomM(2,iflip,k)**2+emomM(3,iflip,k)**2) &
            + 2*ham%kaniso(2,iflip)+emomM(1,iflip,k)*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2

            totfield(2) = totfield(2) &
               + 2*ham%kaniso(1,iflip)*emomM(2,iflip,k)*(emomM(3,iflip,k)**2+emomM(1,iflip,k)**2) &
               + 2*ham%kaniso(2,iflip)+emomM(2,iflip,k)*emomM(3,iflip,k)**2*emomM(1,iflip,k)**2

               totfield(3) = totfield(3) &
               + 2*ham%kaniso(1,iflip)*emomM(3,iflip,k)*(emomM(1,iflip,k)**2+emomM(2,iflip,k)**2) &
               + 2*ham%kaniso(2,iflip)+emomM(3,iflip,k)*emomM(1,iflip,k)**2*emomM(2,iflip,k)**2
         endif
         ! When both Cubic and Uniaxial are switched on
         if (ham%taniso(iflip)==7) then

            ! Uniaxial anisotropy
            tta=sum(emomM(:,iflip,k)*ham%eaniso(:,iflip))

            ! K1*(sin theta)^2
            totfield(1:3) = totfield(1:3)  &
               - 2.0_dblprec*ham%kaniso(1,iflip)*tta*ham%eaniso(1:3,iflip) &
               ! K2*(sin theta)^4
               - 4.0_dblprec*ham%kaniso(2,iflip)*(tta**2)*tta*ham%eaniso(1:3,iflip)
            ! Cubic anisotropy
            ! K1*(sin theta)^2
            aw1=ham%kaniso(1,iflip)*ham%sb(iflip)
            aw2=ham%kaniso(2,iflip)*ham%sb(iflip)
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
            if (ham%taniso_diff(iflip)==1) then
               tta=sum(emomM(:,iflip,k)*ham%eaniso_diff(:,iflip))

               ! K1*(sin theta)^2
               totfield(1:3) = totfield(1:3)  &
                  - 2.0_dblprec*ham%kaniso_diff(1,iflip)*tta*ham%eaniso_diff(1:3,iflip) &
                  ! K2*(sin theta)^4
                  - 4.0_dblprec*ham%kaniso_diff(2,iflip)*(tta**2)*tta*ham%eaniso_diff(1:3,iflip)
            ! Cubic anisotropy
            elseif (ham%taniso_diff(iflip)==2) then
               ! K1*(sin theta)^2
               totfield(1) = totfield(1) &
                  + 2*ham%kaniso_diff(1,iflip)*emomM(1,iflip,k)*(emomM(2,iflip,k)**2+emomM(3,iflip,k)**2) &
                  + 2*ham%kaniso_diff(2,iflip)+emomM(1,iflip,k)*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2

               totfield(2) = totfield(2) &
               + 2*ham%kaniso_diff(1,iflip)*emomM(2,iflip,k)*(emomM(3,iflip,k)**2+emomM(1,iflip,k)**2) &
               + 2*ham%kaniso_diff(2,iflip)+emomM(2,iflip,k)*emomM(3,iflip,k)**2*emomM(1,iflip,k)**2

               totfield(3) = totfield(3) &
                  + 2*ham%kaniso_diff(1,iflip)*emomM(3,iflip,k)*(emomM(1,iflip,k)**2+emomM(2,iflip,k)**2) &
                  + 2*ham%kaniso_diff(2,iflip)+emomM(3,iflip,k)*emomM(1,iflip,k)**2*emomM(2,iflip,k)**2

            endif
            ! When both Cubic and Uniaxial are switched on
            if (ham%taniso_diff(iflip)==7) then
               ! Uniaxial anisotropy
               tta=sum(emomM(:,iflip,k)*ham%eaniso_diff(:,iflip))

               ! K1*(sin theta)^2
               totfield(1:3) = totfield(1:3)  &
                  - 2.0_dblprec*ham%kaniso_diff(1,iflip)*tta*ham%eaniso_diff(1:3,iflip) &
                  ! K2*(sin theta)^4
                  - 4.0_dblprec*ham%kaniso_diff(2,iflip)*(tta**2)*tta*ham%eaniso_diff(1:3,iflip)
               ! Cubic anisotropy
               ! K1*(sin theta)^2
               aw1=ham%kaniso_diff(1,iflip)*ham%sb_diff(iflip)
               aw2=ham%kaniso_diff(2,iflip)*ham%sb_diff(iflip)
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
      endif
      ! DM interaction
      if(do_dm==1) then
         do j=1,ham%dmlistsize(iflip_ham)
            totfield(1) = totfield(1) + ham%dm_vect(3,j,iflip_ham)*emomM(2,ham%dmlist(j,iflip),k) -&
               ham%dm_vect(2,j,iflip_ham)*emomM(3,ham%dmlist(j,iflip),k)
            totfield(2) = totfield(2) + ham%dm_vect(1,j,iflip_ham)*emomM(3,ham%dmlist(j,iflip),k) -&
               ham%dm_vect(3,j,iflip_ham)*emomM(1,ham%dmlist(j,iflip),k)
            totfield(3) = totfield(3) + ham%dm_vect(2,j,iflip_ham)*emomM(1,ham%dmlist(j,iflip),k) -&
               ham%dm_vect(1,j,iflip_ham)*emomM(2,ham%dmlist(j,iflip),k)
         end do
      end if

      ! SA interaction
      if(do_sa==1) then
         do j=1,ham%salistsize(iflip_ham)
            totfield(1) = totfield(1) + ham%sa_vect(3,j,iflip_ham)*emomM(2,ham%salist(j,iflip),k) +&
               ham%sa_vect(2,j,iflip_ham)*emomM(3,ham%salist(j,iflip),k)
            totfield(2) = totfield(2) + ham%sa_vect(1,j,iflip_ham)*emomM(3,ham%salist(j,iflip),k) +&
               ham%sa_vect(3,j,iflip_ham)*emomM(1,ham%salist(j,iflip),k)
            totfield(3) = totfield(3) + ham%sa_vect(2,j,iflip_ham)*emomM(1,ham%salist(j,iflip),k) +&
               ham%sa_vect(1,j,iflip_ham)*emomM(2,ham%salist(j,iflip),k)
         end do
      end if

      ! PD interaction
      if(do_pd==1) then
        ! do j=1,ham%pdlistsize(iflip_ham)
        !    totfield(1) = totfield(1) + ham%pd_vect(1,j,iflip_ham)*emomM(1,ham%pdlist(j,iflip),k) +&
        !       ham%pd_vect(4,j,iflip_ham)*emomM(2,ham%pdlist(j,iflip),k) +&
        !       ham%pd_vect(5,j,iflip_ham)*emomM(3,ham%pdlist(j,iflip),k)
        !    totfield(2) = totfield(2) + ham%pd_vect(4,j,iflip_ham)*emomM(1,ham%pdlist(j,iflip),k) +&
        !       ham%pd_vect(2,j,iflip_ham)*emomM(2,ham%pdlist(j,iflip),k) +&
        !       ham%pd_vect(6,j,iflip_ham)*emomM(3,ham%pdlist(j,iflip),k)
        !    totfield(3) = totfield(3) + ham%pd_vect(5,j,iflip_ham)*emomM(1,ham%pdlist(j,iflip),k) +&
        !       ham%pd_vect(6,j,iflip_ham)*emomM(2,ham%pdlist(j,iflip),k) +&
        !       ham%pd_vect(3,j,iflip_ham)*emomM(3,ham%pdlist(j,iflip),k)
        ! end do

          do j=1,ham%pdlistsize(iflip_ham)
            totfield(1) = totfield(1) + ham%pd_vect(1,j,iflip_ham)*emomM(1,ham%pdlist(j,iflip),k) +&
               ham%pd_vect(2,j,iflip_ham)*emomM(2,ham%pdlist(j,iflip),k) +&
               ham%pd_vect(3,j,iflip_ham)*emomM(3,ham%pdlist(j,iflip),k)
            totfield(2) = totfield(2) + ham%pd_vect(4,j,iflip_ham)*emomM(1,ham%pdlist(j,iflip),k) +&
               ham%pd_vect(5,j,iflip_ham)*emomM(2,ham%pdlist(j,iflip),k) +&
               ham%pd_vect(6,j,iflip_ham)*emomM(3,ham%pdlist(j,iflip),k)
            totfield(3) = totfield(3) + ham%pd_vect(7,j,iflip_ham)*emomM(1,ham%pdlist(j,iflip),k) +&
               ham%pd_vect(8,j,iflip_ham)*emomM(2,ham%pdlist(j,iflip),k) +&
               ham%pd_vect(9,j,iflip_ham)*emomM(3,ham%pdlist(j,iflip),k)
         end do

      end if

      ! BIQDM interaction
      if(do_biqdm==1) then
         do j=1,ham%biqdmlistsize(iflip_ham)
            sxy = emomM(1,iflip,k)*emomM(2,ham%biqdmlist(j,iflip),k)-&
               emomM(2,iflip,k)*emomM(1,ham%biqdmlist(j,iflip),k)
            syz = emomM(2,iflip,k)*emomM(3,ham%biqdmlist(j,iflip),k)-&
               emomM(3,iflip,k)*emomM(2,ham%biqdmlist(j,iflip),k)
            szx = emomM(3,iflip,k)*emomM(1,ham%biqdmlist(j,iflip),k)-&
               emomM(1,iflip,k)*emomM(3,ham%biqdmlist(j,iflip),k)
            totfield(1) = totfield(1) + 2.0_dblprec*ham%biqdm_vect(1,j,iflip_ham)*(&
               szx*emomM(3,ham%biqdmlist(j,iflip),k)-&
               sxy*emomM(2,ham%biqdmlist(j,iflip),k))
            totfield(2) = totfield(2) + 2.0_dblprec*ham%biqdm_vect(1,j,iflip_ham)*(&
               sxy*emomM(1,ham%biqdmlist(j,iflip),k)-&
               syz*emomM(3,ham%biqdmlist(j,iflip),k))
            totfield(3) = totfield(3) + 2.0_dblprec*ham%biqdm_vect(1,j,iflip_ham)*(&
               syz*emomM(2,ham%biqdmlist(j,iflip),k)-&
               szx*emomM(1,ham%biqdmlist(j,iflip),k))
         end do
      end if

      ! BQ interaction scaled down to match Heatbath
      if(do_bq==1) then
         do j=1,ham%bqlistsize(iflip_ham)
            ! current spin
            !!!!!++ Lars changed the last emom to emomM
            bqmdot=emomM(1,ham%bqlist(j,iflip),k)*emomM(1,iflip,k)+&
                   emomM(2,ham%bqlist(j,iflip),k)*emomM(2,iflip,k)+&
                   emomM(3,ham%bqlist(j,iflip),k)*emomM(3,iflip,k)

!           bqmdot=emomM(1,ham%bqlist(j,iflip_ham),k)*emomM(1,iflip,k)+&
!              emomM(2,ham%bqlist(j,iflip_ham),k)*emomM(2,iflip,k)+&
!              emomM(3,ham%bqlist(j,iflip_ham),k)*emomM(3,iflip,k)
            totfield(:) = totfield(:) - 2.0_dblprec*ham%j_bq(j,iflip_ham)*bqmdot*emomM(:,ham%bqlist(j,iflip),k)
         end do
      end if
      
      fbef=totfield
      if(do_ring==1) then
       do j=1,ham%ringlistsize(iflip_ham)
          ! current spin
         ringmdotkl=emomM(1,ham%ringlist(iflip,j,2),k)*emomM(1,ham%ringlist(iflip,j,3),k)+&
                    emomM(2,ham%ringlist(iflip,j,2),k)*emomM(2,ham%ringlist(iflip,j,3),k)+&
                    emomM(3,ham%ringlist(iflip,j,2),k)*emomM(3,ham%ringlist(iflip,j,3),k)

         ringmdotkj=emomM(1,ham%ringlist(iflip,j,2),k)*emomM(1,ham%ringlist(iflip,j,1),k)+&
                    emomM(2,ham%ringlist(iflip,j,2),k)*emomM(2,ham%ringlist(iflip,j,1),k)+&
                    emomM(3,ham%ringlist(iflip,j,2),k)*emomM(3,ham%ringlist(iflip,j,1),k)

         ringmdotjl=emomM(1,ham%ringlist(iflip,j,1),k)*emomM(1,ham%ringlist(iflip,j,3),k)+&
                    emomM(2,ham%ringlist(iflip,j,1),k)*emomM(2,ham%ringlist(iflip,j,3),k)+&
                    emomM(3,ham%ringlist(iflip,j,1),k)*emomM(3,ham%ringlist(iflip,j,3),k)

         totfield(1) = totfield(1)-ham%j_ring(iflip_ham,j)*ringmdotkl*emomM(1,ham%ringlist(iflip,j,1),k)-&
                       ham%j_ring(iflip_ham,j)*ringmdotkj*emomM(1,ham%ringlist(iflip,j,3),k)+&
                       ham%j_ring(iflip_ham,j)*ringmdotjl*emomM(1,ham%ringlist(iflip,j,2),k)

         totfield(2) = totfield(2)-ham%j_ring(iflip_ham,j)*ringmdotkl*emomM(2,ham%ringlist(iflip,j,1),k)-&
                       ham%j_ring(iflip_ham,j)*ringmdotkj*emomM(2,ham%ringlist(iflip,j,3),k)+&
                       ham%j_ring(iflip_ham,j)*ringmdotjl*emomM(2,ham%ringlist(iflip,j,2),k)

         totfield(3) = totfield(3)-ham%j_ring(iflip_ham,j)*ringmdotkl*emomM(3,ham%ringlist(iflip,j,1),k)-&
                       ham%j_ring(iflip,j)*ringmdotkj*emomM(3,ham%ringlist(iflip,j,3),k)+&
                       ham%j_ring(iflip,j)*ringmdotjl*emomM(3,ham%ringlist(iflip,j,2),k)
       end do
    end if
    faft=totfield-fbef
    !   print '(i4,9f12.6)',iflip,faft-fbef,totfield


      ! "Chiral interaction term"
      if(do_chir==1) then
         do j=1,ham%chirlistsize(iflip_ham)
            im1=ham%chirlist(2,j,iflip)
            ip1=ham%chirlist(1,j,iflip)
            im2=ham%chirlist(2,j,im1)
            if(im2==0) im2=im1
            ip2=ham%chirlist(1,j,ip1)
            if(ip2==0) ip2=ip1
            !print '(2x,a,2i4,5x,5i4)','->  ', i,j,im2,im1,i,ip1,ip2
            totfield(1) = totfield(1)  &
               - ham%chir_coup(j,iflip_ham)*emomM(2,ip1,k)*emomM(3,im1,k) + ham%chir_coup(j,iflip_ham)*emomM(3,ip1,k)*emomM(2,im1,k) &
               - ham%chir_coup(j,iflip_ham)*emomM(2,im2,k)*emomM(3,im1,k) + ham%chir_coup(j,iflip_ham)*emomM(3,im2,k)*emomM(2,im1,k) &
               - ham%chir_coup(j,iflip_ham)*emomM(2,ip1,k)*emomM(3,ip2,k) + ham%chir_coup(j,iflip_ham)*emomM(3,ip1,k)*emomM(2,ip2,k)

            totfield(2) = totfield(2)  &
               - ham%chir_coup(j,iflip_ham)*emomM(3,ip1,k)*emomM(1,im1,k) + ham%chir_coup(j,iflip_ham)*emomM(1,ip1,k)*emomM(3,im1,k) &
               - ham%chir_coup(j,iflip_ham)*emomM(3,im2,k)*emomM(1,im1,k) + ham%chir_coup(j,iflip_ham)*emomM(1,im2,k)*emomM(3,im1,k) &
               - ham%chir_coup(j,iflip_ham)*emomM(3,ip1,k)*emomM(1,ip2,k) + ham%chir_coup(j,iflip_ham)*emomM(1,ip1,k)*emomM(3,ip2,k)

            totfield(3) = totfield(3)  &
               - ham%chir_coup(j,iflip_ham)*emomM(1,ip1,k)*emomM(2,im1,k) + ham%chir_coup(j,iflip_ham)*emomM(2,ip1,k)*emomM(1,im1,k) &
               - ham%chir_coup(j,iflip_ham)*emomM(1,im2,k)*emomM(2,im1,k) + ham%chir_coup(j,iflip_ham)*emomM(2,im2,k)*emomM(1,im1,k) &
               - ham%chir_coup(j,iflip_ham)*emomM(1,ip1,k)*emomM(2,ip2,k) + ham%chir_coup(j,iflip_ham)*emomM(2,ip1,k)*emomM(1,ip2,k)
         end do
      end if

      ! Add static field
      totfield(1:3) = totfield(1:3)+extfield(1:3)
      !
   end subroutine calculate_efield

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calculate_efield_tensor
   !> @brief
   !> Calculate total field of a single spin of tensorial exchange
   !---------------------------------------------------------------------------------
   subroutine calculate_efield_tensor(Natom, Mensemble, &
         do_biqdm, do_bq, emomM,mult_axis,iflip,extfield,k, totfield, do_dip,do_anisotropy)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom      !< Number of atoms in system
      integer, intent(in) :: Mensemble  !< Number of ensembles
      integer, intent(in) :: do_biqdm     !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_bq     !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      integer, intent(in) :: iflip !< Atom to flip spin for
      real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
      integer, intent(in) :: k !< Current ensemble
      real(dblprec), dimension(3), intent(inout) :: totfield !< Total effective field
      real(dblprec) :: aw1,aw2
      integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
      integer, intent(in) :: do_anisotropy
      character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
      !.. Local scalars
      integer :: j
      real(dblprec) :: tta
      real(dblprec) :: bqmdot
      real(dblprec) :: sxy, syz, szx

      !.. Local arrays
      real(dblprec), dimension(3) :: befftemp

      !.. Executable statements

      ! First calculate effective field
      totfield(:) = 0.0_dblprec
      befftemp(1) = 0.0_dblprec
      befftemp(2) = 0.0_dblprec
      befftemp(3) = 0.0_dblprec

      ! Exchange
      do j=1,ham%nlistsize(iflip)
         !befftemp = befftemp +matmul(ham%j_tens(:,:,j,iflip),emomM(:,ham%nlist(j,iflip),k))
         befftemp = befftemp +  ham%j_tens(1,:,j,iflip)*emomM(1,ham%nlist(j,iflip),k)  &
                             +  ham%j_tens(2,:,j,iflip)*emomM(2,ham%nlist(j,iflip),k)  &
                             +  ham%j_tens(3,:,j,iflip)*emomM(3,ham%nlist(j,iflip),k)
         !!! befftemp(1) = befftemp(1)+ ham%j_tens(1,1,j,iflip)*emomM(1,ham%nlist(j,iflip),k) + &
         !!!    ham%j_tens(1,2,j,iflip)*emomM(2,ham%nlist(j,iflip),k) + ham%j_tens(1,3,j,iflip)*emomM(3,ham%nlist(j,iflip),k)
         !!! befftemp(2) = befftemp(2)+ ham%j_tens(2,1,j,iflip)*emomM(1,ham%nlist(j,iflip),k) + &
         !!!    ham%j_tens(2,2,j,iflip)*emomM(2,ham%nlist(j,iflip),k) + ham%j_tens(2,3,j,iflip)*emomM(3,ham%nlist(j,iflip),k)
         !!! befftemp(3) = befftemp(3)+ ham%j_tens(3,1,j,iflip)*emomM(1,ham%nlist(j,iflip),k) + &
         !!!    ham%j_tens(3,2,j,iflip)*emomM(2,ham%nlist(j,iflip),k) + ham%j_tens(3,3,j,iflip)*emomM(3,ham%nlist(j,iflip),k)
      end do

      ! BIQDM interaction
      if(do_biqdm==1) then
         do j=1,ham%biqdmlistsize(iflip)
            sxy = emomM(1,iflip,k)*emomM(2,ham%biqdmlist(j,iflip),k)-&
               emomM(2,iflip,k)*emomM(1,ham%biqdmlist(j,iflip),k)
            syz = emomM(2,iflip,k)*emomM(3,ham%biqdmlist(j,iflip),k)-&
               emomM(3,iflip,k)*emomM(2,ham%biqdmlist(j,iflip),k)
            szx = emomM(3,iflip,k)*emomM(1,ham%biqdmlist(j,iflip),k)-&
               emomM(1,iflip,k)*emomM(3,ham%biqdmlist(j,iflip),k)
            befftemp(1) = befftemp(1) + 2.0_dblprec*ham%biqdm_vect(1,j,iflip)*(&
               szx*emomM(3,ham%biqdmlist(j,iflip),k)-&
               sxy*emomM(2,ham%biqdmlist(j,iflip),k))
            befftemp(2) = befftemp(2) + 2.0_dblprec*ham%biqdm_vect(1,j,iflip)*(&
               sxy*emomM(1,ham%biqdmlist(j,iflip),k)-&
               syz*emomM(3,ham%biqdmlist(j,iflip),k))
            befftemp(3) = befftemp(3) + 2.0_dblprec*ham%biqdm_vect(1,j,iflip)*(&
               syz*emomM(2,ham%biqdmlist(j,iflip),k)-&
               szx*emomM(1,ham%biqdmlist(j,iflip),k))
         end do
      end if

      ! BQ interaction
      if(do_bq==1) then
         do j=1,ham%bqlistsize(iflip)
            bqmdot=emomM(1,ham%bqlist(j,iflip),k)*emomM(1,iflip,k)+&
               emomM(2,ham%bqlist(j,iflip),k)*emomM(2,iflip,k)+&
               emomM(3,ham%bqlist(j,iflip),k)*emomM(3,iflip,k)

            !!!!!++ Jonathan used + instead of original -
            !!!!!++ Lars rescaled the 2.0 to 1.0
            befftemp(1) = befftemp(1)- 1.0_dblprec*ham%j_bq(j,iflip)*bqmdot*emomM(1,ham%bqlist(j,iflip),k)
            befftemp(2) = befftemp(2)- 1.0_dblprec*ham%j_bq(j,iflip)*bqmdot*emomM(2,ham%bqlist(j,iflip),k)
            befftemp(3) = befftemp(3)- 1.0_dblprec*ham%j_bq(j,iflip)*bqmdot*emomM(3,ham%bqlist(j,iflip),k)
         end do
      end if

      ! Anisotropy
      ! Nikos  added  also multi axis  in tensor form
      if (do_anisotropy==1) then
         ! Uniaxial anisotropy
         if (ham%taniso(iflip)==1) then
            tta=(emomM(1,iflip,k)*ham%eaniso(1,iflip)+emomM(2,iflip,k)*ham%eaniso(2,iflip)+emomM(3,iflip,k)*ham%eaniso(3,iflip))

            !!!!!++ Lars rescaled 2.0 and 4.0 to 1.0
            ! K1*(sin theta)^2
            ! Nikos updates Coeffs. These give correct fields
            befftemp(1:3) = befftemp(1:3) - 2.0_dblprec*ham%kaniso(1,iflip)*&
               tta*ham%eaniso(1:3,iflip) &

               ! K2*(sin theta)^4
               - 4.0_dblprec*ham%kaniso(2,iflip)*(tta**2)*tta*ham%eaniso(1:3,iflip)

            ! Cubic anisotropy
         elseif (ham%taniso(iflip)==2) then
            ! K1*(sin theta)^2
            befftemp(1) = befftemp(1) &
               + 2*ham%kaniso(1,iflip)*emomM(1,iflip,k)*(emomM(2,iflip,k)**2+emomM(3,iflip,k)**2) &
               + 2*ham%kaniso(2,iflip)+emomM(1,iflip,k)*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2

            befftemp(2) = befftemp(2) &
               + 2*ham%kaniso(1,iflip)*emomM(2,iflip,k)*(emomM(3,iflip,k)**2+emomM(1,iflip,k)**2) &
               + 2*ham%kaniso(2,iflip)+emomM(2,iflip,k)*emomM(3,iflip,k)**2*emomM(1,iflip,k)**2

            befftemp(3) = befftemp(3) &
               + 2*ham%kaniso(1,iflip)*emomM(3,iflip,k)*(emomM(1,iflip,k)**2+emomM(2,iflip,k)**2) &
               + 2*ham%kaniso(2,iflip)+emomM(3,iflip,k)*emomM(1,iflip,k)**2*emomM(2,iflip,k)**2

         endif

         ! When both Cubic and Uniaxial are switched on
         if (ham%taniso(iflip)==7) then

            ! Uniaxial anisotropy
            tta=(emomM(1,iflip,k)*ham%eaniso(1,iflip)+emomM(2,iflip,k)*ham%eaniso(2,iflip)+emomM(3,iflip,k)*ham%eaniso(3,iflip))
            ! K1*(sin theta)^2
            befftemp(1:3) = befftemp(1:3)  &
               - 2.00_dblprec*ham%kaniso(1,iflip)*&
               tta&
               *ham%eaniso(1:3,iflip) &
               ! K2*(sin theta)^4
            -4.00_dblprec*ham%kaniso(2,iflip)*&
               (tta**2)*tta*ham%eaniso(1:3,iflip)

            ! Cubic anisotropy
            ! K1*(sin theta)^2
            aw1=ham%kaniso(1,iflip)*ham%sb(iflip)
            aw2=ham%kaniso(2,iflip)*ham%sb(iflip)
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
         
         !!!   if (mult_axis=='Y') then
         !!!    ! Uniaxial anisotropy
         !!!    if (ham%taniso_diff(iflip)==1) then
         !!!       tta=sum(emomM(:,iflip,k)*ham%eaniso_diff(:,iflip))

         !!!       ! K1*(sin theta)^2
         !!!       totfield(1:3) = totfield(1:3)  &
         !!!          - 2.0_dblprec*ham%kaniso_diff(1,iflip)*tta*ham%eaniso_diff(1:3,iflip) &
         !!!          ! K2*(sin theta)^4
         !!!          - 4.0_dblprec*ham%kaniso_diff(2,iflip)*(tta**2)*tta*ham%eaniso_diff(1:3,iflip)
         !!!    ! Cubic anisotropy
         !!!    elseif (ham%taniso_diff(iflip)==2) then
         !!!       ! K1*(sin theta)^2
         !!!       totfield(1) = totfield(1) &
         !!!          + 2*ham%kaniso_diff(1,iflip)*emomM(1,iflip,k)*(emomM(2,iflip,k)**2+emomM(3,iflip,k)**2) &
         !!!          + 2*ham%kaniso_diff(2,iflip)+emomM(1,iflip,k)*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2

         !!!       totfield(2) = totfield(2) &
         !!!       + 2*ham%kaniso_diff(1,iflip)*emomM(2,iflip,k)*(emomM(3,iflip,k)**2+emomM(1,iflip,k)**2) &
         !!!       + 2*ham%kaniso_diff(2,iflip)+emomM(2,iflip,k)*emomM(3,iflip,k)**2*emomM(1,iflip,k)**2

         !!!       totfield(3) = totfield(3) &
         !!!          + 2*ham%kaniso_diff(1,iflip)*emomM(3,iflip,k)*(emomM(1,iflip,k)**2+emomM(2,iflip,k)**2) &
         !!!          + 2*ham%kaniso_diff(2,iflip)+emomM(3,iflip,k)*emomM(1,iflip,k)**2*emomM(2,iflip,k)**2

         !!!    endif
         !!!    ! When both Cubic and Uniaxial are switched on
         !!!    if (ham%taniso_diff(iflip)==7) then
         !!!       ! Uniaxial anisotropy
         !!!       tta=sum(emomM(:,iflip,k)*ham%eaniso_diff(:,iflip))

         !!!       ! K1*(sin theta)^2
         !!!       totfield(1:3) = totfield(1:3)  &
         !!!          - 2.0_dblprec*ham%kaniso_diff(1,iflip)*tta*ham%eaniso_diff(1:3,iflip) &
         !!!          ! K2*(sin theta)^4
         !!!          - 4.0_dblprec*ham%kaniso_diff(2,iflip)*(tta**2)*tta*ham%eaniso_diff(1:3,iflip)
         !!!       ! Cubic anisotropy
         !!!       ! K1*(sin theta)^2
         !!!       aw1=ham%kaniso_diff(1,iflip)*ham%sb_diff(iflip)
         !!!       aw2=ham%kaniso_diff(2,iflip)*ham%sb_diff(iflip)
         !!!       totfield(1) = totfield(1) &
         !!!       + 2*aw1*emomM(1,iflip,k)*(emomM(2,iflip,k)**2+emomM(3,iflip,k)**2) &
         !!!       + 2*aw2+emomM(1,iflip,k)*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2

         !!!       totfield(2) = totfield(2) &
         !!!          + 2*aw1*emomM(2,iflip,k)*(emomM(3,iflip,k)**2+emomM(1,iflip,k)**2) &
         !!!          + 2*aw2+emomM(2,iflip,k)*emomM(3,iflip,k)**2*emomM(1,iflip,k)**2

         !!!       totfield(3) = totfield(3) &
         !!!          + 2*aw1*emomM(3,iflip,k)*(emomM(1,iflip,k)**2+emomM(2,iflip,k)**2) &
         !!!          + 2*aw2+emomM(3,iflip,k)*emomM(1,iflip,k)**2*emomM(2,iflip,k)**2
         !!!    endif
         !!! endif
         
      endif
         ! Dipolar Interaction Jonathan 19-07-2012
!        if(present(ham%Qdip)) then
            if(do_dip==1) then
               do j=1,Natom
                  !          do k=1,Mensemble
                  befftemp(1) = befftemp(1) + ham%Qdip(1,1,j,iflip)*emomM(1,j,k) + ham%Qdip(2,1,j,iflip)*emomM(2,j,k) + ham%Qdip(3,1,j,iflip)*emomM(3,j,k)
                  befftemp(2) = befftemp(2) + ham%Qdip(1,2,j,iflip)*emomM(1,j,k) + ham%Qdip(2,2,j,iflip)*emomM(2,j,k) + ham%Qdip(3,2,j,iflip)*emomM(3,j,k)
                  befftemp(3) = befftemp(3) + ham%Qdip(1,3,j,iflip)*emomM(1,j,k) + ham%Qdip(2,3,j,iflip)*emomM(2,j,k) + ham%Qdip(3,3,j,iflip)*emomM(3,j,k)
                  !          end do
               end do
            end if
!        end if

         !-----------------------------------------------------------------------------------------------------------------------

         ! Add static field
         totfield(1:3) = totfield(1:3) + befftemp(1:3)+extfield(1:3)

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
         write (filn,'(''moment.'',a,''.out'')') trim(simid)
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
            if (.not. allocated(iflip_a)) then
               allocate(iflip_a(natom),stat=i_stat)
               call memocc(i_stat,product(shape(iflip_a))*kind(iflip_a),'iflip_a','allocate_mcdata')
            end if
         else
            i_all=-product(shape(iflip_a))*kind(iflip_a)
            deallocate(iflip_a,stat=i_stat)
            call memocc(i_stat,i_all,'iflip_a','allocate_mcdata')
         endif
      end subroutine allocate_mcdata

   end module MonteCarlo
