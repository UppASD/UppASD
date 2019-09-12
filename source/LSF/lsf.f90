!-------------------------------------------------------------------------------
! MODULE: LSF
!> @brief
!> Collection of LSF routines including Interpolation to obtain the moment size dependent LSF energy and exchange parameters
!> @author
!> Fan Pan and Lars Bergqvist
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module LSF
   use Parameters
   use Profiling
   use inputdata, only : ammom_inp,nch
   use ChemicalData, only : achtype
   use HamiltonianData, only : ham

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
   subroutine mc_update_LSF(Natom,Nchmax,Mensemble,nHam, conf_num,do_lsf,emomM, emom, mmom, temperature, temprescale,  &
         extfield,mult_axis,mode,lsf_interpolate,lsf_field,lsf_window,lsf_metric,exc_inter,iflip_a,&
         ind_mom_flag,do_dip,Num_macro,mmom_macro,emom_macro,emomM_macro,do_anisotropy)

      !
      use RandomNumbers, only: rng_uniform,rng_uniformP, rng_gaussian, rng_gaussianP, use_vsl
      use montecarlo_common
      use Constants,only : mub,k_bolt

      !.. Implicit declarations
      implicit none
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: nHam  !< Number of atoms in Hamiltonian
      integer, intent(in) :: Nchmax  !< Number of chemical type
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: conf_num  !< number of configurations for LSF
      character(len=1), intent(in)  ::  do_lsf     !< Including LSF energy
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
      real(dblprec), intent(in) :: temperature !< Temperature
      real(dblprec), intent(in) :: temprescale !< Temperature rescaling according to QHB
      real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
      character(len=1), intent(in) :: mult_axis !< Multiple uniaxial anisotropies
      character(len=1) :: mode !< Simulation mode (M=MC, H=MC Heat Bath)
      character(len=1), intent(in)  ::  lsf_interpolate     !< Interpolate LSF or not
      character(len=1), intent(in)  ::  lsf_field           !< LSF field contribution (Local/Total)
      real(dblprec),intent(in) :: lsf_window                !< Range of moment variation in LSF
      integer,intent(in) :: lsf_metric                !< LSF metric in phase space integration (1=Murata-Doniach,2=Jacobian)
      character(len=1), intent(in)  ::  exc_inter           !< Exchange interpolation between FM/DLM (Y/N)
      integer,dimension(natom),intent(in) :: iflip_a !< Flipping pattern
      character(len=1), intent(in) :: ind_mom_flag
      integer, intent(in) :: do_dip
      integer, intent(in) :: Num_macro
      integer, intent(in) :: do_anisotropy
      real(dblprec), dimension(Num_macro,Mensemble), intent(inout) :: mmom_macro
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emom_macro
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emomM_macro

      !.. Local scalars
      !
      integer :: i,ip,k,icell
      real(dblprec) :: de, newmmom !< New trial magnitude of moment
      real(dblprec) :: macro_mag_trial,delta
      !
      !.. Local arrays
      !
      integer, dimension(Natom,mensemble) :: iin
      real(dblprec), dimension(3) :: newmom !< New trial moment
      real(dblprec), dimension(3) :: totfield  !<Total effective field acting on each moment
      real(dblprec), dimension(3) :: macro_trial

      real(dblprec), dimension(Natom,mensemble) :: newmmom_a,flipprob_a,mflip
      real(dblprec), dimension(Natom,mensemble,2) :: rn
      real(dblprec),dimension(3,natom,mensemble) :: newmom_a,flipprob_g,flipprob_m

      delta=(2.0/25.0)*(k_bolt*temperature/mub)**(0.20_dblprec)

      call rng_uniformP(flipprob_m,3*natom*mensemble)
      call rng_uniformP(rn,natom*mensemble*2)
      call rng_gaussianP(flipprob_g,3*natom*mensemble,1.0_dblprec)
      iin=floor(3.0_dblprec*rn(:,:,1))
#ifdef VSL
      !$omp parallel do default(shared),private(i,k,newmom,newmmom),schedule(auto),collapse(2)
#endif
      do i=1,Natom
         do k=1,mensemble
            if (iin(i,k)==0) then
               call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,delta,flipprob_m(:,i,k),flipprob_g(:,i,k))
               newmom_a(1:3,i,k)=newmom(1:3)
               newmmom_a(i,k) = mmom(i,k)
            elseif(iin(i,k)==1) then
               newmom_a(1:3,i,k)=emom(1:3,i,k)
               call vary_moment_magnitude(Nchmax,i,mmom(i,k),newmmom, &
                  ammom_hlim,ammom_llim,lsf_window,delta,rn(i,k,2))
               newmmom_a(i,k) = newmmom
            elseif(iin(i,k)==2) then
               call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,delta,flipprob_m(:,i,k),flipprob_g(:,i,k))
               newmom_a(1:3,i,k)=newmom(1:3)
               call vary_moment_magnitude(Nchmax,i,mmom(i,k),newmmom, &
                  ammom_hlim,ammom_llim,lsf_window,delta,rn(i,k,2))
               newmmom_a(i,k) = newmmom
            else
               write(*,*) 'ERROR'
            endif
         enddo
      enddo
#ifdef VSL
      !$omp end parallel do
#endif
      if(mode=='H') call rng_uniformP(mflip,natom*mensemble)
      call rng_uniformP(flipprob_a,natom*mensemble)
      ! Calculate energy and flip spin if preferred
      !$omp parallel do default(shared), private(i,k,de,totfield), schedule(auto),collapse(2)
      do i=1, Natom
         do k=1, Mensemble
            if (mode=='H') then
               call calculate_field_wLSF(Natom,Nchmax,Mensemble, nHam, conf_num,emomM, emom, mmom, iflip_a(i),&
                  newmmom_a(iflip_a(i),k),extfield, k,lsf_interpolate,lsf_field,totfield,exc_inter,do_anisotropy)
               call flip_h(Natom, Mensemble, emom, emomM, newmmom_a(iflip_a(i),k), mmom(iflip_a(i),k), &
                  iflip_a(i),temperature,temprescale, k,flipprob_a(i,k),totfield,mflip(i,k))
            else
               call calculate_energy_wLSF(Natom,Nchmax,Mensemble,nHam, conf_num,emomM,emom,mmom,&
                  iflip_a(i),newmom_a(1:3,iflip_a(i),k),newmmom_a(iflip_a(i),k),extfield,de,&
                  k,lsf_interpolate,lsf_field,exc_inter,do_anisotropy)
            endif
            if(mode=='D') then
               call flip_g(Natom,Mensemble,emom,emomM,mmom,iflip_a(i),newmom_a(1:3,iflip_a(i),k),&
                  newmmom_a(iflip_a(i),k),de,temperature,temprescale,do_lsf,k,flipprob_a(i,k),&
                  lsf_metric,ham%ind_nlistsize,ham%ind_nlist,ind_mom_flag,ham%max_no_neigh,ham%sus_ind,&
                  do_dip,Num_macro,icell,macro_mag_trial,macro_trial,mmom_macro,emom_macro,emomM_macro)
            elseif (mode=='M') then
               call flip_a(Natom, Mensemble,emom,emomM, mmom,iflip_a(i),newmom_a(1:3,iflip_a(i),k),&
                  newmmom_a(iflip_a(i),k),de,temperature,temprescale,do_lsf,k,flipprob_a(i,k),&
                  lsf_metric,ham%ind_nlistsize,ham%ind_nlist,ind_mom_flag,ham%max_no_neigh,ham%sus_ind,ham%ind_list_full,&
                  do_dip,Num_macro,icell,macro_mag_trial,macro_trial,mmom_macro,emom_macro,emomM_macro)
            endif
         enddo
      enddo
      !$omp end parallel do

   end subroutine mc_update_LSF

   !---------------------------------------------------------------------------
   !> @brief
   !! Metropolis Monte Carlo inclsuing LSF
   !> @author
   !! Lars Bergqvist and Fan Pan
   !--------------------------------------------------------------------------
   subroutine calculate_energy_wLSF(Natom,Nchmax,Mensemble,nHam,conf_num,emomM,emom,mmom,iflip,&
         newmom,newmmom, extfield,de,k,lsf_interpolate,lsf_field,exc_inter,do_anisotropy)

      use Constants, only : mub,mry
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Nchmax  !< Number of chemical type
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: nHam  !< Number of atoms in Hamiltonian
      integer, intent(in) :: conf_num   !< Number of configurations for LSF
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
      integer, intent(in) :: do_anisotropy
      real(dblprec) :: aw1,aw2
      real(dblprec),dimension(2,2) :: lsf_i   !< lsf energy for a single site
      real(dblprec),dimension(Nchmax) :: nbsum     !< sum of the moments in all neighbours
      real(dblprec),dimension(3) :: nbsumfield     !< sum of the moments in all neighbours
      real(dblprec), dimension(Nchmax,2) :: ammom_inter  !< moments on grids of  each of the ch_type
      real(dblprec), dimension(ham%max_no_neigh,2) :: ncoup_i  !<interpolated ncoup only at flipping site
      real(dblprec), dimension(ham%max_no_neigh,2) :: fs_ncoup_i !< first-shell  ncoup  on  flipping site
      !.. Local scalars
      integer :: inttype              !< (0/1) pick the nearest grids or do interpolation
      integer :: i, j, iflip_h
      real(dblprec) :: tt, tta, ttb, e_c, e_t, fc

      !.. Local arrays
      integer(dblprec),dimension(nchmax) :: counter     ! counting for every Nchtype

      real(dblprec), dimension(3) :: beff_c,beff_t, trialmom
      real(dblprec), dimension(3) :: fs_beff_c, fs_beff_t  !< beff generated from first shell

      iflip_h=ham%aham(iflip)
      e_c=0.0_dblprec; e_t=0.0_dblprec; tt=0.0_dblprec; lsf_i=0.0_dblprec; beff_c=0.0_dblprec; beff_t=0.0_dblprec
      fc=mry/mub
      fs_beff_c=0.0_dblprec; fs_beff_t=0.0_dblprec
      trialmom(:)=newmom(:)*newmmom
      inttype=1
      if (lsf_interpolate=='L') then
         inttype=0
      elseif(lsf_interpolate=='H') then
         inttype=2
      elseif(lsf_interpolate=='I') then
         inttype=3
      endif
      !-----------------------------
      !calculate B_i field
      !find the average of the another ch_type surrounded to the center m
      nbsum(:) = 0.0_dblprec ; nbsumfield(:)=0.0_dblprec
      counter(:) = 0 ; ammom_inter(:,:) = 0.0_dblprec 
#if _OPENMP && ( defined __INTEL_COMPILER )
      !DIR$ LOOP COUNT min(8)
#endif
      do j = 1, ham%nlistsize(iflip_h)
         nbsum(achtype(ham%nlist(j,iflip))) = nbsum(achtype(ham%nlist(j,iflip))) &
            + mmom(ham%nlist(j,iflip),k)
         counter(achtype(ham%nlist(j,iflip))) = counter(achtype(ham%nlist(j,iflip))) + 1
      enddo
      do i = 1,Nchmax
         if (counter(i) == 0) then
            ammom_inter(i,1) = minval(ammom_int(i,:))
            ammom_inter(i,2) = minval(ammom_int(i,:))
         else
            ammom_inter(i,1) = nbsum(i)/counter(i)
            ammom_inter(i,2) = nbsum(i)/counter(i)
         endif
      enddo
      ammom_inter(achtype(iflip),1) = mmom(iflip,k)
      ammom_inter(achtype(iflip),2) = newmmom
      call do_interpolation_ncoup_and_lsf(Natom,Mensemble,nHam,Nchmax,conf_num, &
         ammom_inter,iflip,ncoup_i,lsf_i,k,inttype,exc_inter,achtype(iflip),emom,2)
      ! effective field of the nearest neighbours surrounding site i
      if (lsf_field == 'L') then
         do j=1,ham%fs_nlistsize(iflip_h)
            fs_ncoup_i(j,:) = ncoup_i(ham%nind(j,iflip),:)
            fs_beff_c(:) = fs_beff_c(:) + fs_ncoup_i(j,1)*emomM(:,ham%fs_nlist(j,iflip),k)
            fs_beff_t(:) = fs_beff_t(:) + fs_ncoup_i(j,2)*emomM(:,ham%fs_nlist(j,iflip),k)
         enddo
      endif
      ! Find effective field surrounding site i

!      if (exc_inter=='N') then
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
      !DIR$ LOOP COUNT min(8)
      !$omp simd reduction(+:beff_c,beff_t)
#endif
         do j=1,ham%nlistsize(iflip_h)
!            beff_c(:) = beff_c(:)+ ham%ncoup(j,iflip_h,11 )*emomM(:,ham%nlist(j,iflip),k)
!            beff_t(:) = beff_t(:)+ ham%ncoup(j,iflip_h,11 )*emomM(:,ham%nlist(j,iflip),k)
            beff_c(:) = beff_c(:)+ ncoup_i(j,1)*emomM(:,ham%nlist(j,iflip),k)
            beff_t(:) = beff_t(:)+ ncoup_i(j,2)*emomM(:,ham%nlist(j,iflip),k)
         end do
!     else
!#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
!      !DIR$ LOOP COUNT min(8)
!      !$omp simd
!#endif
!         do j=1,ham%nlistsize(iflip_h)
!           beff_c(:)=beff_c(:)+((excscale*ham%ncoup(j,iflip_h,11 )+(1._dblprec-excscale)*ham%ncoupD(j,iflip_h,11 )))* &
!               emomM(:,ham%nlist(j,iflip),k)
!           beff_t(:)=beff_t(:)+((excscale*ham%ncoup(j,iflip_h,11 )+(1._dblprec-excscale)*ham%ncoupD(j,iflip_h,11 )))* &
!               emomM(:,ham%nlist(j,iflip),k)
!         enddo
!      endif

      ! LSF energy + exchange
      if(lsf_field=='L') then
         e_c=e_c+lsf_i(1,1)*fc-lsf_i(2,1)*sum(emomM(:,iflip,k)*fs_beff_c(:))
         e_t=e_t+lsf_i(1,2)*fc-lsf_i(2,2)*sum(trialmom(:)*fs_beff_t(:))
      else
         e_c=e_c+lsf_i(1,1)*fc-lsf_i(2,1)*sum(emomM(:,iflip,k)*beff_c(:))
         e_t=e_t+lsf_i(1,2)*fc-lsf_i(2,2)*sum(trialmom(:)*beff_t(:))
      endif
      ! Anisotropy
      if (do_anisotropy==1) then
         ! Uniaxial anisotropy
         if (ham%taniso(iflip)==1) then
            tta=sum(emomM(:,iflip,k)*ham%eaniso(:,iflip))
            ttb=sum(trialmom(:)*ham%eaniso(:,iflip))
            e_c=e_c+ham%kaniso(1,iflip)*(tta**2)+ham%kaniso(2,iflip)*(tta**4)
            e_t=e_t+ham%kaniso(1,iflip)*(ttb**2)+ham%kaniso(2,iflip)*(ttb**4)
            ! Cubic anisotropy
         elseif (ham%taniso(iflip)==2) then
            e_c=e_c+ham%kaniso(1,iflip)*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2+ &
               emomM(2,iflip,k)**2*emomM(3,iflip,k)**2+emomM(3,iflip,k)**2*emomM(1,iflip,k)**2)+&
               ham%kaniso(2,iflip)*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2)
               e_t=e_t+ham%kaniso(1,iflip)*(trialmom(1)**2*trialmom(2)**2+ &
               trialmom(2)**2*trialmom(3)**2+trialmom(3)**2*trialmom(1)**2)+ &
               ham%kaniso(2,iflip)*(trialmom(1)**2*trialmom(2)**2*trialmom(3)**2)
         endif
         ! When both Cubic and Uniaxial are switched on
         if (ham%taniso(iflip)==7) then
            ! Uniaxial anisotropy
            tta=sum(emomM(:,iflip,k)*ham%eaniso(:,iflip))
            ttb=sum(trialmom(:)*ham%eaniso(1,iflip))
            e_c=e_c+ham%kaniso(1,iflip)*(tta**2)+ham%kaniso(2,iflip)*(tta**4)
            e_t=e_t+ham%kaniso(1,iflip)*(ttb**2)+ham%kaniso(2,iflip)*(ttb**4)

            ! Cubic anisotropy
            aw1=ham%kaniso(1,iflip)*ham%sb(iflip)
            aw2=ham%kaniso(2,iflip)*ham%sb(iflip)

            e_c=e_c+aw1*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2+ &
            emomM(2,iflip,k)**2*emomM(3,iflip,k)**2+emomM(3,iflip,k)**2*emomM(1,iflip,k)**2)+&
            aw2*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2)
            e_t=e_t+aw1*(trialmom(1)**2*trialmom(2)**2+ &
            trialmom(2)**2*trialmom(3)**2+trialmom(3)**2*trialmom(1)**2)+ &
            aw2*(trialmom(1)**2*trialmom(2)**2*trialmom(3)**2)
         endif
      endif
      !Energy difference
      tt=(e_t-e_c)
      de=mub*tt

   end subroutine calculate_energy_wLSF


   !---------------------------------------------------------------------------
   !> @brief
   !! Metropolis Monte Carlo inclsuing LSF
   !> @author
   !! Lars Bergqvist and Fan Pan
   !--------------------------------------------------------------------------
   subroutine calculate_field_wLSF(Natom,Nchmax,Mensemble,nHam,conf_num,&
         emomM,emom,mmom,iflip,newmmom,extfield,k,lsf_interpolate,lsf_field,totfield,&
         exc_inter,do_anisotropy)

      use Constants, only : mub,mry
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Nchmax  !< Number of chemical type
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: nHam  !< Number of atoms in system
      integer, intent(in) :: conf_num   !< Number of configurations for LSF
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      integer, intent(in) :: iflip !< Atom to flip spin for
      real(dblprec), intent(in) :: newmmom !< New trial magnitude of moment
      real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
      integer, intent(in) :: k !< Current ensemble
      character(len=1),intent(in) :: lsf_interpolate !< Interpolate LSF or not
      character(len=1),intent(in) :: lsf_field       !< LSF field contribution (Local/Total)
      real(dblprec), dimension(3), intent(out) :: totfield !< Total field
      character(len=1),intent(in) :: exc_inter    !< Rescaling of exchange interactions (Y/N)
      integer, intent(in) :: do_anisotropy
      real(dblprec) :: aw1,aw2
      real(dblprec),dimension(2) :: lsf_t   !< lsf energy for a single site
      real(dblprec),dimension(Nchmax) :: nbsum     !< sum of the moments in all neighbours
      real(dblprec),dimension(3) :: nbsumfield     !< sum of the moments in all neighbours
      real(dblprec), dimension(Nchmax) :: ammom_inter  !< moments on grids of  each of the ch_type
      real(dblprec), dimension(ham%max_no_neigh) :: ncoup_t, ncoup_tg  !<interpolated ncoup only at flipping site
      !.. Local scalars
      integer :: inttype              !< (0/1) pick the nearest grids or do interpolation
      integer :: i, j, iflip_h
      real(dblprec) :: tt, tta, ttb, fc, excscale, lsf_tf

      !.. Local arrays
      integer(dblprec),dimension(nchmax) :: counter     ! counting for every Nchtype
      !    real(dblprec),dimension(3) :: excfield1,excfield2,lsffield

      iflip_h=ham%aham(iflip)
      tt=0.0_dblprec; totfield=0.0_dblprec
      fc=mry/mub
      inttype=1
      if (lsf_interpolate=='L') then
         inttype=0
      elseif(lsf_interpolate=='H') then
         inttype=2
      elseif(lsf_interpolate=='I') then
         inttype=3
      endif
      !-----------------------------
      !calculate B_i field
      !find the average of the another ch_type surrounded to the center m
      nbsum(:) = 0.0_dblprec; nbsumfield(:)=0.0_dblprec
      counter(:) = 0 ; ammom_inter(:) = 0.0_dblprec ; excscale=1.0_dblprec
#if _OPENMP && ( defined __INTEL_COMPILER )
      !DIR$ LOOP COUNT min(8)
#endif
      do j = 1, ham%nlistsize(iflip_h)
         nbsum(achtype(ham%nlist(j,iflip))) = nbsum(achtype(ham%nlist(j,iflip))) &
            + mmom(ham%nlist(j,iflip),k)
         counter(achtype(ham%nlist(j,iflip))) = counter(achtype(ham%nlist(j,iflip))) + 1
      enddo
      do i = 1,Nchmax
         if (counter(i) == 0) then
            ammom_inter(i) = minval(ammom_int(i,:))
         else
            ammom_inter(i) = nbsum(i)/counter(i)
         endif
      enddo
      ammom_inter(achtype(iflip)) = newmmom
      call do_interpolation_ncoup_and_lsf_gradient(Natom,Mensemble,nHam,Nchmax,conf_num, &
         ammom_inter,iflip,ncoup_t,ncoup_tg,lsf_t,lsf_tf,k,inttype,exc_inter,achtype(iflip),emom)

      ! Find effective field surrounding site i
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
      !DIR$ LOOP COUNT min(8)
      !!$omp simd reduction(+:totfield)
#endif
      do j=1,ham%nlistsize(iflip_h)
!         totfield(:) = totfield(:)+ lsf_t(2)*ham%ncoup(j,iflip_h,11)*emomM(:,ham%nlist(j,iflip),k)
         totfield(:) = totfield(:)+ lsf_t(2)*ncoup_t(j)*emomM(:,ham%nlist(j,iflip),k)
         totfield(:) = totfield(:)+ ncoup_tg(j)*sum(emomM(:,iflip,k)*emomM(:,ham%nlist(j,iflip),k))*emom(:,iflip,k)
      end do
      ! LSF energy
      totfield(:)=totfield(:)+lsf_tf*emom(:,iflip,k)*fc

      ! Anisotropy
      if (do_anisotropy==1) then
         ! Uniaxial anisotropy  scaled down to match heatbath
         if (ham%taniso(iflip)==1) then
            tta=sum(emomM(:,iflip,k)*ham%eaniso(:,iflip))

            ! K1*(sin theta)^2
            totfield(1:3) = totfield(1:3)  &
               - ham%kaniso(1,iflip)*tta*ham%eaniso(1:3,iflip) &
               ! K2*(sin theta)^4
               - ham%kaniso(2,iflip)*(tta**2)*tta*ham%eaniso(1:3,iflip)
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
               - ham%kaniso(1,iflip)*tta*ham%eaniso(1:3,iflip) &
               ! K2*(sin theta)^4
               - ham%kaniso(2,iflip)*(tta**2)*tta*ham%eaniso(1:3,iflip)
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
      endif
      ! Add static field
      totfield(1:3) = totfield(1:3)+extfield(1:3)

   end subroutine calculate_field_wLSF

   !---------------------------------------------------------------------------
   !> @brief
   !> Vary the magnitude of the moment
   !> @author
   !! Lars Bergqvist and Fan Pan
   !---------------------------------------------------------------------------
   subroutine vary_moment_magnitude(Nchmax,ip,oldmmom,newmmom,ammom_hlim,ammom_llim,wind,delta,rn)

      use RandomNumbers, only : rng_gaussian,rng_uniform
      use Constants

      implicit none

      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, intent(in) :: ip  !< site number
      real(dblprec), intent(in) :: oldmmom !< chosen magnitude of moment to be varied
      real(dblprec), intent(in), dimension(Nchmax) :: ammom_hlim,ammom_llim
      real(dblprec), intent(out) :: newmmom !< New trial magnitude of moment
      real(dblprec), intent(in) :: wind !< window for changing moment
      real(dblprec), intent(in):: delta !< broadening 
      real(dblprec), intent(in):: rn !< Random number 

      real(dblprec) :: locwind
      real(dblprec), dimension(1) :: gasrun,tempmom !< Gaussian RN
      integer :: check

      if (wind >= 1.0_dblprec) then
         locwind=(ammom_hlim(achtype(ip))-ammom_llim(achtype(ip)))*delta
         check=1
         do while (check==1)
            call rng_gaussian(gasrun,1,locwind)
            check=0
            tempmom = oldmmom + gasrun
            if (tempmom(1)>ammom_hlim(achtype(ip))) check=1
            if (tempmom(1)<ammom_llim(achtype(ip))) check=1
         end do
         newmmom=tempmom(1)
      elseif(wind >= 0.0_dblprec .and. wind < 1.0_dblprec) then
         locwind=(ammom_hlim(achtype(ip))-ammom_llim(achtype(ip)))*wind
         check=1
         do while (check==1)
            call rng_gaussian(gasrun,1,locwind+1e-15_dblprec)
            check=0
            tempmom = oldmmom + gasrun
            if (tempmom(1)>ammom_hlim(achtype(ip))) check=1
            if (tempmom(1)<ammom_llim(achtype(ip))) check=1
         end do
         newmmom=tempmom(1)
      elseif (wind < 0.0_dblprec) then
         newmmom = (ammom_hlim(achtype(ip))-ammom_llim(achtype(ip))-2.0_dblprec*dbl_tolerance)*rn+ammom_llim(achtype(ip))+dbl_tolerance
      else
         write(*,*) 'LSF_window must be specified'
      endif
   end subroutine vary_moment_magnitude


   !---------------------------------------------------------------------------
   !> @brief
   !> Extract total energy in LSF
   !> @author
   !! Lars Bergqvist
   !> @todo
   !> Reinstate site resolved energies
   !---------------------------------------------------------------------------
   subroutine totalenergy_LSF(Natom, Nchmax, Mensemble, nHam, conf_num, emom, emomM, mmom,simid, &
         plotenergy, mstep, extfield,eenergy, aenergy, fenergy, lsfenergy, exc_inter, &
         do_lsf,inttype,lsf_field,Temp,site_energy)
      use Constants

      implicit none


      integer, intent(in) :: mstep   !< Current simulation step
      integer, intent(in) :: Natom   !< Number of atoms in system
      integer, intent(in) :: Nchmax  !< Number of chemical type
      integer, intent(in) :: Mensemble    !< Number of ensembles
      integer, intent(in) :: nHam  !< Number of atoms in Hamiltonian
      integer, intent(in) :: conf_num  !< number of configurations for LSF
      integer, intent(in) :: plotenergy   !< Calculate and plot energy (0/1)
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom     !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom    !< Current magnetic moment
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: extfield !< External magnetic field
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in)  ::  do_lsf     !< Including LSF energy
      character(len=1),intent(in) :: exc_inter !< Interpolation of Jij between FM/DLM
      integer, intent(in) :: inttype !< Interpolatation type in LSF
      character(len=1), intent(in) :: lsf_field       !< LSF field contribution (Local/Total)
      real(dblprec), intent(in) :: Temp               !< Temperature
      real(dblprec),dimension(10,natom,Mensemble),intent(inout) :: site_energy

      !.. Subroutine output
      real(dblprec), dimension(Mensemble), intent(out) :: aenergy !< Anisotropy energy
      real(dblprec), dimension(Mensemble), intent(out) :: eenergy !< Total exchange (transversal) energy
      real(dblprec), dimension(Mensemble), intent(out) :: fenergy !< Total external energy
      real(dblprec), dimension(Mensemble), intent(out) :: lsfenergy !< Total LSF (longitudinal) energy

      !...Local variables
      integer :: i, j, k,i_h
      real(dblprec) :: fcinv,fc
      real(dblprec) :: excscale !< Interpolation parameter FM/DLM
      real(dblprec) :: ieenergy, exptemp
      real(dblprec) :: xu1,xu2

      !...Local arrays
      real(dblprec), dimension(Mensemble) :: tt
      real(dblprec), dimension(3,Mensemble) :: ttv
      real(dblprec), dimension(Mensemble) :: tempk1, tempk2
      real(dblprec), dimension(Mensemble) :: aeatom
      real(dblprec),dimension(2,1) :: lsfE
      real(dblprec),dimension(Nchmax) :: nbsum     !< sum of the moments in all neighbours
      real(dblprec), dimension(Nchmax,1) :: ammom_inter  !< moments on grids of  each of the ch_type
      real(dblprec), dimension(ham%max_no_neigh,1) :: ncoup_i  !<interpolated ncoup only at flipping site
      real(dblprec), dimension(ham%max_no_neigh,1) :: fs_ncoup  !< first-shell  ncoup  on  flipping site
      real(dblprec), dimension(3) :: fs_beff
      integer(dblprec),dimension(nchmax) :: counter     ! counting for every Nchtype

      ! Factor for energy scale
      fcinv = mub/mry
      fc = mry/mub

      !$omp parallel do default(shared),private(i,j,k,i_h,nbsum,counter,ammom_inter,lsfE, &
      !$omp tempk1,tempk2,xu1,xu2,ncoup_i,ttv,tt,aeatom,ieenergy,exptemp),&
      !$omp reduction(+:lsfenergy,eenergy,aenergy,fenergy),schedule(auto),collapse(2)
      do i=1,natom
         do k=1,mensemble
            i_h=ham%aham(i)
            nbsum(:) = 0.0_dblprec ; counter(:) = 0 ; ttv(:,:)=0.0_dblprec
            do j = 1, ham%nlistsize(i_h)
               nbsum(achtype(ham%nlist(j,i))) = nbsum(achtype(ham%nlist(j,i))) &
                  + mmom(ham%nlist(j,i),k)
               counter(achtype(ham%nlist(j,i))) = counter(achtype(ham%nlist(j,i))) + 1
            enddo
            do j = 1,Nchmax
               if (counter(j) == 0) then
                  ammom_inter(j,1) = minval(ammom_int(j,:))
               else
                  ammom_inter(j,1) = nbsum(j)/counter(j)
               endif
            enddo
            ammom_inter(achtype(i),1) = mmom(i,k)
            call do_interpolation_ncoup_and_lsf(Natom,Mensemble,nHam,Nchmax,conf_num, &
               ammom_inter,i,ncoup_i,lsfE,k,inttype,exc_inter,achtype(i),emom,1)
            !> effective field of the nearest neighbours surrounding site i
            if (lsf_field=='L') then
               fs_beff=0.0_dblprec
               do j=1,ham%fs_nlistsize(i_h)
                  fs_ncoup(j,1) = ncoup_i(ham%nind(j,i),1)
                  fs_beff(:) = fs_beff(:) + fs_ncoup(j,1)*emomM(:,ham%fs_nlist(j,i),k)
               enddo
            end if
!#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422) && __INTEL_COMPILER >= 1800
!            !DIR$ LOOP COUNT min(8)
!            !$omp simd reduction(+:ttv)
!#endif
            do j = 1, ham%nlistsize(i_h)
               ttv(:,k) = ttv(:,k)+ ncoup_i(j,1)*emomM(:,ham%nlist(j,i),k)
            enddo

            lsfenergy(k)=lsfenergy(k)+2.0_dblprec*lsfE(1,1)*fc
            if (plotenergy==2) site_energy(9,i,k)=lsfE(1,1)*fc
            if (lsf_field=='T') then
               ieenergy=-lsfE(2,1)*sum(emomM(:,i,k)*ttv(:,k))
               eenergy(k)=eenergy(k)+ieenergy
               if(plotenergy==2) site_energy(1,i,k)=0.5_dblprec*ieenergy
            else
               ieenergy=-sum(emomM(:,i,k)*ttv(:,k))-lsfE(2,1)*sum(emomM(:,i,k)*fs_beff(:))
               eenergy(k)=eenergy(k)+ieenergy
               if(plotenergy==2) site_energy(1,i,k)=0.5_dblprec*ieenergy
            endif

            ! External field energy
            fenergy(k)=fenergy(k)-sum(extfield(:,i,k)*emomM(:,i,k))
            if(plotenergy==2) site_energy(8,i,k)=-sum(extfield(:,i,k)*emomM(:,i,k))

            if (ham%taniso(i)==1) then
               ! Calculate uniaxial anisotropy energy
               tt(k) = ham%eaniso(1,i)*emomM(1,i,k)+ham%eaniso(2,i)*emomM(2,i,k)+ham%eaniso(3,i)*emomM(3,i,k)
               aeatom(k) = (ham%kaniso(1,i)*tt(k)**2) + ham%kaniso(2,i)*(tt(k)**2)**2
               aenergy(k) = aenergy(k)+aeatom(k)
               if(plotenergy==2) site_energy(7,i,k)=(ham%kaniso(1,i)*tt(k)**2) + ham%kaniso(2,i)*(tt(k)**2)**2
            elseif (ham%taniso(i)==2) then
               ! Calculate cubic anisotropy energy
               tempk1(k) = emomM(1,i,k)**2*emomM(2,i,k)**2 + emomM(2,i,k)**2*emomM(3,i,k)**2 +&
                  emomM(3,i,k)**2*emomM(1,i,k)**2
               tempk2(k) = emomM(1,i,k)**2 * emomM(2,i,k)**2 * emomM(3,i,k)**2
               aeatom(k) = -(ham%kaniso(1,i)*tempk1(k) + ham%kaniso(2,i)*tempk2(k))
               aenergy(k) = aenergy(k)+aeatom(k)
               if(plotenergy==2) site_energy(7,i,k)=-(ham%kaniso(1,i)*tempk1(k) + ham%kaniso(2,i)*tempk2(k))
            endif
         enddo  ! End loop(Mensemble)
      enddo
      !$omp end parallel do

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

      do i = 1,conf_num
         num=1
         do j=1,na
            do k=num,nch(j)
               ammom_int(k,i) = ammom_inp(j,k,i)
            enddo
            num=nch(j)+1
         enddo
      end do

      do j=1,nchmax
         mask=.false.
         do i=1,conf_num
            num=count(ammom_int(j,i)==ammom_int(j,:))
            if (num==1) then
               mask(i)=.true.
            else
               if(.not. any(ammom_int(j,i)==ammom_int(j,:) .and. mask) ) mask(i)=.true.
            endif
         enddo

         mint_tmp(:,j)=pack(ammom_int(j,:), mask)
         ngrids(j)=count(mask)

      enddo

      !< find max/min value in ammom_int
      do i = 1, Nchmax
         ammom_llim(i) = minval(ammom_int(i,:))
         ammom_hlim(i) = maxval(ammom_int(i,:))
      end do

      allocate(mom_int(maxval(ngrids),nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(mom_int))*kind(mom_int),'mom_int','LSF_datareshape')

      mom_int(:,:)=mint_tmp(1:maxval(ngrids),:)

   end subroutine LSF_datareshape

   !> Interpolation of LSF energy and exchange interactions to given moment size from existing grid
   subroutine do_interpolation_ncoup_and_lsf(Natom,Mensemble,nHam, Nchmax,conf_num,ammom_inter, &
                   iflip,itp_noup,obj,k,inttype,exc_inter,chemtype,emom,nstep)

      !
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles in system
      integer, intent(in) :: nHam !< Number of atoms in Hamiltonian
      integer, intent(in) :: Nchmax  !< Number of chemical type
      integer, intent(in) :: k !< Current ensemble
      integer, intent(in) :: conf_num   !< Number of configurations for LSF
      real(dblprec), intent(in), dimension(nchmax,nstep) :: ammom_inter !< moments at grids
      integer, intent(in) :: iflip !< Atom to flip spin for
      real(dblprec), dimension(ham%max_no_neigh,nstep),intent(out) :: itp_noup  !< Interpolated ncoup
      real(dblprec),dimension(2,nstep), intent(out) :: obj !< the thing to be interpolated (possibly LSF_energy)
      integer, intent(in) :: inttype !< (0/1) pick the nearest grids or do interpolation
      character(len=1),intent(in) :: exc_inter !< Exchange rescaling
      integer,intent(in) :: chemtype !< Chemical type of the trial moment
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      integer,intent(in) :: nstep
      !local variables
      integer, dimension(nchmax) :: ind  !< the lower-limit to the moment on grids
      integer, dimension(2**nchmax) :: ac  !< array of conf-num for the interpolating mash
      integer :: i, j, iconf, ii, jj,istep, iflip_h
      real(dblprec),dimension(4) :: temp
      real(dblprec) :: invtemp
      real(dblprec), dimension(ham%max_no_neigh,nstep) :: tmp_noup  !< Interpolated ncoup
      real(dblprec) :: excscale !< Interpolaton parameter FM/DLM

      iflip_h=ham%aham(iflip)
      do istep=1,nstep

         if(inttype==1) then          ! do interpolation
            if (Nchmax== 1) then
               iconf = maxloc (mom_int(:,1),1,mom_int(:,1) .LT. ammom_inter(1,istep))
               ac(1) = iconf
               ac(2) = iconf +1
               invtemp =(ammom_inter(1,istep)-mom_int(ac(1),1))/ (mom_int(ac(2),1)-mom_int(ac(1),1))
               if (exc_inter=='N') then
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
                  !dir$ loop count min(8)
                  !$omp simd
#endif
                  do j=1,ham%nlistsize(iflip_h)
                     itp_noup(j,istep)=ham%ncoup(j,iflip_h,ac(1))+invtemp*(ham%ncoup(j,iflip_h,ac(2))-ham%ncoup(j,iflip_h,ac(1)))
                  enddo
                  obj(:,istep)=LSF_energy(ac(1),:)+invtemp*(LSF_energy(ac(2),:)-LSF_energy(ac(1),:))

               elseif (exc_inter=='Y') then
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
                  !dir$ loop count min(8)
                  !$omp simd private(excscale)
#endif
                  do j=1,ham%nlistsize(iflip_h)
                     excscale=abs(sum(emom(:,ham%nlist(j,iflip),k)*emom(:,iflip,k)))
                     itp_noup(j,istep)=(ham%ncoupD(j,iflip_h,ac(1))+excscale*(ham%ncoup(j,iflip_h,ac(1))-ham%ncoupD(j,iflip_h,ac(1))))+ &
                        invtemp*((ham%ncoupD(j,iflip_h,ac(2))+excscale*(ham%ncoup(j,iflip_h,ac(2))-ham%ncoupD(j,iflip_h,ac(2))))- &
                        (ham%ncoupD(j,iflip_h,ac(1))+excscale*(ham%ncoup(j,iflip_h,ac(1))-ham%ncoupD(j,iflip_h,ac(1)))))
                  enddo
                  obj(:,istep)=LSF_energy(ac(1),:)+invtemp*(LSF_energy(ac(2),:)-LSF_energy(ac(1),:))
               endif
            elseif (Nchmax == 2) then
               do i=1,nchmax
                  ind(i) = maxloc (mom_int(:,i),1,mom_int(:,i) .LT. ammom_inter(i,istep))
               enddo
               iconf=(ind(2)-1)*ngrids(1)+ind(1)
               ii=ind(1) ; jj=ind(2)
               if(chemtype==1) then
                  temp(1)=(mom_int(jj+1,2)-ammom_inter(2,istep))/(mom_int(jj+1,2)-mom_int(jj,2))
                  temp(2)=(ammom_inter(2,istep)-mom_int(jj,2))/(mom_int(jj+1,2)-mom_int(jj,2))
                  temp(3)=(mom_int(ii+1,1)-ammom_inter(1,istep))/(mom_int(ii+1,1)-mom_int(ii,1))
                  temp(4)=(ammom_inter(1,istep)-mom_int(ii,1))/(mom_int(ii+1,1)-mom_int(ii,1))
                  ac(1) = iconf
                  ac(2) = iconf + ngrids(1)
                  ac(3) = iconf +1
                  ac(4) = iconf + ngrids(1) + 1
               else
                  temp(1)=(mom_int(ii+1,1)-ammom_inter(1,istep))/(mom_int(ii+1,1)-mom_int(ii,1))
                  temp(2)=(ammom_inter(1,istep)-mom_int(ii,1))/(mom_int(ii+1,1)-mom_int(ii,1))
                  temp(3)=(mom_int(jj+1,2)-ammom_inter(2,istep))/(mom_int(jj+1,2)-mom_int(jj,2))
                  temp(4)=(ammom_inter(2,istep)-mom_int(jj,2))/(mom_int(jj+1,2)-mom_int(jj,2))
                  ac(1) = iconf
                  ac(2) = iconf +1
                  ac(3) = iconf + ngrids(1)
                  ac(4) = iconf + ngrids(1) + 1
               endif
               if (exc_inter=='N' ) then
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
                  !dir$ loop count min(8)
                  !$omp simd
#endif
                  do j=1,ham%nlistsize(iflip_h)
                     itp_noup(j,istep)= (ham%ncoup(j,iflip_h,ac(1))*temp(1)+ham%ncoup(j,iflip_h,ac(2))*temp(2))*temp(3)+ &
                        (ham%ncoup(j,iflip_h,ac(3))*temp(1)+ham%ncoup(j,iflip_h,ac(4))*temp(2))*temp(4)
                  end do
                  obj(:,istep)=(LSF_energy(ac(1),:)*temp(1)+LSF_energy(ac(2),:)*temp(2))*temp(3)+&
                     (LSF_energy(ac(3),:)*temp(1)+LSF_energy(ac(4),:)*temp(2))*temp(4)
               elseif (exc_inter=='Y') then
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
                  !dir$ loop count min(8)
                  !$omp simd private(excscale)
#endif
                  do j=1,ham%nlistsize(iflip_h)
                     excscale=abs(sum(emom(:,ham%nlist(j,iflip),k)*emom(:,iflip,k)))
                     itp_noup(j,istep)=(ham%ncoup(j,iflip_h,ac(1))*temp(1)+ham%ncoup(j,iflip_h,ac(2))*temp(2))*temp(3)+ &
                        (ham%ncoup(j,iflip_h,ac(3))*temp(1)+ham%ncoup(j,iflip_h,ac(4))*temp(2))*temp(4)


                     tmp_noup(j,istep)= (ham%ncoupD(j,iflip_h,ac(1))*temp(1)+ham%ncoupD(j,iflip_h,ac(2))*temp(2))*temp(3)+ &
                        (ham%ncoupD(j,iflip_h,ac(3))*temp(1)+ham%ncoupD(j,iflip_h,ac(4))*temp(2))*temp(4)
                     itp_noup(j,istep)=tmp_noup(j,istep)+excscale*(itp_noup(j,istep)-tmp_noup(j,istep))
                  end do
                  obj(:,istep)=(LSF_energy(ac(1),:)*temp(1)+LSF_energy(ac(2),:)*temp(2))*temp(3)+&
                     (LSF_energy(ac(3),:)*temp(1)+LSF_energy(ac(4),:)*temp(2))*temp(4)
               endif
            else
               write(*,*) "three component or higher dimension interpolation not being implemented"
            end if
         elseif(inttype==3) then !poor mans interpolation, average between low and high grid point
            if (nchmax == 1) then
               ac(1)  = maxloc (mom_int(:,1),1,mom_int(:,1) .LT. ammom_inter(1,istep))
               ac(2)  = ac(1)+1
               temp(1)= (mom_int(ac(1),1)/ammom_inter(1,istep))**2
               temp(2)= (ammom_inter(1,istep)/mom_int(ac(2),1))**2
            else
               do i=1,nchmax
                  ind(i) = maxloc (mom_int(:,i),1,mom_int(:,i) .LT. ammom_inter(i,istep))
               enddo
               ac(1)  = (ind(2)-1)*ngrids(1)+ind(1)
               ac(2)  = ac(1)+ngrids(1)+1
               temp(1)= ((mom_int(ind(1),1)*mom_int(ind(2),2))/(ammom_inter(1,istep)*ammom_inter(2,istep)))**1
               temp(2)=((ammom_inter(1,istep)*ammom_inter(2,istep))/(mom_int(ind(1)+1,1)*mom_int(ind(2)+1,2)))**1
            endif
            if (exc_inter=='N') then
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
               !dir$ loop count min(8)
               !$omp simd
#endif
               do j=1,ham%nlistsize(iflip_h)
                  itp_noup(j,istep)=ham%ncoup(j,iflip_h,ac(1))*temp(1)
                  tmp_noup(j,istep)=ham%ncoup(j,iflip_h,ac(2))*temp(2)
               enddo
               itp_noup(:,istep)=0.5_dblprec*(itp_noup(:,istep)+tmp_noup(:,istep))
               obj(:,istep) = 0.5_dblprec*(LSF_energy(ac(1),:)*temp(1)+LSF_energy(ac(2),:)*temp(2))
            elseif (exc_inter=='Y') then
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
               !dir$ loop count min(8)
               !$omp simd private(excscale)
#endif
               do j=1,ham%nlistsize(iflip_h)
                  excscale=abs(sum(emom(:,ham%nlist(j,iflip),k)*emom(:,iflip,k)))
                  itp_noup(j,istep)=(ham%ncoupD(j,iflip_h,ac(1))+excscale*(ham%ncoup(j,iflip_h,ac(1))-ham%ncoupD(j,iflip_h,ac(1))))*temp(1)
                  tmp_noup(j,istep)=(ham%ncoupD(j,iflip_h,ac(2))+excscale*(ham%ncoup(j,iflip_h,ac(2))-ham%ncoupD(j,iflip_h,ac(2))))*temp(2)
               enddo
               itp_noup(:,istep)=0.5_dblprec*(itp_noup(:,istep)+tmp_noup(:,istep))
               obj(:,istep) = 0.5_dblprec*(LSF_energy(ac(1),:)*temp(1)+LSF_energy(ac(2),:)*temp(2))
            endif
         else          ! no interpolation, pick nearest lower grid point
            if (nchmax == 1) then
               iconf = maxloc (mom_int(:,1),1,mom_int(:,1) .LT. ammom_inter(1,istep))
               if(inttype==0)  then        ! no interpolation, pick nearest lower grid point
                  invtemp=(mom_int(iconf,1)/ammom_inter(1,istep))**2
               else
                  iconf=iconf+1
                  invtemp=(ammom_inter(1,istep)/mom_int(iconf,1))**2
               endif
            else
               do i=1,nchmax
                  ind(i) = maxloc (mom_int(:,i),1,mom_int(:,i) .LT. ammom_inter(i,istep))
               enddo
               iconf=(ind(2)-1)*ngrids(1)+ind(1)
               if(inttype==0)  then        ! no interpolation, pick nearest lower grid point
                  invtemp=((mom_int(ind(1),1)*mom_int(ind(2),2))/(ammom_inter(1,istep)*ammom_inter(2,istep)))**1
               else
                  ind(1)=ind(1)+1
                  ind(2)=ind(2)+1
                  iconf=iconf+ngrids(1)+1
                  invtemp=((ammom_inter(1,istep)*ammom_inter(2,istep))/(mom_int(ind(1),1)*mom_int(ind(2),2)))**1
                endif
            endif
            if (exc_inter=='N') then
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
               !dir$ loop count min(8)
               !$omp simd
#endif
               do j=1,ham%nlistsize(iflip_h)
                  itp_noup(j,istep)=ham%ncoup(j,iflip_h,iconf)*invtemp
               enddo
               obj(:,istep) = LSF_energy(iconf,:)*invtemp
            elseif (exc_inter=='Y') then
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
               !dir$ loop count min(8)
               !$omp simd private(excscale)
#endif
               do j=1,ham%nlistsize(iflip_h)
                  excscale=abs(sum(emom(:,ham%nlist(j,iflip),k)*emom(:,iflip,k)))
                  itp_noup(j,istep)=(ham%ncoupD(j,iflip_h,iconf)+excscale*(ham%ncoup(j,iflip_h,iconf)-ham%ncoupD(j,iflip_h,iconf)))*invtemp
               enddo
               obj(:,istep) = LSF_energy(iconf,:)*invtemp
            endif
         endif
      enddo
   end subroutine do_interpolation_ncoup_and_lsf


   !> Interpolation of LSF energy and exchange interactions to given moment size from existing grid
   subroutine do_interpolation_ncoup_and_lsf_gradient(Natom,Mensemble,nHam,Nchmax,conf_num,ammom_inter, &
         iflip,itp_noup,itp_noupg,obj,lsff,k,inttype,exc_inter,chemtype,emom)

      !
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles in system
      integer, intent(in) :: nHam!< Number of atoms in Hamiltonian
      integer, intent(in) :: Nchmax  !< Number of chemical type
      integer, intent(in) :: k !< Current ensemble
      integer, intent(in) :: conf_num   !< Number of configurations for LSF
      real(dblprec), intent(in), dimension(:) :: ammom_inter !< moments at grids
      integer, intent(in) :: iflip !< Atom to flip spin for
      real(dblprec), dimension(ham%max_no_neigh),intent(out) :: itp_noup  !< Interpolated ncoup
      real(dblprec), dimension(ham%max_no_neigh),intent(out) :: itp_noupg  !< Interpolated ncoup
      real(dblprec),dimension(2), intent(out) :: obj !< the thing to be interpolated (possibly LSF_energy)
      real(dblprec), intent(out) :: lsff !< lsf field
      integer, intent(in) :: inttype !< (0/1) pick the nearest grids or do interpolation
      character(len=1),intent(in) :: exc_inter !< Exchange rescaling
      integer,intent(in) :: chemtype !< Chemical type of the trial moment
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      !local variables
      integer, dimension(nchmax) :: ind  !< the lower-limit to the moment on grids
      integer, dimension(nchmax) :: indp  !< gradient interpolation index
      integer, dimension(2**nchmax) :: ac  !< array of conf-num for the interpolating mash
      integer :: i, j, iconf, ii, jj, iflip_h
      real(dblprec),dimension(4) :: temp
      real(dblprec) :: invtemp,invtemp2
      real(dblprec), dimension(ham%max_no_neigh,2) :: ncouptmp  !< Interpolated ncoup on boundary
      real(dblprec), dimension(2,2) :: lsftmp
      real(dblprec) :: excscale !< Interpolaton parameter FM/DLM

      iflip_h=ham%aham(iflip)

      if(inttype==1) then          ! do interpolation
         if (Nchmax== 1) then
            iconf = maxloc (mom_int(:,1),1,mom_int(:,1) .LT. ammom_inter(1))
            ac(1) = iconf
            ac(2) = iconf +1
            invtemp =(ammom_inter(1)-mom_int(ac(1),1))/ (mom_int(ac(2),1)-mom_int(ac(1),1))
            invtemp2 =1.0_dblprec / (mom_int(ac(2),1)-mom_int(ac(1),1))
            if (exc_inter=='N') then

#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
               !dir$ loop count min(8)
               !$omp simd
#endif
               do j=1,ham%nlistsize(iflip_h)
                  itp_noup(j)=ham%ncoup(j,iflip_h,ac(1))+invtemp*(ham%ncoup(j,iflip_h,ac(2))-ham%ncoup(j,iflip_h,ac(1)))
                  itp_noupg(j)=invtemp2*(ham%ncoup(j,iflip_h,ac(2))-ham%ncoup(j,iflip_h,ac(1)))
               enddo
               obj(:)=LSF_energy(ac(1),:)+invtemp*(LSF_energy(ac(2),:)-LSF_energy(ac(1),:))
               lsff=-invtemp2*(LSF_energy(ac(2),1)-LSF_energy(ac(1),1))
            elseif (exc_inter=='Y') then
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
               !dir$ loop count min(8)
               !$omp simd private(excscale)
#endif
               do j=1,ham%nlistsize(iflip_h)
                  excscale=abs(sum(emom(:,ham%nlist(j,iflip),k)*emom(:,iflip,k)))
                  itp_noup(j)= (ham%ncoupD(j,iflip_h,ac(1))+excscale*(ham%ncoup(j,iflip_h,ac(1))-ham%ncoupD(j,iflip_h,ac(1))))+ &
                     invtemp*(  (ham%ncoupD(j,iflip_h,ac(2))+excscale*(ham%ncoup(j,iflip_h,ac(2))-ham%ncoupD(j,iflip_h,ac(2))))- &
                     (ham%ncoupD(j,iflip_h,ac(1))+excscale*(ham%ncoup(j,iflip_h,ac(1))-ham%ncoupD(j,iflip_h,ac(1)))))
                  itp_noupg(j)=invtemp2*( (ham%ncoupD(j,iflip_h,ac(2))+excscale*(ham%ncoup(j,iflip_h,ac(2))-ham%ncoupD(j,iflip_h,ac(2))))- &
                     ((ham%ncoupD(j,iflip_h,ac(1))+excscale*(ham%ncoup(j,iflip_h,ac(1))-ham%ncoupD(j,iflip_h,ac(1))))))
               enddo
               obj(:)=LSF_energy(ac(1),:)+invtemp*(LSF_energy(ac(2),:)-LSF_energy(ac(1),:))
               lsff=-invtemp2*(LSF_energy(ac(2),1)-LSF_energy(ac(1),1))
            endif

         elseif (Nchmax == 2) then
            do i=1,nchmax
               ind(i) = maxloc (mom_int(:,i),1,mom_int(:,i) .LT. ammom_inter(i))
            enddo
            iconf=(ind(2)-1)*ngrids(1)+ind(1)
            ii=ind(1) ; jj=ind(2)
            if(chemtype==1) then
               temp(1)=(mom_int(jj+1,2)-ammom_inter(2))/(mom_int(jj+1,2)-mom_int(jj,2))
               temp(2)=(ammom_inter(2)-mom_int(jj,2))/(mom_int(jj+1,2)-mom_int(jj,2))
               temp(3)=(mom_int(ii+1,1)-ammom_inter(1))/(mom_int(ii+1,1)-mom_int(ii,1))
               temp(4)=(ammom_inter(1)-mom_int(ii,1))/(mom_int(ii+1,1)-mom_int(ii,1))
               ac(1) = iconf
               ac(2) = iconf + ngrids(1)
               ac(3) = iconf +1
               ac(4) = iconf + ngrids(1) + 1
            else
               temp(1)=(mom_int(ii+1,1)-ammom_inter(1))/(mom_int(ii+1,1)-mom_int(ii,1))
               temp(2)=(ammom_inter(1)-mom_int(ii,1))/(mom_int(ii+1,1)-mom_int(ii,1))
               temp(3)=(mom_int(jj+1,2)-ammom_inter(2))/(mom_int(jj+1,2)-mom_int(jj,2))
               temp(4)=(ammom_inter(2)-mom_int(jj,2))/(mom_int(jj+1,2)-mom_int(jj,2))
               ac(1) = iconf
               ac(2) = iconf +1
               ac(3) = iconf + ngrids(1)
               ac(4) = iconf + ngrids(1) + 1
            endif
            invtemp =1.0_dblprec / (mom_int(ind(chemtype)+1,chemtype)-mom_int(ind(chemtype),chemtype))
            if (exc_inter=='N' ) then
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
               !dir$ loop count min(8)
               !$omp simd
#endif
               do j=1,ham%nlistsize(iflip_h)
                  itp_noup(j)= (ham%ncoup(j,iflip_h,ac(1))*temp(1)+ham%ncoup(j,iflip_h,ac(2))*temp(2))*temp(3)+ &
                     (ham%ncoup(j,iflip_h,ac(3))*temp(1)+ham%ncoup(j,iflip_h,ac(4))*temp(2))*temp(4)
                  itp_noupg(j)=invtemp*((ham%ncoup(j,iflip_h,ac(3))*temp(1)+ham%ncoup(j,iflip_h,ac(4))*temp(2))- &
                     (ham%ncoup(j,iflip_h,ac(1))*temp(1)+ham%ncoup(j,iflip_h,ac(2))*temp(2)))
               end do

               obj(:)=(LSF_energy(ac(1),:)*temp(1)+LSF_energy(ac(2),:)*temp(2))*temp(3)+&
                  (LSF_energy(ac(3),:)*temp(1)+LSF_energy(ac(4),:)*temp(2))*temp(4)
               lsff=-invtemp*((LSF_energy(ac(3),1)*temp(1)+LSF_energy(ac(4),1)*temp(2))-&
                  (LSF_energy(ac(1),1)*temp(1)+LSF_energy(ac(2),1)*temp(2)))
            elseif (exc_inter=='Y') then

#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
               !dir$ loop count min(8)
               !$omp simd private(excscale)
#endif
               do j=1,ham%nlistsize(iflip_h)
                  excscale=abs(sum(emom(:,ham%nlist(j,iflip),k)*emom(:,iflip,k)))
                  itp_noup(j)= ((ham%ncoupD(j,iflip_h,ac(1))+excscale*(ham%ncoup(j,iflip_h,ac(1))-ham%ncoupD(j,iflip_h,ac(1))))  *temp(1)+ &
                     (ham%ncoupD(j,iflip_h,ac(2))+excscale*(ham%ncoup(j,iflip_h,ac(2))-ham%ncoupD(j,iflip_h,ac(2))))  *temp(2))*temp(3)+ &
                     ((ham%ncoupD(j,iflip_h,ac(3))+excscale*(ham%ncoup(j,iflip_h,ac(3))-ham%ncoupD(j,iflip_h,ac(3))))  *temp(1)+ &
                     (ham%ncoupD(j,iflip_h,ac(4))+excscale*(ham%ncoup(j,iflip_h,ac(4))-ham%ncoupD(j,iflip_h,ac(4))))  *temp(2))*temp(4)

                  itp_noupg(j)=invtemp*( &
                     ((ham%ncoupD(j,iflip_h,ac(3))+excscale*(ham%ncoup(j,iflip_h,ac(3))-ham%ncoupD(j,iflip_h,ac(3))))  *temp(1)+ &
                     (ham%ncoupD(j,iflip_h,ac(4))+excscale*(ham%ncoup(j,iflip_h,ac(4))-ham%ncoupD(j,iflip_h,ac(4))))  *temp(2)) - &
                     ((ham%ncoupD(j,iflip_h,ac(1))+excscale*(ham%ncoup(j,iflip_h,ac(1))-ham%ncoupD(j,iflip_h,ac(1))))  *temp(1)+ &
                     (ham%ncoupD(j,iflip_h,ac(2))+excscale*(ham%ncoup(j,iflip_h,ac(2))-ham%ncoupD(j,iflip_h,ac(2))))  *temp(2)))
               enddo
               obj(:)=(LSF_energy(ac(1),:)*temp(1)+LSF_energy(ac(2),:)*temp(2))*temp(3)+&
                  (LSF_energy(ac(3),:)*temp(1)+LSF_energy(ac(4),:)*temp(2))*temp(4)
               lsff=-invtemp*((LSF_energy(ac(3),1)*temp(1)+LSF_energy(ac(4),1)*temp(2))-&
                  (LSF_energy(ac(1),1)*temp(1)+LSF_energy(ac(2),1)*temp(2)))
            endif
         else
            write(*,*) "three component or higher dimension interpolation not being implemented"
         end if
      else                            ! no interpolation, pick nearest lower grid point
         do i=1,nchmax
            ind(i) = maxloc (mom_int(:,i),1,mom_int(:,i) .LT. ammom_inter(i))
         enddo
         if (nchmax==1) then
            iconf = ind(1)
            if(inttype==0) then
               invtemp=(mom_int(iconf,1)/ammom_inter(1))**2
               ac(1) = iconf
               ac(2) = iconf +1
            else
               iconf = iconf+1
               invtemp=(ammom_inter(1)/mom_int(iconf,1))**2
               ac(1) = iconf -1
               ac(2) = iconf
            endif
         else
            iconf=(ind(2)-1)*ngrids(1)+ind(1)
            if(inttype==0) then
               invtemp=(ammom_int(1,iconf)*ammom_int(2,iconf))/(ammom_inter(1)*ammom_inter(2))
               ac(1) = iconf
               ac(2) = iconf +1
               ac(3) = iconf + ngrids(1)
               ac(4) = iconf + ngrids(1) + 1
            else
               ind(:)=ind(:)+1
               iconf=iconf*ngrids(1)+1
               invtemp=(ammom_inter(1)*ammom_inter(2))/(ammom_int(1,iconf)*ammom_int(2,iconf))
               ac(1) = iconf - ngrids(1) -1
               ac(2) = iconf - ngrids(1)
               ac(3) = iconf - 1
               ac(4) = iconf
            endif
         endif
         indp(1)=ac(1) ; indp(2)=ac(chemtype+1)
         invtemp2 =1.0_dblprec / (mom_int(ind(chemtype)+1,chemtype)-mom_int(ind(chemtype),chemtype))
         obj(:) = LSF_energy(iconf,:)*invtemp
         ! approx gradient
         lsff=-invtemp2*(LSF_energy(indp(2),1)-LSF_energy(indp(1),1))

         if (exc_inter=='N') then
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
            !dir$ loop count min(8)
            !$omp simd
#endif
            do j=1,ham%nlistsize(iflip_h)
               itp_noup(j)=ham%ncoup(j,iflip_h,iconf)*invtemp
               itp_noupg(j)=invtemp2*(ham%ncoup(j,iflip_h,indp(2))-ham%ncoup(j,iflip_h,indp(1)))
            enddo
         elseif (exc_inter=='Y') then
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
            !dir$ loop count min(8)
            !$omp simd private(excscale)
#endif
            do j=1,ham%nlistsize(iflip_h)
               excscale=abs(sum(emom(:,ham%nlist(j,iflip),k)*emom(:,iflip,k)))
               itp_noup(j)=(ham%ncoupD(j,iflip_h,iconf)+excscale*(ham%ncoup(j,iflip_h,iconf)-ham%ncoupD(j,iflip_h,iconf)))*invtemp
               ! approx gradient on the moment square boundary for now
               itp_noupg(j)=invtemp2*( (ham%ncoupD(j,iflip_h,indp(2))+excscale*(ham%ncoup(j,iflip_h,indp(2))- &
                  ham%ncoupD(j,iflip_h,indp(2))))-(ham%ncoupD(j,iflip_h,indp(1))+excscale*(ham%ncoup(j,iflip_h,indp(1))- &
                  ham%ncoupD(j,iflip_h,indp(1)))))
            enddo
         endif
      endif
   end subroutine do_interpolation_ncoup_and_lsf_gradient


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

      open(ifileno, file=adjustl(lsffile))

      ! Reads the configurations dependent LSF energy
      do m =1, conf_num
         read (ifileno,*) iconf,LSF_energy(iconf,1),LSF_energy(iconf,2)
      enddo

      if (mconf==0) then
         mconftmp=minloc(LSF_energy,1)
         mconf=mconftmp(1)
      endif

      close (ifileno)

   end subroutine read_LSF


end module LSF
