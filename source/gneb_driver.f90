!-------------------------------------------------------------------------------
! MODULE: gneb_driver
!> @brief Driver for the calculation of the GNEB modes
!> @details This contains the capacity of use the VPO minimization scheme for the
!> inital phase as well as the GNEB alogrithms to find the Minimum energy path (MEP)
!> that allows one to calculate the energy barrier between the states, as well as the
!> lifetimes and prefactors via the Harmonic Transition State Theory (HTST).
!> @author Pavel Bessarab
!> @todo Remove the local im_ene() arrays to use the ene% used in the energy routine
!> @todo Create a local reading function to modularize the code even more
!> @todo Incorporate the measurement rotuines from the main code into GNEB
!> @updated by Nikos Ntallis
!> @im_ene array has been removed. After the minimazation the ene structure holds the energy per image. 
!> @if plotenergy flag is used energies/atom are printed as in the main code.
!-------------------------------------------------------------------------------
module gneb_driver
   use Energy,   only: ene,allocate_energies,calc_energy
   implicit none

   public

contains

   !---------------------------------------------------------------------------
   !> @brief
   !> Energy minimization initial phase
   !
   !> @author
   !> Pavel Bessarab
   !---------------------------------------------------------------------------
   subroutine em_iphase()
      !
      use InputData
      use FieldData,            only : sitefld, external_field
      use FieldPulse
      use SystemData
      use DemagField
      use MomentData
      use Temperature
      use Measurements
      use UpdateMoments
      use MicroWaveField
      use SimulationData,        only : bn
      use HamiltonianData
      use CalculateFields
      use optimizationRoutines
      use HamiltonianActions
      use EnergyMinima
      use macrocells
     
      !
      implicit none

      real(dblprec) :: energy

      if(plotenergy.lt.1) then
         call allocate_energies(1, Mensemble)
      endif
      write (*,'(1x,a)') "Performing initial phase: energy minimization"
      write (*,'(1x,a,i3)') "Minimization algorithm: ", minalgo

      ! Calculate demagnetization field
      if(demag=='Y') then
         call calc_demag(Natom, Mensemble, demag1, demag2, demag3, demagvol, emomM)
      endif

      ! Calculate external field
      call calc_external_fields(Natom,Mensemble,iphfield,anumb,external_field,   &
         do_bpulse,sitefld,sitenatomfld)

      if (mode=='G') call save_ifmom(Natom,Mensemble,simid,0,mmom,emom,mode,do_mom_legacy)

      ! --Optimization Region-- !
      ! Allocation and initial opt build
      if(OPT_ON) then
         call timing(0,'Hamiltonian   ','OF')
         call timing(0,'BuildOptReg   ','ON')
         call buildOptimizationRegions(na,natom,nHam,Mensemble,ham%max_no_neigh,    &
            atype,emom,mmom,ham%ncoup,ham%nlist,ham%nlistsize,ham%aham,             &
            cellPosNumNeigh,cos_thr)
         call timing(0,'BuildOptReg   ','OF')
         call timing(0,'Hamiltonian   ','ON')
      end if

      write (*,'(1x,a)') "Energy minimization in progress"
      if (minalgo==1) then

         !!! call vpo_min(nHam,mintraj_step,Natom,do_dm,do_pd,do_bq,do_ring,do_chir,do_dip, &
         !!!    Nchmax,minitrmax,OPT_flag,do_biqdm,conf_num,2,Num_macro,do_jtensor,         &
         !!!    plotenergy,do_anisotropy,max_no_constellations,minftol,vpomass,vpodt,       &
         !!!    simid,do_lsf,lsf_field,exc_inter,mult_axis,lsf_interpolate,cell_index,      &
         !!!    macro_nlistsize,mmom(:,1:Mensemble:(Mensemble-1)),                          &
         !!!    emom(:,:,1:Mensemble:(Mensemble-1)),emomM_macro,                            &
         !!!    external_field(:,:,1:Mensemble:(Mensemble-1)),maxNoConstl,unitCellType,     &
         !!!    constlNCoup,constellations,constellationsNeighType,energy,                  &
         !!!    emomM(:,:,1:Mensemble:(Mensemble-1)),NA,N1,N2,N3,mode,do_mom_legacy)
         if (mode=='G') then
            call vpo_min(nHam,mintraj_step,Natom,Nchmax,minitrmax,OPT_flag,conf_num,2,     &
               Num_macro,plotenergy,max_no_constellations,minftol,vpomass,vpodt,           &
               simid,do_lsf,lsf_field,lsf_interpolate,cell_index,                          &
               macro_nlistsize,mmom(:,1:Mensemble:(Mensemble-1)),                          &
               emom(:,:,1:Mensemble:(Mensemble-1)),emomM_macro,                            &
               external_field(:,:,1:Mensemble:(Mensemble-1)),maxNoConstl,unitCellType,     &
               constlNCoup,constellations,constellationsNeighType,energy,                  &
               emomM(:,:,1:Mensemble:(Mensemble-1)),NA,N1,N2,N3,mode,do_mom_legacy)
         else
            call vpo_min_single(nHam,mintraj_step,Natom,Nchmax,minitrmax,OPT_flag,conf_num,Mensemble, &
               Num_macro,plotenergy,max_no_constellations,minftol,vpomass,vpodt,           &
               simid,do_lsf,lsf_field,lsf_interpolate,cell_index,                          &
               macro_nlistsize,mmom,                                                       &
               emom,emomM_macro,                                                           &
               external_field,maxNoConstl,unitCellType,                                    &
               constlNCoup,constellations,constellationsNeighType,energy,                  &
               emomM,NA,N1,N2,N3,mode,do_mom_legacy)
         end if

      end if
      write (*,'(1x,a)') "Done"
      write (*,'(1x,a)') "------------------------------------------"

      if (mode=='G') call save_ifmom(Natom,Mensemble,simid,1,mmom,emom,mode,do_mom_legacy)

   end subroutine em_iphase

   !---------------------------------------------------------------------------
   !> @brief
   !> Calculation of minimum energy path (MEP) using geodesic nudged elastic band (GNEB) method
   !
   !> @author
   !> Pavel Bessarab
   !---------------------------------------------------------------------------
   subroutine mep_mphase()
      !
      use InputData
      use FieldData,            only : beff, beff1, beff2, sitefld, external_field
      use FieldPulse
      use SystemData
      use DemagField
      use MomentData
      use Temperature
      use Measurements
      use UpdateMoments
      use MicroWaveField
      use SimulationData,        only : bn
      use HamiltonianData
      use CalculateFields
      use optimizationRoutines
      use HamiltonianActions
      use ErrorHandling
      use PathInit
      use GNEB
      use HermiteLib,           only : hermite_fit,spline_hermite_val
      use Hessian
      use macrocells
      use Restart, only: GNEB_read_wrapper
      !
      implicit none

      integer :: ci !< Index of the climbing image
      integer :: ii, jj,i1, rstep
      integer :: i_stat,i_all
      real(dblprec) :: denergy,fcinv
      character(len=35) :: num
      character(len=35) :: beff_filn
      logical :: exists, selected
      real(dblprec), dimension(3,Natom,Mensemble) :: tef !< External time-dependent magnetic field, zero here
      ! .. Local allocatable variables
      real(dblprec), dimension(:), allocatable :: rx
      real(dblprec), dimension(:), allocatable :: xx
      real(dblprec), dimension(:), allocatable :: yy
      real(dblprec), dimension(:), allocatable :: dyy
      real(dblprec), dimension(:), allocatable :: rx0
      real(dblprec), dimension(:), allocatable :: ene0
      real(dblprec), dimension(:), allocatable :: dene
      real(dblprec), dimension(:), allocatable :: dene0
     ! real(dblprec), dimension(:), allocatable :: im_ene !! Energy of the images
      real(dblprec), dimension(:,:), allocatable :: coef !! Coefficients of the piecewise Hermite polynomials
      real(dblprec), dimension(:,:,:), allocatable :: emomsp !! Magnetic configuration of the saddle point
      real(dblprec), dimension(:,:,:), allocatable :: emom_tmp !! Temporal storage of the magnetic moments

      ! Allocate the needed arrays for the GNEB measurement phase
      allocate(coef(4,Mensemble),stat=i_stat)
      call memocc(i_stat,product(shape(coef))*kind(coef),'coef','mep_mphase')
      coef=0.0_dblprec
      allocate(rx(Mensemble),stat=i_stat)
      call memocc(i_stat,product(shape(rx))*kind(rx),'rx','mep_mphase')
      rx=0.0_dblprec
      allocate(xx(sample_num),stat=i_stat)
      call memocc(i_stat,product(shape(xx))*kind(xx),'xx','mep_mphase')
      xx=0.0_dblprec
      allocate(yy(sample_num),stat=i_stat)
      call memocc(i_stat,product(shape(yy))*kind(yy),'yy','mep_mphase')
      yy=0.0_dblprec
      allocate(rx0(1),stat=i_stat)
      call memocc(i_stat,product(shape(rx0))*kind(rx0),'rx0','mep_mphase')
      rx0=0.0_dblprec
      allocate(dyy(sample_num),stat=i_stat)
      call memocc(i_stat,product(shape(dyy))*kind(dyy),'dyy','mep_mphase')
      dyy=0.0_dblprec
      allocate(ene0(1),stat=i_stat)
      call memocc(i_stat,product(shape(ene0))*kind(ene0),'ene0','mep_mphase')
      ene0=0.0_dblprec
      allocate(dene(Mensemble),stat=i_stat)
      call memocc(i_stat,product(shape(dene))*kind(dene),'dene','mep_mphase')
      dene=0.0_dblprec
      allocate(dene0(1),stat=i_stat)
      call memocc(i_stat,product(shape(dene0))*kind(dene0),'dene0','mep_mphase')
      dene0=0.0_dblprec
     ! allocate(im_ene(Mensemble),stat=i_stat)
     ! call memocc(i_stat,product(shape(im_ene))*kind(im_ene),'im_ene','mep_mphase')
     ! im_ene=0.0_dblprec
      allocate(emom_tmp(3,Natom,Mensemble),stat=i_stat)
      call memocc(i_stat,product(shape(emom_tmp))*kind(emom_tmp),'emom_tmp','mep_mphase')
      emom_tmp=0.0_dblprec

      fcinv=mub/mry
      selected = .false.

      if ((do_gneb=='N').and.(do_gneb_ci=='N')) then
         call ErrorHandling_ERROR("Both 'do_gneb' and 'do_gneb_ci' can not be equal to 'N' at the same time!")
      end if

      ! Calculate demagnetization field
      if(demag=='Y') then
         call calc_demag(Natom,Mensemble,demag1,demag2,demag3,demagvol,emomM)
      endif

      ! Calculate external field
      call calc_external_fields(Natom,Mensemble,iphfield,anumb,external_field,   &
         do_bpulse,sitefld,sitenatomfld)

      ! Generate the path
      if (initpath==2) then
         write (*,'(1x,a,a)') "Initial guess for the path from file ",trim(adjustl(restartfile))
         call GNEB_read_wrapper(Natom,Mensemble,amp_rnd,'path',relaxed_if,          &
            do_mom_legacy,rstep,exists,restartfile,mmom,emom,emomM)
         if (.not.exists) then
            write (*,'(1x,a,a)') "File ",trim(adjustl(restartfile))," does not exist"
            call geodesic_path(Natom,Mensemble,emom(:,:,1),emom(:,:,Mensemble),     &
               amp_rnd_path,emom_tmp)
            write (*,'(1x,a,a)') "Geodesic path generated"
            !$omp parallel do default(shared), private(ii,jj)
            do ii=2,(Mensemble-1)
               do jj = 1,Natom
                  emom(:,jj,ii) = emom_tmp(:,jj,ii)
                  emomM(:,jj,ii) =emom(:,jj,ii)*mmom(jj,1)
               end do
            end do
            !$omp end parallel do
         end if
      else
         call geodesic_path(Natom,Mensemble,emom(:,:,1),emom(:,:,Mensemble),        &
            amp_rnd_path,emom_tmp)
         write (*,'(1x,a,a)') "Geodesic path generated"
         !$omp parallel do default(shared), private(ii,jj)
         do ii=2,(Mensemble-1)
            do jj = 1,Natom
               emom(:,jj,ii) = emom_tmp(:,jj,ii)
               emomM(:,jj,ii) =emom(:,jj,ii)*mmom(jj,1)
            end do
         end do
         !$omp end parallel do
      end if

      ! Deallocate the temporal magnetic moment array
      i_all=-product(shape(emom_tmp))*kind(emom_tmp)
      deallocate(emom_tmp,stat=i_stat)
      call memocc(i_stat,i_all,'emom_tmp','mep_mphase')
      
      call save_path(Natom,Mensemble,simid,0,emom,mmom,mode,do_mom_legacy)

      call print_GNEB_siminfo()

      ! Perform GNEB calculation
      if (do_gneb=='Y') then
         write (*,'(1x,a)') "       GNEB calculation in progress       "
         if (minalgo==1) then
            call gneb_mep(Mensemble,nHam,Natom,meptraj_step,ham_inp%do_dm,ham_inp%do_pd,ham_inp%do_bq,      &
               ham_inp%do_ring,ham_inp%do_chir,ham_inp%do_dip,mepitrmax,Nchmax,conf_num,ham_inp%do_biqdm,Num_macro, &
               ham_inp%do_jtensor,plotenergy,ham_inp%do_anisotropy,max_no_constellations,vpomass,   &
               mepftol,spring,vpodt,simid,do_lsf,en_zero,fixed_if,ham_inp%mult_axis,        &
               ham_inp%exc_inter,lsf_field,lsf_interpolate,OPT_flag,cell_index,             &
               macro_nlistsize,mmom,emom,emomM_macro,external_field,maxNoConstl,    &
               unitCellType,constellationsNeighType,constlNCoup,constellations,     &
               denergy,emomM,rx,dene,NA,N1,N2,N3,mode,do_mom_legacy)
         end if
         write (*,'(1x,a)') " Done"
         write (*,'(1x,a)') "------------------------------------------"
      end if

      ! Perform climbing image GNEB calculation
      if (do_gneb_ci=='Y') then
         write (*,'(1x,a)') "     CI-GNEB calculation in progress      "
         if (minalgo==1) then
            call gneb_ci_mep(Mensemble,nHam,Natom,meptraj_step,ham_inp%do_dm,ham_inp%do_pd,ham_inp%do_bq,   &
               ham_inp%do_ring,ham_inp%do_chir,ham_inp%do_dip,mepitrmax,Nchmax,conf_num,ham_inp%do_biqdm,Num_macro, &
               ham_inp%do_jtensor,plotenergy,ham_inp%do_anisotropy,max_no_constellations,vpomass,   &
               mepftol_ci,spring,vpodt,simid,do_lsf,en_zero,fixed_if,ham_inp%mult_axis,     &
               ham_inp%exc_inter,lsf_field,lsf_interpolate,OPT_flag,cell_index,             &
               macro_nlistsize,mmom,emom,emomM_macro,external_field,maxNoConstl,    &
               unitCellType,constellationsNeighType,constlNCoup,constellations,     &
               denergy,emomM,ci,rx,dene,NA,N1,N2,N3,mode,do_mom_legacy)

            write(*,'(2x,a,2x,i8)') 'Climbing Image:',ci

         end if
         write (*,'(1x,a)') "Done"
         write (*,'(1x,a)') "------------------------------------------"
      end if

      ! Interpolate the energy along the Minimum Energy Path
      call hermite_fit(Mensemble,sample_num,rx,ene%energy(:),dene,xx,yy,dyy,coef)

      do ii=1,Mensemble
         do jj=1,Natom
            emom(1:3,jj,ii) = emomM(1:3,jj,ii)/mmom(jj,ii)
         end do
      end do

      call save_path(Natom,Mensemble,simid,1,emom,mmom,mode,do_mom_legacy)

      ! Save the energy along the Minimum Energy Path
      call save_en(Mensemble,rx,ene%energy(:),dene,rx(Mensemble),                          &
         'en_path.'//trim(adjustl(simid))//'.out',do_norm_rx)

      ! Save the interpolated energy along the Minimum Energy Path
      call save_en(sample_num,xx,yy,dyy,rx(Mensemble),                              &
         'enfit_path.'//trim(adjustl(simid))//'.out',do_norm_rx)

      allocate(emomsp(3,Natom,1),stat=i_stat)
      call memocc(i_stat,product(shape(emomsp))*kind(emomsp),'emomsp','mep_mphase')
      emomsp=0.0_dblprec

      ! Identify the Saddle Point
      if (do_gneb_ci=='Y') then ! CI gives the Saddle Point directly
         rx0(1) = rx(ci)
         do ii=1,Natom
            emomsp(1,ii,1) = emom(1,ii,ci)
            emomsp(2,ii,1) = emom(2,ii,ci)
            emomsp(3,ii,1) = emom(3,ii,ci)
         end do
      else
         ! Otherwise, find a maximum along the Minimum Energy Path
         call find_SP(Mensemble,rx,coef,rx0(1))
         i1=1
         ii=1
         do while(ii<=Mensemble-1.and. .not. selected)
            if ((rx(ii)<=rx0(1)).and.(rx(ii+1)>=rx0(1))) then
               i1 = ii
               selected=.true.
            end if
            ii=ii+1
         end do

         ! Find the Saddle Point configuration
         call find_SP_conf(Natom,emom(:,:,i1),emom(:,:,i1+1),rx(i1),rx(i1+1),rx0(1),&
            emomsp(:,:,1))
      end if

      ! Calculate the Saddle Point energy
      call spline_hermite_val(Mensemble,rx,coef,rx0(1),ene0(1),dene0(1))

      ! Save the Saddle Point energy
      call save_en(1,rx0(1),ene0(1),dene0(1),rx(Mensemble), 'en_sp.'//trim(adjustl(simid))//'.out',do_norm_rx)

      ! Save the Saddle Point configuration
      call save_path(Natom,1,simid,4,emomsp(:,:,:),mmom(:,:),mode,do_mom_legacy)

      call timing(0,'Hamiltonian   ','ON')
      ! Calculate effective fields at each image
      call effective_field(Natom,Mensemble,1,Natom,emomM,mmom,    &
         external_field,tef,beff,beff1,beff2,OPT_flag,max_no_constellations,        &
         maxNoConstl,unitCellType,constlNCoup,constellations,                       &
         constellationsNeighType,denergy,Num_macro,cell_index,emomM_macro,&
         macro_nlistsize,NA,N1,N2,N3)
      call timing(0,'Hamiltonian   ','OF')

      if (prn_gneb_fields=='Y') then
         ! Save the effective fields in separate files
         do ii=1,Mensemble
            write(num,'(i8)') ii
            beff_filn='beff_'//trim(adjustl(num))//'.'//trim(adjustl(simid))//'.out'
            open(ofileno,file=beff_filn,access='sequential',action='write',status='replace')
            write(ofileno,'(a8,2x,a16,2x,a16,2x,a16)') "# Site","B_x","B_y","B_z"
            do jj=1,Natom
               write(ofileno,'(i8,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3)',advance='yes') &
                  jj,beff(1,jj,ii)*fcinv,beff(2,jj,ii)*fcinv,beff(3,jj,ii)*fcinv
            end do
            close(ofileno)
         end do
      endif

      !-------------------------------------------------------------------------
      ! Calculation of the Hessian
      !-------------------------------------------------------------------------
      if ((do_hess_ini=='Y').or.(do_hess_fin=='Y').or.(do_hess_sp=='Y')) then
         call hessian_wrapper(N1,N2,N3,Mensemble,nHam,Natom,ham_inp%do_dm,ham_inp%do_pd,ham_inp%do_bq,      &
            ham_inp%do_ring,ham_inp%do_chir,ham_inp%do_dip,ham_inp%do_biqdm,Num_macro,ham_inp%do_jtensor,plotenergy,        &
            ham%max_no_neigh,ham_inp%do_anisotropy,ham%max_no_dmneigh,max_no_constellations,&
            eig_0,BC1,BC2,BC3,simid,do_lsf,is_afm,ham_inp%mult_axis,ham_inp%exc_inter,lsf_field,    &
            do_hess_sp,do_hess_ini,do_hess_fin,lsf_interpolate,OPT_flag,ham%aHam,   &
            ham%taniso,ham%nlistsize,ham%dmlistsize,cell_index,macro_nlistsize,     &
            ham%nlist,ham%dmlist,C1,C2,C3,ene0,ene%energy(:),coord,ham%eaniso,ham%kaniso,  &
            emomsp,ham%ncoup,mmom,emomM_macro,ham%dm_vect,external_field,           &
            maxNoConstl,unitCellType,constellationsNeighType,constlNCoup,           &
            constellations,emomM,beff,beff1,beff2,NA,emom)
      endif

      if(plotenergy.lt.1) then
      call allocate_energies(-1, Mensemble)
      endif
      ! Deallocation and cleanup of the local arrays used in the GNEB measurement phase
      i_all=-product(shape(coef))*kind(coef)
      deallocate(coef,stat=i_stat)
      call memocc(i_stat,i_all,'coef','mep_mphase')
      i_all=-product(shape(rx))*kind(rx)
      deallocate(rx,stat=i_stat)
      call memocc(i_stat,i_all,'rx','mep_mphase')
      i_all=-product(shape(xx))*kind(xx)
      deallocate(xx,stat=i_stat)
      call memocc(i_stat,i_all,'xx','mep_mphase')
      i_all=-product(shape(yy))*kind(yy)
      deallocate(yy,stat=i_stat)
      call memocc(i_stat,i_all,'yy','mep_mphase')
      i_all=-product(shape(dyy))*kind(dyy)
      deallocate(dyy,stat=i_stat)
      call memocc(i_stat,i_all,'dyy','mep_mphase')
      i_all=-product(shape(rx0))*kind(rx0)
      deallocate(rx0,stat=i_stat)
      call memocc(i_stat,i_all,'rx0','mep_mphase')
      i_all=-product(shape(ene0))*kind(ene0)
      deallocate(ene0,stat=i_stat)
      call memocc(i_stat,i_all,'ene0','mep_mphase')
      i_all=-product(shape(dene))*kind(dene)
      deallocate(dene,stat=i_stat)
      call memocc(i_stat,i_all,'dene','mep_mphase')
      i_all=-product(shape(dene0))*kind(dene0)
      deallocate(dene0,stat=i_stat)
      call memocc(i_stat,i_all,'dene0','mep_mphase')
     ! i_all=-product(shape(im_ene))*kind(im_ene)
     ! deallocate(im_ene,stat=i_stat)
      !call memocc(i_stat,i_all,'im_ene','mep_mphase')
      i_all=-product(shape(emomsp))*kind(emomsp)
      deallocate(emomsp,stat=i_stat)
      call memocc(i_stat,i_all,'emomsp','mep_mphase')

   end subroutine mep_mphase

   !----------------------------------------------------------------------------
   !> @brief
   !> Print brief information about the simulation to stdout for GNEB applications
   !
   !> @author
   !> Jonathan Chico based in the routine by Anders Bergman
   !----------------------------------------------------------------------------
   subroutine print_GNEB_siminfo()

      use SystemData
      use InputData

      implicit none

      write (*,'(1x,a)')          "------------------------------------------"
      write (*,'(1x,a)')          "           Simulation data                "
      write (*,'(1x,a)')          "------------------------------------------"
      write (*,'(1x,a,x,i8)')     " Number of atoms  :", Natom
      write (*,'(1x,a,x,i8)')     " Number of images :", Mensemble
      write (*,'(1x,a,x,a)')      " GNEB    :",do_gneb
      write (*,'(1x,a,x,a)')      " GNEB_CI :",do_gneb_ci
      write (*,'(1x,a,x,i3)')     " Minimization algorithm :", minalgo
      write (*,'(1x,a,x,es16.8)') " Spring constant        :",spring
      if (minalgo.eq.1) then
         write (*,'(1x,a,x,G11.4)')  " VPO step :", vpodt
         write (*,'(1x,a,x,G11.4)')  " VPO mass :", vpomass
         write (*,'(1x,a,x,G11.4)')  " Force tolerance    :", mepftol
         write (*,'(1x,a,x,G11.4)')  " CI Force tolerance :", mepftol_ci
         write (*,'(1x,a,x,i10)')    " Maximum number of iterations :", mepitrmax
         write(*,'(1x,a,x,i10,x,a)') " Path saved after every       :",meptraj_step, "steps"
      endif
      if(ipmode.ne.'G') then
         write(*,'(1x,a)') 'WARNING: initial and/or final states have not been relaxed!'
      endif
      write (*,'(1x,a)')          "------------------------------------------"
      write (*,'(1x,a)')          "     Energy interpolation parameters      "
      write (*,'(1x,a)')          "------------------------------------------"
      write (*,'(1x,a,x,a)')      " Normalization of reaction coordinates  : ", do_norm_rx
      write (*,'(1x,a,x,i10)')    " Num. points in interpolation of energy : ", sample_num
      write (*,'(1x,a)')          "------------------------------------------"
 

   end subroutine print_GNEB_siminfo

end module gneb_driver
