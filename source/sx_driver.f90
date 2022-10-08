!******************************************************************
!*                     *** UppASD ***                             *
!*                                                                *
!*               Uppsala Atomic Spin Dynamics                     *
!*                                                                *
!*                   Version 5.0 Mar 2017                         *
!*                                                                *
!*       Anders Bergman                                           *
!*       Johan Hellsvik             Lars Bergqvist                *
!*       Jonathan Chico             Thomas Nystrand               *
!*       Bjorn Skubic               Olle Eriksson                 *
!*       Andrea Taroni              Lars Nordstrom                *
!*       + other contributors.                                    *
!*                                                                *
!*       Division of Materials Theory                             *
!*       Department of Physics and Materials Science              *
!*       Uppsala University                                       *
!*       Sweden                                                   *
!*                                                                *
!*       and                                                      *
!*                                                                *
!*       Department of Materials and Nano Science                 *
!*       KTH Royal Institute of Technology                        *
!*       Sweden                                                   *
!*                                                                *
!******************************************************************
module sx_driver
   use Constants
   use Parameters

   implicit none

   integer :: sx_step = 50
   integer :: sx_numrep

   real(dblprec),dimension(:,:,:,:), allocatable :: sx_emom   !< Unit moments for PT replicas
   real(dblprec),dimension(:,:,:), allocatable :: sx_emom_1q   !< Unit moments for PT replicas
   real(dblprec),dimension(:,:,:), allocatable :: sx_emom_3q   !< Unit moments for PT replicas
   real(dblprec),dimension(:,:,:), allocatable :: sx_emom_rn   !< Unit moments for PT replicas
   real(dblprec),dimension(:,:,:), allocatable :: sx_emom_fm   !< Unit moments for PT replicas
   real(dblprec),dimension(:,:,:,:), allocatable :: sx_emomM  !< Moments for PT replicas
   real(dblprec),dimension(:,:,:), allocatable :: sx_mmom   !< Moment magnitudes for PT replicas

   real(dblprec),dimension(:,:,:,:), allocatable :: sx_emom_macro   !< Unit macrospin moments for PT replicas
   real(dblprec),dimension(:,:,:,:), allocatable :: sx_emomM_macro  !< Macrospin moments for PT replicas
   real(dblprec),dimension(:,:,:), allocatable :: sx_mmom_macro   !< Macrospin moment magnitudes for PT replicas

   character(len=2) :: sx_mode = 'H'

   public

contains

   !---------------------------------------------------------------------------
   !> @brief
   !> Replica Exchange / Parallel Tempering Monte Carlo initial phase.
   !> This phase of the program is used for equillibration.
   !> No measurements are performed.
   !> Currently only Heat Bath MC is supported
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine sx_iphase()

      use Energy,          only : calc_energy
      use InputData
      use FieldData,       only : external_field, sitefld, time_external_field, beff, beff1, beff2
      use SystemData
      use MonteCarlo
      use DemagField
      use CalculateFields, only : calc_external_fields
      use MomentData
      use HamiltonianData
      use HamiltonianActions, only : effective_field
      use Measurements, only : calc_mavrg
      use AMS
      use QHB
      use Restart
      use macrocells
      use optimizationRoutines
      use RandomNumbers, only: rng_uniform
      use sd_driver, only : sd_minimal
      use mc_driver, only : mc_minimal

      integer :: irep, ipmcstep
      integer :: i_all, i_stat, ia, ik
      real(dblprec) :: temprescale
      character(len=30) :: filn
      integer :: mc_nmeas,mc_nsamp
      real(dblprec) :: mavg
      real(dblprec) :: n_hit, n_miss, n_trial
      real(dblprec), dimension(:), allocatable :: repene  !< Array of replica energies
      real(dblprec), dimension(:), allocatable :: repene_sum  !< Array of replica energies averages
      integer, dimension(:), allocatable :: replist   !< List to denominate which temperature each replica has
      integer, dimension(:), allocatable :: templist  !< List to denominate which replica has which temperature
      ! I.e. if replica 3 has the highest temperature, then replist(3)=1, and templist(1)=3
      ! At start, replist(i)=templist(i)=i.
      real(dblprec), dimension(:), allocatable :: sx_flipprob !< Array for RNG flipping probabilities.

      sx_numrep=4 
      ! Write header for moment file
      write (filn,'(''sxinitial.'',a,''.out'')') trim(simid)
      open(11, file=filn, position="append")
      write(11,'(a)') "#  Iter.   M_avg.    U_Bind.    Susc."

      ! Allocate work arrays for Metropolis algorithm
      call allocate_mcdata(Natom,1)

      n_hit=0.0_dblprec;n_miss=0.0_dblprec;n_trial=0.0_dblprec

      ! Allocate local moment arrays
      allocate(sx_flipprob(sx_numrep),stat=i_stat)
      call memocc(i_stat,product(shape(sx_flipprob))*kind(sx_flipprob),'sx_flipprob','sx_iphase')
      allocate(templist(sx_numrep),stat=i_stat)
      call memocc(i_stat,product(shape(templist))*kind(templist),'templist','sx_iphase')
      allocate(replist(sx_numrep),stat=i_stat)
      call memocc(i_stat,product(shape(replist))*kind(replist),'replist','sx_iphase')
      allocate(repene(sx_numrep),stat=i_stat)
      call memocc(i_stat,product(shape(repene))*kind(repene),'repene','sx_iphase')
      allocate(repene_sum(sx_numrep),stat=i_stat)
      call memocc(i_stat,product(shape(repene_sum))*kind(repene_sum),'repene_sum','sx_iphase')
      repene_sum=0.0_dblprec
      allocate(sx_emom(3,Natom,Mensemble,sx_numrep),stat=i_stat)
      call memocc(i_stat,product(shape(sx_emom))*kind(sx_emom),'sx_emom','sx_iphase')
      allocate(sx_emomM(3,Natom,Mensemble,sx_numrep),stat=i_stat)
      call memocc(i_stat,product(shape(sx_emomM))*kind(sx_emomM),'sx_emomM','sx_iphase')
      allocate(sx_mmom(Natom,Mensemble,sx_numrep),stat=i_stat)
      call memocc(i_stat,product(shape(sx_mmom))*kind(sx_mmom),'sx_mmom','sx_iphase')

      if(Num_macro>0) then
         allocate(sx_emom_macro(3,Num_macro,Mensemble,sx_numrep),stat=i_stat)
         call memocc(i_stat,product(shape(sx_emom_macro))*kind(sx_emom_macro),'sx_emom','sx_iphase')
         allocate(sx_emomM_macro(3,Num_macro,Mensemble,sx_numrep),stat=i_stat)
         call memocc(i_stat,product(shape(sx_emomM_macro))*kind(sx_emomM_macro),'sx_emomM_macro','sx_iphase')
         allocate(sx_mmom_macro(Num_macro,Mensemble,sx_numrep),stat=i_stat)
         call memocc(i_stat,product(shape(sx_mmom_macro))*kind(sx_mmom_macro),'sx_mmom_macro','sx_iphase')
      end if


      ! Write output to stdout
      write (*,'(a28,i3,a10,G11.4,a10,i10,a10,a10)') &
         "Performing SXPT,  Replica:", irep ," Temp: ", ipTemp(1)

      sx_emom(:,:,:,1)=sx_emom_1q
      sx_emom(:,:,:,2)=sx_emom_3q
      sx_emom(:,:,:,3)=sx_emom_rn
      sx_emom(:,:,:,4)=sx_emom_fm

      if(Num_macro>0) then
         sx_mmom_macro(:,:,irep)=mmom_macro
         sx_emom_macro(:,:,:,irep)=emom_macro
         sx_emomM_macro(:,:,:,irep)=emomM_macro
      end if

      do irep=1,sx_numrep

         sx_mmom(:,:,irep)=mmom
         replist(irep)=irep
         templist(irep)=irep

         do ik=1, Mensemble
            do ia=1,Natom
               sx_emomM(:,ia,ik,irep)=sx_emom(:,ia,ik,irep)*sx_mmom(ia,ik,irep)
            end do
         end do

      end do


      ! Calculate the static magnetic fields which will be calculated only once as they are not time dependent
      call calc_external_fields(Natom, Mensemble, hfield, anumb, external_field, do_bpulse,sitefld,sitenatomfld)

      ! Set up order for Metropolis sweeps
      call choose_random_atom_x(Natom,iflip_a)

      ! Zero data for averaging
      mc_nmeas=0
      mc_nsamp=0
      !
      !do irep=1,sx_numrep
      !   call mc_minimal(sx_emomM(:,:,:,irep),sx_emom(:,:,:,irep),sx_mmom(:,:,irep),sx_step,'M',ipTemp(irep))
      !   call sd_minimal(sx_emomM(:,:,:,irep),sx_emom(:,:,:,irep),sx_mmom(:,:,irep),sx_step,1,ipTemp(irep))
      !end do
      ! Loop over steps of sweeps
      do ipmcstep=1,ipmcnstep(1),sx_step
         do irep=1,sx_numrep
            !print *,'-----------------------',ipTemp(irep)
            !print '(3f10.4)',emom(:,:,1)
            !Calculate energy for PT exchange
!           if( mod(ipmcstep, sx_step)==0 ) then
               repene(irep)=0.0_dblprec

               call effective_field(Natom,Mensemble,1,Natom,sx_emomM(:,:,:,irep),   &
                  sx_mmom(:,:,irep),            &
                  external_field,time_external_field,beff,beff1,beff2,OPT_flag,     &
                  max_no_constellations,maxNoConstl,unitCellType,constlNCoup,       &
                  constellations,constellationsNeighType,repene(irep),    &
                  Num_macro,cell_index,sx_emomM_macro(:,:,:,irep),macro_nlistsize,  &
                  NA,N1,N2,N3)

               !repene_sum(irep) = repene(irep)/Natom
               if (ipmcstep>ipmcnstep(1)/2.0_dblprec) then
                  repene_sum(irep) = repene_sum(irep) + repene(irep)/Natom
                  mc_nsamp = mc_nsamp + 1
               end if


!           end if

            ! Rescaling of temperature according to Quantum Heat bath
            temprescale=1.0_dblprec
            !if (do_qhb=='Q' .or. do_qhb=='R' .or. do_qhb=='T') then
            !   if(qhb_mode=='MT') then
            !      call calc_mavrg(Natom,Mensemble,emomM,mavg)
            !      call qhb_rescale(iptemp(irep),temprescale,do_qhb,qhb_mode,mavg)
            !   else
            !      call qhb_rescale(iptemp(irep),temprescale,do_qhb,qhb_mode,dummy)
            !   endif
            !endif

            !! Perform metropolis algorithm
            !call mc_evolve(Natom,Nchmax,Mensemble,nHam,ipTemp(irep),temprescale,ipmode,     &
            !   conf_num,lsf_metric,lsf_window,do_lsf,lsf_field,exc_inter,                &
            !   lsf_interpolate,do_jtensor,do_dm, do_pd, do_biqdm,do_bq,do_chir,mult_axis,&
            !   iflip_a,emomM,emom,mmom,ind_mom_flag,iphfield(1:3),do_dip,Num_macro,      &
            !   max_num_atom_macro_cell,cell_index,macro_nlistsize,macro_atom_nlist,      &
            !   emomM_macro,emom_macro,mmom_macro,do_anisotropy)

            call choose_random_atom_x(Natom,iflip_a)
            if (sx_mode.ne.'S') then
               !!! call mc_evolve(Natom,Nchmax,Mensemble,nHam,ipTemp(irep),temprescale,'M',&
               !!!    conf_num,lsf_metric,lsf_window,do_lsf,lsf_field,ham_inp%exc_inter,           &
               !!!    lsf_interpolate,ham_inp%do_jtensor,ham_inp%do_dm, ham_inp%do_pd, ham_inp%do_biqdm,ham_inp%do_bq, ham_inp%do_ring,    &
               !!!    ham_inp%do_chir, ham_inp%mult_axis,iflip_a,sx_emomM(:,:,:,irep),sx_emom(:,:,:,irep), &
               !!!    sx_mmom(:,:,irep),ind_mom_flag,iphfield(1:3),ham_inp%do_dip,Num_macro,       &
               !!!    max_num_atom_macro_cell,cell_index,macro_nlistsize,macro_atom_nlist, &
               !!!    sx_emomM_macro(:,:,:,irep),sx_emom_macro(:,:,:,irep),                &
               !!!    sx_mmom_macro(:,:,irep),ham_inp%do_anisotropy)
               !call sd_minimal(sx_emomM(:,:,:,irep),sx_emom(:,:,:,irep),sx_mmom(:,:,irep),sx_step,66,ipTemp(irep))
               call mc_minimal(sx_emomM(:,:,:,irep),sx_emom(:,:,:,irep),sx_mmom(:,:,irep),sx_step,sx_mode,ipTemp(irep))
            else if (sx_mode=='S') then
               !call mc_minimal(sx_emomM(:,:,:,irep),sx_emom(:,:,:,irep),sx_mmom(:,:,irep),sx_step,'M',ipTemp(irep))
               !print *,'temp:',ipTemp(irep)
               call sd_minimal(sx_emomM(:,:,:,irep),sx_emom(:,:,:,irep),sx_mmom(:,:,irep),sx_step,1,ipTemp(irep))
            end if



         end do

         !end if

         !!! ! Print m_avg
!        !!! if (mod(ipmcstep,ipmcnstep(1)/10)==0) then

         !!!    !print *,n_hit, n_trial
         !!!    !print '(4f14.6)',repene(:)/Natom
         !!!    print '(2x,a,i5,a,f10.0,a,f6.2,a)', "IP SX ", ipmcstep*100/ipmcnstep(1),"% done.  No. trials: ",n_trial, " Hitrate: ",n_hit/n_trial*100, "%"
         if (ipmcstep>ipmcnstep(1)/2.0_dblprec) then
            print '(2x,a,i5,a,f14.6,a,i3)', "IP SX ", ipmcstep*100/ipmcnstep(1),"% done.  E_min: ", &
               minval(repene_sum/mc_nsamp*sx_numrep) , " mRy/atom for sample ", minloc(repene_sum)
         end if
         !!!    !call calc_mavrg(Natom,Mensemble,emomM,mavg)
         !!!    !if(ipmcstep<=3*ipmcnstep(i)/4) then
         !!!    !   write(*,'(2x,a,i5,a,f10.6)') "IP MC ",ipmcstep*100/ipmcnstep(i),"% done. Mbar:", mavg
         !!!    !else
         !!!    !   write(*,'(2x,a,i5,a,f10.6,a,f10.6,a,f10.6)') "IP MC ",ipmcstep*100/ipmcnstep(i),"% done. Mbar:", mavg,&
         !!!    !      "  U(L):",mc_avcum/mc_nsamp, &
         !!!    !      "  X:",mc_avsus/mc_nsamp/(8.617343d-5)/ (ipTemp(i)+1.0d-15)*Mensemble
         !!!    !end if
         !!!    call sx_prn_mag_conf_clone(Natom,0,sx_numrep,'X',simid,sx_mmom(:,1,:),sx_emom(:,:,1,:),'',mode,iptemp)

         !!!    call choose_random_atom_x(Natom,iflip_a)
!        !!! end if

         !!! !Adjust QHB
         !!! if(do_qhb=='Q' .or. do_qhb=='R' .or. do_qhb=='T') then
         !!!    if (qhb_mode=='MT') then
         !!!       if (mod(ipmcstep,ipmcnstep(i)/20)==0) then
         !!!          call calc_mavrg(Natom,Mensemble,emomM,mavg)
         !!!          call qhb_rescale(iptemp(i),temprescale,do_qhb,qhb_mode,mavg)
         !!!       endif
         !!!    endif
         !!! endif
      enddo
      !!! !PT exchange here (Not for the last 10% of the run (to equilibrate
      !!! !           !temperature)
      !!! !           if( mod(ipmcstep, sx_step)==0 .and.  ipmcstep/(1.0d0*ipmcnstep(1))<0.90d0)  then

      !!! call rng_uniform(sx_flipprob(:),sx_numrep)
      repene_sum = repene_sum / mc_nsamp * sx_numrep

      !!! do ishift=0,1
      !!!    do irep=1+ishift,sx_numrep-1,2

      !!!       n_trial=n_trial+1.0_dblprec


      !!!       temp_1=ipTemp(irep)
      !!!       temp_2=ipTemp(irep+1)
      !!!       beta_1=1.0_dblprec/k_bolt/(temprescale*temp_1+1.0d-15)
      !!!       beta_2=1.0_dblprec/k_bolt/(temprescale*temp_2+1.0d-15)
      !!!       ene_1=repene_sum(irep)*mry/Natom
      !!!       ene_2=repene_sum(irep+1)*mry/Natom
      !!!       if(beta_1==beta_2) then
      !!!          try_pow=(ene_2-ene_1)
      !!!          if(ene_2<ene_1) then
      !!!             trial_prob=1.0_dblprec
      !!!          else
      !!!             trial_prob=0.0_dblprec
      !!!          end if
      !!!       else
      !!!          try_pow=(ene_2-ene_1)*(beta_2-beta_1)
      !!!          trial_prob=min(1.0_dblprec,exp(try_pow))
      !!!       end if

      !!!       if(sx_flipprob(irep)<trial_prob) then
      !!!          temp_tmp=ipTemp(irep)
      !!!          ipTemp(irep)=ipTemp(irep+1)
      !!!          ipTemp(irep+1)=temp_tmp
      !!!          n_hit=n_hit+1.0_dblprec
      !!!       else
      !!!          n_miss=n_miss+1.0_dblprec
      !!!       end if

      !!!    end do
      !!! end do


      do irep=1,sx_numrep
         !!! !Calculate energy for PT exchange
         !!! repene(irep)=0.0_dblprec

         !!! call effective_field(Natom,Mensemble,1,Natom, &
         !!!    sx_emomM(:,:,:,irep),sx_mmom(:,:,irep),external_field,                  &
         !!!    time_external_field,beff,beff1,beff2,OPT_flag,max_no_constellations,    &
         !!!    maxNoConstl,unitCellType,constlNCoup,constellations,                    &
         !!!    constellationsNeighType,repene(irep),Num_macro,cell_index,    &
         !!!    sx_emomM_macro(:,:,:,irep),macro_nlistsize,NA,N1,N2,N3)

         call calc_mavrg(Natom,Mensemble,sx_emomM(:,:,:,irep),mavg)

         print '(a,i4,a,i4,a,f12.6a,f12.6,a,f12.6)', 'Replica #',irep, " List position: ", templist(irep), &
            " Average magnetization:", mavg, " Temperature: ", iptemp(irep),"  Energy: " , repene_sum(irep)!/Natom

      end do

      ! Copy equillibrated moments for use in measurement phase
      ! TODO: replace with minloc()
      irep=1
      !do while (templist(irep)<sx_numrep)
      do while (repene_sum(irep)>minval(repene_sum))
         irep=irep+1
      end do
      !!!irep = minloc(repene_sum)
      print *,' Selected replica:',irep
      emomM=sx_emomM(:,:,:,irep)
      emom=sx_emom(:,:,:,irep)
      mmom=sx_mmom(:,:,irep)
      emom2=sx_emom(:,:,:,irep)
      mmom2=sx_mmom(:,:,irep)
      !print '(1x,a,f10.0,a,f5.2)',"No. PT trials: ",n_trial, " Hitrate: ",n_hit/n_trial*100, "%"

      call timing(0,'Initial       ','OF')
      call timing(0,'PrintRestart  ','ON')
      if (do_mom_legacy.ne.'Y') then 
         call prn_mag_conf(Natom,0,Mensemble,'R',simid,mmom,emom,'',mode)
      else
         call prnrestart(Natom, Mensemble, simid, 0, emom, mmom)
      endif
      call timing(0,'PrintRestart  ','OF')
      call timing(0,'Initial       ','ON')

      ! Deallocate work arrays
      call allocate_mcdata(Natom,-1)

      close(11)


      i_all=-product(shape(sx_emomM))*kind(sx_emomM)
      deallocate(sx_emomM,stat=i_stat)
      call memocc(i_stat,i_all,'sx_emomM','sx_iphase')
      i_all=-product(shape(sx_emomM))*kind(sx_emomM)
      deallocate(sx_emom,stat=i_stat)
      call memocc(i_stat,i_all,'sx_emom','sx_iphase')
      i_all=-product(shape(sx_mmom))*kind(sx_mmom)
      deallocate(sx_mmom,stat=i_stat)
      call memocc(i_stat,i_all,'sx_mmom','sx_iphase')
      !
      i_all=-product(shape(sx_emom_1q))*kind(sx_emom_1q)
      deallocate(sx_emom_1q,stat=i_stat)
      call memocc(i_stat,i_all,'sx_emom_1q','sx_iphase')
      i_all=-product(shape(sx_emom_3q))*kind(sx_emom_3q)
      deallocate(sx_emom_3q,stat=i_stat)
      call memocc(i_stat,i_all,'sx_emom_3q','sx_iphase')
      i_all=-product(shape(sx_emom_rn))*kind(sx_emom_rn)
      deallocate(sx_emom_rn,stat=i_stat)
      call memocc(i_stat,i_all,'sx_emom_rn','sx_iphase')
      i_all=-product(shape(sx_emom_fm))*kind(sx_emom_fm)
      deallocate(sx_emom_fm,stat=i_stat)
      call memocc(i_stat,i_all,'sx_emom_fm','sx_iphase')
      !

      if(Num_macro>0) then
         i_all=-product(shape(sx_emomM_macro))*kind(sx_emomM_macro)
         deallocate(sx_emomM_macro,stat=i_stat)
         call memocc(i_stat,i_all,'sx_emomM_macro','sx_iphase')
         i_all=-product(shape(sx_emomM_macro))*kind(sx_emomM_macro)
         deallocate(sx_emom_macro,stat=i_stat)
         call memocc(i_stat,i_all,'sx_emom_macro','sx_iphase')
         i_all=-product(shape(sx_mmom_macro))*kind(sx_mmom_macro)
         deallocate(sx_mmom_macro,stat=i_stat)
         call memocc(i_stat,i_all,'sx_mmom_macro','sx_iphase')
      end if

      i_all=-product(shape(sx_flipprob))*kind(sx_flipprob)
      deallocate(sx_flipprob,stat=i_stat)
      call memocc(i_stat,i_all,'sx_flipprob','sx_iphase')
      i_all=-product(shape(templist))*kind(templist)
      deallocate(templist,stat=i_stat)
      call memocc(i_stat,i_all,'templist','sx_iphase')
      i_all=-product(shape(replist))*kind(replist)
      deallocate(replist,stat=i_stat)
      call memocc(i_stat,i_all,'replist','sx_iphase')
      i_all=-product(shape(repene))*kind(repene)
      deallocate(repene,stat=i_stat)
      call memocc(i_stat,i_all,'repene','sx_iphase')
      !
      !
      return
      !
   end subroutine sx_iphase

   subroutine sx_copy_ens(irep,hfield)
      !
      use Profiling
      use RandomNumbers, only : rng_uniform
      use InputData, only : Mensemble,Natom
      use MomentData, only : emom
      !
      implicit none
      !
      integer, intent(in) :: irep
      real(dblprec), dimension(3), intent(in) :: hfield
      !
      integer :: i_stat
      integer :: ia,ik
      real(dblprec), dimension(3) :: rn, mom_hat
      real(dblprec) :: sx, sy, sz
      !
      if (irep==1) then
         allocate(sx_emom_1q(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(sx_emom_1q))*kind(sx_emom_1q),'sx_emom_1q','sx_copy_ens')
         sx_emom_1q = emom
      else if (irep==2) then
         allocate(sx_emom_3q(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(sx_emom_3q))*kind(sx_emom_3q),'sx_emom_3q','sx_copy_ens')
         sx_emom_3q = emom
      else if (irep==3) then
         allocate(sx_emom_rn(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(sx_emom_rn))*kind(sx_emom_rn),'sx_emom_rn','sx_copy_ens')
         do ik=1,Mensemble
            do ia=1,Natom
                     ! Call the random number generator to generate the spin directions
                     call rng_uniform(rn,3)
                     sx=2.0_dblprec*(rn(1)-0.50_dblprec)
                     sy=2.0_dblprec*(rn(2)-0.50_dblprec)
                     sz=2.0_dblprec*(rn(3)-0.50_dblprec)
                     do while (sx**2+sy**2+sz**2>1)
                        call rng_uniform(rn,3)
                        sx=rn(1)-0.50_dblprec
                        sy=rn(2)-0.50_dblprec
                        sz=rn(3)-0.50_dblprec
                     end do
                     ! Normalize the spins directions
                     sx_emom_rn(1,ia,ik) = sx/sqrt(sx**2+sy**2+sz**2)
                     sx_emom_rn(2,ia,ik) = sy/sqrt(sx**2+sy**2+sz**2)
                     sx_emom_rn(3,ia,ik) = sz/sqrt(sx**2+sy**2+sz**2)
            end do
         end do
         !sx_emom_rn = emom
      else if (irep==4) then
         allocate(sx_emom_fm(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(sx_emom_fm))*kind(sx_emom_fm),'sx_emom_fm','sx_copy_ens')
         if (norm2(hfield)>0.0_dblprec) then
            mom_hat=hfield/norm2(hfield)
         else
            mom_hat(1)=0.0_dblprec;mom_hat(2)=0.0_dblprec;mom_hat(3)=1.0_dblprec
         end if
         do ik=1,Mensemble
            do ia=1,Natom
               sx_emom_fm(1:3,ia,ik)=mom_hat
            end do
         end do
      end if

   end subroutine sx_copy_ens



   !---------------------------------------------------------------------------------
   ! SUBROUTINE: prn_mag_conf_iter
   !> @brief Prints a given magnetic configuration for either a restartfile or a momentfile
   !> @details Prints a magnetic configuration, the objective is to make all the types
   !> of printing honogeneous, of that way restartfiles, momentfiles and GNEB files
   !> would have all the same structure.
   !> Cloned from prn_mag_conf_iter() 
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine sx_prn_mag_conf_clone(Natom,mstep,Mensemble,type,simid,mmom,emom,suffix, mode,Mtemp)

      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: mstep     !< Current simulation step
      integer, intent(in) :: Mensemble !< Number of ensembles 
      character(len=1), intent(in) :: type   !< type to see whether the file is a restartfile or a momentfile
      character(len=1), intent(in) :: mode   !< Simulation mode (S=SD, M=MC, H=MC Heat Bath, P=LD, C=SLD, G=GNEB)
      character(len=8), intent(in) :: simid  !< Name of simulation
      character(len=*), intent(in) :: suffix !< Suffix to be appended to the files
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(Mensemble), intent(in) :: Mtemp !< Temperature array

      !.. Local variables
      integer :: ii,jj
      character(len=40) filn,fil_pos

      ! Write the name of the restartfile and the position of the writing of file
         write (filn,'(''sx_restart'',a,''.'',a,''.out'')') suffix,trim(simid)
         fil_pos="rewind"
         open(ofileno,file=filn,position=trim(fil_pos))
         write(ofileno,'(a)') repeat("#",80)
         write(ofileno,'(a,1x,a)') "# File type:", type
         write(ofileno,'(a,1x,a)') "# Simulation type:", mode
         write(ofileno,'(a,1x,i8)')"# Number of atoms: ", Natom
         write(ofileno,'(a,1x,i8)')"# Number of ensembles: ", Mensemble
         write(ofileno,'(a)') repeat("#",80)
         write(ofileno,'(a8,a,a8,a16,a16,a16,a16)') "#iter","Temp","iatom","|Mom|","M_x","M_y","M_z"

      do ii=1, Mensemble
         do jj=1, Natom
            write(ofileno,10003) mstep,Mtemp(ii),jj, mmom(jj,ii), emom(1,jj,ii),emom(2,jj,ii),emom(3,jj,ii)
         enddo
      enddo
      close(ofileno)
      10003 format(i8,f8.3,i8,2x,es16.8,es16.8,es16.8,es16.8)

   end subroutine sx_prn_mag_conf_clone

end module sx_driver
