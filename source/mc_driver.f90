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
module mc_driver

   implicit none

   public

contains

   !---------------------------------------------------------------------------------
   !> @brief
   !> Monte Carlo initial phase.
   !> This phase of the program is used for equillibration.
   !> No measurements are performed.
   !> Choose between Monte Carlo, Heat Baths or Spin Dynamics.
   !
   !> @author
   !> Anders Bergman, Lars Bergqvist, Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine mc_iphase()

      use InputData
      use SystemData
      use MonteCarlo
      use DemagField
      use MomentData
      use HamiltonianData
      use Measurements, only : calc_mavrg
      !use InducedMoments,        only : renorm_ncoup_ind
      use AMS
      use QHB
      use Restart
      use macrocells

      integer :: i, k, ipmcstep
      real(dblprec) :: energy,temprescale,temprescalegrad, dummy
      character(len=30) :: filn
      integer :: mc_nmeas,mc_nsamp
      real(dblprec) :: mc_mavg,mc_mavg2,mc_mavg4,mc_minst
      real(dblprec) :: mc_avcum,mc_avsus
      real(dblprec) :: mavg

      ! Write header for moment file
      write (filn,'(''mcinitial.'',a,''.out'')') trim(simid)
      open(11, file=filn, position="append")
      write(11,'(a)') "#  Iter.   M_avg.    U_Bind.    Susc."

      ! Allocate work arrays for Metropolis algorithm
      call allocate_mcdata(Natom,1)

      do i=1,ipmcnphase

         ! Write output to stdout
         write (*,'(a28,i3,a10,G11.4,a10,i10,a10,a10)') &
            "Performing MC initial phase: ", i ," Temp: ", ipTemp(i), "MCS/s: ", ipmcnstep(i), " Mode: ",ipmode

         ! Calculate demagnetization field
         if(demag=='Y') then
            call calc_demag(Natom,Mensemble,demag1,demag2,demag3,demagvol,emomM)
         endif

         ! Rescaling of temperature according to Quantum Heat bath
         temprescale=1.0_dblprec
         temprescalegrad=0.0_dblprec
         if (do_qhb=='Q' .or. do_qhb=='R' .or. do_qhb=='P' .or. do_qhb=='T') then
            if(qhb_mode=='MT') then
               call calc_mavrg(Natom,Mensemble,emomM,mavg)
               call qhb_rescale(iptemp(i),temprescale,temprescalegrad,do_qhb,qhb_mode,mavg)
            else
               call qhb_rescale(iptemp(i),temprescale,temprescalegrad,do_qhb,qhb_mode,dummy)
            endif
         endif

         ! Set up order for Metropolis sweeps
         call choose_random_atom_x(Natom,iflip_a)

         ! Zero data for averaging
         mc_nmeas=0
         mc_nsamp=0
         mc_mavg=0.0_dblprec
         mc_mavg2=0.0_dblprec
         mc_mavg4=0.0_dblprec
         mc_avcum=0.0_dblprec
         mc_avsus=0.0_dblprec
         !
         ! Set energy to zero
         energy = 0.0_dblprec
         ! Loop over steps of sweeps
         do ipmcstep=1,ipmcnstep(i)
            ! Perform metropolis algorithm
            call mc_evolve(Natom,Nchmax,Mensemble,nHam,ipTemp(i),temprescale,ipmode,   &
               conf_num,lsf_metric,lsf_window,do_lsf,lsf_field,ham_inp%exc_inter,              &
               lsf_interpolate,ham_inp%do_jtensor,ham_inp%do_dm, ham_inp%do_pd, &
               ham_inp%do_biqdm,ham_inp%do_bq,ham_inp%do_ring,ham_inp%do_chir,ham_inp%do_sa,&
               ham_inp%mult_axis,iflip_a,emomM,emom,mmom,ind_mom_flag,iphfield(1:3),ham_inp%do_dip,    &
               Num_macro,max_num_atom_macro_cell,cell_index,macro_nlistsize,           &
               macro_atom_nlist,emomM_macro,emom_macro,mmom_macro,ham_inp%do_anisotropy)

            ! Sample m, m2, m4 for second half of run (every tenth step)
            if(ipmcstep>ipmcnstep(i)/2.and.mod(ipmcstep-1,10)==0) then

               ! Calculate m_avg
               mc_minst=0.0_dblprec
               do k=1,Mensemble
                  mc_minst=sqrt(sum(emomM(1,:,k))**2+sum(emomM(2,:,k))**2+sum(emomM(3,:,k))**2)/natom/Mensemble
                  mc_mavg=mc_mavg+mc_minst
                  mc_mavg2=mc_mavg2+mc_minst*mc_minst
                  mc_mavg4=mc_mavg4+mc_minst*mc_minst*mc_minst*mc_minst
                  mc_nmeas=mc_nmeas+1
               end do

               ! Calculate Binder cumulant for final quarter of run
               if(ipmcstep>3*ipmcnstep(i)/4) then
                  mc_avcum=mc_avcum+1.0_dblprec-(mc_mavg4/mc_nmeas)/((mc_mavg2/mc_nmeas)**2)/3.0_dblprec
                  mc_avsus=mc_avsus+((mc_mavg2/mc_nmeas)-(mc_mavg/mc_nmeas)**2)
                  mc_nsamp=mc_nsamp+1

                  ! Print Binder cumulant and susceptibility for final tenth of run
                  if(ipmcstep>ipmcnstep(i)*0.9_dblprec) then
                     write(11,'(i8,4f10.6,3f18.6)') ipmcstep,mc_mavg/mc_nmeas,&
                        mc_avcum/mc_nsamp, mc_avsus/mc_nsamp/(8.617343d-5)/ (ipTemp(i)+1.0d-15) *Mensemble
                  end if
               end if
            end if

            ! Print m_avg
            if (mod(ipmcstep,ipmcnstep(i)/10)==0) then

               call calc_mavrg(Natom,Mensemble,emomM,mavg)
               if(ipmcstep<=3*ipmcnstep(i)/4) then
                  write(*,'(2x,a,i5,a,f10.6)') "IP MC ",ipmcstep*100/ipmcnstep(i),"% done. Mbar:", mavg
               else
                  write(*,'(2x,a,i5,a,f10.6,a,f10.6,a,f10.6)') "IP MC ",ipmcstep*100/ipmcnstep(i),"% done. Mbar:", mavg,&
                     "  U(L):",mc_avcum/mc_nsamp, &
                     "  X:",mc_avsus/mc_nsamp/(8.617343d-5)/ (ipTemp(i)+1.0d-15)*Mensemble
               end if

               call choose_random_atom_x(Natom,iflip_a)
            end if
            !Adjust QHB
            if(do_qhb=='Q' .or. do_qhb=='R' .or. do_qhb=='P' .or. do_qhb=='T') then
               if (qhb_mode=='MT') then
                  if (mod(ipmcstep,ipmcnstep(i)/20)==0) then
                     call calc_mavrg(Natom,Mensemble,emomM,mavg)
                     call qhb_rescale(iptemp(i),temprescale,temprescalegrad,do_qhb,qhb_mode,mavg)
                  endif
               endif
            endif
         !------------------------------------------------------------------------
         ! Induced moments treatment
         !------------------------------------------------------------------------
         ! If there are induced moments in the system one must renormalize the
         ! the Heisenberg exchange and Dzyaloshinskii-Moriya vectors, as well as
         ! the magnitude and direction of the induced moments
         !if (ind_mom_flag=='Y') then
         !   call renorm_ncoup_ind(do_dm,Natom,conf_num,Mensemble, &
         !      ham%max_no_neigh,ham%max_no_dmneigh,ham%max_no_neigh_ind, &
         !      ham%nlistsize,ham%dmlistsize,ham%ind_list_full,ham%ind_nlistsize, &
         !      ham%nlist,ham%dmlist,ham%ind_nlist,ham%sus_ind,mmom,emom,emomM, &
         !      ham%ncoup,ham%dm_vect,ham%fix_list,ham%fix_num)
         !endif
         !------------------------------------------------------------------------
         ! End ! of ! the ! induced ! moments ! treatment
         !------------------------------------------------------------------------

         enddo

         call timing(0,'Initial       ','OF')
         call timing(0,'PrintRestart  ','ON')
         if (do_mom_legacy.ne.'Y') then
            call prn_mag_conf(Natom,0,Mensemble,'R',simid,mmom,emom,'',mode)
         else
            call prnrestart(Natom, Mensemble, simid, 0, emom, mmom)
         endif
         call timing(0,'PrintRestart  ','OF')
         call timing(0,'Initial       ','ON')
      enddo

      ! Copy equillibrated moments for use in measurement phase
      emom2=emom
      mmom2=mmom

      ! Deallocate work arrays
      call allocate_mcdata(Natom,-1)

      close(11)
      !
   end subroutine mc_iphase



   !---------------------------------------------------------------------------
   !> @brief
   !> Monte Carlo measurement phase
   !
   !> @author
   !> Anders Bergman, Lars Bergqvist
   !---------------------------------------------------------------------------
   subroutine mc_mphase()
      !
      use Energy,          only : calc_energy
      use InputData
      use FieldData,       only : external_field, sitefld, time_external_field,&
         allocation_field_time, thermal_field, beff, beff1, beff3
      use SystemData
      use MonteCarlo
      use DemagField
      use MomentData
      use FieldPulse
      use Correlation
      use Polarization
      use prn_averages
      use Measurements
      use ChemicalData,    only : achem_ch, asite_ch
      use MicroWaveField
      use CalculateFields
      use HamiltonianData
      use AutoCorrelation,       only : autocorr_sample, do_autocorr
      use ChemicalData, only : achtype
      use QHB, only : qhb_rescale, do_qhb, qhb_mode, qhb_Tmix_prn,&
             do_qhb_mix, qhb_Tmix_cumu 
      !use InducedMoments,        only : renorm_ncoup_ind
      use macrocells
      use optimizationRoutines

      !
      implicit none
      !
      integer :: cgk_flag, scount_pulse, bcgk_flag, cgk_flag_pc, mcmstep
      real(dblprec) :: temprescale,temprescalegrad,dummy,totene

      call timing(0,'MonteCarlo    ','ON')

      ! Flag for G(q) sampling and G(r) calculation
      cgk_flag=0 ; scount_pulse=1 ; bcgk_flag=0 ; cgk_flag_pc=0

      ! Allocate work arrays for MC
      call allocate_mcdata(Natom,1)

      ! Calculation of cumulant
      if(do_cumu .ne. 'A') do_cumu = 'Y'

      ! Setup order for Metropolis sweeps
      call choose_random_atom_x(Natom,iflip_a)


      ! Calculate demagnetizing field
      if(demag=='Y') then
         call calc_demag(Natom, Mensemble, demag1, demag2, demag3, demagvol, emomM)
      endif

      ! Rescaling of temperature according to Quantum Heat bath
      temprescale=1.0_dblprec
      temprescalegrad=0.0_dblprec
      if (do_qhb=="Q" .or. do_qhb=='R' .or. do_qhb=='P' .or. do_qhb=='T') then
         if(qhb_mode=='MT') then
            call calc_mavrg(Natom,Mensemble,emomM,mavg)
            call qhb_rescale(Temp,temprescale,temprescalegrad,do_qhb,qhb_mode,mavg)
         else
            call qhb_rescale(Temp,temprescale,temprescalegrad,do_qhb,qhb_mode,dummy)
         endif
      endif

      if (do_pol=='Y') then
         call init_polarization(Natom,Mensemble,ham%max_no_neigh,ham%nlist,coord,1)
      end if

      ! Calculate the static magnetic fields which will be calculated only once as they are not time dependent
      call calc_external_fields(Natom,Mensemble,hfield,anumb,external_field,     &
         do_bpulse,sitefld,sitenatomfld)

      ! Perform MC sweeps
      do mcmstep=1,mcnstep

         call timing(0,'MonteCarlo    ','OF')
         call timing(0,'Measurement   ','ON')

         ! Sample autocorrelation
         if(do_autocorr=='Y') then
            call autocorr_sample(Natom, Mensemble, simid, mcmstep, emom)
         end if

         ! Measure averages and trajectories
         call measure(Natom,Mensemble,NT,NA,nHam,N1,N2,N3,simid,mcmstep,emom,emomM, &
            mmom,Nchmax,do_ralloy,Natom_full,asite_ch,achem_ch,atype,plotenergy,    &
            Temp,temprescale,temprescalegrad,'N',1.0_dblprec,logsamp,               &
            ham%max_no_neigh,ham%nlist,ham%ncoup,                                   &
            ham%nlistsize,ham%aham,thermal_field,beff,beff1,beff3,coord,            &
            ham%ind_list_full,ham%ind_nlistsize,ham%ind_nlist,ham%max_no_neigh_ind, &
            ham%sus_ind,do_mom_legacy,mode)
         call timing(0,'Measurement   ','OF')

         ! Calculate total and term resolved energies
         if(plotenergy>0.and.mod(mcmstep-1,cumu_step)==0) then

            call timing(0,'Energy        ','ON')
            totene=0.0_dblprec

            call calc_energy(nHam,mcmstep,Natom,Nchmax,       &
               conf_num,Mensemble,Natom,Num_macro,1, &
               plotenergy,Temp,1.0_dblprec,do_lsf,    &
               lsf_field,lsf_interpolate,'N',simid,cell_index,macro_nlistsize,mmom,     &
               emom,emomM,emomM_macro,external_field,time_external_field,               &
               max_no_constellations,maxNoConstl,unitCellType,constlNCoup,              &
               constellations,OPT_flag,constellationsNeighType,totene,NA,N1,N2,N3)
            call timing(0,'Energy        ','OF')
         endif

         call timing(0,'SpinCorr      ','ON')
         ! Spin correlation
         ! Sample magnetic moments for correlation functions
         ! Sample S(r) through S(q)
         call correlation_wrapper(Natom,Mensemble,coord,simid,emomM,mcmstep,delta_t,&
            NT,atype,Nchmax,achtype,sc,do_sc,do_sr,cgk_flag)

         call timing(0,'SpinCorr      ','OF')
         call timing(0,'MonteCarlo    ','ON')

         ! Metropolis sweeps

            call mc_evolve(Natom,Nchmax,Mensemble,nHam,Temp,temprescale,mode,   &
               conf_num,lsf_metric,lsf_window,do_lsf,lsf_field,ham_inp%exc_inter,              &
               lsf_interpolate,ham_inp%do_jtensor,ham_inp%do_dm, ham_inp%do_pd, &
               ham_inp%do_biqdm,ham_inp%do_bq,ham_inp%do_ring,ham_inp%do_chir,ham_inp%do_sa,&
               ham_inp%mult_axis,iflip_a,emomM,emom,mmom,ind_mom_flag,hfield,ham_inp%do_dip,    &
               Num_macro,max_num_atom_macro_cell,cell_index,macro_nlistsize,           &
               macro_atom_nlist,emomM_macro,emom_macro,mmom_macro,ham_inp%do_anisotropy)

         !call mc_evolve(Natom,Nchmax,Mensemble,nHam,Temp,temprescale,mode,conf_num,       &
         !   lsf_metric,lsf_window,do_lsf,lsf_field,exc_inter,lsf_interpolate,             &
         !   do_jtensor,do_dm,do_pd,do_biqdm,do_bq,do_ring,do_chir,mult_axis,iflip_a,      &
         !   emomM,emom,mmom,ind_mom_flag,hfield,do_dip,Num_macro,max_num_atom_macro_cell, &
         !   cell_index,macro_nlistsize,macro_atom_nlist,emomM_macro,emom_macro,           &
         !   mmom_macro,do_anisotropy)
         
         ! Calculate Tmix (classic+quantum Tsim)
         call qhb_Tmix_cumu(mcmstep,cumu_step)

         ! Calculate and print m_avg
         if(mcnstep>20) then
            if(mod(mcmstep,mcnstep/20)==0) then

               call calc_mavrg(Natom,Mensemble,emomM,mavg)
               write(*,'(2x,a,i3,a,f10.6)',advance='no') &
                  "MP MC ",mcmstep*100/(mcnstep),"% done. Mbar:",mavg
               if(do_qhb_mix=='Y') then
                  write(*,'(a,f7.2)',advance='no') ". Tmix:",qhb_Tmix_prn
               end if 
               if(plotenergy>0) then
                  write(*,'(a,f12.6,a,f8.5,a)') ". Ebar:", totene,". U:",binderc,"."
               else
                  write(*,'(a,f8.5,a)') ". U:",binderc,"."
               end if
            end if
         else
            write(*,'(2x,a,i3,a,G13.6)')   "Iteration",mcmstep," Mbar ",mavg
         end if

         !Adjust QHB
         if(do_qhb=='Q' .or. do_qhb=='R' .or. do_qhb=='P' .or. do_qhb=='T') then
            if(qhb_mode=='MT') then
               if (mod(mcmstep,mcnstep/40)==0) then
                  call calc_mavrg(Natom,Mensemble,emomM,mavg)
                  call qhb_rescale(Temp,temprescale,temprescalegrad,do_qhb,qhb_mode,mavg)
               endif
            endif
         endif

         if(mod(mcmstep,mcnstep/10)==0) then
            ! Change order for Metropolis sweeps
            call choose_random_atom_x(Natom,iflip_a)
         end if

         !------------------------------------------------------------------------
         ! Induced moments treatment
         !------------------------------------------------------------------------
         ! If there are induced moments in the system one must renormalize the
         ! the Heisenberg exchange and Dzyaloshinskii-Moriya vectors, as well as
         ! the magnitude and direction of the induced moments
         !if (ind_mom_flag=='Y') then
         !   call renorm_ncoup_ind(do_dm,Natom,conf_num,Mensemble, &
         !      ham%max_no_neigh,ham%max_no_dmneigh,ham%max_no_neigh_ind, &
         !      ham%nlistsize,ham%dmlistsize,ham%ind_list_full,ham%ind_nlistsize, &
         !      ham%nlist,ham%dmlist,ham%ind_nlist,ham%sus_ind,mmom,emom,emomM, &
         !      ham%ncoup,ham%dm_vect,ham%fix_list,ham%fix_num)
         !endif
         !------------------------------------------------------------------------
         ! End ! of ! the ! induced ! moments ! treatment
         !------------------------------------------------------------------------


      enddo

      call timing(0,'MonteCarlo    ','OF')
      call timing(0,'Measurement   ','ON')

      ! Print remaining measurements
      call flush_measurements(Natom,Mensemble,NT,NA,N1,N2,N3,simid,mcmstep,emom,    &
         mmom,Nchmax,atype,'N',mcnstep,ham%ind_list_full,do_mom_legacy,mode)

      if (do_pol=='Y') then
         call init_polarization(Natom,Mensemble,ham%max_no_neigh,ham%nlist,coord,-1)
      end if

      ! Deallocate work arrays
      call allocate_mcdata(Natom,-1)
      close(11)
      call timing(0,'Measurement   ','OF')
   end subroutine mc_mphase

   !---------------------------------------------------------------------------------
   !> @brief
   !> Monte Carlo minimal driver
   !> This phase of the program is used for equillibration.
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   subroutine mc_minimal(emomM_io,emom_io,mmom_io,niter,mcmode,mctemp)

      use InputData
      use SystemData
      use MonteCarlo
      use DemagField
      use MomentData
      use HamiltonianData
      use Measurements, only : calc_mavrg
      !use InducedMoments,        only : renorm_ncoup_ind
      use AMS
      use QHB
      use Restart
      use macrocells

      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM_io
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom_io
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom_io
      integer, intent(in) :: niter
      character(len=2), intent(in) :: mcmode
      real(dblprec), intent(in) :: mctemp

      integer :: ipmcstep, ia,ik
      real(dblprec) :: energy,temprescale,temprescalegrad, dummy
      real(dblprec) :: mavg

      ! Copy inmoments to working array
      do ik=1,Mensemble
         do ia=1,Natom
            emom(:,ia,ik)=emom_io(:,ia,ik)
            mmom(ia,ik)=mmom_io(ia,ik)
            emomM(:,ia,ik)=emomM_io(:,ia,ik)
            mmomi(ia,ik) = 1.0_dblprec/mmom(ia,ik)
         end do
      end do

      ! Calculate demagnetization field
      if(demag=='Y') then
         call calc_demag(Natom,Mensemble,demag1,demag2,demag3,demagvol,emomM)
      endif

      ! Rescaling of temperature according to Quantum Heat bath
      temprescale=1.0_dblprec
      temprescalegrad=0.0_dblprec
      if (do_qhb=='Q' .or. do_qhb=='R' .or. do_qhb=='P' .or. do_qhb=='T') then
         if(qhb_mode=='MT') then
            call calc_mavrg(Natom,Mensemble,emomM,mavg)
            call qhb_rescale(iptemp(1),temprescale,temprescalegrad,do_qhb,qhb_mode,mavg)
         else
            call qhb_rescale(iptemp(1),temprescale,temprescalegrad,do_qhb,qhb_mode,dummy)
         endif
      endif

      ! Set up order for Metropolis sweeps
      call choose_random_atom_x(Natom,iflip_a)
      !
      ! Set energy to zero
      energy = 0.0_dblprec
      ! Loop over steps of sweeps
      do ipmcstep=1,niter
         ! Perform metropolis algorithm
         call mc_evolve(Natom,Nchmax,Mensemble,nHam,mctemp,temprescale,mcmode,   &
            conf_num,lsf_metric,lsf_window,do_lsf,lsf_field,ham_inp%exc_inter,              &
            lsf_interpolate,ham_inp%do_jtensor,ham_inp%do_dm, ham_inp%do_pd, &
            ham_inp%do_biqdm,ham_inp%do_bq,ham_inp%do_ring,ham_inp%do_chir,ham_inp%do_sa,&
            ham_inp%mult_axis,iflip_a,emomM,emom,mmom,ind_mom_flag,iphfield(1:3),ham_inp%do_dip,    &
            Num_macro,max_num_atom_macro_cell,cell_index,macro_nlistsize,           &
            macro_atom_nlist,emomM_macro,emom_macro,mmom_macro,ham_inp%do_anisotropy)
      enddo

      ! Copy equillibrated moments for use in measurement phase
      emom2=emom
      mmom2=mmom

      ! Copy working moments to outdata
      do ik=1,Mensemble
         do ia=1,Natom
            emom_io(:,ia,ik)=emom(:,ia,ik)
            mmom_io(ia,ik)=mmom(ia,ik)
            emomM_io(:,ia,ik)=emomM(:,ia,ik)
         end do
      end do

      !
   end subroutine mc_minimal


end module mc_driver
