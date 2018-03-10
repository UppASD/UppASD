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

   !---------------------------------------------------------------------------
   !> @brief
   !> Monte Carlo initial phase.
   !> This phase of the program is used for equillibration.
   !> No measurements are performed.
   !> Choose between Monte Carlo, Heat Baths or Spin Dynamics.
   !
   !> @author
   !> Anders Bergman, Lars Bergqvist, Jonathan Chico
   !---------------------------------------------------------------------------
   subroutine mc_iphase()

      !      use Ewaldmom
      use InputData
      use SystemData
      use MonteCarlo
      use DemagField
      use MomentData
      use HamiltonianData
      use ChemicalData, only : atype_ch
      use Measurements, only : calc_mavrg
      use AMS
      use QHB
      use Restart

      integer :: i, k, ipmcstep
      real(dblprec) :: energy,temprescale,dummy
      character(len=30) :: filn
      integer :: mc_nmeas,mc_nsamp
      real(dblprec) :: mc_mavg,mc_mavg2,mc_mavg4,mc_minst
      real(dblprec) :: mc_avcum,mc_avsus
      real(dblprec) :: mavg

      ! Write header for moment file
      write (filn,'(''mcinitial.'',a8,''.out'')') simid
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
            call calc_demag(Natom, Mensemble, demag1, demag2, demag3, demagvol, emomM)
         endif

         ! Rescaling of temperature according to Quantum Heat bath
         temprescale=1.d0
         if (do_qhb=='Q' .or. do_qhb=='R' .or. do_qhb=='T') then
            if(qhb_mode=='MT') then
               call calc_mavrg(Natom,Mensemble,emomM,mavg)
               call qhb_rescale(iptemp(i),temprescale,do_qhb,qhb_mode,mavg)
            else
               call qhb_rescale(iptemp(i),temprescale,do_qhb,qhb_mode,dummy)
            endif
         endif

         !! Initialize counter for MC spin flips

         ! Set up order for Metropolis sweeps
         call choose_random_atom_x(Natom,iflip_a)

         ! Zero data for averaging
         mc_nmeas=0
         mc_nsamp=0
         mc_mavg=0.0d0
         mc_mavg2=0.0d0
         mc_mavg4=0.0d0
         mc_avcum=0.0d0
         mc_avsus=0.0d0
         !
         ! Set energy to zero
         energy = 0.0d0
         ! Loop over steps of sweeps
         do ipmcstep=1,ipmcnstep(i)
            ! Perform metropolis algorithm
            call mc_evolve(NA,Natom,Nchmax,Mensemble,do_ralloy,Natom_full,atype,atype_ch,ipTemp(i),temprescale,&
               ipmode,conf_num,lsf_metric,fs_nlistsize,fs_nlist,nind,lsf_window,do_lsf,lsf_field,exc_inter,&
               lsf_interpolate,do_jtensor,max_no_neigh,nlistsize,nlist,ncoup,ncoupD,j_tens,do_dm,max_no_dmneigh,&
               dmlistsize,dmlist,dm_vect,do_pd,nn_pd_tot,pdlistsize,pdlist,pd_vect,do_biqdm,nn_biqdm_tot,&
               biqdmlistsize,biqdmlist,biqdm_vect,do_bq,nn_bq_tot,bqlistsize,bqlist,j_bq,taniso,taniso_diff,&
               eaniso,eaniso_diff,kaniso,kaniso_diff,sb,sb_diff,mult_axis,iflip_a,emomM,emom,mmom,ind_nlistsize,&
               ind_nlist,ind_mom,sus_ind,ind_mom_flag,iphfield(1:3),do_dip,Qdip)

            ! Sample m, m2, m4 for second half of run (every tenth step)
            if(ipmcstep>ipmcnstep(i)/2.and.mod(ipmcstep-1,10)==0) then

               ! Calculate m_avg
               mc_minst=0.0d0
               do k=1,Mensemble
                  mc_minst=sqrt(sum(emomM(1,:,k))**2+sum(emomM(2,:,k))**2+sum(emomM(3,:,k))**2)/natom/Mensemble
                  mc_mavg=mc_mavg+mc_minst
                  mc_mavg2=mc_mavg2+mc_minst*mc_minst
                  mc_mavg4=mc_mavg4+mc_minst*mc_minst*mc_minst*mc_minst
                  mc_nmeas=mc_nmeas+1
               end do

               ! Calculate Binder cumulant for final quarter of run
               if(ipmcstep>3*ipmcnstep(i)/4) then
                  mc_avcum=mc_avcum+1.0d0-(mc_mavg4/mc_nmeas)/((mc_mavg2/mc_nmeas)**2)/3.0d0
                  mc_avsus=mc_avsus+((mc_mavg2/mc_nmeas)-(mc_mavg/mc_nmeas)**2)
                  mc_nsamp=mc_nsamp+1

                  ! Print Binder cumulant and susceptibility for final tenth of run
                  if(ipmcstep>ipmcnstep(i)*0.9d0) then
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
            if(do_qhb=='Q' .or. do_qhb=='R' .or. do_qhb=='T') then
               if (qhb_mode=='MT') then
                  if (mod(ipmcstep,ipmcnstep(i)/20)==0) then
                     call calc_mavrg(Natom,Mensemble,emomM,mavg)
                     call qhb_rescale(iptemp(i),temprescale,do_qhb,qhb_mode,mavg)
                  endif
               endif
            endif

         enddo

         call timing(0,'Initial       ','OF')
         call timing(0,'PrintRestart  ','ON')
         call prnrestart(Natom, Mensemble, simid, 0, emom, mmom)
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
      use Correlation,          only : correlation_wrapper
      use Polarization
      use prn_averages
      use Measurements
      use ChemicalData,    only : achem_ch, asite_ch,atype_ch
      use MicroWaveField
      use CalculateFields
      use HamiltonianData
      use AutoCorrelation, only : autocorr_sample
      use ChemicalData, only : achtype
      use QHB, only : qhb_rescale, do_qhb, qhb_mode

      !
      implicit none

      !
      integer :: cgk_flag, scount_pulse, bcgk_flag, cgk_flag_pc, mcmstep
      real(dblprec) :: temprescale,dummy

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
      temprescale=1.d0
      if (do_qhb=="Q" .or. do_qhb=='R' .or. do_qhb=='T') then
         if(qhb_mode=='MT') then
            call calc_mavrg(Natom,Mensemble,emomM,mavg)
            call qhb_rescale(Temp,temprescale,do_qhb,qhb_mode,mavg)
         else
            call qhb_rescale(Temp,temprescale,do_qhb,qhb_mode,dummy)
         endif
      endif

      if (do_pol=='Y') then
         call init_polarization(Natom, Mensemble, max_no_neigh, &
            nlist, coord, 1)
      end if

      ! Calculate the static magnetic fields which will be calculated only once as they are not time dependent
      call calc_external_fields(Natom, Mensemble, NA, hfield, anumb, external_field, &
         do_bpulse,sitefld,sitenatomfld)

      ! Perform MC sweeps
      do mcmstep=1,mcnstep

         call timing(0,'MonteCarlo    ','OF')
         call timing(0,'Measurement   ','ON')

         ! Measure averages and trajectories
         call measure(Natom, Mensemble, NT, NA, N1, N2, N3, simid, mcmstep, emom, emomM, mmom, &
            Nchmax, do_ralloy, Natom_full, asite_ch, achem_ch, atype,  plotenergy, Temp, &
            'N',1.0d0, logsamp, max_no_neigh, nlist, ncoup,nlistsize,thermal_field, &
            beff,beff1,beff3,coord,ind_mom,ind_nlistsize,ind_nlist,atype_ch)

         ! Calculate total and term resolved energies
         if(plotenergy>0.and.mod(mcmstep-1,cumu_step)==0) then

            call timing(0,'Measurement   ','OF')
            call timing(0,'Energy        ','ON')
            totene=0.0d0
            call calc_energy(Natom, Nchmax,Mensemble, conf_num, emom, emomM, mmom,simid, plotenergy, mcmstep, external_field,time_external_field,tenergy, eenergy, lsfenergy, totene, &
               max_no_neigh, nlistsize, nlist, ncoup, ncoupD,exc_inter,do_dm, max_no_dmneigh, dmlistsize, dmlist, dm_vect, &
               do_pd, nn_pd_tot, pdlistsize, pdlist, pd_vect, &
               do_biqdm, nn_biqdm_tot, biqdmlistsize, biqdmlist, biqdm_vect, &
               do_bq, nn_bq_tot, bqlistsize, bqlist, j_bq, &
               do_dip, qdip,taniso, eaniso, kaniso,sb,&
               mult_axis,taniso_diff, eaniso_diff, kaniso_diff,sb_diff,&
               1.0d0,'N',do_lsf,fs_nlist,fs_nlistsize,nind,lsf_interpolate,lsf_field,Temp)
            call timing(0,'Energy        ','OF')
            call timing(0,'Measurement   ','ON')
         endif

         ! Spin correlation
         ! Sample magnetic moments for correlation functions
         ! Sample S(r) through S(q)
         call correlation_wrapper(Natom, Mensemble, coord, simid,emomM, mcmstep, delta_t, NT, atype, Nchmax,achtype,cgk_flag,cgk_flag_pc)

         call timing(0,'Measurement   ','OF')
         call timing(0,'MonteCarlo    ','ON')

         ! Metropolis sweeps

         call mc_evolve(NA,Natom,Nchmax,Mensemble,do_ralloy,Natom_full,atype,atype_ch,Temp,temprescale,&
            mode,conf_num,lsf_metric,fs_nlistsize,fs_nlist,nind,lsf_window,do_lsf,lsf_field,exc_inter,&
            lsf_interpolate,do_jtensor,max_no_neigh,nlistsize,nlist,ncoup,ncoupD,j_tens,do_dm,max_no_dmneigh,&
            dmlistsize,dmlist,dm_vect,do_pd,nn_pd_tot,pdlistsize,pdlist,pd_vect,do_biqdm,nn_biqdm_tot,&
            biqdmlistsize,biqdmlist,biqdm_vect,do_bq,nn_bq_tot,bqlistsize,bqlist,j_bq,taniso,taniso_diff,&
            eaniso,eaniso_diff,kaniso,kaniso_diff,sb,sb_diff,mult_axis,iflip_a,emomM,emom,mmom,ind_nlistsize,&
            ind_nlist,ind_mom,sus_ind,ind_mom_flag,hfield,do_dip,Qdip)

         ! Calculate and print m_avg
         if(mcnstep>20) then
            if(mod(mcmstep,mcnstep/20)==0) then

               call calc_mavrg(Natom,Mensemble,emomM,mavg)
               write(*,'(2x,a,i3,a,f10.6)',advance='no') &
                  "MP MC ",mcmstep*100/(mcnstep),"% done. Mbar:",mavg
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
         if(do_qhb=='Q' .or. do_qhb=='R' .or. do_qhb=='T') then
            if(qhb_mode=='MT') then
               if (mod(mcmstep,mcnstep/40)==0) then
                  call calc_mavrg(Natom,Mensemble,emomM,mavg)
                  call qhb_rescale(Temp,temprescale,do_qhb,qhb_mode,mavg)
               endif
            endif
         endif

         if(mod(mcmstep,mcnstep/10)==0) then
            ! Change order for Metropolis sweeps
            call choose_random_atom_x(Natom,iflip_a)
         end if

      enddo

      call timing(0,'MonteCarlo    ','OF')
      call timing(0,'Measurement   ','ON')

      ! Print remaining measurements
      call flush_measurements(Natom, Mensemble, NT, NA, N1, N2, N3, simid, mcmstep, emom, mmom, &
         Nchmax,atype,'N',mcnstep,do_ralloy,Natom_full,atype_ch,ind_mom)

      if (do_pol=='Y') then
         call init_polarization(Natom, Mensemble, max_no_neigh, &
            nlist, coord, -1)
      end if

      ! Deallocate work arrays
      call allocate_mcdata(Natom,-1)
      close(11)
      call timing(0,'Measurement   ','OF')
   end subroutine mc_mphase


end module mc_driver
