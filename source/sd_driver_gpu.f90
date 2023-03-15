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
module sd_driver_gpu

   implicit none

   public

contains


   !---------------------------------------------------------------------------
   !> @brief
   !> Spin Dynamics measurement phase
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine sd_mphase_gpu()
      !
      use KMC
      use QHB,                   only : qhb_rescale, do_qhb, qhb_mode
      use BLS,                   only : calc_bls
      use Energy,                only : calc_energy
      use Sparse
      use Damping
      use KMCData
      use FixedMom
      use Gradients
      use Evolution
      use Evolution_gpu
      use InputData
      use FieldData,             only : beff,beff1,beff2,beff3,b2eff,sitefld,       &
         external_field,field1,field2,time_external_field,allocation_field_time,    &
         thermal_field
      use macrocells
      use DemagField
      use MomentData
      use FieldPulse
      use SystemData,            only: coord
      use SystemData,            only : atype, anumb, Landeg
      use Correlation
      use Temperature
      use SpinTorques
      use ChemicalData
      use prn_averages
      use Measurements
      use Polarization
      use UpdateMoments
      !use InducedMoments,        only : renorm_ncoup_ind
      use MicroWaveField
      use SimulationData,        only : bn, rstep, mstep
      use Math_functions, only : f_logstep
      use HamiltonianData
      use CalculateFields
      use AutoCorrelation,       only : autocorr_sample, do_autocorr
      use prn_trajectories
      use HamiltonianActions
      use HamiltonianActions_gpu
      use OptimizationRoutines
      use AdaptiveTimeStepping
      use MetaTypes

      implicit none
      logical :: time_dept_flag, deltat_correction_flag
      integer :: cgk_flag,cgk_flag_p,adapt_step,cr_flag,spt_flag,ntmp
      integer :: bcgk_flag,cgk_flag_pc
      integer :: scount_pulse, sstep
      ! Phase flag indicator (true for sd initial phase; false for sd measurement phase)
      ! Used in order to separate between phases in an adaptive time step environment with spin correlation
      logical :: sd_phaseflag

      real(dblprec) :: temprescale, temprescalegrad, totene, totenergy,dummy

      ! Spin correlation measurements allowed
      adapt_step = 0
      deltat_correction_flag = .true.
      ! Measurement phase indicator (false here)
      sd_phaseflag = .false.
      call timing(0,'Measurement   ','ON')

      ! Flag for G(q) sampling and G(r) calculation
      cgk_flag=0 ; cgk_flag_p=0 ; cr_flag=0 ; bcgk_flag=0 ; cgk_flag_pc=0
      spt_flag=0

      scount_pulse = 1

      time_dept_flag=allocation_field_time(mwf,do_gauss,mwf_gauss,mov_gauss,        &
         mwf_mov_gauss,mwf_gauss_spatial,mov_circle,mwf_mov_circle,mov_square,      &
         mwf_mov_square)

      time_dept_flag = time_dept_flag .or. (do_bpulse.gt.0.and.do_bpulse.lt.5)

      ! Calculate demagnetization field
      if(demag=='Y') then
         call calc_demag(Natom,Mensemble,demag1,demag2,demag3,demagvol,emomM)
      endif

      call calc_bls(N1,N2,N3,C1,C2,C3,Natom,Mensemble,simid,coord,emomM,sstep,      &
         delta_t,0)

      !  Iniitialize magnetic field pulse
      if (do_bpulse.gt.0.and.do_bpulse.lt.5) then
         bpulse_time = delta_t*rstep
         call bpulse(do_bpulse, bpulse_time)
      end if

      ! Initialize polarization including local chirality and polarity
      if (do_spintemp=='Y') then
         ntmp=nstep/spintemp_step
         call spintemperature(Natom,Mensemble,mstep,ntmp,simid,emomM,beff,0)
      end if

      if (do_pol=='Y') then
         call init_polarization(Natom,Mensemble,ham%max_no_neigh,ham%nlist,coord,1)
      end if

      ! Calculate the external static fields, the can be calculated once and it does not needs to be calculated again
      call calc_external_fields(Natom,Mensemble,hfield,anumb,external_field,     &
         do_bpulse,sitefld,sitenatomfld)
      !
      ! Adaptive Time Stepping Region
      if(adaptive_time_flag) then
         call calculate_omegainit(omega_max, larmor_numrev, delta_t)
      end if

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

      ! --Optimization Region-- !
      ! Allocation and initial opt build

      if (OPT_ON) then
         call timing(0,'Hamiltonian   ','OF')
         call timing(0,'BuildOptReg   ','ON')
         call buildOptimizationRegions(na,natom,nHam, Mensemble,ham%max_no_neigh,   &
            atype,emom,mmom,ham%ncoup,ham%nlist,ham%nlistsize,ham%aham,             &
            cellPosNumNeigh,cos_thr)
         call timing(0,'BuildOptReg   ','OF')
         call timing(0,'Hamiltonian   ','ON')
      end if

      ! Initialize cumulant counter
      Navrgcum=0
      ! Perform nstep simulation steps
      ! Observe that Fortran cannot handle a regular do loop if nstep changes dynamically
      ! Hence, a do while loop has been introduced to allow for adaptive time stepping
      mstep=rstep+1
      !------------------------------------------------------------------------------
      ! If KMC is turned on calculate the initial rate
      !------------------------------------------------------------------------------
      if (do_kmc=='Y') then
         ! Calling the wrapper KMC routine for all the different methods
         call initial_wrapper_kmc(mstep,Natom,NA_KMC,do_ralloy,Mensemble,kmc_method,&
            rate0,delta_t,Temp_array,coord,atype,kmc_index_prv,kmc_index_aft,       &
            kmc_time_steps,time_efield)
      endif

      !------------------------------------------------------------------------------
      ! End of initial KMC wrapper
      !------------------------------------------------------------------------------

      do while (mstep.LE.rstep+nstep) !+1

         if (time_dept_flag) then
            ! Calculate Microwave fields (time dependent fields)
            call calculate_mwf_fields(Natom,mstep,rstep+nstep-1,delta_t,coord,0)

            ! Calculate the total time dependent fields
            call calc_external_time_fields(Natom, Mensemble, time_external_field,   &
               do_bpulse, demag, mwf, mwf_gauss_spatial, do_gauss, mwf_gauss,       &
               mov_gauss, mwf_mov_gauss, bpulsefield, demagfld, mwffield,           &
               gauss_mwffield,site_mwffield,gauss_spatial_site_mwffield,            &
               gauss_spatial,gauss_site_mwffield,mov_gauss_spatial,                 &
               mwf_mov_gauss_spatial,mov_circle,mwf_mov_circle,mov_circle_spatial,  &
               mwf_mov_circle_spatial,mov_square,mwf_mov_square,mov_square_spatial, &
               mwf_mov_square_spatial)
         else
            time_external_field=0.0_dblprec
         endif

         ! Measure averages and trajectories
         call measure(Natom,Mensemble,NT,NA,nHam,N1,N2,N3,simid,mstep,emom,emomM,   &
            mmom,Nchmax,do_ralloy,Natom_full,asite_ch,achem_ch,atype,plotenergy,    &
            Temp,temprescale,temprescalegrad,real_time_measure,delta_t,logsamp,     &
            ham%max_no_neigh,                                                       &
            ham%nlist,ham%ncoup,ham%nlistsize,ham%aham,thermal_field,beff,beff1,    &
            beff3,coord,ham%ind_list_full,ham%ind_nlistsize,ham%ind_nlist,          &
            ham%max_no_neigh_ind,ham%sus_ind,do_mom_legacy,mode)

         ! Calculate total and term resolved energies
         if(plotenergy>0) then
            call timing(0,'Measurement   ','OF')
            call timing(0,'Energy        ','ON')
            if (mod(mstep-1,avrg_step)==0) then
               call calc_energy(nHam,mstep,Natom,Nchmax, &
                  conf_num,Mensemble,Natom,Num_macro,1,         &
                  plotenergy,Temp,delta_t,do_lsf,        &
                  lsf_field,lsf_interpolate,real_time_measure,simid,cell_index,            &
                  macro_nlistsize,mmom,emom,emomM,emomM_macro,external_field,              &
                  time_external_field,max_no_constellations,maxNoConstl,                   &
                  unitCellType,constlNCoup,constellations,OPT_flag,                        &
                  constellationsNeighType,totene,NA,N1,N2,N3)
            end if
            call timing(0,'Energy        ','OF')
            call timing(0,'Measurement   ','ON')

         endif

         call timing(0,'Measurement   ','OF')
         call timing(0,'SpinCorr      ','ON')

         ! Spin correlation
         ! Sample magnetic moments for correlation functions
         call correlation_wrapper(Natom,Mensemble,coord,simid,emomM,mstep,delta_t,  &
            NT_meta,atype_meta,Nchmax,achtype,sc,do_sc,do_sr,cgk_flag)

         call calc_bls(N1,N2,N3,C1,C2,C3,Natom,Mensemble,simid,coord,emomM,mstep,   &
            delta_t,1)

         call timing(0,'SpinCorr      ','OF')
         call timing(0,'Measurement   ','ON')

         ! Write simulation status for each 5% of the simulation length
         if(nstep>20) then
            if(mod(mstep,(rstep+nstep)/20)==0) then
               call calc_mavrg(Natom,Mensemble,emomM,mavg)
               ! Only calculate binderc if it is not calculated elsewhere
               ! If calculated elsewhere, display binderc will not be real binderc
               if(do_cumu=='N') then
                  call calc_and_print_cumulant(Natom,Mensemble,emomM,simid,Temp,    &
                     temprescale,temprescalegrad,plotenergy,cumu_buff,.false.)
               endif
               write(*,'(2x,a,i4,a,f10.6)',advance='no') &
                  "MP",mstep*100/(rstep+nstep),"% done. Mbar:",mavg
               if(plotenergy>0) then
                  write(*,'(a,f12.6)',advance='no') ". Ebar:", totene
               endif
               sstep = f_logstep(mstep,logsamp)
               if(mod(sstep,cumu_step)==0)then
                  write(*,'(a,f8.5,a)',advance='no') ". U:",binderc,"."
               endif
               write(*,*) ''
            endif
         else
            write(*,'(2x,a,i3,a,G13.6)')   "Iteration",mstep," Mbar ",mavg
         endif

         !Adjust QHB
         if(do_qhb=='Q' .or. do_qhb=='R' .or. do_qhb=='P' .or. do_qhb=='T') then
            if(qhb_mode=='MT') then
               if (mod(mstep,(rstep+nstep)/40)==0) then
                  call calc_mavrg(Natom,Mensemble,emomM,mavg)
                  call qhb_rescale(Temp,temprescale,temprescalegrad,do_qhb,qhb_mode,mavg)
               endif
            endif
         endif
         call timing(0,'Measurement   ','OF')
         call timing(0,'Hamiltonian   ','ON')

         ! Calculate spin transfer torque contributions to the local field
         if(stt=='A'.or.stt=='F'.or.do_she=='Y'.or.do_sot=='Y') then
            call calculate_spintorques(Natom, Mensemble,lambda1_array,emom,mmom)
         end if

         ! --Optimization Region-- !
         if (OPT_ON) then
            if (mod(mstep-rstep,OPT_rebuild_time)==0) then
               call timing(0,'Hamiltonian   ','OF')
               call timing(0,'BuildOptReg   ','ON')
               call buildOptimizationRegions(na,natom,nHam, Mensemble,              &
                  ham%max_no_neigh,atype,emom,mmom,ham%ncoup,ham%nlist,             &
                  ham%nlistsize,ham%aham,cellPosNumNeigh,cos_thr)
               call timing(0,'BuildOptReg   ','OF')
               call timing(0,'Hamiltonian   ','ON')
            end if
         end if

         ! Apply Hamiltonian to obtain effective field
         if(do_sparse=='Y') then
            if(ham_inp%do_dm==1.or.ham_inp%do_jtensor==1) then
               call effective_field_SparseBlock(Natom,Mensemble,emomM,beff)
            else
               call effective_field_SparseScalar(Natom,Mensemble,emomM,beff)
            end if
            beff1=beff
            beff2=external_field+time_external_field
            beff=beff1+beff2
         else

            call effective_field_gpu(Natom,Mensemble,1,Natom,emomM,   &
               mmom,external_field,time_external_field,beff,beff1,beff2,OPT_flag,   &
               max_no_constellations, maxNoConstl,unitCellType, constlNCoup,        &
               constellations, constellationsNeighType, totenergy,        &
               Num_macro,cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)
         end if

         call timing(0,'Hamiltonian   ','OF')
         call timing(0,'Evolution     ','ON')

            ! Calculate average field for use with multiple heat baths
            if (llg==2.or.llg==3) then
               call calc_averagefield(Natom,Mensemble,beff1,beff2,field1,field2)
            endif
            ! Check if this changes
            thermal_field=0.0_dblprec
            call evolve_first_gpu(Natom,Mensemble,Landeg,llg,SDEalgh,bn,lambda1_array,  &
               lambda2_array,NA,compensate_drift,delta_t,relaxtime,Temp_array,      &
               temprescale,beff,b2eff,thermal_field,beff2,btorque,field1,field2,    &
               emom,emom2,emomM,mmom,mmomi,stt,do_site_damping,ham%nlist,           &
               ham%nlistsize,constellationsUnitVec,constellationsUnitVec2,          &
               constellationsMag,constellations,unitCellType,OPT_flag,cos_thr,      &
               max_no_constellations,do_she,she_btorque,Nred,red_atom_list,         &
               do_sot,sot_btorque)

         call timing(0,'Evolution     ','OF')
         call timing(0,'Hamiltonian   ','ON')

         if(SDEalgh==1.or.SDEalgh==4.or.SDEalgh==5.or.SDEalgh==6.or.SDEalgh==7.or.SDEalgh==11) then

            if (time_dept_flag) then
               ! Calculate Microwave fields (time dependent fields)
               call calculate_mwf_fields(Natom,mstep,rstep+nstep-1,delta_t,coord,1)
               ! Calculate the total time dependent fields
               call calc_external_time_fields(Natom,Mensemble,time_external_field,  &
                  do_bpulse,demag,mwf,mwf_gauss_spatial,do_gauss,mwf_gauss,         &
                  mov_gauss,mwf_mov_gauss,bpulsefield,demagfld,mwffield,            &
                  gauss_mwffield,site_mwffield,gauss_spatial_site_mwffield,         &
                  gauss_spatial,gauss_site_mwffield,mov_gauss_spatial,              &
                  mwf_mov_gauss_spatial,mov_circle,mwf_mov_circle,                  &
                  mov_circle_spatial,mwf_mov_circle_spatial,mov_square,             &
                  mwf_mov_square,mov_square_spatial,mwf_mov_square_spatial)
            else
               time_external_field=0.0_dblprec
            endif

            ! Apply Hamiltonian to obtain effective field
            if(do_sparse=='Y') then
               if(ham_inp%do_dm==1.or.ham_inp%do_jtensor==1) then
                  call effective_field_SparseBlock(Natom,Mensemble,emomM,beff)
               else
                  call effective_field_SparseScalar(Natom,Mensemble,emomM,beff)
               end if
               beff1=beff
               beff2=external_field+time_external_field
               beff=beff1+beff2
            else
               !---------------------------------------------------------------------
               !! End of the induced moments treatment
               !---------------------------------------------------------------------
               call effective_field_gpu(Natom,Mensemble,1,Natom,  &
                  emomM,mmom,external_field,time_external_field,     &
                  beff,beff1,beff2,OPT_flag,max_no_constellations,maxNoConstl,      &
                  unitCellType,constlNCoup,constellations,constellationsNeighType,  & 
                  totenergy,Num_macro,cell_index,emomM_macro,             &
                  macro_nlistsize,NA,N1,N2,N3)
            end if

         end if
         call timing(0,'Hamiltonian   ','OF')
         call timing(0,'Evolution     ','ON')

         ! Re-Calculate spin transfer torque
         if(stt=='A'.or.stt=='F'.or.do_she=='Y'.or.do_sot=='Y') then
            call calculate_spintorques(Natom, Mensemble,lambda1_array,emom2,mmom)
         end if

         ! Perform second (corrector) step of SDE solver
         call evolve_second_gpu(Natom,Mensemble,Landeg,llg,SDEalgh,bn,lambda1_array,    &
            delta_t,relaxtime,beff,beff2,b2eff,btorque,emom,emom2,stt,              &
            ham%nlist,ham%nlistsize,constellationsUnitVec,constellationsUnitVec2,   &
            constellationsMag,constellations,unitCellType,OPT_flag,cos_thr,         &
            max_no_constellations,do_she,she_btorque,Nred,red_atom_list,            &
            do_sot,sot_btorque)

         call timing(0,'Evolution     ','OF')
         call timing(0,'Moments       ','ON')

         ! Update magnetic moments after time evolution step
         call moment_update(Natom,Mensemble,mmom,mmom0,mmom2,emom,emom2,emomM,mmomi,&
            mompar,initexc)

         call timing(0,'Moments       ','OF')
         call timing(0,'Measurement   ','ON')

         ! Update magnetic field pulse
         if (do_bpulse.gt.0.and.do_bpulse.lt.5) then
            if (scount_pulse == bpulse_step) then
               bpulse_time = delta_t*mstep
               call bpulse(do_bpulse, bpulse_time)
               scount_pulse=1
            else
               scount_pulse=scount_pulse+1
            end if
         end if
         ! Sample autocorrelation
         if(do_autocorr=='Y') then
            call autocorr_sample(Natom, Mensemble, simid, mstep, emom)
         end if

         ! Sample spin temperature
         if (do_spintemp=='Y') then
            if(mod(mstep,spintemp_step)==0) then
               call spintemperature(Natom,Mensemble,mstep,1,simid,emomM,beff,1)
            end if
         endif

         mstep = mstep + 1
         !----------------------------------------------------------------------------
         ! If the macrospin dipole-dipole interaction is used one needs to re-calculate
         ! the macrospins
         !----------------------------------------------------------------------------
         if (ham_inp%do_dip==2) then
            call calc_macro_mom(Natom,Num_macro,Mensemble,max_num_atom_macro_cell,  &
               macro_nlistsize,macro_atom_nlist,mmom,emom,emomM,mmom_macro,         &
               emom_macro,emomM_macro)
         endif
         !----------------------------------------------------------------------------
         ! End of the macrospin re-calculation
         !----------------------------------------------------------------------------

         !----------------------------------------------------------------------------
         ! Start of KMC wrapper
         !----------------------------------------------------------------------------
         if (do_kmc=='Y') then
            call updater_wrapper_kmc(mstep,Natom,NA_KMC,do_ralloy,Mensemble,        &
               Natom_full,kmc_method,ham%max_no_neigh,ham%nlistsize,ham%nlist,rate0,&
               delta_t,Temp_array,coord,atype,kmc_index_prv,kmc_index_aft,          &
               kmc_time_steps,achem_ch,mmom,ham%ncoup,emom,emomM,time_efield)
         endif
         !----------------------------------------------------------------------------
         ! End of KMC wrapper
         !----------------------------------------------------------------------------

         !----------------------------------------------------------------------------
         ! Induced moments treatment
         !----------------------------------------------------------------------------
         ! If there are induced moments in the system one must renormalize the
         ! the Heisenberg exchange and Dzyaloshinskii-Moriya vectors, as well as
         ! the magnitude and direction of the induced moments
         !if (ind_mom_flag=='Y') then
         !   call renorm_ncoup_ind(do_dm,Natom,conf_num,Mensemble,ham%max_no_neigh,  &
         !      ham%max_no_dmneigh,ham%max_no_neigh_ind,ham%nlistsize,ham%dmlistsize,&
         !      ham%ind_list_full,ham%ind_nlistsize,ham%nlist,ham%dmlist,            &
         !      ham%ind_nlist,ham%sus_ind,mmom,emom,emomM,ham%ncoup,ham%dm_vect,     &
         !      ham%fix_list,ham%fix_num)
         !endif
         !----------------------------------------------------------------------------
         ! End of the induced moments treatment
         !----------------------------------------------------------------------------

      end do ! End loop over simulation steps

      ! Measure averages and trajectories
      call measure(Natom,Mensemble,NT,NA,nHam,N1,N2,N3,simid,mstep,emom,emomM,mmom, &
         Nchmax,do_ralloy,Natom_full,asite_ch,achem_ch,atype,plotenergy,Temp,temprescale,&
         temprescalegrad,real_time_measure,delta_t,logsamp,ham%max_no_neigh,ham%nlist,ham%ncoup,&
         ham%nlistsize,ham%aham,thermal_field,beff,beff1,beff3,coord,               &
         ham%ind_list_full,ham%ind_nlistsize,ham%ind_nlist,ham%max_no_neigh_ind,    &
         ham%sus_ind,do_mom_legacy,mode)

      ! Print remaining measurements
      call flush_measurements(Natom,Mensemble,NT,NA,N1,N2,N3,simid,mstep,emom,mmom, &
         Nchmax,atype,real_time_measure,rstep+nstep,ham%ind_list_full,do_mom_legacy,&
         mode)

      ! Print final polarization, chirality and local polarization
      if (do_pol=='Y') then
         call init_polarization(Natom,Mensemble,ham%max_no_neigh,ham%nlist,coord,-1)
      end if

      if (do_spintemp=='Y') then
         call spintemperature(Natom,Mensemble,mstep,1,simid,emomM,beff,2)
      endif
      call timing(0,'Measurement   ','OF')
   end subroutine sd_mphase_gpu

end module sd_driver_gpu
