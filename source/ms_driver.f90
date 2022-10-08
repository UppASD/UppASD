!******************************************************************
!*                     *** UppASD ***                             *
!*                                                                *
!*               Uppsala Atomic Spin Dynamics                     *
!*                                                                *
!*                   Version 5.0 Mar 2020                         *
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
module ms_driver

   implicit none

   public

contains

   !---------------------------------------------------------------------------------
   !> @brief
   !> Multiscale Spin Dynamics initial phase
   !
   !> @author
   !> !> authors
   !> Edgar Mendez
   !> Nikos  Ntallis
   !> Manuel Pereiro
   !---------------------------------------------------------------------------------
   subroutine ms_iphase()
      !
      use KMC
      use QHB,                   only : qhb_rescale
      use BLS,                   only : calc_bls
      use Energy,                only : calc_energy
      use Sparse
      use Damping
      use KMCData
      use FixedMom
      use Gradients
      use Evolution_ms
      use InputData
      use FieldData,             only : beff,beff1,beff2,b2eff,sitefld,       &
         external_field,field1,field2,time_external_field,allocation_field_time,    &
         thermal_field
      use macrocells
      use DemagField
      use MomentData
      use FieldPulse
      use SystemData,            only : atype, anumb, Landeg
      use Correlation,           only : correlation_wrapper
      use Temperature
      use SpinTorques
      use ChemicalData
      !use prn_averages
      use Measurements
      use Polarization
      use UpdateMoments
      !use InducedMoments,        only : renorm_ncoup_ind
      use MicroWaveField
      use SimulationData,        only : bn
      use HamiltonianData
      use CalculateFields
      use AutoCorrelation,       only : autocorr_sample
      use prn_trajectories
      use HamiltonianActions
      use OptimizationRoutines
      use AdaptiveTimeStepping
      use MultiscaleDampingBand
      use MultiscaleInterpolation



      !
      implicit none

      ! Adaptive time stepping
      integer :: i, iprstep, adapt_step, ipstep
      logical :: deltat_correction_flag, sd_phaseflag
      real(dblprec) :: mavg
      real(dblprec) :: temprescale, temprescalegrad, totenergy
      ! Adaptive time stepping with spin correlation

      ! No spin correlation in the initial phase
      adapt_step = 0
      deltat_correction_flag = .true.
      ! Measurement phase indicator (true here)
      sd_phaseflag = .true.

      ! Loop over initial phases
      do i=1,ipnphase

         ! Adaptive Time Stepping Region
            if(ip_adaptive_time_flag) then
            call calculate_omegainit(omega_max, larmor_numrev, ipdelta_t(i))
           end if

         ! Write output to stdout
            if (do_site_ip_damping=='Y') then
            write (*,'(2x,a25,i3,a10,G11.4,a10,i10,a10,G14.4)') &
               "Performing initial phase: ", i, " Temp. ", ipTemp(i), "nsteps ",    &
               ipnstep(i), "dt ", ipdelta_t(i)
            else
            write (*,'(2x,a25,i3,a10,G11.4,a10,i10,a10,G14.4,a10,2G14.4)') &
               "Performing initial phase: ", i, " Temp. ", ipTemp(i), "nsteps ",    &
               ipnstep(i), "dt ", ipdelta_t(i), "lambda ", iplambda1(i), iplambda2(i)
           endif

         ! Calculate demagnetization field
         

         ! Rescaling of temperature according to Quantum Heat bath
         temprescale=1.0_dblprec
         temprescalegrad=0.0_dblprec
         
         ! Calculate external fields
           call calc_external_fields(Natom,Mensemble,iphfield,anumb,external_field,&
            do_bpulse,sitefld,sitenatomfld)


         ! --Optimization Region-- !
         ! Allocation and initial opt build
            if(OPT_ON) then
               call timing(0,'Initial       ','OF')
               call timing(0,'BuildOptReg   ','ON')
               call buildOptimizationRegions(na,natom,nHam,Mensemble,ham%max_no_neigh, &
               atype,emom,mmom,ham%ncoup,ham%nlist,ham%nlistsize,ham%aham,          &
               cellPosNumNeigh,cos_thr)
               call timing(0,'BuildOptReg   ','OF')
               call timing(0,'Initial       ','ON')
           end if

         ! Adaptive time stepping
           iprstep = 0
           ipstep=1

         ! Main initial loop
         !--------------------------------------!
         do while(ipstep.LE.ipnstep(i))

            ! Write output to stdout
                  if(ipnstep(i)<10.or.mod(ipstep,ipnstep(i)/10)==0) then
                    call calc_mavrg(Natom,Mensemble,emomM,mavg)
                    write(*,'(2x,a,i2,i5,a,f10.6)') "IP ",i,ipstep*100/ipnstep(i),"% done. Mbar:",mavg
                  end if
            !Adjust QHB
            
            ! --Optimization Region-- !
                  if(OPT_ON) then
                     if (mod(ipstep,OPT_rebuild_time)==0) then
                     call timing(0,'Initial       ','OF')
                     call timing(0,'BuildOptReg   ','ON')
                     call buildOptimizationRegions(na,natom,nHam,Mensemble,            &
                          ham%max_no_neigh,atype,emom,mmom,ham%ncoup,ham%nlist,          &
                          ham%nlistsize,ham%aham, cellPosNumNeigh,cos_thr)
                     call timing(0,'BuildOptReg   ','OF')
                     call timing(0,'Initial       ','ON')
                    end if
                  end if

            ! Apply Hamiltonian to obtain effective field
            call timing(0,'Initial       ','OF')
            call timing(0,'Hamiltonian   ','ON')
            
              call effective_field(Natom,Mensemble,1,Natom,emomM,   &
               mmom,external_field,time_external_field,beff,beff1,beff2,OPT_flag,   &
               max_no_constellations, maxNoConstl,unitCellType, constlNCoup,        &
               constellations, constellationsNeighType, totenergy,        &
               Num_macro,cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)
            
            call timing(0,'Hamiltonian   ','OF')
            call timing(0,'Initial       ','ON')

            ! Calculate average field for use with multiple heat baths
            !if(llg>=2.and.llg<=3) call calc_averagefield(Natom,Mensemble,beff1,beff2,field1,field2)

            ! Try to see if this fixes the energy bouncing around. It did. This is need to avoid incremental temp.
            !thermal_field=0.0_dblprec

             call dampingBandPreinterpolation(dampingBand,emom)

               ! Perform first (predictor) step of SDE solver
               call timing(0,'Initial       ','OF')
               call timing(0,'Evolution     ','ON')

              ! call evolve_first_ms(Natom,Mensemble,Landeg,llg,ipSDEalgh,bn,             &
              !    iplambda1_array(i,:),iplambda2_array(i,:),NA,compensate_drift,    &
              !    ipdelta_t(i),relaxtime,ipTemp_array(:,i),temprescale,beff,b2eff,  &
              !    thermal_field,beff2,btorque,field1,field2,emom,emom2,emomM,mmom,  &
              !    mmomi,'N',do_site_ip_damping,ham%nlist,ham%nlistsize,             &
              !    constellationsUnitVec,constellationsUnitVec2,constellationsMag,   &
              !    constellations,unitCellType,OPT_flag,cos_thr,                     &
              !    max_no_constellations,'N',she_btorque,Nred,red_atom_list,         &
              !    do_fixed_mom,'N',sot_btorque,dampingBand)

               call evolve_first_ms(Natom,Mensemble,Landeg,llg,SDEalgh,bn,iplambda1_array(i,:),  &
               iplambda2_array(i,:),NA,compensate_drift,ipdelta_t(i),relaxtime,            & 
               ipTemp_array(:,i),temprescale,beff,b2eff,thermal_field,beff2,btorque,field1,field2,    &
               emom,emom2,emomM,mmom,mmomi,stt,do_site_damping,ham%nlist,           &
               ham%nlistsize,constellationsUnitVec,constellationsUnitVec2,          &
               constellationsMag,constellations,unitCellType,OPT_flag,cos_thr,      &
               max_no_constellations,do_she,she_btorque,Nred,red_atom_list,         &
               do_fixed_mom,do_sot,sot_btorque,dampingBand)

               !print*, thermal_field

               call dampingBandCorrPreinterpolation(dampingBand,emom2)
                 ! if (do_multiscale) then            
                   call multiscaleInterpolateInterfaces()
                      if(SDEalgh==1) then
                      call dampingBandCorrPreinterpolationAvg(dampingBand,emom2)
                      elseif(SDEalgh==5) then
                      call dampingBandCorrPreinterpolation(dampingBand,emom)
                     end if
                ! end if

               

                  ! Apply Hamiltonian to obtain effective field
                  call timing(0,'Evolution     ','OF')
                  call timing(0,'Hamiltonian   ','ON')
                  

                     call effective_field(Natom,Mensemble,1,Natom,emomM,   &
               mmom,external_field,time_external_field,beff,beff1,beff2,OPT_flag,   &
               max_no_constellations, maxNoConstl,unitCellType, constlNCoup,        &
               constellations, constellationsNeighType, totenergy,        &
               Num_macro,cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)
                 
                  call timing(0,'Hamiltonian   ','OF')
                  call timing(0,'Evolution     ','ON')
              

               ! Perform second (corrector) step of SDE solver
              ! call evolve_second_ms(Natom,Mensemble,Landeg,llg,ipSDEalgh,bn,            &
              !    iplambda1_array(i,:),ipdelta_t(i),relaxtime,beff,beff2,b2eff,     &
              !    btorque,emom,emom2,'N',ham%nlist,ham%nlistsize,                   &
              !    constellationsUnitVec,constellationsUnitVec2,constellationsMag,   &
              !    constellations,unitCellType,OPT_flag,cos_thr,                     &
              !    max_no_constellations,'N',she_btorque,Nred,red_atom_list,         &
              !    do_fixed_mom,'N',sot_btorque,dampingBand)

               call evolve_second_ms(Natom,Mensemble,Landeg,llg,SDEalgh,bn,iplambda1_array(i,:),    &
            ipdelta_t(i),relaxtime,beff,beff2,b2eff,btorque,emom,emom2,stt,              &
            ham%nlist,ham%nlistsize,constellationsUnitVec,constellationsUnitVec2,   &
            constellationsMag,constellations,unitCellType,OPT_flag,cos_thr,         &
            max_no_constellations,do_she,she_btorque,Nred,red_atom_list,            &
            do_fixed_mom,do_sot,sot_btorque,dampingBand)

             !if (do_multiscale) then
                call multiscaleInterpolateInterfaces()
            ! end if


            ! Update magnetic moments after time evolution step
            call moment_update(Natom,Mensemble,mmom,mmom0,mmom2,emom,emom2,emomM,   &
               mmomi,mompar,initexc)

           
            
            ipstep = ipstep + 1
            !------------------------------------------------------------------------
            !------------------------------------------------------------------------
            ! If the macrospin dipole-dipole interaction is used one needs to re-calculate
            ! the macrospins
            !------------------------------------------------------------------------
            
            !------------------------------------------------------------------------
            ! End of the macrospin re-calculation
            !------------------------------------------------------------------------
            call timing(0,'Evolution     ','OF')
            call timing(0,'Initial       ','ON')

         enddo
      enddo

   end subroutine ms_iphase

   !---------------------------------------------------------------------------
   !> @brief
   !> Multiscale Spin Dynamics measurement phase
   !
   !> authors
   !> Edgar Mendez
   !> Nikos  Ntallis
   !> Manuel Pereiro
   !---------------------------------------------------------------------------
   subroutine ms_mphase()
      !
      use KMC
      use QHB,                   only : qhb_rescale
      use BLS,                   only : calc_bls
      use Energy,                only : calc_energy
      use Sparse
      use Damping
      use KMCData
      use FixedMom
      use Gradients
      use Evolution_ms
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
      use Correlation,           only : correlation_wrapper
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
      use OptimizationRoutines
      use AdaptiveTimeStepping
      use MultiscaleDampingBand
      use MultiscaleInterpolation

      implicit none
      logical :: time_dept_flag, deltat_correction_flag
      integer :: cgk_flag, cgk_flag_p, adapt_step, cr_flag, spt_flag
      integer :: bcgk_flag,cgk_flag_pc
      integer :: scount_pulse, sstep
      ! Phase flag indicator (true for sd initial phase; false for sd measurement phase)
      ! Used in order to separate between phases in an adaptive time step environment with spin correlation
      logical :: sd_phaseflag

      real(dblprec) :: temprescale, temprescalegrad, totene, totenergy

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
      

      call calc_bls(N1,N2,N3,C1,C2,C3,Natom,Mensemble,simid,coord,emomM,sstep,      &
         delta_t,0)

      

      !  Iniitialize magnetic field pulse
      if (do_bpulse.gt.0.and.do_bpulse.lt.5) then
         bpulse_time = delta_t*rstep
         call bpulse(do_bpulse, bpulse_time)
      end if

      ! Initialize polarization including local chirality and polarity
      

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
         !call correlation_wrapper(Natom,Mensemble,coord,simid,emomM,mstep,delta_t,  &
         !   NT,atype,Nchmax,achtype,cgk_flag,cgk_flag_p)

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
         

             call effective_field(Natom,Mensemble,1,Natom,emomM,   &
               mmom,external_field,time_external_field,beff,beff1,beff2,OPT_flag,   &
               max_no_constellations, maxNoConstl,unitCellType, constlNCoup,        &
               constellations, constellationsNeighType, totenergy,        &
               Num_macro,cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)
        
            call dampingBandPreinterpolation(dampingBand,emom)

            call timing(0,'Hamiltonian   ','OF')
            call timing(0,'Evolution     ','ON')

            ! Calculate average field for use with multiple heat baths
           ! if (llg==2.or.llg==3) then
           !    call calc_averagefield(Natom,Mensemble,beff1,beff2,field1,field2)
           ! endif
            ! Check if this changes

            thermal_field=0.0_dblprec

            call evolve_first_ms(Natom,Mensemble,Landeg,llg,SDEalgh,bn,lambda1_array,  &
               lambda2_array,NA,compensate_drift,delta_t,relaxtime,Temp_array,      &
               temprescale,beff,b2eff,thermal_field,beff2,btorque,field1,field2,    &
               emom,emom2,emomM,mmom,mmomi,stt,do_site_damping,ham%nlist,           &
               ham%nlistsize,constellationsUnitVec,constellationsUnitVec2,          &
               constellationsMag,constellations,unitCellType,OPT_flag,cos_thr,      &
               max_no_constellations,do_she,she_btorque,Nred,red_atom_list,         &
               do_fixed_mom,do_sot,sot_btorque,dampingBand)

                     
             call multiscaleInterpolateInterfaces()
             if(SDEalgh==1) then
                call dampingBandCorrPreinterpolationAvg(dampingBand,emom2)
             elseif(SDEalgh==5) then
                call dampingBandCorrPreinterpolation(dampingBand,emom)
             end if
          


         call timing(0,'Evolution     ','OF')
         call timing(0,'Hamiltonian   ','ON')

        

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
            
               !---------------------------------------------------------------------
               !! End of the induced moments treatment
               !---------------------------------------------------------------------
               call effective_field(Natom,Mensemble,1,Natom,emomM,   &
               mmom,external_field,time_external_field,beff,beff1,beff2,OPT_flag,   &
               max_no_constellations, maxNoConstl,unitCellType, constlNCoup,        &
               constellations, constellationsNeighType, totenergy,        &
               Num_macro,cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)
            
           
         
         call timing(0,'Hamiltonian   ','OF')
         call timing(0,'Evolution     ','ON')

         ! Re-Calculate spin transfer torque
         if(stt=='A'.or.stt=='F'.or.do_she=='Y'.or.do_sot=='Y') then
            call calculate_spintorques(Natom, Mensemble,lambda1_array,emom2,mmom)
         end if

         ! Perform second (corrector) step of SDE solver
         call evolve_second_ms(Natom,Mensemble,Landeg,llg,SDEalgh,bn,lambda1_array,    &
            delta_t,relaxtime,beff,beff2,b2eff,btorque,emom,emom2,stt,              &
            ham%nlist,ham%nlistsize,constellationsUnitVec,constellationsUnitVec2,   &
            constellationsMag,constellations,unitCellType,OPT_flag,cos_thr,         &
            max_no_constellations,do_she,she_btorque,Nred,red_atom_list,            &
            do_fixed_mom,do_sot,sot_btorque,dampingBand)

              call multiscaleInterpolateInterfaces()
         


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
         

         mstep = mstep + 1
         

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
      
      call timing(0,'Measurement   ','OF')
   end subroutine ms_mphase

   

end module ms_driver
