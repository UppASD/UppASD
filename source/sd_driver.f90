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
module sd_driver

   implicit none

   public

contains



   !---------------------------------------------------------------------------
   !> @brief
   !> Spin Dynamics initial phase
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine sd_iphase()
      !
      use Damping,              only : iplambda1_array,iplambda2_array
      !      use Ewaldmom
      use Evolution
      use InputData
      use FieldData,            only : beff, beff1,beff2, b2eff, sitefld, external_field, &
         field1, field2, thermal_field, time_external_field
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
      use AdaptiveTimestepping
      use HamiltonianActions
      use SpinTorques
      use QHB, only : qhb_rescale, do_qhb, qhb_mode
      !
      implicit none

      ! Implicit midpoint methods
      integer :: imp_count
      logical :: converged
      integer :: j, k

      ! Adaptive time stepping
      integer :: i, iprstep, adapt_step, ipstep
      logical :: deltat_correction_flag, sd_phaseflag
      real(dblprec) :: temprescale,dummy, denergy
      real(dblprec) :: mavg

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
               "Performing initial phase: ", i, " Temp. ", ipTemp(i), "nsteps ",&
               ipnstep(i), "dt ", ipdelta_t(i)
         else
            write (*,'(2x,a25,i3,a10,G11.4,a10,i10,a10,G14.4,a10,2G14.4)') &
               "Performing initial phase: ", i, " Temp. ", ipTemp(i), "nsteps ",&
               ipnstep(i), "dt ", ipdelta_t(i), "lambda ", iplambda1(i), iplambda2(i)
         endif

         ! Calculate demagnetization field
         if(demag=='Y') then
            call calc_demag(Natom, Mensemble, demag1, demag2, demag3, demagvol, emomM)
         endif

         ! Rescaling of temperature according to Quantum Heat bath
         temprescale=1.d0
         if (do_qhb=="Q" .or. do_qhb=='R' .or. do_qhb=='T') then
            if(qhb_mode=='MT') then
               call calc_mavrg(Natom,Mensemble,emomM,mavg)
               call qhb_rescale(iptemp(i),temprescale,do_qhb,qhb_mode,mavg)
            else
               call qhb_rescale(iptemp(i),temprescale,do_qhb,qhb_mode,dummy)
            endif
         endif

         call calc_external_fields(Natom, Mensemble, NA, iphfield, anumb, external_field, &
            do_bpulse,sitefld,sitenatomfld)


         ! --Optimization Region-- !
         ! Allocation and initial opt build
         if(OPT_ON) then
            call timing(0,'Initial       ','OF')
            call timing(0,'BuildOptReg   ','ON')
            call buildOptimizationRegions(na,natom,Mensemble,max_no_neigh,atype,&
               emom,mmom,ncoup,nlist,nlistsize, cellPosNumNeigh,cos_thr)
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
            if(do_qhb=='Q' .or. do_qhb=='R' .or. do_qhb=='T') then
               if(qhb_mode=='MT') then
                  if (mod(ipstep,ipnstep(i)/20)==0) then
                     call calc_mavrg(Natom,Mensemble,emomM,mavg)
                     call qhb_rescale(iptemp(i),temprescale,do_qhb,qhb_mode,mavg)
                  endif
               endif
            endif
            ! --Optimization Region-- !
            if(OPT_ON) then
               if (mod(ipstep,OPT_rebuild_time)==0) then
                  call timing(0,'Initial       ','OF')
                  call timing(0,'BuildOptReg   ','ON')
                  call buildOptimizationRegions(na,natom,Mensemble,max_no_neigh,atype, &
                     emom,mmom,ncoup,nlist,nlistsize,cellPosNumNeigh,cos_thr)
                  call timing(0,'BuildOptReg   ','OF')
                  call timing(0,'Initial       ','ON')
               end if
            end if

            ! Apply Hamiltonian to obtain effective field
            call timing(0,'Initial       ','OF')
            call timing(0,'Hamiltonian   ','ON')
            call effective_field(Natom, Mensemble, 1, Natom, &
               do_jtensor, exc_inter, do_dm, do_pd, do_biqdm, do_bq, &
               do_dip,emomM, mmom, external_field,time_external_field, beff, beff1, beff2, &
               OPT_flag, max_no_constellations, maxNoConstl, &
               unitCellType, constlNCoup, constellations, constellationsNeighType, mult_axis,denergy)
            call timing(0,'Hamiltonian   ','OF')
            call timing(0,'Initial       ','ON')

            ! Calculate average field for use with multiple heat baths
            if(llg>=2.and.llg<=3) call calc_averagefield(Natom,Mensemble,beff1,beff2,field1,field2)

            ! Try to see if this fixes the energy bouncing around. It did. This is need to avoid incremental temp.
            thermal_field=0.0D0

            ! One- and two-step solvers
            if(SDEalgh<=20) then

               ! Perform first (predictor) step of SDE solver
               call timing(0,'Initial       ','OF')
               call timing(0,'Evolution     ','ON')
               call evolve_first(Natom, Mensemble, Landeg, llg, SDEalgh, bn, iplambda1_array(i,:), iplambda2_array(i,:), NA, &
                  compensate_drift, ipdelta_t(i), relaxtime,ipTemp_array(:,i),temprescale, beff, b2eff, thermal_field,beff2, &
                  btorque, field1, field2, emom, emom2, emomM, mmom, mmomi ,'N', do_site_ip_damping,&
                  nlist,nlistsize,constellationsUnitVec,constellationsUnitVec2,constellationsMag, &
                  constellations,unitCellType,OPT_flag,cos_thr,max_no_constellations,'N',she_btorque)

               if(SDEalgh==1.or.SDEalgh==4.or.SDEalgh==5.or.SDEalgh==6.or.SDEalgh==7.or.SDEalgh==11) then

                  ! Apply Hamiltonian to obtain effective field
                  call timing(0,'Evolution     ','OF')
                  call timing(0,'Hamiltonian   ','ON')
                  call effective_field(Natom, Mensemble, 1, Natom, &
                     do_jtensor, exc_inter, do_dm, do_pd, do_biqdm, do_bq, &
                     do_dip,emomM, mmom, external_field,time_external_field, beff, beff1, beff2, &
                     OPT_flag, max_no_constellations, maxNoConstl, &
                     unitCellType, constlNCoup, constellations, constellationsNeighType, mult_axis,denergy)
                  call timing(0,'Hamiltonian   ','OF')
                  call timing(0,'Evolution     ','ON')

               endif

               ! Perform second (corrector) step of SDE solver
               call evolve_second(Natom, Mensemble, Landeg, llg, SDEalgh, bn, iplambda1_array(i,:), &
                  ipdelta_t(i), relaxtime, beff, beff2,b2eff, btorque, emom, emom2,'N', &
                  nlist,nlistsize,constellationsUnitVec,constellationsUnitVec2,constellationsMag, &
                  constellations,unitCellType,OPT_flag,cos_thr,max_no_constellations,'N',she_btorque) !!!,&

               ! Implicit solvers
            else

               ! Fix-point iteration until convergence.
               imp_count = 1
               converged = .false.

               ! Before first iteration, copy emom to emom2
               !$omp parallel do default(shared) private(j,k) collapse(2)
               do j=1, Mensemble
                  do k=1, Natom
                     emom2(:,k,j)  = emom(:,k,j)
                  end do
               end do
               !$omp end parallel do

               do
                  call evolve_imp(Natom, Mensemble, Landeg, llg, SDEalgh, bn, iplambda1_array(i,:), iplambda2_array(i,:), NA, &
                     compensate_drift, ipdelta_t(i), relaxtime,ipTemp_array(:,i),temprescale, beff, b2eff, thermal_field,beff2, &
                     btorque, field1, field2, emom, emom2, emomM, mmom, mmomi ,'N', do_site_ip_damping,&
                     nlist,nlistsize,constellationsUnitVec,constellationsUnitVec2,constellationsMag, &
                     constellations,unitCellType,OPT_flag,cos_thr,max_no_constellations,'N',she_btorque, converged)
                  write(*,*) 'Fix point iteration nr: ', imp_count
                  if (converged) exit
                  imp_count = imp_count + 1
                  call timing(0,'Evolution     ','OF')
                  call timing(0,'Hamiltonian   ','ON')
                  call effective_field(Natom, Mensemble, 1, Natom, &
                     do_jtensor, exc_inter, do_dm, do_pd, do_biqdm, do_bq, &
                     do_dip,emomM, mmom, external_field,time_external_field, beff, beff1, beff2, &
                     OPT_flag, max_no_constellations, maxNoConstl, &
                     unitCellType, constlNCoup, constellations, constellationsNeighType, mult_axis,denergy)
                  call timing(0,'Hamiltonian   ','OF')
                  call timing(0,'Evolution     ','ON')
               end do
               write(*,*) 'Fix point iteration nr: ', imp_count, ' converged.'

            end if

            ! Update magnetic moments after time evolution step
            call moment_update(Natom, Mensemble, mmom, mmom0, mmom2, emom, emom2, emomM, mmomi, mompar, initexc)

            ipstep = ipstep + 1
            call timing(0,'Evolution     ','OF')
            call timing(0,'Initial       ','ON')

         enddo

      enddo
   end subroutine sd_iphase


   !---------------------------------------------------------------------------
   !> @brief
   !> Spin Dynamics measurement phase
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine sd_mphase()
      !
      use BLS,                 only : calc_bls
      use Energy,              only : calc_energy
      use Damping
      use Evolution
      use InputData
      use FieldData,           only : beff, beff1, beff2, beff3, b2eff, sitefld, external_field, &
         field1, field2, time_external_field, allocation_field_time, thermal_field
      use DemagField
      use MomentData
      use FieldPulse
      use SystemData,           only : atype, anumb, Landeg
      use Correlation,          only : correlation_wrapper
      use Temperature
      use ChemicalData
      use Measurements
      use Polarization
      use UpdateMoments
      use MicroWaveField
      use SimulationData,       only : bn, rstep, mstep
      use HamiltonianData
      use CalculateFields
      use AutoCorrelation,      only : autocorr_sample, do_autocorr
      use SystemData,          only: coord
      use OptimizationRoutines
      use AdaptiveTimeStepping
      use HamiltonianActions
      use diamag
      use prn_averages
      use prn_trajectories
      use Gradients
      use SpinTorques
      use QHB, only : qhb_rescale, do_qhb, qhb_mode

      implicit none
      logical :: time_dept_flag, deltat_correction_flag
      integer :: cgk_flag, cgk_flag_p, adapt_step, cr_flag, spt_flag, ntmp, bcgk_flag, cgk_flag_pc
      integer :: j, k, scount_pulse, sstep
      ! Phase flag indicator (true for sd initial phase; false for sd measurement phase)
      ! Used in order to separate between phases in an adaptive time step environment with spin correlation
      logical :: sd_phaseflag

      real(dblprec) :: temprescale,dummy, denergy

      ! Implicit midpoint methods
      integer :: imp_count
      logical :: converged

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

      time_dept_flag=allocation_field_time(mwf,do_gauss,mwf_gauss,mov_gauss,mwf_mov_gauss,mwf_gauss_spatial,&
         mov_circle,mwf_mov_circle,mov_square,mwf_mov_square)

      time_dept_flag = time_dept_flag .or. (do_bpulse.gt.0.and.do_bpulse.lt.5)

      ! Calculate demagnetization field
      if(demag=='Y') then
         call calc_demag(Natom, Mensemble, demag1, demag2, demag3, demagvol, emomM)
      endif

      call calc_bls(N1,N2,N3,C1,C2,C3,Natom, Mensemble, simid, coord,emomM, sstep, delta_t, 0)

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
         call init_polarization(Natom, Mensemble, max_no_neigh, &
            nlist, coord, 1)
      end if

      ! Calculate the external static fields, the can be calculated once and it does not needs to be calculated again
      call calc_external_fields(Natom, Mensemble, NA, hfield, anumb, external_field, &
         do_bpulse,sitefld,sitenatomfld)
      !
      ! Adaptive Time Stepping Region
      if(adaptive_time_flag) then
         call calculate_omegainit(omega_max, larmor_numrev, delta_t)
      end if

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

      ! --Optimization Region-- !
      ! Allocation and initial opt build

      if (OPT_ON) then
         call timing(0,'Hamiltonian   ','OF')
         call timing(0,'BuildOptReg   ','ON')
         call buildOptimizationRegions(na,natom,Mensemble,max_no_neigh,&
            atype,emom,mmom,ncoup,nlist,nlistsize,cellPosNumNeigh,cos_thr)
         call timing(0,'BuildOptReg   ','OF')
         call timing(0,'Hamiltonian   ','ON')
      end if

      ! Initialize cumulant counter
      Navrgcum=0
      ! Perform nstep simulation steps
      ! Observe that Fortran cannot handle a regular do loop if nstep changes dynamically
      ! Hence, a do while loop has been introduced to allow for adaptive time stepping
      mstep=rstep+1

      do while (mstep.LE.rstep+nstep)

         if (time_dept_flag) then
            ! Calculate Microwave fields (time dependent fields)
            call calculate_mwf_fields(Natom,mstep,rstep+nstep-1,delta_t,coord,0)

            ! Calculate the total time dependent fields
            call calc_external_time_fields(Natom, Mensemble, time_external_field, &
               do_bpulse, demag, mwf, mwf_gauss_spatial, do_gauss, mwf_gauss, mov_gauss, mwf_mov_gauss, &
               bpulsefield, demagfld, mwffield,gauss_mwffield,site_mwffield,&
               gauss_spatial_site_mwffield,gauss_spatial,gauss_site_mwffield,mov_gauss_spatial,mwf_mov_gauss_spatial,&
               mov_circle,mwf_mov_circle,mov_circle_spatial,mwf_mov_circle_spatial,mov_square,&
               mwf_mov_square,mov_square_spatial,mwf_mov_square_spatial)
         else
            time_external_field=0.0D0
         endif

         ! Measure averages and trajectories
         call measure(Natom, Mensemble, NT, NA, N1, N2, N3, simid, mstep, emom, emomM, mmom, &
            Nchmax, do_ralloy, Natom_full, asite_ch, achem_ch, atype, plotenergy, Temp, &
            real_time_measure,delta_t,logsamp, max_no_neigh, nlist, ncoup,nlistsize,&
            thermal_field,beff,beff1,beff3,coord,ind_mom,ind_nlistsize,ind_nlist,atype_ch)

         ! Calculate total and term resolved energies
         if(plotenergy>0) then
            call timing(0,'Measurement   ','OF')
            call timing(0,'Energy        ','ON')
            if (mod(mstep-1,avrg_step)==0) then
               tenergy=0.0d0
               call calc_energy(Natom, Nchmax,Mensemble, conf_num,emom, emomM, mmom,simid, plotenergy, &
                  mstep, external_field,time_external_field,tenergy, eenergy, lsfenergy, totene, &
                  max_no_neigh, nlistsize, nlist, ncoup, ncoupD,exc_inter, &
                  do_dm, max_no_dmneigh, dmlistsize, dmlist, dm_vect, &
                  do_pd, nn_pd_tot, pdlistsize, pdlist, pd_vect, &
                  do_biqdm, nn_biqdm_tot, biqdmlistsize, biqdmlist, biqdm_vect, &
                  do_bq, nn_bq_tot, bqlistsize, bqlist, j_bq, &
                  do_dip, qdip,taniso, eaniso, kaniso,sb,&
                  mult_axis,taniso_diff, eaniso_diff, kaniso_diff,sb_diff,&
                  delta_t,real_time_measure,&
                  do_lsf,fs_nlist,fs_nlistsize,nind,lsf_interpolate,lsf_field,Temp)
            end if
            call timing(0,'Energy        ','OF')
            call timing(0,'Measurement   ','ON')

         endif

         call timing(0,'Measurement   ','OF')
         call timing(0,'SpinCorr      ','ON')

         ! Spin correlation
         ! Sample magnetic moments for correlation functions
         call correlation_wrapper(Natom, Mensemble, coord, simid,emomM, mstep, delta_t, NT, atype,Nchmax,achtype,cgk_flag,cgk_flag_p)

         call calc_bls(N1,N2,N3,C1,C2,C3,Natom, Mensemble, simid, coord,emomM, mstep, delta_t, 1)

         call timing(0,'SpinCorr      ','OF')
         call timing(0,'Measurement   ','ON')

         ! Write simulation status for each 5% of the simulation length
         if(nstep>20) then
            if(mod(mstep,(rstep+nstep)/20)==0) then
               call calc_mavrg(Natom,Mensemble,emomM,mavg)
               ! Only calculate binderc if it is not calculated elsewhere
               ! If calculated elsewhere, display binderc will not be real binderc
               if(do_cumu=='N') then
                  call calc_and_print_cumulant(Natom, Mensemble, emomM, simid, Temp, plotenergy, cumu_buff, .false.)
               endif
               write(*,'(2x,a,i4,a,f10.6)',advance='no') &
                  "MP",mstep*100/(rstep+nstep),"% done. Mbar:",mavg
               if(plotenergy>0) then
                  write(*,'(a,f12.6)',advance='no') ". Ebar:", totene
               endif
               sstep = logstep(mstep,logsamp)
               if(mod(sstep,cumu_step)==0)then
                  write(*,'(a,f8.5,a)',advance='no') ". U:",binderc,"."
               endif
               write(*,*) ''
            endif
         else
            write(*,'(2x,a,i3,a,G13.6)')   "Iteration",mstep," Mbar ",mavg
         endif

         !Adjust QHB
         if(do_qhb=='Q' .or. do_qhb=='R' .or. do_qhb=='T') then
            if(qhb_mode=='MT') then
               if (mod(mstep,(rstep+nstep)/40)==0) then
                  call calc_mavrg(Natom,Mensemble,emomM,mavg)
                  call qhb_rescale(Temp,temprescale,do_qhb,qhb_mode,mavg)
               endif
            endif
         endif
         call timing(0,'Measurement   ','OF')
         call timing(0,'Hamiltonian   ','ON')

         ! Calculate spin transfer torque contributions to the local field
         if(stt=='A'.or.stt=='F'.or.do_she=='Y') then
            call calculate_spintorques(Natom, Mensemble,lambda1_array,emom,mmom)
         end if

         ! --Optimization Region-- !
         if (OPT_ON) then
            if (mod(mstep-rstep,OPT_rebuild_time)==0) then
               call timing(0,'Hamiltonian   ','OF')
               call timing(0,'BuildOptReg   ','ON')
               call buildOptimizationRegions(na,natom,Mensemble,max_no_neigh,atype,&
                  emom,mmom,ncoup,nlist,nlistsize,cellPosNumNeigh,cos_thr)
               call timing(0,'BuildOptReg   ','OF')
               call timing(0,'Hamiltonian   ','ON')
            end if
         end if

         ! Apply Hamiltonian to obtain effective field
         call effective_field(Natom, Mensemble, 1, Natom, &
            do_jtensor, exc_inter, do_dm, do_pd, do_biqdm, do_bq, &
            do_dip,emomM, mmom, external_field,time_external_field, beff, beff1, beff2, &
            OPT_flag, max_no_constellations, maxNoConstl, &
            unitCellType, constlNCoup, constellations, constellationsNeighType, mult_axis,denergy)

         call timing(0,'Hamiltonian   ','OF')
         call timing(0,'Evolution     ','ON')

         ! One- and two-step solvers
         if(SDEalgh<=20) then

            ! Calculate average field for use with multiple heat baths
            if (llg==2.or.llg==3) call calc_averagefield(Natom,Mensemble,beff1,beff2,field1,field2)
            ! Check if this changes
            thermal_field=0.0D0
            call evolve_first(Natom, Mensemble, Landeg, llg, SDEalgh, bn, lambda1_array, lambda2_array, NA, &
               compensate_drift, delta_t, relaxtime, Temp_array, temprescale, beff, b2eff, thermal_field, &
               beff2, btorque, field1, field2, emom, emom2, emomM, mmom, mmomi, stt, do_site_damping,&
               nlist,nlistsize,constellationsUnitVec,constellationsUnitVec2,constellationsMag, &
               constellations,unitCellType,OPT_flag,cos_thr,max_no_constellations,do_she,she_btorque) !!!,&
            call timing(0,'Evolution     ','OF')
            call timing(0,'Hamiltonian   ','ON')

            if(SDEalgh==1.or.SDEalgh==4.or.SDEalgh==5.or.SDEalgh==6.or.SDEalgh==7.or.SDEalgh==11) then

               if (time_dept_flag) then
                  ! Calculate Microwave fields (time dependent fields)
                  call calculate_mwf_fields(Natom,mstep,rstep+nstep-1,delta_t,coord,1)

                  ! Calculate the total time dependent fields
                  call calc_external_time_fields(Natom, Mensemble, time_external_field, &
                     do_bpulse, demag, mwf, mwf_gauss_spatial, do_gauss, mwf_gauss, mov_gauss, mwf_mov_gauss, &
                     bpulsefield, demagfld, mwffield,gauss_mwffield,site_mwffield,&
                     gauss_spatial_site_mwffield,gauss_spatial,gauss_site_mwffield,mov_gauss_spatial,mwf_mov_gauss_spatial,&
                     mov_circle,mwf_mov_circle,mov_circle_spatial,mwf_mov_circle_spatial,mov_square,&
                     mwf_mov_square,mov_square_spatial,mwf_mov_square_spatial)
               else
                  time_external_field=0.0D0
               endif

               ! Apply Hamiltonian to obtain effective field
               call effective_field(Natom, Mensemble, 1, Natom, &
                  do_jtensor, exc_inter, do_dm, do_pd, do_biqdm, do_bq, &
                  do_dip,emomM, mmom, external_field,time_external_field, beff, beff1, beff2, &
                  OPT_flag, max_no_constellations, maxNoConstl, &
                  unitCellType, constlNCoup, constellations, constellationsNeighType, mult_axis, denergy)

            end if
            call timing(0,'Hamiltonian   ','OF')
            call timing(0,'Evolution     ','ON')

            ! Re-Calculate spin transfer torque
            if(stt=='A'.or.stt=='F'.or.do_she=='Y') then
               call calculate_spintorques(Natom, Mensemble,lambda1_array,emom,mmom)
            end if


            ! Perform second (corrector) step of SDE solver
            call evolve_second(Natom, Mensemble, Landeg, llg, SDEalgh, bn, lambda1_array, &
               delta_t, relaxtime, beff, beff2,b2eff, btorque, emom, emom2, stt, &
               nlist,nlistsize,constellationsUnitVec,constellationsUnitVec2,constellationsMag, &
               constellations,unitCellType,OPT_flag,cos_thr,max_no_constellations,do_she,she_btorque)

            ! Implicit solvers
         else

            ! Fix-point iteration until convergence.
            imp_count = 1
            converged = .false.

            ! Before first iteration, copy emom to emom2
            !$omp parallel do default(shared) private(j,k) collapse(2)
            do j=1, Mensemble
               do k=1, Natom
                  emom2(:,k,j)  = emom(:,k,j)
               end do
            end do
            !$omp end parallel do

            do
               call evolve_imp(Natom, Mensemble, Landeg, llg, SDEalgh, bn, lambda1_array, lambda2_array, NA, &
                  compensate_drift, delta_t, relaxtime, Temp_array, temprescale, beff, b2eff, thermal_field, &
                  beff2, btorque, field1, field2, emom, emom2, emomM, mmom, mmomi, stt, do_site_damping,&
                  nlist,nlistsize,constellationsUnitVec,constellationsUnitVec2,constellationsMag, &
                  constellations,unitCellType,OPT_flag,cos_thr,max_no_constellations,do_she,she_btorque, converged)
               write(*,*) 'Fix point iteration nr: ', imp_count
               if (converged) exit
               imp_count = imp_count + 1
               call effective_field(Natom, Mensemble, 1, Natom, &
                  do_jtensor, exc_inter, do_dm, do_pd, do_biqdm, do_bq, &
                  do_dip,emomM, mmom, external_field,time_external_field, beff, beff1, beff2, &
                  OPT_flag, max_no_constellations, maxNoConstl, &
                  unitCellType, constlNCoup, constellations, constellationsNeighType, mult_axis, denergy)
            end do
            write(*,*) 'Fix point iteration nr: ', imp_count, ' converged.'
         end if

         call timing(0,'Evolution     ','OF')
         call timing(0,'Moments       ','ON')

         ! Update magnetic moments after time evolution step
         call moment_update(Natom, Mensemble, mmom, mmom0, mmom2, emom, emom2, emomM, mmomi, mompar, initexc)

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

      end do ! End loop over simulation steps

      ! Measure averages and trajectories
      call measure(Natom, Mensemble, NT, NA, N1, N2, N3, simid, mstep,emom, emomM, mmom, &
         Nchmax, do_ralloy, Natom_full, asite_ch, achem_ch, atype, plotenergy, Temp, &
         real_time_measure,delta_t,logsamp,max_no_neigh, nlist,ncoup,nlistsize, &
         thermal_field,beff,beff1,beff3,coord,ind_mom,ind_nlistsize,ind_nlist,atype_ch)

      ! Print remaining measurements
      call flush_measurements(Natom, Mensemble, NT, NA, N1, N2, N3, simid, mstep,emom, mmom, &
         Nchmax,atype,real_time_measure,rstep+nstep,do_ralloy,Natom_full,atype_ch,ind_mom)


      ! Print final polarization, chirality and local polarization
      if (do_pol=='Y') then
         call init_polarization(Natom, Mensemble, max_no_neigh, &
            nlist, coord, -1)
      end if

      if (do_spintemp=='Y') then
         call spintemperature(Natom,Mensemble,mstep,1,simid,emomM,beff,2)
      endif
      call timing(0,'Measurement   ','OF')
   end subroutine sd_mphase

   !---------------------------------------------------------------------------
   !> @brief
   !> CUDA implemented sd measurement phase
   !
   !> @author
   !> Niklas Fejes
   !
   !> @date 2014/08/22 : Thomas Nystrand
   !> - Moved to separate routine
   !---------------------------------------------------------------------------
   subroutine sd_mphaseCUDA()
#ifdef CUDA
      use Chelper
#else
      use NoCuda
#endif
      use Damping
      use SpinTorques, only : btorque, stt
      use InputData, only : gpu_mode

      ! Common stuff
      ! Copy core fortran data needed by CPP and CUDA solver to local cpp class
      !!! TEMPORARY COMMENTED OUT
      call FortranData_Initiate(stt,btorque)
      !!! TEMPORARY COMMENTED OUT

      ! Let the fortran timing think we are in Measurement
      call timing(0,'Measurement   ','ON')

      ! Start simulation
      if (gpu_mode==1) then  !CUDA
         call cudaMdSim_initiateConstants()
         call cudaMdSim_initiateMatrices()
         call cudaMdSim_measurementPhase()

      else if (gpu_mode==2) then     !C/C++
         call cMdSim_initiateConstants() ! calls mdSimulation.cpp to copy initial constants from fortrandata.hpp
         call cMdSim_initiateFortran()   ! calls mdSimulation.cpp to copy and initialize matrices from fortrandata.hpp
         call cMdSim_measurementPhase()

      else
         stop "Invalid gpu_mode"
      endif
      call timing(0,'Measurement   ','OF')
   end subroutine sd_mphaseCUDA

end module sd_driver
