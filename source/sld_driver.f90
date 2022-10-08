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
module sld_driver

   implicit none

   public

contains

   !---------------------------------------------------------------------------------
   !> @brief
   !> Lattice Dynamics inital phase
   !
   !> @author
   !> Johan Hellsvik
   !---------------------------------------------------------------------------------
   subroutine ld_iphase()
      !
      use InputData
      use MomentData
      use ChemicalData
      use OptimizationRoutines
      !use LatticeMeasurements
      use LatticeInputData
      use LatticeFieldData
      use SimulationData
      use LatticeData
      use LatticeHamiltonianData
      use LatticeHamiltonianActions
      !use LatticeCorrelation
      use Verlet
      use Temperature
      use RandomNumbers

      implicit none

      !logical :: deltat_correction_flag
      !integer :: adapt_step
      !integer :: ucgk_flag, ucgk_flag_p, ucr_flag
      integer :: i, ipstep
      real(dblprec) :: temprescale,temprescalegrad
      real(dblprec) :: lattdamptol

      real(dblprec), dimension(ipnphase,Natom) :: lattdampvec

      integer :: kk,l

      do kk=1,ipnphase
         do l=1,Natom
            lattdampvec(kk,l) = iplattdamp(kk)
         end do
      end do

      temprescale = 1.0_dblprec
      temprescalegrad = 0.0_dblprec

      ! Used to switch between regular velocity Verlet and the GJF Verlet.
      ! AB - note : Not used currently
      lattdamptol =  1e-17

      ! Loop over initial phases
      do i=1,ipnphase

         ipstep=1

         ! Main initial loop
         !--------------------------------------!
         do while(ipstep.LE.ipnstep(i))

            ! Perform nstep simulation steps
            ! Observe that Fortran cannot handle a regular do loop if nstep changes dynamically
            ! Hence, a do while loop has been introduced to allow for adaptive time stepping
            !mstep=rstep+1
            !do while (mstep.LE.rstep+nstep) !+1

            call timing(0,'Hamiltonian   ','ON')

            ! Apply lattice Hamiltonian to obtain effective field: First call
            call effective_latticefield(Natom,Mensemble,1,Natom,do_ll,       &
               do_ml,do_mml,mode,uvec,emomM,latt_external_field,    &
               latt_time_external_field,eeff,eeff1,eeff2,eeff3)

            call timing(0,'Hamiltonian   ','OF')
            call timing(0,'Evolution     ','ON')

            !If needed, generate stochastic forces
            !if(iplattdamp(i) .gt. lattdamptol) then
            if(iptemp(i)>0.0_dblprec) then
               call lattrannum(Natom,Mensemble,NA,delta_t,lattdampvec,              &
                  compensate_drift,Temp_array,temprescale)
            else
               lattranv = 0_dblprec
            end if

            ! Update ionic displacements
            !if(iplattdamp(i) .lt. lattdamptol) then
            !   call u_vverlet(Natom,Mensemble,mioninv,eeff,uvec,vvec,delta_t)
            !else
               call u_gjfverlet(Natom,Mensemble,mioninv,eeff,uvec,vvec,delta_t,     &
                  iplattdamp(i))
            !end if

            call timing(0,'Evolution     ','OF')
            call timing(0,'Hamiltonian   ','ON')

            ! Apply lattice Hamiltonian to obtain effective field: Second call
            call effective_latticefield(Natom,Mensemble,1,Natom,do_ll,       &
               do_ml,do_mml,mode,uvec,emomM,latt_external_field,    &
               latt_time_external_field,e2eff,eeff1,eeff2,eeff3)

            call timing(0,'Hamiltonian   ','OF')
            call timing(0,'Evolution     ','ON')

            ! Update ionic velocities
            !if(iplattdamp(i) .lt. lattdamptol) then
            !   call v_vverlet(Natom,Mensemble,mioninv,eeff,e2eff,vvec,delta_t)
            !else
               call v_gjfverlet(Natom,Mensemble,mioninv,eeff,e2eff,uvec,uvec2,vvec, &
                  delta_t,iplattdamp(i))
               call set_avrgu0(Natom, Mensemble, uvec, vvec)
               call set_avrgp0(Natom, Mensemble, uvec, vvec)
            !end if

            call timing(0,'Evolution     ','OF')

            ipstep = ipstep + 1
            !mstep = mstep + 1

         end do ! End loop over simulation steps
      end do ! End loop over initial phases

   end subroutine ld_iphase

   !---------------------------------------------------------------------------
   !> @brief
   !> Lattice Dynamics measurement phase
   !
   !> @author
   !> Johan Hellsvik
   !---------------------------------------------------------------------------
   subroutine ld_mphase()
      !
      use InputData
      use SystemData,           only : atype, coord
      use MomentData
      use ChemicalData
      use OptimizationRoutines
      use LatticeMeasurements
      use LatticeInputData
      use LatticeFieldData
      use SimulationData
      use LatticeData
      use LatticeHamiltonianData
      use LatticeHamiltonianActions
!     use LatticeCorrelation
      use Correlation
      use Verlet
      use Temperature
      use RandomNumbers

      implicit none

      integer :: ucgk_flag
      real(dblprec) :: temprescale,temprescalegrad
      real(dblprec) :: lattdamptol

      real(dblprec), dimension(Natom) :: lattdampvec
      real(dblprec), dimension(Mensemble) :: ll_energy
      real(dblprec), dimension(Mensemble) :: ml_energy
      real(dblprec), dimension(Mensemble) :: mml_energy
      real(dblprec), dimension(Mensemble) :: ldpot_energy
      real(dblprec), dimension(Mensemble) :: sdpot_energy
      real(dblprec), dimension(Mensemble) :: sldpot_energy
      real(dblprec), dimension(Mensemble) :: totpot_energy
      real(dblprec), dimension(Mensemble) :: sld_single_energy
      real(dblprec), dimension(Mensemble) :: mm_energy0

      lattdampvec = lattdamp
      temprescale = 1.0_dblprec
      temprescalegrad = 0.0_dblprec

      ! Used to switch between regular velocity Verlet and the GJF Verlet.
      lattdamptol = 1e-17

      ucgk_flag = 0

      ! Perform nstep simulation steps
      ! Observe that Fortran cannot handle a regular do loop if nstep changes dynamically
      ! Hence, a do while loop has been introduced to allow for adaptive time stepping
      mstep=rstep+1
      do while (mstep.LE.rstep+nstep) !+1

         call timing(0,'Hamiltonian   ','ON')

         ! Apply lattice Hamiltonian to obtain effective field: First call
         call effective_latticefield(Natom,Mensemble,1,Natom,do_ll,  &
            do_ml,do_mml,mode,uvec,emomM,latt_external_field,               &
            latt_time_external_field,eeff,eeff1,eeff2,eeff3)

         ! Calculate energies of the lattice (and the spin-lattice) Hamiltonians
         ! The energy of the magnetic system is zero for pure lattice dynamics
         sdpot_energy = 0.0_dblprec
         call calc_lattenergies(Natom,Mensemble,1,Mensemble,1,Natom,do_ll,   &
            do_ml,do_mml,mode,uvec,emomM,latt_external_field,       &
            latt_time_external_field,ll_energy,ml_energy,    &
            mml_energy,ldpot_energy,sdpot_energy,sldpot_energy,         &
            totpot_energy,sld_single_energy,mm_energy0,ammom_inp,aemom_inp,NA)
         !ldpot_energy, sdpot_energy, sldpot_energy, totpot_energy, sld_single_energy)
         !ldpot_energy, sdpot_energy, sldpot_energy, totpot_energy)

         call timing(0,'Hamiltonian   ','OF')
         call timing(0,'Measurement   ','ON')

         call local_angular_momentum(Natom, Mensemble, mion,  uvec, vvec, lvec)

         call lattmeasure(simid,Natom,Mensemble,NT,NA,N1,N2,N3,mstep,logsamp,coord, &
            mion,uvec,vvec,lvec,eeff,eeff1,eeff3,ethermal_field,Nchmax,do_ralloy,        &
            Natom_full,asite_ch,achem_ch,atype,Temp,ldpot_energy,sdpot_energy,      &
            sldpot_energy,totpot_energy)

         call timing(0,'Measurement   ','OF')
         call timing(0,'LattCorr      ','ON')

         ! Lattice correlation
         ! Sample displacements and velocities for correlation functions
         call correlation_wrapper(Natom,Mensemble,coord,simid,uvec,mstep,delta_t,  &
            NT,atype,Nchmax,achtype,uc,do_uc,do_ur,ucgk_flag)

         call correlation_wrapper(Natom,Mensemble,coord,simid,vvec,mstep,delta_t,  &
            NT,atype,Nchmax,achtype,vc,do_vc,do_vr,ucgk_flag)

         call correlation_wrapper(Natom,Mensemble,coord,simid,lvec,mstep,delta_t,  &
            NT,atype,Nchmax,achtype,lc,do_lc,do_lr,ucgk_flag)

         !call lattcorrelation_wrapper(Natom,Mensemble,coord,simid,uvec,vvec,mstep,  &
         !   delta_t,NT,atype,Nchmax,achtype,ucgk_flag,ucgk_flag_p)

         call timing(0,'LattCorr      ','OF')
         call timing(0,'Measurement   ','ON')

         ! Write simulation status for each 5% of the simulation length
         !       if(nstep>20) then
         !          if(mod(mstep,(rstep+nstep)/20)==0) then
         !             call calc_mavrg(Natom,Mensemble,emomM,mavg)
         ! Only calculate binderc if it is not calculated elsewhere
         ! If calculated elsewhere, display binderc will not be real binderc
         !             if(do_cumu=='N') then
         !                call calc_and_print_cumulant(Natom, Mensemble, emomM, simid, Temp, plotenergy, cumu_buff, .false.)
         !             endif
         !             write(*,'(2x,a,i4,a,f10.6)',advance='no') &
         !                  "MP",mstep*100/(rstep+nstep),"% done. Mbar:",mavg
         !             if(plotenergy>0) then
         !                write(*,'(a,f12.6)',advance='no') ". Ebar:", totene/Mensemble
         !             endif
         !             sstep = f_logstep(mstep,logsamp)
         !             if(mod(sstep,cumu_step)==0)then
         !                write(*,'(a,f8.5,a)',advance='no') ". U:",binderc,"."
         !             endif
         !             write(*,*) ''
         !          endif
         !       else
         !          write(*,'(2x,a,i3,a,G13.6)')   "Iteration",mstep," Mbar ",mavg
         !       endif
         !Adjust QHB
         !       if(do_qhb=='Q' .or. do_qhb=='R') then
         !          if(qhb_mode=='MT') then
         !             if (mod(mstep,(rstep+nstep)/40)==0) then
         !                call calc_mavrg(Natom,Mensemble,emomM,mavg)
         !                call qhb_rescale(Temp,temprescale,do_qhb,qhb_mode,mavg)
         !              endif
         !          endif
         !       endif

         call timing(0,'Measurement   ','OF')
         call timing(0,'Evolution     ','ON')

         !If needed, generate stochastic forces
         !if(lattdamp .gt. lattdamptol) then
         if(Temp > 0.0_dblprec) then
            call lattrannum(Natom,Mensemble,NA,delta_t,lattdampvec,compensate_drift,&
               Temp_array,temprescale)
         else
            lattranv = 0_dblprec
         end if

         ! Update ionic displacements
         !if(lattdamp .lt. lattdamptol) then
         !   call u_vverlet(Natom, Mensemble, mioninv, eeff, uvec, vvec, delta_t)
         !else
            call u_gjfverlet(Natom,Mensemble,mioninv,eeff,uvec,vvec,delta_t,lattdamp)
         !end if

         call timing(0,'Evolution     ','OF')
         call timing(0,'Hamiltonian   ','ON')

         ! Apply lattice Hamiltonian to obtain effective field: Second call
         call effective_latticefield(Natom,Mensemble,1,Natom,do_ll,  &
            do_ml,do_mml,mode,uvec,emomM,latt_external_field,               &
            latt_time_external_field,e2eff,eeff1,eeff2,eeff3)

         call timing(0,'Hamiltonian   ','OF')
         call timing(0,'Evolution     ','ON')

         ! Update ionic velocities
         !if(lattdamp .lt. lattdamptol) then
         !   call v_vverlet(Natom, Mensemble,mioninv,eeff,e2eff,vvec,delta_t)
         !else
            call v_gjfverlet(Natom,Mensemble,mioninv,eeff,e2eff,uvec,uvec2,vvec,    &
               delta_t,lattdamp)
            call set_avrgu0(Natom, Mensemble, uvec, vvec)
            call set_avrgp0(Natom, Mensemble, uvec, vvec)
         !end if

         call timing(0,'Evolution     ','OF')

         mstep = mstep + 1

      end do ! End loop over simulation steps

      ! Apply lattice Hamiltonian to obtain effective field
      call effective_latticefield(Natom,Mensemble,1,Natom,do_ll,     &
         do_ml,do_mml,mode,uvec,emomM,latt_external_field,                  &
         latt_time_external_field,eeff,eeff1,eeff2,eeff3)

      ! Calculate energies of the lattice and the spin-lattice Hamiltonians
      ! The energy of the magnetic system is zero for pure lattice dynamics
      sdpot_energy = 0.0_dblprec
      call calc_lattenergies(Natom,Mensemble,1,Mensemble,1,Natom,do_ll,      &
         do_ml,do_mml,mode,uvec,emomM,latt_external_field,          &
         latt_time_external_field,ll_energy,ml_energy,       &
         mml_energy,ldpot_energy,sdpot_energy,sldpot_energy,            &
         totpot_energy,sld_single_energy,mm_energy0,ammom_inp,aemom_inp,NA)
      !ldpot_energy, sdpot_energy, sldpot_energy, totpot_energy, sld_single_energy)
      !ldpot_energy, sdpot_energy, sldpot_energy, totpot_energy)

      call local_angular_momentum(Natom, Mensemble, mion,  uvec, vvec, lvec)

      ! Measure averages and trajectories
      call lattmeasure(simid,Natom,Mensemble,NT,NA,N1,N2,N3,mstep,logsamp,coord,    &
         mion,uvec,vvec,lvec,eeff,eeff1,eeff3,ethermal_field,Nchmax,do_ralloy,Natom_full,&
         asite_ch,achem_ch,atype,Temp,ldpot_energy,sdpot_energy,sldpot_energy,      &
         totpot_energy)

      ! Print remaining measurements
      call flush_lattmeasurements(simid,Natom,Mensemble,NT,NA,N1,N2,N3,mstep,mion,  &
         uvec,vvec,coord,Nchmax,atype,Temp)

   end subroutine ld_mphase


   !---------------------------------------------------------------------------
   !> @brief
   !> Spin Lattice Dynamics initial phase
   !
   !> @author
   !> Johan Hellsvik
   !---------------------------------------------------------------------------
   subroutine sld_iphase()
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
      use SystemData,           only : atype, anumb, Landeg,coord
      use SimulationData
      !use Correlation,          only : correlation_wrapper
      use Temperature
      use ChemicalData
      use Measurements
      use Polarization
      use UpdateMoments
      use MicroWaveField
      use AdaptiveTimeStepping
      use HamiltonianData
      use CalculateFields
      !use AutoCorrelation,      only : autocorr_sample, do_autocorr
      use optimizationRoutines
      use HamiltonianActions
      use Gradients
      use prn_averages
      use LatticeMeasurements
      use LatticeInputData
      use LatticeFieldData
      use LatticeData
      use LatticeHamiltonianData
      use LatticeHamiltonianActions
      !use LatticeCorrelation
      use Verlet
      use SpinTorques
      use QHB, only : qhb_rescale, do_qhb, qhb_mode
      use Temperature
      use RandomNumbers
      use VelocityRescale
      use FixedMom
      use macrocells

      implicit none

      logical :: time_dept_flag
      integer :: ucgk_flag, ucgk_flag_p, ucr_flag
      integer :: k, scount_pulse
      integer :: ii, ipstep

      real(dblprec), dimension(ipnphase,Natom) :: lattdampvec

      ! Phase flag indicator (true for sd initial phase; false for sd measurement phase)
      ! Used in order to separate between phases in an adaptive time step environment with spin correlation
      logical :: sd_phaseflag


      real(dblprec) :: temprescale, temprescalegrad, dummy, totenergy

      ! Implicit midpoint methods

      ! Canonical velocity rescaling
      integer :: ndeg, i!, j
      real(dblprec) :: v2, kinenrg0, kinenrg1, kinrescaled, alpha

      integer :: kk,l
      real(dblprec) :: lattdamptol

      do kk=1,ipnphase
         do l=1,Natom
            lattdampvec(kk,l) = iplattdamp(kk)
         end do
      end do

      temprescale = 1.0_dblprec
      temprescalegrad = 0.0_dblprec

      ! Used to switch between regular velocity Verlet and the GJF Verlet.
      lattdamptol = 1e-17

      ! Spin correlation measurements allowed
      !adapt_step = 0
      !deltat_correction_flag = .true.

      ! Measurement phase indicator (false here)
      sd_phaseflag = .false.
      !call timing(0,'Measurement   ','ON')

      ! Flag for G(q) sampling and G(r) calculation
      !cgk_flag=0 ; cgk_flag_p=0 ; cr_flag=0 ; cgk_flag_pc=0
      !spt_flag=0

      ! Flag for G(q) sampling and G(r) calculation
      ucgk_flag=0 ; ucgk_flag_p=0 ; ucr_flag=0
      !uc_tidx=0  ; upt_flag=0

      scount_pulse = 1

      time_dept_flag=allocation_field_time(mwf,do_gauss,mwf_gauss,mov_gauss,mwf_mov_gauss,mwf_gauss_spatial,&
         mov_circle,mwf_mov_circle,mov_square,mwf_mov_square)

      !call print_siminfo()

      ! Calculate demagnetization field
      if(demag=='Y') then
         call calc_demag(Natom, Mensemble, demag1, demag2, demag3, demagvol, emomM)
      endif

      !call calc_bls(N1,N2,N3,C1,C2,C3,Natom, Mensemble, simid, coord,emomM, sstep,delta_t, 0)

      !  Iniitialize magnetic field pulse
      !if (do_bpulse.gt.0.and.do_bpulse.lt.5) then
      !   bpulse_time = delta_t*rstep
      !   call bpulse(do_bpulse, bpulse_time)
      !end if

      ! Initialize polarization including local chirality and polarity
      !if (do_spintemp=='Y') then
      !   ntmp=nstep/spintemp_step
      !   call spintemperature(Natom,Mensemble,mstep,ntmp,simid,emomM,beff,0)
      !end if

      !if (do_pol=='Y') then
      !   call init_polarization(Natom, Mensemble, max_no_neigh, &
      !        nlist, coord, 1)
      !end if

      ! Calculate the external static fields, then can be calculated once and do not need to be calculated again
      call calc_external_fields(Natom,Mensemble,hfield,anumb,external_field,     &
         do_bpulse,sitefld,sitenatomfld)
      !
      ! Adaptive Time Stepping Region
      !if(adaptive_time_flag) then
      !   call calculate_omegainit(omega_max, larmor_numrev, delta_t)
      !end if

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
      !if (OPT_ON) then
      !call timing(0,'Hamiltonian   ','OF')
      !   call timing(0,'BuildOptReg   ','ON')
      !   call buildOptimizationRegions(na,natom,Mensemble,max_no_neigh,&
      !        atype,emom,mmom,ncoup,nlist,nlistsize,cellPosNumNeigh,cos_thr)
      !   call timing(0,'BuildOptReg   ','OF')
      !   !call timing(0,'Hamiltonian   ','ON')
      !end if

      ! Initialize cumulant counter
      !Navrgcum=0
      ! Perform nstep simulation steps
      ! Observe that Fortran cannot handle a regular do loop if nstep changes dynamically
      ! Hence, a do while loop has been introduced to allow for adaptive time stepping
      !mstep=rstep+1

      ! Loop over initial phases

      write(*,*) "Starts initial phases ipnphase ", ipnphase

      do ii=1,ipnphase

         write(*,*) "Starts initial phase ", ii, " ipnstep(ii) ", ipnstep(ii)
         ipstep=1

         ! Write output to stdout
         if (do_site_ip_damping=='Y') then
            write (*,'(2x,a25,i3,a10,G11.4,a10,i10,a10,G14.4)') &
            "Performing initial phase: ", ii, " Temp. ", ipTemp(ii), "nsteps ",     &
            ipnstep(ii), "dt ", ipdelta_t(ii)
         else
            write (*,'(2x,a25,i3,a10,G11.4,a10,i10,a10,G14.4,a10,2G14.4,a12,G14.4)')&
            !write (*,'(2x,a25,i3,a10,G11.4,a10,i10,a10,G14.4,a10,2G14.4)') &
            "Performing initial phase: ", ii, " Temp. ", ipTemp(ii), "nsteps ",     &
            ipnstep(ii), "dt ",ipdelta_t(ii),"lambda ",iplambda1(ii),iplambda2(ii), &
            'Lattdamp(ii) ', iplattdamp(ii)
            !ipnstep(ii), "dt ", ipdelta_t(ii), "lambda ", iplambda1(ii), iplambda2(ii)
         endif

         ! Main initial loop
         !--------------------------------------!
         do while(ipstep.LE.ipnstep(ii))

            !call timing(0,'Hamiltonian   ','ON')

            if (time_dept_flag) then
               ! Calculate Microwave fields (time dependent fields)
               call calculate_mwf_fields(Natom,mstep,rstep+nstep-1,ipdelta_t(ii),   &
                  coord,0)

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

            ! --Optimization Region-- !
            if (OPT_ON) then
               if (mod(mstep-rstep,OPT_rebuild_time)==0) then
                  !call timing(0,'Hamiltonian   ','OF')
                  !call timing(0,'BuildOptReg   ','ON')
                  call buildOptimizationRegions(na,natom,nHam,Mensemble,            &
                     ham%max_no_neigh,atype,emom,mmom,ham%ncoup,ham%nlist,          &
                     ham%nlistsize,ham%aham,cellPosNumNeigh,cos_thr)
                  !call timing(0,'BuildOptReg   ','OF')
                  !call timing(0,'Hamiltonian   ','ON')
               end if
            end if

            ! Apply Hamiltonian to obtain effective field
            call effective_field(Natom,Mensemble,1,Natom,emomM,mmom, &
               external_field,time_external_field,beff,beff1,beff2,OPT_flag,           &
               max_no_constellations,maxNoConstl,unitCellType,constlNCoup,             &
               constellations, constellationsNeighType,totenergy,           &
               Num_macro,cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)

            call effective_bmixedfield(Natom,Mensemble,1,Natom,do_ml,do_mml,&
               uvec,emomM,beff,beff3)

            ! Apply lattice Hamiltonian to obtain effective field: First call
            call effective_latticefield(Natom,Mensemble,1,Natom,do_ll,       &
               do_ml,do_mml,mode,uvec,emomM,latt_external_field,    &
               latt_time_external_field,eeff,eeff1,eeff2,eeff3)

            !call timing(0,'Hamiltonian   ','OF')
            !call timing(0,'Measurement   ','ON')


            !Adjust QHB
            if(do_qhb=='Q' .or. do_qhb=='R' .or. do_qhb=='P' .or. do_qhb=='T') then
               if(qhb_mode=='MT') then
                  if (mod(mstep,(rstep+nstep)/40)==0) then
                     call calc_mavrg(Natom,Mensemble,emomM,mavg)
                     call qhb_rescale(Temp,temprescale,temprescalegrad,do_qhb,qhb_mode,mavg)
                  endif
               endif
            endif
            !call timing(0,'Measurement   ','OF')
            !call timing(0,'Hamiltonian   ','ON')

            ! Here both the B-field and the E-field have been calculated for the first time of the time-step
            ! Note that the E-field was calculated early in the loop over mstep to have fields
            ! available for lattmeasurements

            ! One- and two-step solvers
            if(ipSDEalgh<=20) then

               ! Calculate average field for use with multiple heat baths
               if (llg==2.or.llg==3) then
                  call calc_averagefield(Natom,Mensemble,beff1,beff2,field1,field2)
               endif
               ! Check if this changes
               thermal_field=0.0_dblprec
               call evolve_first(Natom,Mensemble,Landeg,llg,ipSDEalgh,bn,             &
                  iplambda1_array(ii,:),iplambda2_array(ii,:),NA,compensate_drift,  &
                  ipdelta_t(ii),relaxtime,ipTemp_array(:,ii),temprescale,beff,b2eff,&
                  thermal_field,beff2,btorque,field1,field2,emom,emom2,emomM,mmom,  &
                  mmomi,stt,do_site_damping,ham%nlist,ham%nlistsize,                &
                  constellationsUnitVec,constellationsUnitVec2,constellationsMag,   &
                  constellations,unitCellType,OPT_flag,cos_thr,                     &
                  max_no_constellations,do_she,she_btorque,Nred,red_atom_list,      &
                  do_sot,sot_btorque)

               !If needed, generate stochastic forces
               !if(iplattdamp(ii) .gt. lattdamptol) then
               if(iptemp(ii)>0.0_dblprec) then
                  !if(iplattdamp .gt. lattdamptol) then
                  call lattrannum(Natom,Mensemble,NA,ipdelta_t(ii),                 &
                     lattdampvec(ii,:),compensate_drift,ipTemp_array(:,ii),temprescale)
                     !lattdampvec(ii,:),compensate_drift,Temp_array,temprescale)
                  !call lattrannum(Natom, Mensemble, NA, ipdelta_t(ii), lattdampvec, compensate_drift, Temp_array, temprescale)
               else
                  lattranv = 0_dblprec
               end if

               ! Update ionic displacements
                  call u_gjfverlet(Natom,Mensemble,mioninv,eeff,uvec,vvec,          &
                     ipdelta_t(ii),iplattdamp(ii))

               if (time_dept_flag) then
                  ! Calculate Microwave fields (time dependent fields)
                  call calculate_mwf_fields(Natom,mstep,rstep+nstep-1,ipdelta_t(ii),&
                     coord,1)

                  ! Calculate the total time dependent fields
                  call calc_external_time_fields(Natom,Mensemble,                   &
                     time_external_field,do_bpulse,demag,mwf,mwf_gauss_spatial,     &
                     do_gauss,mwf_gauss,mov_gauss,mwf_mov_gauss,bpulsefield,        &
                     demagfld, mwffield,gauss_mwffield,site_mwffield,               &
                     gauss_spatial_site_mwffield,gauss_spatial,gauss_site_mwffield, &
                     mov_gauss_spatial,mwf_mov_gauss_spatial,mov_circle,            &
                     mwf_mov_circle,mov_circle_spatial,mwf_mov_circle_spatial,      &
                     mov_square,mwf_mov_square,mov_square_spatial,                  &
                     mwf_mov_square_spatial)
               else
                  time_external_field=0.0_dblprec
               endif

               ! Apply Hamiltonian to obtain effective field
               call effective_field(Natom,Mensemble,1,Natom, &
                  emomM,mmom,external_field,time_external_field,beff,beff1,beff2,   &
                  OPT_flag,max_no_constellations,maxNoConstl,unitCellType,          &
                  constlNCoup,constellations,constellationsNeighType,               &
                  totenergy,Num_macro,cell_index,emomM_macro,macro_nlistsize,NA,N1, &
                  N2,N3)

               call effective_bmixedfield(Natom,Mensemble,1,Natom,do_ml,do_mml,     &
                  uvec,emomM,beff,beff3)

            end if
            !call timing(0,'Hamiltonian   ','OF')
            !call timing(0,'Evolution     ','ON')

            ! Re-Calculate spin transfer torque
            if(stt=='A'.or.stt=='F'.or.do_she=='Y'.or.do_sot=='Y') then
               call calculate_spintorques(Natom, Mensemble,iplambda1_array,emom,mmom)
            end if

            ! Perform second (corrector) step of SDE solver
            ! SLDTODO check how this choice of having the second step of spin evolution occuring
            ! SLDTODO before the second call to effective_latticefield performs in practice
            call evolve_second(Natom,Mensemble,Landeg,llg,ipSDEalgh,bn,             &
               iplambda1_array(ii,:),ipdelta_t(II),relaxtime,beff,beff2,b2eff,      &
               btorque,emom,emom2,stt,ham%nlist,ham%nlistsize,constellationsUnitVec,&
               constellationsUnitVec2,constellationsMag,constellations,unitCellType,&
               OPT_flag,cos_thr,max_no_constellations,do_she,she_btorque,Nred,      &
               red_atom_list,do_sot,sot_btorque)

            !call timing(0,'Evolution     ','OF')
            !call timing(0,'Hamiltonian   ','ON')

            ! Apply lattice Hamiltonian to obtain effective field: Second call
            call effective_latticefield(Natom,Mensemble,1,Natom,do_ll,       &
               do_ml,do_mml,mode,uvec,emomM,latt_external_field,    &
               latt_time_external_field,e2eff,eeff1,eeff2,eeff3)

            !call timing(0,'Hamiltonian   ','OF')
            !call timing(0,'Evolution     ','ON')

            ! Update ionic velocities
               call v_gjfverlet(Natom,Mensemble,mioninv,eeff,e2eff,uvec,uvec2,vvec, &
                  ipdelta_t(ii),iplattdamp(ii))
               call set_avrgu0(Natom, Mensemble, uvec, vvec)
               call set_avrgp0(Natom, Mensemble, uvec, vvec)

            !call timing(0,'Evolution     ','OF')
            !call timing(0,'Moments       ','ON')

            ! Update magnetic moments after time evolution step
            call moment_update(Natom,Mensemble,mmom,mmom0,mmom2,emom,emom2,emomM,   &
               mmomi,mompar,initexc)

            !call timing(0,'Moments       ','OF')
            !call timing(0,'Measurement   ','ON')

            ! Update magnetic field pulse
            if (do_bpulse.gt.0.and.do_bpulse.lt.5) then
               if (scount_pulse == bpulse_step) then
                  bpulse_time = ipdelta_t(ii)*mstep
                  call bpulse(do_bpulse, bpulse_time)
                  scount_pulse=1
               else
                  scount_pulse=scount_pulse+1
               end if
            end if

            ipstep = ipstep + 1

            ! Canonical velocity rescaling
            if ( do_velrsc=='Y' .and. mod(mstep-1,velrsc_step)==0 ) then
               ndeg = 3 * Natom - 3
               do k=1,Mensemble
                  kinenrg0 = 0.0_dblprec
                  do i=1, Natom
                     v2 = vvec(1,i,k)*vvec(1,i,k) + vvec(2,i,k)*vvec(2,i,k) + vvec(3,i,k)*vvec(3,i,k)
                     kinenrg0 = kinenrg0 + 0.5_dblprec * amu * mion(i,1) * v2 * angstrom**2 / ev
                  end do
                  kinenrg1 = ndeg * k_bolt * temp / 2 / ev
                  velrsc_taut=10_dblprec

                  kinrescaled = resamplekin(kinenrg0, kinenrg1, ndeg, velrsc_taut)
                  alpha = sqrt(kinrescaled/kinenrg0)

                  ! Rescale velocities
                  call rescale_vel(Natom, Mensemble, uvec, vvec, alpha)

                  ! Set average linear momenta to zero
                  call set_avrgp0(Natom, Mensemble, uvec, vvec)

                  ! Set average displacement to zero
                  call set_avrgu0(Natom, Mensemble, uvec, vvec)

               end do
            end if

            if(do_set_avrgp0 == 'Y') then
               ! Set average linear momenta to zero
               call set_avrgp0(Natom, Mensemble, uvec, vvec)
            end if

            if(do_set_avrgu0 == 'Y') then
               ! Set average displacement to zero
               call set_avrgu0(Natom, Mensemble, uvec, vvec)
            end if

         end do ! End loop over simulation steps
      end do ! End loop over initial phases

      ! Apply lattice Hamiltonian to obtain effective field
      call effective_latticefield(Natom,Mensemble,1,Natom,do_ll,     &
         do_ml,do_mml,mode,uvec,emomM,latt_external_field,                  &
         latt_time_external_field,eeff,eeff1,eeff2,eeff3)

   end subroutine sld_iphase


   !---------------------------------------------------------------------------
   !> @brief
   !> Spin Lattice Dynamics measurement phase
   !
   !> @author
   !> Johan Hellsvik
   !---------------------------------------------------------------------------
   subroutine sld_mphase()
      !
      use BLS,                 only : calc_bls
      use Energy,              only : calc_energy
      use Damping
      use Evolution
      use InputData
      use FieldData,           only : beff, beff1, beff2, beff3, b2eff, sitefld,    &
      external_field,field1, field2, time_external_field, allocation_field_time,    &
      thermal_field
      use DemagField
      use MomentData
      use FieldPulse
      use SystemData,           only : atype, anumb, Landeg,coord
      use SimulationData
      use Correlation
      use Temperature
      use ChemicalData
      use Measurements
      use Polarization
      use UpdateMoments
      use MicroWaveField
      use AdaptiveTimeStepping
      use HamiltonianData
      use CalculateFields
      use AutoCorrelation,      only : autocorr_sample, do_autocorr
      use optimizationRoutines
      use HamiltonianActions
      use Gradients
      use prn_averages
      use prn_latticeaverages, only : f_iontemp
      use LatticeMeasurements
      use LatticeInputData
      use LatticeFieldData
      use LatticeData
      use LatticeHamiltonianData
      use LatticeHamiltonianActions
!     use LatticeCorrelation
      use Verlet
      use SpinTorques
      use QHB, only : qhb_rescale, do_qhb, qhb_mode
      use Temperature
      use RandomNumbers
      use VelocityRescale
      use FixedMom
      use Math_functions, only : f_logstep
      use Temperature_3TM
      use macrocells

      implicit none

      logical :: time_dept_flag, deltat_correction_flag
      integer :: cgk_flag, cgk_flag_p, adapt_step, cr_flag, cgk_flag_pc, spt_flag, ntmp
      integer :: ucgk_flag, ucr_flag
      integer :: vcgk_flag, vcr_flag
      integer :: lcgk_flag, lcr_flag
      integer :: k, scount_pulse, sstep

      real(dblprec), dimension(Natom) :: lattdampvec
      real(dblprec), dimension(Mensemble) :: ll_energy
      real(dblprec), dimension(Mensemble) :: ml_energy
      real(dblprec), dimension(Mensemble) :: mml_energy
      real(dblprec), dimension(Mensemble) :: ldpot_energy
      real(dblprec), dimension(Mensemble) :: sdpot_energy
      real(dblprec), dimension(Mensemble) :: sldpot_energy
      real(dblprec), dimension(Mensemble) :: totpot_energy
      real(dblprec), dimension(Mensemble) :: sld_single_energy
      real(dblprec), dimension(Mensemble) :: mm_energy0

      real(dblprec) :: Temp_s !< Spin temperature (for 3TM)
      real(dblprec) :: Temp_l !< Lattice temperature (for 3TM)
      real(dblprec) :: Temp_e !< Electron temperature (for 3TM)
      real(dblprec) :: t_in !< Time (for 3TM)

      ! Phase flag indicator (true for sd initial phase; false for sd measurement phase)
      ! Used in order to separate between phases in an adaptive time step environment with spin correlation
      logical :: sd_phaseflag

      real(dblprec) :: temprescale, temprescalegrad, dummy, totenergy

      ! Implicit midpoint methods

      ! Canonical velocity rescaling
      integer :: ndeg, i!, j
      real(dblprec) :: v2, kinenrg0, kinenrg1, kinrescaled, alpha

      real(dblprec) :: lattdamptol

      lattdampvec = lattdamp
      temprescale = 1.0_dblprec
      temprescalegrad = 0.0_dblprec

      ! Used to switch between regular velocity Verlet and the GJF Verlet.
      lattdamptol = 1e-17

      ! Spin correlation measurements allowed
      adapt_step = 0
      deltat_correction_flag = .true.

      ! Measurement phase indicator (false here)
      sd_phaseflag = .false.
      !call timing(0,'Measurement   ','ON')

      ! Flag for G(q) sampling and G(r) calculation
      cgk_flag=0 ; cgk_flag_p=0 ; cr_flag=0 ; cgk_flag_pc=0
      spt_flag=0

      ! Flag for G(q) sampling and G(r) calculation
      ucgk_flag=0 ; ucr_flag=0
      vcgk_flag=0 ; vcr_flag=0
      lcgk_flag=0 ; lcr_flag=0
      !uc_tidx=0  ; upt_flag=0

      scount_pulse = 1

      time_dept_flag=allocation_field_time(mwf,do_gauss,mwf_gauss,mov_gauss,        &
         mwf_mov_gauss,mwf_gauss_spatial,mov_circle,mwf_mov_circle,mov_square,      &
         mwf_mov_square)

      !call print_siminfo()
      if(mml_diag.and.do_mml>0) then
         print *,' Optimized diagonal MML evaluation enabled'
         do_mml=2
      end if

      ! Calculate demagnetization field
      if(demag=='Y') then
         call calc_demag(Natom, Mensemble, demag1, demag2, demag3, demagvol, emomM)
      endif

      call calc_bls(N1,N2,N3,C1,C2,C3,Natom,Mensemble,simid,coord,emomM,sstep,      &
         delta_t,0)

      !  Iniitialize magnetic field pulse
      if (do_bpulse.gt.0.and.do_bpulse.lt.5) then
         bpulse_time = delta_t*rstep
         call bpulse(do_bpulse, bpulse_time)
      end if

      ! Initialize spin temperature measurement
      if (do_spintemp=='Y') then
         ntmp=nstep/spintemp_step
         call spintemperature(Natom,Mensemble,mstep,ntmp,simid,emomM,beff,0)
      end if

      ! Initialize polarization including local chirality and polarity
      if (do_pol=='Y') then
         call init_polarization(Natom, Mensemble,ham%max_no_neigh,ham%nlist,coord,1)
      end if

      ! Calculate the external static fields, then can be calculated once and do not need to be calculated again
      call calc_external_fields(Natom,Mensemble,hfield,anumb,external_field,     &
         do_bpulse,sitefld,sitenatomfld)
      !
      ! Adaptive Time Stepping Region
      if(adaptive_time_flag) then
         call calculate_omegainit(omega_max, larmor_numrev, delta_t)
      end if

      ! Initialize 3TM functionality if enabled 
      if(do_3tm=='Y'.or.do_3tm=='E') then
!        call allocate_3tm(Natom,0,1)
         call init_3tm_cv(1)
         call unify_3tm_params(C1,C2,C3,alat,NA)
         if (do_3tm=='E') then
            call effective_field(Natom,Mensemble,1,Natom,emomM,mmom, &
               external_field,time_external_field,beff,beff1,beff2,OPT_flag,           &
               max_no_constellations,maxNoConstl,unitCellType,constlNCoup,             &
               constellations,constellationsNeighType,totenergy,Num_macro,   &
               cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)
            Temp_s = f_spintemp(Natom,Mensemble,emomM,beff)
            Temp_l = f_iontemp(Natom, Mensemble, mion, vvec)
            call threetemp_elec_init(simid,Temp,Temp_s,Temp_l,Temp_e)
         else 
            call set_initial_temp_3tm(Temp,Temp_s,Temp_l,Temp_e)
            call threetemp_print(rstep,nstep,delta_t,simid)
         end if
         print '(1x,a,3f14.5)', 'Initial 3TM temperatures:',Temp_s,Temp_l,Temp_e
      else
         Temp_s=Temp
         Temp_l=Temp
         Temp_e=Temp
      end if

      ! Rescaling of temperature according to Quantum Heat bath
      temprescale=1.0_dblprec
      temprescalegrad=0.0_dblprec
      if (do_qhb=="Q" .or. do_qhb=='R' .or. do_qhb=='P' .or. do_qhb=='T') then
         if(qhb_mode=='MT') then
            call calc_mavrg(Natom,Mensemble,emomM,mavg)
            call qhb_rescale(Temp_e,temprescale,temprescalegrad,do_qhb,qhb_mode,mavg)
         else
            call qhb_rescale(Temp_e,temprescale,temprescalegrad,do_qhb,qhb_mode,dummy)
         endif
      endif

      ! --Optimization Region-- !
      ! Allocation and initial opt build
      !if (OPT_ON) then
      !call timing(0,'Hamiltonian   ','OF')
      !   call timing(0,'BuildOptReg   ','ON')
      !   call buildOptimizationRegions(na,natom,Mensemble,max_no_neigh,&
      !        atype,emom,mmom,ncoup,nlist,nlistsize,cellPosNumNeigh,cos_thr)
      !   call timing(0,'BuildOptReg   ','OF')
      !   !call timing(0,'Hamiltonian   ','ON')
      !end if

      ! Initialize cumulant counter
      Navrgcum=0
      ! Perform nstep simulation steps
      ! Observe that Fortran cannot handle a regular do loop if nstep changes dynamically
      ! Hence, a do while loop has been introduced to allow for adaptive time stepping
      mstep=rstep+1

      do while (mstep.LE.rstep+nstep) !+1

         !call timing(0,'Hamiltonian   ','ON')
         if (do_3tm=='Y') then
            t_in=delta_t*mstep
            call threetemp_single(t_in,delta_t,Temp_s,Temp_l,Temp_e)
         else if (do_3tm=='E') then
            t_in=delta_t*mstep
            Temp_s = f_spintemp(Natom,Mensemble,emomM,beff)
            Temp_l = f_iontemp(Natom, Mensemble, mion, vvec)
            call threetemp_elec(t_in,delta_t,Temp_s,Temp_l,Temp_e)
            call threetemp_elec_print(mstep,simid,Temp_s,Temp_l,Temp_e)
         end if

         if (time_dept_flag) then
            ! Calculate Microwave fields (time dependent fields)
            call calculate_mwf_fields(Natom,mstep,rstep+nstep-1,delta_t,coord,0)

            ! Calculate the total time dependent fields
            call calc_external_time_fields(Natom, Mensemble,time_external_field,    &
               do_bpulse,demag,mwf,mwf_gauss_spatial,do_gauss,mwf_gauss,mov_gauss,  &
               mwf_mov_gauss,bpulsefield,demagfld,mwffield,gauss_mwffield,          &
               site_mwffield,gauss_spatial_site_mwffield,gauss_spatial,             &
               gauss_site_mwffield,mov_gauss_spatial,mwf_mov_gauss_spatial,         &
               mov_circle,mwf_mov_circle,mov_circle_spatial,mwf_mov_circle_spatial, &
               mov_square,mwf_mov_square,mov_square_spatial,mwf_mov_square_spatial)
         else
            time_external_field=0.0_dblprec
         endif

         ! --Optimization Region-- !
         if (OPT_ON) then
            if (mod(mstep-rstep,OPT_rebuild_time)==0) then
               !call timing(0,'Hamiltonian   ','OF')
               !call timing(0,'BuildOptReg   ','ON')
               call buildOptimizationRegions(na,natom,nHam,Mensemble,               &
                  ham%max_no_neigh,atype,emom,mmom,ham%ncoup,ham%nlist,             &
                  ham%nlistsize,ham%aham,cellPosNumNeigh,cos_thr)
               !call timing(0,'BuildOptReg   ','OF')
               !call timing(0,'Hamiltonian   ','ON')
            end if
         end if

         ! Apply Hamiltonian to obtain effective field
         call effective_field(Natom,Mensemble,1,Natom,emomM,mmom, &
            external_field,time_external_field,beff,beff1,beff2,OPT_flag,           &
            max_no_constellations,maxNoConstl,unitCellType,constlNCoup,             &
            constellations,constellationsNeighType,totenergy,Num_macro,   &
            cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)

         call effective_bmixedfield(Natom,Mensemble,1,Natom,do_ml,do_mml,   &
            uvec,emomM,beff,beff3)

         ! Apply lattice Hamiltonian to obtain effective field: First call
         call effective_latticefield(Natom,Mensemble,1,Natom,do_ll,  &
            do_ml,do_mml,mode,uvec,emomM,latt_external_field,               &
            latt_time_external_field,eeff,eeff1,eeff2,eeff3)

         !call timing(0,'Hamiltonian   ','OF')
         !call timing(0,'Measurement   ','ON')

         ! Measure averages and trajectories (only spin part)
         call measure(Natom,Mensemble,NT,NA,nHam,N1,N2,N3,simid,mstep,emom,emomM,   &
            mmom,Nchmax,do_ralloy,Natom_full,asite_ch,achem_ch,atype,plotenergy,    &
            Temp_s,temprescale,temprescalegrad,real_time_measure,delta_t,logsamp,     &
            ham%max_no_neigh,ham%nlist,                                             &
            ham%ncoup,ham%nlistsize,ham%aham,thermal_field,beff,beff1,beff3,coord,  &
            ham%ind_list_full,ham%ind_nlistsize,ham%ind_nlist,ham%max_no_neigh_ind, &
            ham%sus_ind,do_mom_legacy,mode)

         ! Calculate total and term resolved spin energies
         if(plotenergy>0) then
            if (mod(mstep-1,avrg_step)==0) then
               totenergy=0.0_dblprec
               call calc_energy(nHam,mstep,Natom,Nchmax,  &
                  conf_num,Mensemble,Natom,Num_macro,1,  &
                  plotenergy,Temp_s,delta_t,do_lsf, &
                  lsf_field,lsf_interpolate,real_time_measure,simid,cell_index,     &
                  macro_nlistsize,mmom,emom,emomM,emomM_macro,external_field,       &
                  time_external_field,max_no_constellations,maxNoConstl,            &
                  unitCellType,constlNCoup,constellations,OPT_flag,                 &
                  constellationsNeighType,totenergy,NA,N1,N2,N3)
            end if
         endif

         ! Calculate energies of the lattice and the spin-lattice Hamiltonians
         ! The energy of the magnetic system was calculated with the call to calc_energy
         sdpot_energy = Natom*totenergy
         call calc_lattenergies(Natom,Mensemble,1,Mensemble,1,Natom,do_ll,   &
            do_ml,do_mml,mode,uvec,emomM,latt_external_field,       &
            latt_time_external_field,ll_energy,ml_energy,    &
            mml_energy,ldpot_energy,sdpot_energy,sldpot_energy,         &
            totpot_energy,sld_single_energy,mm_energy0,ammom_inp,aemom_inp,NA)

         call local_angular_momentum(Natom, Mensemble, mion,  uvec, vvec, lvec)

         ! Measure lattice averages and trajectories
         call lattmeasure(simid,Natom,Mensemble,NT,NA,N1,N2,N3,mstep,logsamp,coord, &
            mion,uvec,vvec,lvec,eeff,eeff1,eeff3,ethermal_field,Nchmax,do_ralloy,        &
            Natom_full,asite_ch,achem_ch,atype,Temp_l,ldpot_energy,                   &
            sdpot_energy - mm_energy0,sldpot_energy,totpot_energy)

         ! Spin correlation
         ! Sample magnetic moments for correlation functions
         call correlation_wrapper(Natom,Mensemble,coord,simid,emomM,mstep,delta_t,  &
            NT,atype,Nchmax,achtype,sc,do_sc,do_sr,cgk_flag)

         ! Lattice correlation
         ! Sample displacements and velocities for correlation functions
         call correlation_wrapper(Natom,Mensemble,coord,simid,uvec,mstep,delta_t,  &
            NT,atype,Nchmax,achtype,uc,do_uc,do_ur,ucgk_flag)

         call correlation_wrapper(Natom,Mensemble,coord,simid,vvec,mstep,delta_t,  &
            NT,atype,Nchmax,achtype,vc,do_vc,do_vr,vcgk_flag)

         call correlation_wrapper(Natom,Mensemble,coord,simid,lvec,mstep,delta_t,  &
            NT,atype,Nchmax,achtype,lc,do_lc,do_vr,lcgk_flag)

         !call lattcorrelation_wrapper(Natom,Mensemble,coord,simid,uvec,vvec,mstep,  &
         !   delta_t,NT,atype,Nchmax,achtype,ucgk_flag,ucgk_flag_p)


         call calc_bls(N1,N2,N3,C1,C2,C3,Natom,Mensemble,simid,coord,emomM,mstep,   &
            delta_t,1)

         !call timing(0,'SpinCorr      ','OF')
         !call timing(0,'Measurement   ','ON')

         ! Write simulation status for each 5% of the simulation length
         if(nstep>20) then
            if(mod(mstep,(rstep+nstep)/20)==0) then
               call calc_mavrg(Natom,Mensemble,emomM,mavg)
               ! Only calculate binderc if it is not calculated elsewhere
               ! If calculated elsewhere, display binderc will not be real binderc
               if(do_cumu=='N') then
                  call calc_and_print_cumulant(Natom,Mensemble,emomM,simid,Temp_s,    &
                  temprescale,temprescalegrad,plotenergy,cumu_buff,.false.)
               endif
               write(*,'(2x,a,i4,a,f10.6)',advance='no') &
               "MP",mstep*100/(rstep+nstep),"% done. Mbar:",mavg
               if(plotenergy>0) then
                  write(*,'(a,f12.6)',advance='no') ". Ebar:", totenergy
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
                  call qhb_rescale(Temp_s,temprescale,temprescalegrad,do_qhb,qhb_mode,mavg)
               endif
            endif
         endif
         !call timing(0,'Measurement   ','OF')
         !call timing(0,'Hamiltonian   ','ON')

         ! Calculate spin transfer torque
         if(stt=='A'.or.stt=='F'.or.do_she=='Y'.or.do_sot=='Y') then
            call calculate_spintorques(Natom, Mensemble,lambda1_array,emom,mmom)
         end if

         !call timing(0,'Hamiltonian   ','OF')
         call timing(0,'Evolution     ','ON')

         ! Here both the B-field and the E-field have been calculated for the first time of the time-step
         ! Note that the E-field was calculated early in the loop over mstep to have fields
         ! available for lattmeasurements

         ! Calculate average field for use with multiple heat baths
         if (llg==2.or.llg==3) then 
            call calc_averagefield(Natom,Mensemble,beff1,beff2,field1,field2)
         endif
         ! Check if this changes
         thermal_field=0.0_dblprec

         ! Set temperature to spin temperature for 3TM simulationa
         ! TODO: Change to spatially resolved temperatures
         if(do_3tm=='Y') Temp_array=Temp_s
         if(do_3tm=='E') Temp_array=Temp_e

         !
         call evolve_first(Natom,Mensemble,Landeg,llg,SDEalgh,bn,lambda1_array,     &
            lambda2_array,NA,compensate_drift,delta_t,relaxtime,Temp_array,         &
            temprescale,beff,b2eff,thermal_field,beff2,btorque,field1,field2,emom,  &
            emom2,emomM,mmom,mmomi,stt,do_site_damping,ham%nlist,ham%nlistsize,     &
            constellationsUnitVec,constellationsUnitVec2,constellationsMag,         &
            constellations,unitCellType,OPT_flag,cos_thr,max_no_constellations,     &
            do_she,she_btorque,Nred,red_atom_list,do_sot,sot_btorque)

         !If needed, generate stochastic forces
         !if(lattdamp .gt. lattdamptol) then
         if(Temp_l>0.0_dblprec) then
            ! Set temperature to spin temperature for 3TM simulationa
            ! TODO: Change to spatially resolved temperatures
            if(do_3tm=='Y') Temp_array=Temp_l
            if(do_3tm=='E') Temp_array=Temp_e

            call lattrannum(Natom,Mensemble,NA,delta_t,lattdampvec,compensate_drift,&
               Temp_array,temprescale)
         else
            lattranv = 0_dblprec
         end if

         ! Update ionic displacements
         !if(lattdamp .lt. lattdamptol) then
         !   call u_vverlet(Natom, Mensemble, mioninv, eeff, uvec, vvec, delta_t)
         !else
            call u_gjfverlet(Natom,Mensemble,mioninv,eeff,uvec,vvec,delta_t,lattdamp)
         !end if

         call timing(0,'Evolution     ','OF')
         !call timing(0,'Hamiltonian   ','ON')

         ! This if-statement seems unnecessary, the if(SDEalgh<=20) above should suffice
         if(SDEalgh==1.or.SDEalgh==4.or.SDEalgh==5.or.SDEalgh==6.or.SDEalgh==7.or.SDEalgh==11) then

            if (time_dept_flag) then
               ! Calculate Microwave fields (time dependent fields)
               call calculate_mwf_fields(Natom,mstep,rstep+nstep-1,delta_t,coord,1)

               ! Calculate the total time dependent fields
               call calc_external_time_fields(Natom,Mensemble,time_external_field,  &
                  do_bpulse,demag,mwf,mwf_gauss_spatial,do_gauss,mwf_gauss,         &
                  mov_gauss,mwf_mov_gauss,bpulsefield,demagfld,                     &
                  mwffield,gauss_mwffield,site_mwffield,gauss_spatial_site_mwffield,&
                  gauss_spatial,gauss_site_mwffield,mov_gauss_spatial,              &
                  mwf_mov_gauss_spatial,mov_circle,mwf_mov_circle,                  &
                  mov_circle_spatial,mwf_mov_circle_spatial,mov_square,             &
                  mwf_mov_square,mov_square_spatial,mwf_mov_square_spatial)
            else
               time_external_field=0.0_dblprec
            endif

            ! Apply Hamiltonian to obtain effective field
            call effective_field(Natom,Mensemble,1,Natom,emomM,mmom, &
               external_field,time_external_field,beff,beff1,beff2,OPT_flag,        &
               max_no_constellations,maxNoConstl,unitCellType,constlNCoup,          &
               constellations,constellationsNeighType,totenergy,Num_macro,&
               cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)

            call effective_bmixedfield(Natom,Mensemble,1,Natom,do_ml,do_mml,&
               uvec,emomM,beff,beff3)

         end if
         !call timing(0,'Hamiltonian   ','OF')
         !call timing(0,'Evolution     ','ON')

         ! Re-Calculate spin transfer torque
         if(stt=='A'.or.stt=='F'.or.do_she=='Y'.or.do_sot=='Y') then
            call calculate_spintorques(Natom, Mensemble,lambda1_array,emom,mmom)
         end if

         ! Perform second (corrector) step of SDE solver
         ! SLDTODO check how this choice of having the second step of spin evolution occuring
         ! SLDTODO before the second call to effective_latticefield performs in practice
         call evolve_second(Natom, Mensemble,Landeg,llg,SDEalgh,bn,lambda1_array,   &
            delta_t,relaxtime,beff,beff2,b2eff,btorque,emom,emom2,stt,ham%nlist,    &
            ham%nlistsize,constellationsUnitVec,constellationsUnitVec2,             &
            constellationsMag,constellations,unitCellType,OPT_flag,cos_thr,         &
            max_no_constellations,do_she,she_btorque,Nred,red_atom_list,            &
            do_sot,sot_btorque)

         !call timing(0,'Evolution     ','OF')
         !call timing(0,'Hamiltonian   ','ON')

         ! Apply lattice Hamiltonian to obtain effective field: Second call
         call effective_latticefield(Natom, Mensemble, 1, Natom, &
         do_ll, do_ml, do_mml, mode, &
         uvec, emomM, latt_external_field, latt_time_external_field, &
         e2eff, eeff1, eeff2, eeff3)

         !call timing(0,'Hamiltonian   ','OF')
         !call timing(0,'Evolution     ','ON')

         ! Update ionic velocities
         !if(lattdamp .lt. lattdamptol) then
         !  call v_vverlet(Natom, Mensemble, mioninv, eeff, e2eff, vvec, delta_t)
         !else
            call v_gjfverlet(Natom,Mensemble,mioninv,eeff,e2eff,uvec,uvec2,vvec,    &
               delta_t,lattdamp)
            ! Commented out since check for reset of u an v is done further down.
            !call set_avrgu0(Natom, Mensemble, uvec, vvec)
            !call set_avrgp0(Natom, Mensemble, uvec, vvec)
         !end if

         ! Update magnetic moments after time evolution step
         call moment_update(Natom,Mensemble,mmom,mmom0,mmom2,emom,emom2,emomM,mmomi,& 
            mompar,initexc)
         !call timing(0,'Moments       ','OF')
         !call timing(0,'Measurement   ','ON')

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

         !!!AB restruc reinstate later
         !!!if (adaptive_time_flag) then
         !!!   call adapt_time_step(Natom,Mensemble,beff,omega_max,larmor_numrev,larmor_thr,rstep,mstep,nstep,totalsimtime,&
         !!!        thermal_field,do_sc,sc_step,sc_tidx,sd_phaseflag,adapt_time_interval,&
         !!!        adapt_step,adaptive_time_flag,deltat_correction_flag,delta_t)
         !!!endif

         mstep = mstep + 1

         ! Canonical velocity rescaling
         if ( do_velrsc=='Y' .and. mod(mstep-1,velrsc_step)==0 ) then
            ndeg = 3 * Natom - 3
            do k=1,Mensemble
               kinenrg0 = 0.0_dblprec
               do i=1, Natom
                  v2 = vvec(1,i,k)*vvec(1,i,k) + vvec(2,i,k)*vvec(2,i,k) + vvec(3,i,k)*vvec(3,i,k)
                  kinenrg0 = kinenrg0 + 0.5_dblprec * amu * mion(i,1) * v2 * angstrom**2 / ev
               end do
               kinenrg1 = ndeg * k_bolt * Temp_l / 2 / ev
               velrsc_taut=10_dblprec
               kinrescaled = resamplekin(kinenrg0, kinenrg1, ndeg, velrsc_taut)
               alpha = sqrt(kinrescaled/kinenrg0)

               ! Rescale velocities
               call rescale_vel(Natom, Mensemble, uvec, vvec, alpha)

               ! Set average linear momenta to zero
               call set_avrgp0(Natom, Mensemble, uvec, vvec)

               ! Set average displacement to zero
               call set_avrgu0(Natom, Mensemble, uvec, vvec)

            end do
         end if

         if(do_set_avrgp0 == 'Y') then
            ! Set average linear momenta to zero
            call set_avrgp0(Natom, Mensemble, uvec, vvec)
         end if

         if(do_set_avrgu0 == 'Y') then
            ! Set average displacement to zero
            call set_avrgu0(Natom, Mensemble, uvec, vvec)
         end if

      end do ! End loop over simulation steps

      !!!AB restruc reinstate later
      !!!if(do_sc=='C'.and.adaptive_time_flag) then
      !!!   write(*,*) 'Fraction of spin correlation samples relative to input specification:', 100*sc_tidx/sc_nstep, '%'
      !!!end if

      ! Measure averages and trajectories
      call measure(Natom,Mensemble,NT,NA,nHam,N1,N2,N3,simid,mstep,emom,emomM,mmom, &
         Nchmax,do_ralloy,Natom_full,asite_ch,achem_ch,atype,plotenergy,Temp_s,       &
         temprescale,temprescalegrad,real_time_measure,delta_t,logsamp,             &
         ham%max_no_neigh,ham%nlist,                                                &
         ham%ncoup,ham%nlistsize,ham%aham,thermal_field,beff,beff1,beff3,coord,     &
         ham%ind_list_full,ham%ind_nlistsize,ham%ind_nlist,ham%max_no_neigh_ind,    &
         ham%sus_ind,do_mom_legacy,mode)

      ! Print remaining measurements
      call flush_measurements(Natom,Mensemble,NT,NA,N1,N2,N3,simid,mstep,emom,mmom, &
         Nchmax,atype,real_time_measure,rstep+nstep,ham%ind_list_full,do_mom_legacy,&
         mode)

      ! Calculate total and term resolved spin energies
      if(plotenergy>0) then
         if (mod(mstep-1,avrg_step)==0) then
            totenergy=0.0_dblprec
            call calc_energy(nHam,mstep,Natom,Nchmax, &
               conf_num,Mensemble,Natom,Num_macro,1,     &
               plotenergy,Temp_s,delta_t,do_lsf,    &
               lsf_field,lsf_interpolate,real_time_measure,simid,cell_index,        &
               macro_nlistsize,mmom,emom,emomM,emomM_macro,external_field,          &
               time_external_field,max_no_constellations,maxNoConstl,unitCellType,  &
               constlNCoup,constellations,OPT_flag,constellationsNeighType,         &
               totenergy,NA,N1,N2,N3)
         end if
      endif

      ! Apply lattice Hamiltonian to obtain effective field
      call effective_latticefield(Natom,Mensemble,1,Natom,do_ll,     &
         do_ml,do_mml,mode,uvec,emomM,latt_external_field,                  &
         latt_time_external_field,eeff,eeff1,eeff2,eeff3)

      ! Calculate energies of the lattice and the spin-lattice Hamiltonians
      ! The energy of the magnetic system was calculated with the call to calc_energy
      sdpot_energy = Natom*totenergy
      call calc_lattenergies(Natom,Mensemble,1,Mensemble,1,Natom,do_ll,      &
         do_ml,do_mml,mode,uvec,emomM,latt_external_field,          &
         latt_time_external_field,ll_energy,ml_energy,       &
         mml_energy,ldpot_energy,sdpot_energy,sldpot_energy,            &
         totpot_energy,sld_single_energy,mm_energy0,ammom_inp,aemom_inp,NA)

      call local_angular_momentum(Natom, Mensemble, mion,  uvec, vvec, lvec)

      ! Measure averages and trajectories
      call lattmeasure(simid,Natom,Mensemble,NT,NA,N1,N2,N3,mstep,logsamp,coord,    &
         mion,uvec,vvec,lvec,eeff,eeff1,eeff3,ethermal_field,Nchmax,do_ralloy,Natom_full,&
         asite_ch,achem_ch,atype,Temp_l,ldpot_energy,sdpot_energy - mm_energy0,       &
         sldpot_energy, totpot_energy)

      ! Print remaining measurements
      call flush_lattmeasurements(simid, Natom, Mensemble, NT, NA, N1, N2, N3, &
      mstep, mion, uvec, vvec, coord, Nchmax, atype, Temp_l)

      ! Print final polarization, chirality and local polarization
      if (do_pol=='Y') then
         call init_polarization(Natom,Mensemble,ham%max_no_neigh,ham%nlist,coord,-1)
      end if

      if (do_spintemp=='Y') then
         call spintemperature(Natom,Mensemble,mstep,1,simid,emomM,beff,2)
      endif
      !call timing(0,'Measurement   ','OF')

      ! Deallocate 3TM arrays
      if(do_3tm=='Y') then
         call init_3tm_cv(-1)
!        call allocate_3tm(Natom,0,-1)
      end if

   end subroutine sld_mphase


   !  Variant that uses fix point iteration of the implicit midpoint
   !---------------------------------------------------------------------------
   !> @brief
   !> Spin Lattice Dynamics measurement phase
   !
   !> @author
   !> Johan Hellsvik
   !---------------------------------------------------------------------------
   subroutine sld_mphase3()
      !
      use BLS,                only : calc_bls
      use QHB,                only : qhb_rescale,do_qhb,qhb_mode
      use Verlet
      use Energy,             only : calc_energy
      use Damping
      use FixedMom
      use Gradients
      use Evolution
      use InputData
      use FieldData,          only : beff,beff1,beff2,beff3,sitefld,          &
         external_field,field1,field2,time_external_field,allocation_field_time,    &
         thermal_field
      use DemagField
      use MomentData
      use macrocells
      use FieldPulse
      use SystemData,         only : atype, anumb, Landeg,coord
      use Correlation
      use Temperature
      use LatticeData
      use SpinTorques
      use Temperature
      use SphMidpoint
      use ChemicalData
      use Measurements
      use Polarization
      use prn_averages
      use RandomNumbers
      use UpdateMoments
      use MicroWaveField
      use Math_functions, only : f_logstep
      use SimulationData
      use HamiltonianData
      use VelocityRescale
      use CalculateFields
      use AutoCorrelation,    only : autocorr_sample,do_autocorr
      use LatticeInputData
      use LatticeFieldData
!     use LatticeCorrelation
      use HamiltonianActions
      use LatticeMeasurements
      use AdaptiveTimeStepping
      use optimizationRoutines
      use LatticeHamiltonianData
      use LatticeHamiltonianActions

      implicit none

      logical :: time_dept_flag,deltat_correction_flag
      integer :: cgk_flag,cgk_flag_p,adapt_step,cr_flag,cgk_flag_pc,spt_flag,ntmp
      integer :: ucgk_flag, ucr_flag
      integer :: vcgk_flag, vcr_flag, lcgk_flag, lcr_flag
      integer :: k, scount_pulse, sstep

      real(dblprec), dimension(Natom) :: lattdampvec
      real(dblprec), dimension(Mensemble) :: ll_energy
      real(dblprec), dimension(Mensemble) :: ml_energy
      real(dblprec), dimension(Mensemble) :: mml_energy
      real(dblprec), dimension(Mensemble) :: ldpot_energy
      real(dblprec), dimension(Mensemble) :: sdpot_energy
      real(dblprec), dimension(Mensemble) :: sldpot_energy
      real(dblprec), dimension(Mensemble) :: totpot_energy
      real(dblprec), dimension(Mensemble) :: sld_single_energy
      real(dblprec), dimension(Mensemble) :: mm_energy0

      ! Phase flag indicator (true for sd initial phase; false for sd measurement phase)
      ! Used in order to separate between phases in an adaptive time step environment with spin correlation
      logical :: sd_phaseflag

      real(dblprec) :: temprescale, temprescalegrad,dummy, denergy, totenergy

      ! Implicit midpoint methods
      integer :: imp_count
      logical :: converged
      logical :: interrupt

      ! Canonical velocity rescaling
      integer :: ndeg, i!, j
      real(dblprec) :: v2, kinenrg0, kinenrg1, kinrescaled, alpha

      real(dblprec) :: lattdamptol

      real(dblprec) :: lambdatol = 1d-12

      lattdampvec = lattdamp
      temprescale = 1.0_dblprec
      temprescalegrad = 0.0_dblprec

      ! Used to switch between microcanonical LD or Langevin LD
      lattdamptol = 1e-17

      ! Spin correlation measurements allowed
      adapt_step = 0
      deltat_correction_flag = .true.

      ! Measurement phase indicator (false here)
      sd_phaseflag = .false.
      !call timing(0,'Measurement   ','ON')

      ! Flag for G(q) sampling and G(r) calculation
      cgk_flag=0 ; cgk_flag_p=0 ; cr_flag=0 ; cgk_flag_pc=0
      spt_flag=0

      ! Flag for G(q) sampling and G(r) calculation
      ucgk_flag=0 ; ucr_flag=0
      vcgk_flag=0 ; vcr_flag=0
      lcgk_flag=0 ; lcr_flag=0
      !uc_tidx=0  ; upt_flag=0

      scount_pulse = 1

      time_dept_flag=allocation_field_time(mwf,do_gauss,mwf_gauss,mov_gauss,        &
         mwf_mov_gauss,mwf_gauss_spatial,mov_circle,mwf_mov_circle,mov_square,      &
         mwf_mov_square)

      !call print_siminfo()
      if(mml_diag.and.do_mml>0) then
         do_mml=2
         print *,' Optimized diagonal MML evaluation enabled', do_mml
      end if

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

      ! Calculate the external static fields, then can be calculated once and do not need to be calculated again
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

      ! Allocate fix point iteration work arrays
      call allocate_emomfp(Natom,Mensemble,1)

      ! Before time step, copy emom, uvec, and vvec
      ! to fix point iteration work arrays
      emom2=emom
      emoma=emom
      emomb=emom

      uvec2=uvec
      uveca=uvec
      uvecb=uvec

      vvec2=vvec
      vveca=vvec
      vvecb=vvec

      ! Initialize cumulant counter
      Navrgcum=0
      ! Perform nstep simulation steps
      ! Observe that Fortran cannot handle a regular do loop if nstep changes dynamically
      ! Hence, a do while loop has been introduced to allow for adaptive time stepping
      mstep=rstep+1

      do while (mstep.LE.rstep+nstep) !+1

         if (time_dept_flag) then
            ! Calculate Microwave fields (time dependent fields)
            call calculate_mwf_fields(Natom,mstep,rstep+nstep-1,delta_t,coord,0)

            ! Calculate the total time dependent fields
            call calc_external_time_fields(Natom,Mensemble,time_external_field,     &
               do_bpulse,demag,mwf,mwf_gauss_spatial,do_gauss,mwf_gauss,mov_gauss,  &
               mwf_mov_gauss,bpulsefield,demagfld,mwffield,gauss_mwffield,          &
               site_mwffield,gauss_spatial_site_mwffield,gauss_spatial,             &
               gauss_site_mwffield,mov_gauss_spatial,mwf_mov_gauss_spatial,         &
               mov_circle,mwf_mov_circle,mov_circle_spatial,mwf_mov_circle_spatial, &
               mov_square,mwf_mov_square,mov_square_spatial,mwf_mov_square_spatial)
         else
            time_external_field=0.0_dblprec
         endif

         ! Calculate spin transfer torque contributions to the local field
         if(stt=='A'.or.stt=='F'.or.do_she=='Y'.or.do_sot=='Y') then
            call calculate_spintorques(Natom, Mensemble,lambda1_array,emom,mmom)
         end if
         ! Apply Hamiltonian to obtain effective field
         call effective_field(Natom,Mensemble,1,Natom,emomM,mmom, &
            external_field,time_external_field,beff,beff1,beff2,OPT_flag,           &
            max_no_constellations,maxNoConstl,unitCellType,constlNCoup,             &
            constellations,constellationsNeighType,denergy,Num_macro,     &
            cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)

         call effective_bmixedfield(Natom,Mensemble,1,Natom,do_ml,do_mml,   &
            uvec,emomM,beff,beff3)

         ! Apply lattice Hamiltonian to obtain effective field: First call
         call effective_latticefield(Natom,Mensemble,1,Natom,do_ll,  &
            do_ml,do_mml,mode,uvec,emomM,latt_external_field,               &
            latt_time_external_field,eeff,eeff1,eeff2,eeff3)

         ! Measure averages and trajectories
         call measure(Natom,Mensemble,NT,NA,nHam,N1,N2,N3,simid,mstep,emom,emomM,   &
            mmom,Nchmax,do_ralloy,Natom_full,asite_ch,achem_ch,atype,plotenergy,    &
            Temp,temprescale,temprescalegrad,real_time_measure,delta_t,logsamp,     &
            ham%max_no_neigh,ham%nlist,&
            ham%ncoup,ham%nlistsize,ham%aham,thermal_field,beff,beff1,beff3,coord,  &
            ham%ind_list_full,ham%ind_nlistsize,ham%ind_nlist,ham%max_no_neigh_ind, &
            ham%sus_ind,do_mom_legacy,mode)

         ! Calculate total and term resolved spin energies
         if(plotenergy>0) then
            if (mod(mstep-1,avrg_step)==0) then
               totenergy=0.0_dblprec
               call calc_energy(nHam,mstep,Natom,Nchmax,  &
                  conf_num,Mensemble,nHam,Num_macro,1,   &
                  plotenergy,Temp,delta_t,do_lsf, &
                  lsf_field,lsf_interpolate,real_time_measure,simid,cell_index,     &
                  macro_nlistsize,mmom,emom,emomM,emomM_macro,external_field,       &
                  time_external_field,max_no_constellations,maxNoConstl,            &
                  unitCellType,constlNCoup,constellations,OPT_flag,                 &
                  constellationsNeighType,totenergy,NA,N1,N2,N3)
            end if
         endif

         ! Calculate energies of the lattice and the spin-lattice Hamiltonians
         ! The energy of the magnetic system was calculated with the call to calc_energy
         sdpot_energy = Natom*totenergy

         call calc_lattenergies(Natom,Mensemble,1,Mensemble,1,Natom,do_ll,   &
            do_ml,do_mml,mode,uvec,emomM,latt_external_field,       &
            latt_time_external_field,ll_energy,ml_energy,    &
            mml_energy,ldpot_energy,sdpot_energy,sldpot_energy,         &
            totpot_energy,sld_single_energy,mm_energy0,ammom_inp,aemom_inp,NA)

         call local_angular_momentum(Natom, Mensemble, mion,  uvec, vvec, lvec)

         ! Measure lattice averages and trajectories
         call lattmeasure(simid,Natom,Mensemble,NT,NA,N1,N2,N3,mstep,logsamp,coord, &
            mion,uvec,vvec,lvec,eeff,eeff1,eeff3,ethermal_field,Nchmax,do_ralloy,        &
            Natom_full,asite_ch,achem_ch,atype,Temp,ldpot_energy,                   &
            sdpot_energy - mm_energy0,sldpot_energy,totpot_energy)

         ! Spin correlation
         ! Sample magnetic moments for correlation functions
         call correlation_wrapper(Natom,Mensemble,coord,simid,emomM,mstep,delta_t,  &
            NT,atype,Nchmax,achtype,sc,do_sc,do_sr,cgk_flag)

         ! Lattice correlation
         ! Sample displacements and velocities for correlation functions
         call correlation_wrapper(Natom,Mensemble,coord,simid,uvec,mstep,delta_t,  &
            NT,atype,Nchmax,achtype,uc,do_uc,do_ur,ucgk_flag)

         call correlation_wrapper(Natom,Mensemble,coord,simid,vvec,mstep,delta_t,  &
            NT,atype,Nchmax,achtype,vc,do_vc,do_vr,vcgk_flag)

         call correlation_wrapper(Natom,Mensemble,coord,simid,lvec,mstep,delta_t,  &
            NT,atype,Nchmax,achtype,lc,do_lc,do_vr,lcgk_flag)

         !call lattcorrelation_wrapper(Natom,Mensemble,coord,simid,uvec,vvec,mstep,  &
         !   delta_t,NT,atype,Nchmax,achtype,ucgk_flag,ucgk_flag_p)

         call calc_bls(N1,N2,N3,C1,C2,C3,Natom,Mensemble,simid,coord,emomM,mstep,   &
            delta_t,1)

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
               write(*,'(2x,a,i4,a,f10.6)',advance='no')                            &
                  "MP",mstep*100/(rstep+nstep),"% done. Mbar:",mavg
               if(plotenergy>0) then
                  write(*,'(a,f12.6)',advance='no') ". Ebar:", totenergy
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

         if(stt=='A'.or.stt=='F'.or.do_she=='Y') then
            call calculate_spintorques(Natom, Mensemble,lambda1_array,emom,mmom)
         end if

         call timing(0,'Evolution     ','ON')

         ! Here both the B-field and the E-field have been calculated for the first time of the time-step
         ! Note that the E-field was calculated early in the loop over mstep to have fields
         ! available for lattmeasurements

         ! Fix-point iteration until convergence.
         imp_count = 0
         converged = .false.
         interrupt = .false.

         ! If the damping is negliably small ( < lambdatol ), there is
         ! no need to generate random numbers
         if (do_site_damping=='Y') then
            call rannum(Natom,Mensemble,NA,llg,lambda1_array,lambda2_array,         &
               compensate_drift,bn,field1,field2,mmomi,Temp_array,temprescale)
         else
            if(minval(lambda1_array).gt.lambdatol.or.minval(lambda2_array).gt.lambdatol) then
               call rannum(Natom,Mensemble,NA,llg,lambda1_array,lambda2_array,      &
               compensate_drift,bn,field1,field2,mmomi,Temp_array,temprescale)
            else
               !write(*,*) 'Damping zero, no rng'
               ranv=0_dblprec
            end if
         endif

         !If needed, generate stochastic forces
         !if(lattdamp .gt. lattdamptol) then
         if(Temp>0.0_dblprec) then
            call lattrannum(Natom,Mensemble,NA,delta_t*0.5_dblprec,lattdampvec,     &
               compensate_drift,Temp_array,temprescale)
         else
            lattranv = 0_dblprec
         end if

         do

            imp_count = imp_count + 1
            !print *,imp_count,imp_max_count
            if(imp_count .gt. imp_max_count) interrupt = .true.
            !if(imp_count .gt. imp_max_count) exit

            !write(*,*) 'imp_epsilon ', imp_epsilon, ' imp_max_count ', imp_max_count

            call imp_fp(Natom,Mensemble,Landeg,bn,lambda1_array,beff,emom,emom2,    &
               emomM,mmom,delta_t,SDEalgh,imp_epsilon,converged,      &
               interrupt,mioninv,eeff,uvec,uvec2,vvec,vvec2)

            !write(*,*) 'Fix point iteration nr: ', imp_count
            if (converged .or. interrupt) exit
            !if (converged) exit

            call effective_field(Natom,Mensemble,1,Natom,emomM,mmom, &
               external_field,time_external_field,beff,beff1,beff2,OPT_flag,        &
               max_no_constellations,maxNoConstl,unitCellType,constlNCoup,          &
               constellations,constellationsNeighType,denergy,Num_macro,  &
               cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)

            ! Calculate spin transfer torque contributions to the local field
            if(stt=='A'.or.stt=='F'.or.do_she=='Y'.or.do_sot=='Y') then
               call calculate_spintorques(Natom, Mensemble,lambda1_array,emom,mmom)
            end if

            call effective_bmixedfield(Natom,Mensemble,1,Natom,do_ml,do_mml,&
               uvec2,emomM,beff,beff3)

            call effective_latticefield(Natom,Mensemble,1,Natom,do_ll,       &
               do_ml,do_mml,mode,uvec2,emomM,latt_external_field,   &
               latt_time_external_field,eeff,eeff1,eeff2,eeff3)

         end do
         !if (converged) then
         !   write(*,'(a,i8,a,i8,a)') 'mstep ', mstep, ' fix point iteration nr: ', imp_count, ' converged.'
         !else
         !   write(*,'(a,i8,a,i8,a)') 'mstep ', mstep, ' fix point iteration nr: ', imp_count, ' Not converged.'
         !end if

         call timing(0,'Evolution     ','OF')

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

         ! Canonical velocity rescaling
         if ( do_velrsc=='Y' .and. mod(mstep-1,velrsc_step)==0 ) then
            ndeg = 3 * Natom - 3
            do k=1,Mensemble
               kinenrg0 = 0.0_dblprec
               do i=1, Natom
                  v2 = vvec(1,i,k)*vvec(1,i,k)+vvec(2,i,k)*vvec(2,i,k)+vvec(3,i,k)*vvec(3,i,k)
                  kinenrg0 = kinenrg0+0.5_dblprec*amu*mion(i,1)*v2*angstrom**2/ev
               end do
               kinenrg1 = ndeg * k_bolt * temp / 2 / ev
               velrsc_taut=10_dblprec
               kinrescaled = resamplekin(kinenrg0, kinenrg1, ndeg, velrsc_taut)
               alpha = sqrt(kinrescaled/kinenrg0)

               ! Rescale velocities
               call rescale_vel(Natom, Mensemble, uvec, vvec, alpha)

               ! Set average linear momenta to zero
               call set_avrgp0(Natom, Mensemble, uvec, vvec)

               ! Set average displacement to zero
               call set_avrgu0(Natom, Mensemble, uvec, vvec)

            end do

            if(do_set_avrgp0 == 'Y') then
               ! Set average linear momenta to zero
               call set_avrgp0(Natom, Mensemble, uvec, vvec)
            end if

            if(do_set_avrgu0 == 'Y') then
               ! Set average displacement to zero
               call set_avrgu0(Natom, Mensemble, uvec, vvec)
            end if
         end if

      end do ! End loop over simulation steps

      ! Measure averages and trajectories
      call measure(Natom,Mensemble,NT,NA,nHam,N1,N2,N3,simid,mstep,emom,emomM,mmom, &
         Nchmax,do_ralloy,Natom_full,asite_ch,achem_ch,atype,plotenergy,Temp,temprescale,&
         temprescalegrad,real_time_measure,delta_t,logsamp,ham%max_no_neigh,ham%nlist,ham%ncoup,    &
         ham%nlistsize,ham%aham,thermal_field,beff,beff1,beff3,coord,               &
         ham%ind_list_full,ham%ind_nlistsize,ham%ind_nlist,ham%max_no_neigh_ind,    &
         ham%sus_ind,do_mom_legacy,mode)

      ! Print remaining measurements
      call flush_measurements(Natom,Mensemble,NT,NA,N1,N2,N3,simid,mstep,emom,mmom, &
         Nchmax,atype,real_time_measure,rstep+nstep,ham%ind_list_full,do_mom_legacy,&
         mode)

      ! Calculate total and term resolved spin energies
      if(plotenergy>0) then
         if (mod(mstep-1,avrg_step)==0) then
            totenergy=0.0_dblprec
            call calc_energy(nHam,mstep,Natom,Nchmax, &
               conf_num,Mensemble,nHam,Num_macro,1,      &
               plotenergy,Temp,delta_t,do_lsf,    &
               lsf_field,lsf_interpolate,real_time_measure,simid,cell_index,        &
               macro_nlistsize,mmom,emom,emomM,emomM_macro,external_field,          &
               time_external_field,max_no_constellations,maxNoConstl,unitCellType,  &
               constlNCoup,constellations,OPT_flag,constellationsNeighType,         &
               totenergy,NA,N1,N2,N3)
         end if
      endif

      ! Apply lattice Hamiltonian to obtain effective field
      call effective_latticefield(Natom,Mensemble,1,Natom,do_ll,     &
         do_ml,do_mml,mode,uvec,emomM,latt_external_field,                  &
         latt_time_external_field,eeff,eeff1,eeff2,eeff3)

      ! Calculate energies of the lattice and the spin-lattice Hamiltonians
      ! The energy of the magnetic system was calculated with the call to calc_energy
      sdpot_energy = Natom*totenergy
      call calc_lattenergies(Natom,Mensemble,1,Mensemble,1,Natom,do_ll,      &
         do_ml,do_mml,mode,uvec,emomM,latt_external_field,          &
         latt_time_external_field,ll_energy,ml_energy,       &
         mml_energy,ldpot_energy,sdpot_energy,sldpot_energy,            &
         totpot_energy,sld_single_energy,mm_energy0,ammom_inp,aemom_inp,NA)

      call local_angular_momentum(Natom, Mensemble, mion,  uvec, vvec, lvec)

      ! Measure averages and trajectories
      call lattmeasure(simid,Natom,Mensemble,NT,NA,N1,N2,N3,mstep,logsamp,coord,    &
         mion,uvec,vvec,lvec,eeff,eeff1,eeff3,ethermal_field,Nchmax,do_ralloy,Natom_full,&
         asite_ch,achem_ch,atype,Temp,ldpot_energy,sdpot_energy - mm_energy0,       &
         sldpot_energy,totpot_energy)

      ! Print remaining measurements
      call flush_lattmeasurements(simid,Natom,Mensemble,NT,NA,N1,N2,N3,mstep,mion,  &
         uvec,vvec,coord,Nchmax,atype,Temp)

      ! Print final polarization, chirality and local polarization
      if (do_pol=='Y') then
         call init_polarization(Natom,Mensemble,ham%max_no_neigh,ham%nlist,coord,1)
      end if

      if (do_spintemp=='Y') then
         call spintemperature(Natom,Mensemble,mstep,1,simid,emomM,beff,2)
      endif

      ! Deallocate fix point iteration work arrays
      call allocate_emomfp(Natom,Mensemble,-1)

   end subroutine sld_mphase3


end module sld_driver
