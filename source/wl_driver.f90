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
module wl_driver

   implicit none

   public

contains


   !---------------------------------------------------------------------------
   !> @brief
   !> Wang-Landau initial phase
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine wl_iphase()
      !
      use Energy,          only : calc_energy
      use InputData
      use FieldData,       only : external_field, sitefld, time_external_field,&
         allocation_field_time, thermal_field, beff, beff1, beff2
      use SystemData
      use MonteCarlo
      use DemagField
      use MomentData
      use FieldPulse
      !use prn_averages
      use ChemicalData,    only : achem_ch, asite_ch
      use MicroWaveField
      use CalculateFields
      use HamiltonianData
      use ChemicalData, only : achtype
      use macrocells
      use optimizationRoutines
      use HamiltonianActions
      use WangLandau
      use MonteCarlo_common, only : randomize_spins
      use omp_lib

      !
      implicit none
      !
      integer :: cgk_flag, scount_pulse, bcgk_flag, cgk_flag_pc, mcmstep, wl_count
      integer :: ii, ipmcstep, iloop
      character(len=30) :: filn


      !call timing(0,'Initial       ','ON')
      ! Allocate work arrays for MC
      call allocate_mcdata(Natom,1)

      do iloop=1,2

         wl_direction=(-1.0_dblprec)**iloop
         wl_totenergy=0.0_dblprec
         call randomize_spins(Natom,Mensemble,emom,emomM,mmom)

         ! Calculate the starting energy
         call effective_field(Natom,Mensemble,1,Natom,do_jtensor,do_anisotropy,     &
            exc_inter,do_dm,do_pd,do_biqdm,do_bq,do_chir,do_dip,emomM,mmom,         &
            external_field,time_external_field,beff,beff1,beff2,OPT_flag,           &
            max_no_constellations,maxNoConstl,unitCellType,constlNCoup,             &
            constellations,constellationsNeighType,mult_axis,wl_totenergy,Num_macro,&
            cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)

         !!! call calc_energy(nHam,mcmstep,do_dm,do_pd,do_bq,Natom,Nchmax,do_dip,do_biqdm,conf_num,&
         !!!    Mensemble,nHam,Num_macro,1,do_jtensor,plotenergy,Temp,1.0_dblprec,do_lsf,exc_inter,mult_axis,&
         !!!    lsf_field,lsf_interpolate,'N',simid,cell_index,macro_nlistsize,mmom,emom,emomM,&
         !!!    emomM_macro,external_field,time_external_field,max_no_constellations,maxNoConstl,&
         !!!    unitCellType,constlNCoup,constellations,OPT_flag,constellationsNeighType,&
         !!!    wl_totenergy)

         !!! wl_totenergy=wl_totenergy*natom

         !
         write(*,'(1x,a)') 'One-sided minimization in progress'
         print *,'Direction: ',wl_direction

         ! Setup order for Metropolis sweeps
         call choose_random_atom_x(Natom,iflip_a)


         ! Calculate demagnetizing field
         if(demag=='Y') then
            call calc_demag(Natom,Mensemble,demag1,demag2,demag3,demagvol,emomM)
         endif

         ! Calculate the static magnetic fields which will be calculated only once as they are not time dependent
         call calc_external_fields(Natom,Mensemble,NA,hfield,anumb,external_field,  &
            do_bpulse,sitefld,sitenatomfld)

         ! Perform MC sweeps
         print *,'Starting mc'
         do ipmcstep=1,ipmcnstep(1)/2

            ! Calculate total and term resolved energies

            ! Metropolis sweeps
            ! print *,'call wl_warmup'
            call wl_warmup(Natom,Nchmax,Mensemble,nHam,mode,       &
               conf_num,lsf_metric,lsf_window,do_lsf,lsf_field,exc_inter,           &
               lsf_interpolate,do_jtensor,do_dm, do_pd, do_biqdm,do_bq,do_chir,     &
               mult_axis,iflip_a,emomM,emom,mmom,ind_mom_flag,hfield,do_dip,        &
               Num_macro,max_num_atom_macro_cell,cell_index,macro_nlistsize,        &
               macro_atom_nlist,emomM_macro,emom_macro,mmom_macro,do_anisotropy)

            !           ! Calculate and print m_avg
            !           if(mcnstep>20) then
            !              if(mod(mcmstep,mcnstep/20)==0) then

            !              end if
            !           else

            !           end if

            if(mod(mcmstep,mcnstep/10)==0) then
               ! Change order for Metropolis sweeps
               call choose_random_atom_x(Natom,iflip_a)
            end if

            !           write(777,*) wl_totenergy

         enddo

         ! Calculate the starting energy
         call effective_field(Natom,Mensemble,1,Natom,do_jtensor,do_anisotropy,     &
            exc_inter,do_dm,do_pd,do_biqdm,do_bq,do_chir,do_dip,emomM,mmom,         &
            external_field,time_external_field,beff,beff1,beff2,OPT_flag,           &
            max_no_constellations,maxNoConstl,unitCellType,constlNCoup,             &
            constellations,constellationsNeighType,mult_axis,wl_totenergy,Num_macro,&
            cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)

         !!! call calc_energy(nHam,mcmstep,do_dm,do_pd,do_bq,Natom,Nchmax,do_dip,do_biqdm,conf_num,&
         !!!    Mensemble,nHam,Num_macro,1,do_jtensor,plotenergy,Temp,1.0_dblprec,do_lsf,exc_inter,mult_axis,&
         !!!    lsf_field,lsf_interpolate,'N',simid,cell_index,macro_nlistsize,mmom,emom,emomM,&
         !!!    emomM_macro,external_field,time_external_field,max_no_constellations,maxNoConstl,&
         !!!    unitCellType,constlNCoup,constellations,OPT_flag,constellationsNeighType,&
         !!!    wl_totenergy)

         !!! wl_totenergy=wl_totenergy*natom

         if(wl_emax.ne.0.0_dblprec.and.wl_emin.ne.0.0_dblprec) then
            print *,"Energy inverval for Wang-Landau given in input file. No initial phase information used."
         else
            if(iloop==1) then
               wl_emin=wl_totenergy
               print *,' Min energy found as:',wl_emin
               !write(776+iloop,'(3f12.6)') emomM(:,:,1)
            else
               wl_emax=wl_totenergy
               print *,' Max energy found as:',wl_emax
               !write(776+iloop,'(3f12.6)') emomM(:,:,1)
            end if
         end if
      end do

      ! Wang-Landau entity extraction
      !wl_emin=-1.05_dblprec*maxval(abs(cwres))*natom/ry_ev
      wl_emin=wl_lcut*wl_emin
      wl_emax=wl_hcut*wl_emax
      !wl_nhist=0
!      if(wl_nhist==0) wl_nhist=int(natom*abs(wl_emax-wl_emin))
      if(wl_nhist==0) wl_nhist=int(omp_get_max_threads()*abs(wl_emax-wl_emin))
      wl_estep=abs(wl_emax-wl_emin)/(wl_nhist)
      print '(1x,a,f10.2,a,f10.2,a)', ' Wang-Landau energy interval:',wl_emin,'<E<',wl_emax,' (mRy)'
      print '(1x,a,i9)', ' Number of bins in histogram:',wl_nhist
      print '(1x,a,f8.4)', ' Width of histogram bins:',wl_estep
      print '(1x,a,f6.4,a,f8.4,a)',' Broadening:',wl_sigma,' (x energy window) =',wl_sigma*wl_nhist*wl_estep,' (mRy)'
      !     print *,'debug:',wl_emin+wl_estep*wl_nhist

      !wl_maxT=3.0_dblprec*(2.0_dblprec*em/3.0_dblprec/k_bolt_ev/1000)
      ! Deallocate work arrays
      call allocate_mcdata(Natom,-1)

      !call timing(0,'Initial       ','OF')

   end subroutine wl_iphase
   !---------------------------------------------------------------------------
   !> @brief
   !> Wang-Landau measurement phase
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine wl_mphase()
      !
      use Energy,          only : calc_energy
      use InputData
      use FieldData,       only : external_field, sitefld, time_external_field,&
         allocation_field_time, thermal_field, beff, beff1, beff2
      use SystemData
      use MonteCarlo
      use MonteCarlo_common, only : randomize_spins
      use DemagField
      use MomentData
      use FieldPulse
      !use prn_averages
      use ChemicalData,    only : achem_ch, asite_ch
      use MicroWaveField
      use CalculateFields
      use HamiltonianData
      use ChemicalData, only : achtype
      use macrocells
      use optimizationRoutines
      use WangLandau
      use HamiltonianActions
      use omp_lib

      !
      implicit none
      !
      integer :: wl_lhist_min, wl_lhist_max
      integer :: cgk_flag, scount_pulse, bcgk_flag, cgk_flag_pc,wl_count, ii
      integer :: mcmstep, mcmstep_loc, mcmstep_glob
      integer :: num_threads, ithread, i_stat, i_info, printloop
      character(len=30) :: filn
      real(dblprec) :: totenergy,flatness,testene,q_prefac,accrate
      real(dblprec) :: accrate_opt,a,b
      real(dblprec), dimension(3) :: m_avg
      real(dblprec), dimension(:), allocatable :: wl_finaldos, wl_finalhist
      real(dblprec), dimension(:,:), allocatable :: wl_localdos, wl_localhist
      real(dblprec), dimension(:,:), allocatable :: wl_finalmhist
      real(dblprec), dimension(:,:,:), allocatable, save :: wl_emom
      !$omp threadprivate(wl_emom)
      real(dblprec), dimension(:,:,:), allocatable, save :: wl_emomM
      !$omp threadprivate(wl_emomM)
      real(dblprec), dimension(:,:,:), allocatable, save :: ran_w
      !$omp threadprivate(ran_w)

      call timing(0,'MonteCarlo    ','ON')

      ! Flag for G(q) sampling and G(r) calculation
      cgk_flag=0 ; scount_pulse=1 ; bcgk_flag=0 ; cgk_flag_pc=0
      accrate_opt=0.5_dblprec ; a=0.82988_dblprec ; b=0.014625_dblprec
      accrate=0.0_dblprec

      ! Get number of available OpenMP threads
      num_threads=omp_get_max_threads()
      q_prefac=1.0_dblprec/num_threads

      ! Allocate thread-local arrays
      allocate(wl_finalhist(wl_nhist),stat=i_stat)
      call memocc(i_stat,product(shape(wl_finalhist))*kind(wl_finalhist),'wl_finalhist','wl_mphase')
      allocate(wl_finaldos(wl_nhist),stat=i_stat)
      call memocc(i_stat,product(shape(wl_finaldos))*kind(wl_finaldos),'wl_finaldos','wl_mphase')
      allocate(wl_finalmhist(3,wl_nhist),stat=i_stat)
      call memocc(i_stat,product(shape(wl_finalmhist))*kind(wl_finalmhist),'wl_finalmhist','wl_mphase')
      !allocate(wl_localhist(wl_nhist,0:num_threads-1),stat=i_stat)
      !call memocc(i_stat,product(shape(wl_localhist))*kind(wl_localhist),'wl_localhist','wl_mphase')
      !wl_localhist=0.0_dblprec
      !allocate(wl_localdos(wl_nhist,0:num_threads-1),stat=i_stat)
      !call memocc(i_stat,product(shape(wl_localdos))*kind(wl_localdos),'wl_localdos','wl_mphase')
      !wl_localdos=0.0_dblprec
      !$omp parallel
      !allocate(ran_w(3,Natom,Mensemble),stat=i_stat)
      !call memocc(i_stat,product(shape(ran_w))*kind(ran_w),'ran_w','wl_mphase')
      allocate(wl_emom(3,Natom,Mensemble),stat=i_stat)
      allocate(wl_emomM(3,Natom,Mensemble),stat=i_stat)
      !$omp end parallel
      !$omp barrier



      ! Allocate work arrays for MC
      call allocate_mcdata(Natom,1)
      call allocate_wldata(Natom,Mensemble,1)

      ! Start with a fresh spin config.
      !!! emom(:,:,:)=0.0_dblprec
      !!! emom(3,:,:)=1.0_dblprec
      !!! do iii=1,Mensemble
      !!! do ii=1,Natom
      !!!    emomM(:,ii,iii)=emom(:,ii,iii)*mmom(ii,iii)
      !!! end do
      !!! end do
      call randomize_spins(Natom,Mensemble,emom,emomM,mmom)
      !!! open(unit=100,file='fort.777')
      !!! read(100,*) emomm(:,:,1)
      !!! do ii=1,natom
      !!!    mmom(ii,1)=sqrt(sum(emomm(:,ii,1)**2))
      !!!    emom(:,ii,1)=emomm(:,ii,1)/mmom(ii,1)
      !!! end do

      !!! totenergy=0.0_dblprec
      ! Calculate the starting energy
      beff=0.0_dblprec
      beff1=0.0_dblprec
      beff2=0.0_dblprec
      external_field=0.0_dblprec
      time_external_field=0.0_dblprec
      !
      write(*,'(1x,a)')'Wang-Landau sampling in progress.'
      do ii=1,wl_nhist
         wl_ene(ii)=wl_emin+ii*wl_estep
      end do

      ! Initialize Wang-Landau variables
      wl_dos=1.0_dblprec
      wl_hist=0.0_dblprec
      wl_fac=1.0_dblprec
      wl_mhist=0.0_dblprec
      m_avg=0.0_dblprec
      !wl_fac=exp(1.0_dblprec)
      wl_count=1

      ! Not currently used
      !     wl_direction=1.0_dblprec
      ! To be read/set to default
      !wl_nloop=30

      ! Setup order for Metropolis sweeps
      call choose_random_atom_x(Natom,iflip_a)


      ! Calculate demagnetizing field
      if(demag=='Y') then
         call calc_demag(Natom,Mensemble,demag1,demag2,demag3,demagvol,emomM)
      endif

      ! Calculate the static magnetic fields which will be calculated only once as they are not time dependent
      call calc_external_fields(Natom,Mensemble,NA,hfield,anumb,external_field,     &
         do_bpulse,sitefld,sitenatomfld)

      ! Calculate average magnetization for M_hist
      do ii=1,Natom
         m_avg(:)=m_avg(:)+emomM(:,ii,1)
      end do
      ! Not per atom yet
      !m_avg=m_avg/Natom

      wl_lhist_min=1
      wl_lhist_max=wl_nhist

      call omp_set_nested(.false.)
      !$omp parallel default(shared) firstprivate(totenergy,iflip_a,m_avg,wl_lhist_min,wl_lhist_max) private(mcmstep_loc,ithread)
      ithread=omp_get_thread_num()
      call omp_set_nested(.false.)

      wl_emom=emom
      wl_emomM=emomM
      !print *,'randomize_spins', shape(wl_emom), ithread
      call randomize_spins(Natom,Mensemble,wl_emom,wl_emomM,mmom)
      !print *,'effective_field'
      call effective_field(Natom,Mensemble,1,Natom,do_jtensor,do_anisotropy,        &
         exc_inter,do_dm,do_pd,do_biqdm,do_bq,do_chir,do_dip,wl_emomM,mmom,         &
         external_field,time_external_field,beff,beff1,beff2,OPT_flag,              &
         max_no_constellations,maxNoConstl,unitCellType,constlNCoup,constellations, &
         constellationsNeighType,mult_axis,totenergy,Num_macro,cell_index,          &
         emomM_macro,macro_nlistsize,NA,N1,N2,N3)

      !$omp barrier

      print *, 'Starting energy for ensemble',ithread,': ', totenergy
      wl_hit=0.0_dblprec
      wl_miss=0.0_dblprec
      !wl_localdos=1.0_dblprec
      !wl_localhist=0.0_dblprec
      wl_fac=1.0_dblprec
      flatness=0.0_dblprec
      wl_stepsize=1.0_dblprec

      ! Perform MC sweeps
      mcmstep=1
      mcmstep_loc=1

      !wl_localhist(:,ithread)=0.0_dblprec
      !wl_localdos(:,ithread)=wl_dos
      !do while(mcmstep<=mcnstep/num_threads.and.wl_count<=wl_nloop)
      do while(mcmstep<=mcnstep.and.wl_count<=wl_nloop)

         call wl_evolve(Natom,Nchmax,Mensemble,nHam,mode,conf_num, &
            lsf_metric,lsf_window,do_lsf,lsf_field,exc_inter,lsf_interpolate,       &
            do_jtensor,do_dm, do_pd, do_biqdm, do_bq,do_chir,mult_axis,iflip_a,     &
            wl_emomM,wl_emom,mmom,ind_mom_flag,hfield,do_dip,Num_macro,             &
            max_num_atom_macro_cell,cell_index,macro_nlistsize,macro_atom_nlist,    &
            emomM_macro,emom_macro,mmom_macro,wl_hist,wl_dos,totenergy,wl_mhist,    &
            m_avg,wl_lhist_min,wl_lhist_max,do_anisotropy,wl_stepsize)
         !wl_localhist(:,ithread),wl_localdos(:,ithread),totenergy,wl_mhist,m_avg,do_anisotropy)


         if(mod(mcmstep_loc,mcnstep/(num_threads*10))==0) then
            ! Change order for Metropolis sweeps
            call choose_random_atom_x(Natom,iflip_a)
         end if

         !$omp master
         if(mod(mcmstep_loc,10000/num_threads)==0) then

            !-----Adapt WL e-independent stepsize----------------------------------
            accrate=wl_hit/(wl_hit+wl_miss+dbl_tolerance)
            wl_stepsize=min(1.0_dblprec,(log(a*accrate_opt+b)/log(a*accrate+b))*wl_stepsize+0.05_dblprec)
            !--------------------------------------------------------------------
            flatness=minval(wl_hist(2:wl_nhist-1))/(sum(wl_hist(2:wl_nhist-1)+1e-12)/(wl_nhist-2+1e-12))*100.0_dblprec
            !            print '(1x,a,i8,a,f6.2,a,4g12.4)','Step:',mcmstep,', flatness:',flatness, '% ',&
            !            maxval(wl_hist(2:wl_nhist-1)),minval(wl_hist(2:wl_nhist-1)),wl_hit/(1.0_dblprec*Natom*Mensemble*mcmstep*num_threads),sum((wl_hist(1:wl_nhist)-minval(wl_hist(1:wl_nhist)))*wl_estep)
            print '(1x,a,i8,a,f6.2,a,4g12.4)','Step:',mcmstep,', flatness:',flatness, '% ',&
               maxval(wl_hist(2:wl_nhist-1)),minval(wl_hist(2:wl_nhist-1)),accrate,wl_stepsize
            !wl_hist=wl_hist*0.25_dblprec
            wl_hit=0.0_dblprec ; wl_miss=0.0_dblprec

            !if(mod(mcmstep,mcnstep/wl_nloop)==0) then
            if(flatness>70.0_dblprec+wl_count/2.and.maxval(wl_hist)>1.0_dblprec)  then
               ! Reset the WL histogram
               write(filn,'(''wlhistogram.'',a8,''.out'')') simid
               open(ofileno,file=filn, position="append")
               do ii=1,wl_nhist
                  write(ofileno,'(1x,i4,i10,5g20.10)') wl_count,ii,                  &
                     wl_emin+ii*wl_estep,wl_dos(ii),wl_hist(ii),                    &
                     sqrt(sum(wl_mhist(:,ii)**2))/Natom/wl_hist(ii), wl_dos(ii)
               end do
               close(ofileno)
      
               call wl_integration(wl_nhist,wl_ene,wl_dos,wl_hist,simid)
               
               wl_finaldos=wl_dos
               wl_finalhist=wl_hist
               wl_finalmhist=wl_mhist
               wl_fac=wl_fac*0.5_dblprec
               !wl_fac=wl_fac*0.9_dblprec
               !           wl_flatness=minval(wl_hist)/(sum(wl_hist)/wl_nhist)
               !           flatness=minval(wl_hist)/(sum(wl_hist)/wl_nhist)*100.0_dblprec
               write(*,'(1x,a,i3,a,f12.10,a,f12.6)') 'Wang-Landau histogram reset. Sweep: ',wl_count,',  Lambda: ',exp(wl_fac), ', Flatness: ',flatness
               !write(*,'(1x,a,i3,a,f12.10,a,f12.6)') 'Wang-Landau histogram reset. Sweep: ',wl_count,',  Lambda: ',exp(wl_fac), ', Flatness: ',minval(wl_hist)/(sum(wl_hist)/wl_nhist)*100.0_dblprec
               !!wl_dos=wl_dos/maxval(wl_dos)
               !do ii=0,num_threads-1
               !   wl_localhist(:,ii)=0.0_dblprec
               !   wl_localdos(:,ii)=wl_dos
               !end do
               !!$omp atomic
               wl_hist=0.0_dblprec
               wl_mhist=0.0_dblprec
               !              wl_hist=1
               wl_count=wl_count+1
               flatness=0.0_dblprec
               wl_hit=0.0_dblprec
               wl_miss=0.0_dblprec
               !
            end if
            !!! !!$omp barrier
            !!! do ii=0,num_threads-1
            !!!    wl_localdos(:,ii)=wl_dos(:)
            !!!    wl_localhist(:,ii)=wl_hist(:)
            !!! end do
         end if
         !$omp end master

         !print *,totenergy,omp_get_thread_num()

         mcmstep_loc=mcmstep_loc+1
         !!$omp atomic
         mcmstep=mcmstep+1
      enddo

      !$omp end parallel
      call omp_set_nested(.true.)

      call timing(0,'MonteCarlo    ','OF')

      if(wl_count>1) then
         wl_dos=wl_finaldos
         wl_hist=wl_finalhist
         wl_mhist=wl_finalmhist
      end if

      wl_dos=wl_gfac*wl_dos
      write(filn,'(''wlfinal.'',a8,''.out'')') simid
      open(ofileno,file=filn)
      do ii=1,wl_nhist
         !write(ofileno,'(4g20.10)') wl_emin+ii*wl_estep,wl_dos(ii),wl_hist(ii),     &
         write(ofileno,'(4g20.10)') wl_ene(ii),wl_dos(ii),wl_hist(ii),     &
            sqrt(sum(wl_mhist(:,ii)**2))/Natom/(wl_hist(ii)+1.0e-12_dblprec)
      end do
      do ii=1,wl_nhist
         !write(500,'(4g20.10)') wl_emin+ii*wl_estep,wl_mhist(:,ii)/wl_hist(ii)
         write(500,'(4g20.10)') wl_ene(ii),wl_mhist(:,ii)/wl_hist(ii)
      end do
      close(ofileno)

      call wl_integration(wl_nhist,wl_ene,wl_dos,wl_hist,simid)

      ! Deallocate work arrays
      call allocate_mcdata(Natom,-1)
      !     call allocate_wldata(Natom,-1)
      close(11)
   end subroutine wl_mphase


end module wl_driver
