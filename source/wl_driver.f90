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
         allocation_field_time, beff, beff1, beff2
      use SystemData
      use MonteCarlo
      use DemagField
      use MomentData
      use FieldPulse
      use MicroWaveField
      use CalculateFields
      use HamiltonianData
      use macrocells
      use optimizationRoutines
      use HamiltonianActions
      use WangLandau
      use MonteCarlo_common, only : randomize_spins
      use omp_lib

      !
      implicit none
      !
      integer :: ipmcstep, iloop


      ! Allocate work arrays for MC
      call allocate_mcdata(Natom,1)

      do iloop=1,2

         wl_direction=(-1.0_dblprec)**iloop
         wl_totenergy=0.0_dblprec
         call randomize_spins(Natom,Mensemble,emom,emomM,mmom)

         ! Calculate the starting energy
         call effective_field(Natom,Mensemble,1,Natom,emomM,mmom, &
            external_field,time_external_field,beff,beff1,beff2,OPT_flag,           &
            max_no_constellations,maxNoConstl,unitCellType,constlNCoup,             &
            constellations,constellationsNeighType,wl_totenergy,Num_macro,&
            cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)

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
         call calc_external_fields(Natom,Mensemble,hfield,anumb,external_field,  &
            do_bpulse,sitefld,sitenatomfld)

         ! Perform MC sweeps
         print *,'Starting mc'
         do ipmcstep=1,ipmcnstep(1)/2

            ! Calculate total and term resolved energies

            ! Metropolis sweeps
            ! print *,'call wl_warmup'
            call wl_warmup(Natom,Nchmax,Mensemble,nHam,mode,       &
               conf_num,lsf_metric,lsf_window,do_lsf,lsf_field,ham_inp%exc_inter,           &
               lsf_interpolate,ham_inp%do_jtensor,ham_inp%do_dm, ham_inp%do_pd,  &
            ham_inp%do_biqdm,ham_inp%do_bq,ham_inp%do_ring,ham_inp%do_chir, ham_inp%do_sa,&
               ham_inp%mult_axis,iflip_a,emomM,emom,mmom,ind_mom_flag,hfield,ham_inp%do_dip,        &
               Num_macro,max_num_atom_macro_cell,cell_index,macro_nlistsize,        &
               macro_atom_nlist,emomM_macro,emom_macro,mmom_macro,ham_inp%do_anisotropy)

            if(mod(ipmcstep,ipmcnstep(1)/10)==0) then
               ! Change order for Metropolis sweeps
               call choose_random_atom_x(Natom,iflip_a)
            end if

            !           write(777,*) wl_totenergy

         enddo

         ! Calculate the starting energy
         call effective_field(Natom,Mensemble,1,Natom,emomM,mmom, &
            external_field,time_external_field,beff,beff1,beff2,OPT_flag,           &
            max_no_constellations,maxNoConstl,unitCellType,constlNCoup,             &
            constellations,constellationsNeighType,wl_totenergy,Num_macro,&
            cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)

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
      
      if(wl_nhist==0) wl_nhist=int(omp_get_max_threads()*abs(wl_emax-wl_emin))
      wl_estep=abs(wl_emax-wl_emin)/(wl_nhist)
      print '(1x,a,f10.2,a,f10.2,a)', ' Wang-Landau energy interval:',wl_emin,'<E<',wl_emax,' (mRy)'
      print '(1x,a,i9)', ' Number of bins in histogram:',wl_nhist
      print '(1x,a,f8.4)', ' Width of histogram bins:',wl_estep
      print '(1x,a,f6.4,a,f8.4,a)',' Broadening:',wl_sigma,' (x energy window) =',wl_sigma*wl_nhist*wl_estep,' (mRy)'
   
      ! Deallocate work arrays
      call allocate_mcdata(Natom,-1)


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
         allocation_field_time, beff, beff1, beff2
      use SystemData
      use MonteCarlo
      use MonteCarlo_common, only : randomize_spins
      use DemagField
      use MomentData
      use FieldPulse
      use MicroWaveField
      use CalculateFields
      use HamiltonianData
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
      integer :: mcmstep, mcmstep_loc
      integer :: num_threads, ithread, i_stat
      character(len=30) :: filn
      real(dblprec) :: totenergy, flatness, q_prefac, accrate
      real(dblprec) :: accrate_opt,a,b
      real(dblprec), dimension(3) :: m_avg
      real(dblprec), dimension(:), allocatable :: wl_finaldos, wl_finalhist
      real(dblprec), dimension(:,:), allocatable :: wl_finalmhist
      real(dblprec), dimension(:,:,:), allocatable, save :: wl_emom

      !$omp threadprivate(wl_emom)
      real(dblprec), dimension(:,:,:), allocatable, save :: wl_emomM
      !$omp threadprivate(wl_emomM)

      call timing(0,'MonteCarlo    ','ON')

      ! Flag for G(q) sampling and G(r) calculation
      cgk_flag=0 ; scount_pulse=1 ; bcgk_flag=0 ; cgk_flag_pc=0
      accrate_opt=0.5_dblprec ; a=0.82988_dblprec ; b=0.014625_dblprec
      accrate=0.0_dblprec
      totenergy = 0.0_dblprec

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
      !$omp parallel
      allocate(wl_emom(3,Natom,Mensemble),stat=i_stat)
      allocate(wl_emomM(3,Natom,Mensemble),stat=i_stat)
      !$omp end parallel
      !$omp barrier



      ! Allocate work arrays for MC
      call allocate_mcdata(Natom,1)
      call allocate_wldata(Natom,Mensemble,1)

      ! Start with a fresh spin config.
      call randomize_spins(Natom,Mensemble,emom,emomM,mmom)

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
      wl_count=1

      ! Setup order for Metropolis sweeps
      call choose_random_atom_x(Natom,iflip_a)


      ! Calculate demagnetizing field
      if(demag=='Y') then
         call calc_demag(Natom,Mensemble,demag1,demag2,demag3,demagvol,emomM)
      endif

      ! Calculate the static magnetic fields which will be calculated only once as they are not time dependent
      call calc_external_fields(Natom,Mensemble,hfield,anumb,external_field,     &
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
      
      call randomize_spins(Natom,Mensemble,wl_emom,wl_emomM,mmom)
     
      call effective_field(Natom,Mensemble,1,Natom,wl_emomM,mmom,         &
         external_field,time_external_field,beff,beff1,beff2,OPT_flag,              &
         max_no_constellations,maxNoConstl,unitCellType,constlNCoup,constellations, &
         constellationsNeighType,totenergy,Num_macro,cell_index,          &
         emomM_macro,macro_nlistsize,NA,N1,N2,N3)

      !$omp barrier

      print *, 'Starting energy for ensemble',ithread,': ', totenergy
      wl_hit=0.0_dblprec
      wl_miss=0.0_dblprec
      wl_fac=1.0_dblprec
      flatness=0.0_dblprec
      wl_stepsize=1.0_dblprec

      ! Perform MC sweeps
      mcmstep=1
      mcmstep_loc=1

      do while(mcmstep<=mcnstep.and.wl_count<=wl_nloop)

         call wl_evolve(Natom,Nchmax,Mensemble,nHam,mode,conf_num, &
            lsf_metric,lsf_window,do_lsf,lsf_field,ham_inp%exc_inter,lsf_interpolate,       &
            ham_inp%do_jtensor,ham_inp%do_dm, ham_inp%do_pd, &
            ham_inp%do_biqdm,ham_inp%do_bq,ham_inp%do_ring,ham_inp%do_chir,ham_inp%do_sa,&
            ham_inp%mult_axis,iflip_a,wl_emomM,wl_emom,mmom,ind_mom_flag,hfield,ham_inp%do_dip,Num_macro,     &
            max_num_atom_macro_cell,cell_index,macro_nlistsize,macro_atom_nlist,    &
            emomM_macro,emom_macro,mmom_macro,wl_hist,wl_dos,totenergy,wl_mhist,    &
            m_avg,wl_lhist_min,wl_lhist_max,ham_inp%do_anisotropy,wl_stepsize)


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
            print '(1x,a,i8,a,f6.2,a,4g12.4)','Step:',mcmstep,', flatness:',flatness, '% ',&
               maxval(wl_hist(2:wl_nhist-1)),minval(wl_hist(2:wl_nhist-1)),accrate,wl_stepsize
            
            wl_hit=0.0_dblprec ; wl_miss=0.0_dblprec

            if(flatness>70.0_dblprec+wl_count/2.and.maxval(wl_hist)>1.0_dblprec)  then
               ! Reset the WL histogram
               write(filn,'(''wlhistogram.'',a,''.out'')') trim(simid)
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
               write(*,'(1x,a,i3,a,f12.10,a,f12.6)') 'Wang-Landau histogram reset. Sweep: ',wl_count,',  Lambda: ',exp(wl_fac), ', Flatness: ',flatness
               
               !!$omp atomic
               wl_hist=0.0_dblprec
               wl_mhist=0.0_dblprec
               wl_count=wl_count+1
               flatness=0.0_dblprec
               wl_hit=0.0_dblprec
               wl_miss=0.0_dblprec
               !
            end if
         end if
         !$omp end master


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
      write(filn,'(''wlfinal.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn)
      do ii=1,wl_nhist
         write(ofileno,'(4g20.10)') wl_ene(ii),wl_dos(ii),wl_hist(ii),     &
            sqrt(sum(wl_mhist(:,ii)**2))/Natom/(wl_hist(ii)+1.0e-12_dblprec)
      end do
      do ii=1,wl_nhist
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
