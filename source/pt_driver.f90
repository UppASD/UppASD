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
module pt_driver
   use Constants
   use Parameters

   implicit none

   integer :: pt_step = 100
   integer :: pt_numrep

   real(dblprec),dimension(:,:,:,:), allocatable :: pt_emom   !< Unit moments for PT replicas
   real(dblprec),dimension(:,:,:,:), allocatable :: pt_emomM  !< Moments for PT replicas
   real(dblprec),dimension(:,:,:), allocatable :: pt_mmom   !< Moment magnitudes for PT replicas

   real(dblprec),dimension(:,:,:,:), allocatable :: pt_emom_macro   !< Unit macrospin moments for PT replicas
   real(dblprec),dimension(:,:,:,:), allocatable :: pt_emomM_macro  !< Macrospin moments for PT replicas
   real(dblprec),dimension(:,:,:), allocatable :: pt_mmom_macro   !< Macrospin moment magnitudes for PT replicas

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
   subroutine pt_iphase()

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

      integer :: irep, ipmcstep
      integer :: i_all, i_stat
      real(dblprec) :: temprescale
      character(len=30) :: filn
      integer :: mc_nmeas,mc_nsamp, ishift
      real(dblprec) :: mc_mavg, mc_mavg2, mc_mavg4
      real(dblprec) :: mc_avcum,mc_avsus
      real(dblprec) :: mavg
      real(dblprec) :: temp_tmp
      real(dblprec) :: try_pow
      real(dblprec) :: n_hit, n_miss, n_trial
      real(dblprec), dimension(:), allocatable :: repene  !< Array of replica energies
      integer, dimension(:), allocatable :: replist   !< List to denominate which temperature each replica has
      integer, dimension(:), allocatable :: templist  !< List to denominate which replica has which temperature
      ! I.e. if replica 3 has the highest temperature, then replist(3)=1, and templist(1)=3
      ! At start, replist(i)=templist(i)=i.
      real(dblprec), dimension(:), allocatable :: pt_flipprob !< Array for RNG flipping probabilities.

      real(dblprec) :: temp_1,temp_2,beta_1,beta_2, ene_1, ene_2, trial_prob

      pt_numrep=ipmcnphase
      ! Write header for moment file
      write (filn,'(''ptinitial.'',a,''.out'')') trim(simid)
      open(11, file=filn, position="append")
      write(11,'(a)') "#  Iter.   M_avg.    U_Bind.    Susc."

      ! Allocate work arrays for Metropolis algorithm
      call allocate_mcdata(Natom,1)

      n_hit=0.0_dblprec;n_miss=0.0_dblprec;n_trial=0.0_dblprec

      ! Allocate local moment arrays
      allocate(pt_flipprob(pt_numrep),stat=i_stat)
      call memocc(i_stat,product(shape(pt_flipprob))*kind(pt_flipprob),'pt_flipprob','pt_iphase')
      allocate(templist(pt_numrep),stat=i_stat)
      call memocc(i_stat,product(shape(templist))*kind(templist),'templist','pt_iphase')
      allocate(replist(pt_numrep),stat=i_stat)
      call memocc(i_stat,product(shape(replist))*kind(replist),'replist','pt_iphase')
      allocate(repene(pt_numrep),stat=i_stat)
      call memocc(i_stat,product(shape(repene))*kind(repene),'repene','pt_iphase')
      allocate(pt_emom(3,Natom,Mensemble,pt_numrep),stat=i_stat)
      call memocc(i_stat,product(shape(pt_emom))*kind(pt_emom),'pt_emom','pt_iphase')
      allocate(pt_emomM(3,Natom,Mensemble,pt_numrep),stat=i_stat)
      call memocc(i_stat,product(shape(pt_emomM))*kind(pt_emomM),'pt_emomM','pt_iphase')
      allocate(pt_mmom(Natom,Mensemble,pt_numrep),stat=i_stat)
      call memocc(i_stat,product(shape(pt_mmom))*kind(pt_mmom),'pt_mmom','pt_iphase')

      if(Num_macro>0) then
         allocate(pt_emom_macro(3,Num_macro,Mensemble,pt_numrep),stat=i_stat)
         call memocc(i_stat,product(shape(pt_emom_macro))*kind(pt_emom_macro),'pt_emom','pt_iphase')
         allocate(pt_emomM_macro(3,Num_macro,Mensemble,pt_numrep),stat=i_stat)
         call memocc(i_stat,product(shape(pt_emomM_macro))*kind(pt_emomM_macro),'pt_emomM_macro','pt_iphase')
         allocate(pt_mmom_macro(Num_macro,Mensemble,pt_numrep),stat=i_stat)
         call memocc(i_stat,product(shape(pt_mmom_macro))*kind(pt_mmom_macro),'pt_mmom_macro','pt_iphase')
      end if

      do irep=1,pt_numrep

         ! Write output to stdout
         write (*,'(a28,i3,a10,G11.4,a10,i10,a10,a10)') &
            "Performing PT,  Replica:", irep ," Temp: ", ipTemp(irep)

         pt_mmom(:,:,irep)=mmom
         pt_emom(:,:,:,irep)=emom
         pt_emomM(:,:,:,irep)=emomM

         if(Num_macro>0) then
            pt_mmom_macro(:,:,irep)=mmom_macro
            pt_emom_macro(:,:,:,irep)=emom_macro
            pt_emomM_macro(:,:,:,irep)=emomM_macro
         end if

         replist(irep)=irep
         templist(irep)=irep

      end do


      ! Calculate demagnetization field
      if(demag=='Y') then
         call calc_demag(Natom, Mensemble, demag1, demag2, demag3, demagvol, emomM)
      endif

      ! Calculate the static magnetic fields which will be calculated only once as they are not time dependent
      call calc_external_fields(Natom, Mensemble, hfield, anumb, external_field,&
         do_bpulse,sitefld,sitenatomfld)

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
      ! Loop over steps of sweeps
      do ipmcstep=1,ipmcnstep(1)
         do irep=1,pt_numrep
            
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

            call mc_evolve(Natom,Nchmax,Mensemble,nHam,ipTemp(irep),temprescale,'H',&
               conf_num,lsf_metric,lsf_window,do_lsf,lsf_field,ham_inp%exc_inter,           &
               lsf_interpolate,ham_inp%do_jtensor,ham_inp%do_dm, ham_inp%do_pd, ham_inp%do_biqdm,ham_inp%do_bq, ham_inp%do_ring,    &
               ham_inp%do_chir, ham_inp%do_sa, ham_inp%mult_axis,iflip_a,pt_emomM(:,:,:,irep),pt_emom(:,:,:,irep), &
               pt_mmom(:,:,irep),ind_mom_flag,iphfield(1:3),ham_inp%do_dip,Num_macro,       &
               max_num_atom_macro_cell,cell_index,macro_nlistsize,macro_atom_nlist, &
               pt_emomM_macro(:,:,:,irep),pt_emom_macro(:,:,:,irep),                &
               pt_mmom_macro(:,:,irep),ham_inp%do_anisotropy)


            !Calculate energy for PT exchange
            if( mod(ipmcstep, pt_step)==0 ) then
               repene(irep)=0.0_dblprec

               call effective_field(Natom,Mensemble,1,Natom,pt_emomM(:,:,:,irep),   &
                  pt_mmom(:,:,irep),            &
                  external_field,time_external_field,beff,beff1,beff2,OPT_flag,     &
                  max_no_constellations,maxNoConstl,unitCellType,constlNCoup,       &
                  constellations,constellationsNeighType,repene(irep),    &
                  Num_macro,cell_index,pt_emomM_macro(:,:,:,irep),macro_nlistsize,  &
                  NA,N1,N2,N3)


            end if

         end do

            !PT exchange here (Not for the last 10% of the run (to equilibrate
            !temperature)
            if( mod(ipmcstep, pt_step)==0 .and.  ipmcstep/(1.0d0*ipmcnstep(1))<0.90d0)  then

               call rng_uniform(pt_flipprob(:),pt_numrep)

               do ishift=0,1
                  do irep=1+ishift,pt_numrep-1,2

                     n_trial=n_trial+1.0_dblprec
                     
                     temp_1=ipTemp(irep)
                     temp_2=ipTemp(irep+1)
                     beta_1=1.0_dblprec/k_bolt/(temprescale*temp_1+1.0d-15)
                     beta_2=1.0_dblprec/k_bolt/(temprescale*temp_2+1.0d-15)
                     ene_1=repene(irep)*mry/Natom
                     ene_2=repene(irep+1)*mry/Natom
                     try_pow=(ene_2-ene_1)*(beta_2-beta_1)
                     
                     trial_prob=min(1.0_dblprec,exp(try_pow))
                     
                     if(pt_flipprob(irep)<=trial_prob) then
                        temp_tmp=ipTemp(irep)
                        ipTemp(irep)=ipTemp(irep+1)
                        ipTemp(irep+1)=temp_tmp
                        n_hit=n_hit+1.0_dblprec
                     else
                        n_miss=n_miss+1.0_dblprec
                     !  
                  end if

               end do
            end do

         end if

         ! Sample m, m2, m4 for second half of run (every tenth step)
         if(ipmcstep>ipmcnstep(1)/2.and.mod(ipmcstep-1,10)==0) then


         end if

         ! Print m_avg
         if (mod(ipmcstep,ipmcnstep(1)/10)==0) then

            print '(2x,a,i5,a,f10.0,a,f5.2,a)', "IP PT ", ipmcstep*100/ipmcnstep(1),"% done.  No. trials: ",n_trial, " Hitrate: ",n_hit/n_trial*100, "%"
            call pt_prn_mag_conf_clone(Natom,0,pt_numrep,'X',simid,pt_mmom(:,1,:),pt_emom(:,:,1,:),'',mode,iptemp)

            call choose_random_atom_x(Natom,iflip_a)
         end if

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

      do irep=1,pt_numrep
         !Calculate energy for PT exchange
         repene(irep)=0.0_dblprec

         call effective_field(Natom,Mensemble,1,Natom, &
            pt_emomM(:,:,:,irep),pt_mmom(:,:,irep),external_field,                  &
            time_external_field,beff,beff1,beff2,OPT_flag,max_no_constellations,    &
            maxNoConstl,unitCellType,constlNCoup,constellations,                    &
            constellationsNeighType,repene(irep),Num_macro,cell_index,    &
            pt_emomM_macro(:,:,:,irep),macro_nlistsize,NA,N1,N2,N3)

            call calc_mavrg(Natom,Mensemble,pt_emomM(:,:,:,irep),mavg)

            print '(a,i4,a,i4,a,f12.6,a,f12.6,a,f12.6)', 'Replica #',irep, " List position: ", templist(irep), &
               " Average magnetization:", mavg, " Temperature: ", iptemp(irep),"  Energy: " , repene(irep)/Natom

      end do

      ! Copy equillibrated moments for use in measurement phase
      irep=1
      
      do while (repene(irep)>minval(repene))
         irep=irep+1
      end do
      print *,' Selected replica:',irep

      emomM=pt_emomM(:,:,:,irep)
      emom=pt_emom(:,:,:,irep)
      mmom=pt_mmom(:,:,irep)
      emom2=pt_emom(:,:,:,irep)
      mmom2=pt_mmom(:,:,irep)
      print '(1x,a,f10.0,a,f5.2)',"No. PT trials: ",n_trial, " Hitrate: ",n_hit/n_trial*100, "%"

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


      i_all=-product(shape(pt_emomM))*kind(pt_emomM)
      deallocate(pt_emomM,stat=i_stat)
      call memocc(i_stat,i_all,'pt_emomM','pt_iphase')
      i_all=-product(shape(pt_emomM))*kind(pt_emomM)
      deallocate(pt_emom,stat=i_stat)
      call memocc(i_stat,i_all,'pt_emom','pt_iphase')
      i_all=-product(shape(pt_mmom))*kind(pt_mmom)
      deallocate(pt_mmom,stat=i_stat)
      call memocc(i_stat,i_all,'pt_mmom','pt_iphase')

      if(Num_macro>0) then
         i_all=-product(shape(pt_emomM_macro))*kind(pt_emomM_macro)
         deallocate(pt_emomM_macro,stat=i_stat)
         call memocc(i_stat,i_all,'pt_emomM_macro','pt_iphase')
         i_all=-product(shape(pt_emomM_macro))*kind(pt_emomM_macro)
         deallocate(pt_emom_macro,stat=i_stat)
         call memocc(i_stat,i_all,'pt_emom_macro','pt_iphase')
         i_all=-product(shape(pt_mmom_macro))*kind(pt_mmom_macro)
         deallocate(pt_mmom_macro,stat=i_stat)
         call memocc(i_stat,i_all,'pt_mmom_macro','pt_iphase')
      end if

      i_all=-product(shape(pt_flipprob))*kind(pt_flipprob)
      deallocate(pt_flipprob,stat=i_stat)
      call memocc(i_stat,i_all,'pt_flipprob','pt_iphase')
      i_all=-product(shape(templist))*kind(templist)
      deallocate(templist,stat=i_stat)
      call memocc(i_stat,i_all,'templist','pt_iphase')
      i_all=-product(shape(replist))*kind(replist)
      deallocate(replist,stat=i_stat)
      call memocc(i_stat,i_all,'replist','pt_iphase')
      i_all=-product(shape(repene))*kind(repene)
      deallocate(repene,stat=i_stat)
      call memocc(i_stat,i_all,'repene','pt_iphase')
      !
      return
      !
   end subroutine pt_iphase

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: prn_mag_conf_iter
   !> @brief Prints a given magnetic configuration for either a restartfile or a momentfile
   !> @details Prints a magnetic configuration, the objective is to make all the types
   !> of printing honogeneous, of that way restartfiles, momentfiles and GNEB files
   !> would have all the same structure.
   !> Cloned from prn_mag_conf_iter() 
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine pt_prn_mag_conf_clone(Natom,mstep,Mensemble,type,simid,mmom,emom,suffix, mode,Mtemp)

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
         write (filn,'(''pt_restart'',a,''.'',a,''.out'')') suffix,trim(simid)
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

   end subroutine pt_prn_mag_conf_clone

end module pt_driver
