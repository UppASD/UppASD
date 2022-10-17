!> Routine for printing input parameters
!> @copyright
!> GNU Public License.
module PrintInput
   use Parameters
   use Profiling
   use JSON
   implicit none
   public

   real(dblprec), dimension(3) :: fdum
   integer, dimension(3) :: idum
   character, dimension(3) :: cdum

   interface yes_no_print
      module procedure yes_no_print_str, yes_no_print_int, yes_no_print_str_arr
   end interface yes_no_print

   interface yaml_print
      module procedure yaml_print_str, yaml_print_int, yaml_print_real
   end interface yaml_print

contains


   !> Print input parameters
   !! @todo Restructure whole subroutine
   subroutine prninp()
      !
      use InputData
      use prn_averages
      use prn_microwaves
      use prn_trajectories

      !.. Implicit declarations
      implicit none

      !.. Local scalars
      integer :: i, file_id
      character(len=20) :: filn

      !.. Executable statements
      file_id=ofileno

      ! Open outputfile and print starting brackets
      write (filn,'(''inp.'',a,''.json'')') trim(simid)
      open(file_id, file=filn)
      write (file_id,'(a)') '{'

      ! Print header
      call print_header()
      ! Print physical constants 
      ! call print_constants()

      ! Print general system info
      call print_input1()

      ! Print simulation parameters
      call print_simulation()

      ! Print initial magnetization info
      call print_initmagnetization()

      ! Conditional printing of unusual features
      if(roteul > 0 )    call print_rotation()
      if(demag .ne. 'N') call print_demag()

      ! Print initial phase info
      call print_initialphase()

      ! Print measurement phase info
      call print_measurementphase()

      ! Print correlation info
      call print_correlation_info()

      ! print comment and close file
      call json_write(file_id,'Input data from UppASD simulation',key='comment',contflag=.true.)
      write (file_id,'(a)') "  } "
      close(file_id)

      ! Print unported outputs to old text format
      !write (filn,'(''inp.'',a,''.out'')') trim(simid)
      !open(file_id, file=filn)
      !call print_mwf_info()
      !close(file_id)


   CONTAINS


      subroutine print_header()
         call json_write(file_id,simid,key='simid')
      end subroutine print_header

      subroutine print_constants()
             call json_write(file_id,(/k_bolt/),1,key='k_bolt')
             call json_write(file_id,(/mub/),1,key='muB')
             call json_write(file_id,(/mry/),1,key='mRy')
      !!!    write (file_id,'(a)')  " "
      end subroutine print_constants

      subroutine print_input1()
         !
         ! ncell
         idum(1)=N1;idum(2)=N2;idum(3)=N3
         call json_write(file_id,idum,3,key='ncell')

         ! bc
         cdum(1)=BC1;cdum(2)=BC2;cdum(3)=BC3
         call json_write(file_id,cdum,3,key='bc')

         ! cell
         call json_key(file_id,'cell')
         write(file_id,'(a)') '['
         call json_write(file_id,C1,3,indflag=.true.)
         call json_write(file_id,C2,3,indflag=.true.)
         call json_write(file_id,C3,3,indflag=.true.,contflag=.true.)
         write(file_id,'(25x,a)') '],'

         ! positions
         call json_write(file_id,anumb_inp,NA,key='anumb_inp')
         call json_write(file_id,atype_inp,NA,key='atype_inp')
         call json_key(file_id,'bas')
         write(file_id,'(a)') '['
         do i=1,NA-1
            call json_write(file_id,bas(:,i),3,indflag=.true.)
         end do
         call json_write(file_id,bas(:,NA),3,indflag=.true.,contflag=.true.)
         write(file_id,'(25x,a)') '],'

         ! na
         call json_write(file_id,(/NA/),1,key='na')

         ! nt
         call json_write(file_id,(/NT/),1,key='nt')

         ! sym
         call json_write(file_id,(/Sym/),1,key='sym')

         ! natom not used in input
         !write (file_id,'(a, 20x ,i0)') '   "natom" :', Natom

      end subroutine print_input1


      subroutine print_simulation()
         !write (file_id,'(a)') "# Simulation parameters"
         ! llg
         call json_write(file_id,(/llg/),1,key='llg')

         ! sdealgh
         call json_write(file_id,(/sdealgh/),1,key='sdealgh')

         ! mensemble
         call json_write(file_id,(/mensemble/),1,key='mensemble')

      end subroutine print_simulation


      subroutine print_initmagnetization()
         !write (file_id,'(a)') "# Initial magnetization"
         ! initmag
         call json_write(file_id,(/initmag/),1,key='initmag')
         if(initmag==1) then
             call json_write(file_id,(/mseed/),1,key='mseed')
         endif
         if(initmag==2) then
            call json_write(file_id,(/mseed/),1,key='mseed')
            call json_write(file_id,(/theta0/),1,key='theta0')
            call json_write(file_id,(/phi0/),1,key='phi0')
         endif
         if(initmag==3) then
            call json_key(file_id,'moments')
            write(file_id,'(a)') '['
            do i=1,NA-1
               call json_write(file_id,aemom_inp(:,i,1,1),3,indflag=.true.)
            enddo
            call json_write(file_id,aemom_inp(:,i,1,1),3,contflag=.true.,indflag=.true.)
            write(file_id,'(25x,a)') '],'
         endif
         if(initmag==4) then
            call json_write(file_id,restartfile,key='restartfile')
         endif
      end subroutine print_initmagnetization


      subroutine print_rotation()
         call json_write(file_id,(/roteul/),1,key='roteul')
         call json_write(file_id,rotang,3,key='rotang')
      end subroutine print_rotation

      subroutine print_demag()
         call json_write(file_id,(/demag/),1,key='demag')
         call json_write(file_id,(/demag1,demag2,demag3/),3,key='demagfield')
         call json_write(file_id,(/demagvol/),1,key='demagvol')
      end subroutine print_demag

      subroutine print_mwf_info()
         write (file_id,'(a)') "**************** Micro-wave field               ****************"
         if (mwf=='Y'.or.mwf=='P'.or.mwf=='I') then
            write (file_id,'(a,a)') "MWF                            ", mwf
            write (file_id,'(a,3es16.8)') "Direction       ", mwfdir(1), mwfdir(2), mwfdir(3)
            write (file_id,'(a,es16.8)') "Amplitude       ", mwfampl
            write (file_id,'(a,es16.8)') "Frequency       ", mwffreq
            if (mwf=='P') write (file_id,'(a,i8)') "Pulse steps       ", mwf_pulse_time
         else if (mwf=='S'.or.mwf=='W') then
            write (file_id,'(a,a)') "MWF Site dependent             ", mwf
            write (file_id,'(a,3es16.8)') "Direction       ", mwfdir(1), mwfdir(2), mwfdir(3)
            write (file_id,'(a,es16.8)') "Amplitude       ", mwfampl
            write (file_id,'(a,es16.8)') "Frequency       ", mwffreq
            if (mwf=='W') write (file_id,'(a,i8)') "Pulse steps       ", mwf_pulse_time
         endif
         if (mwf_gauss=='Y'.or.mwf_gauss=='P') then
            write (file_id,'(a,a)') "MWF frequency broadened        ", mwf_gauss
            write (file_id,'(a,3es16.8)') "Direction       ", mwf_gauss_dir(1), mwf_gauss_dir(2), mwf_gauss_dir(3)
            write (file_id,'(a,es16.8)') "Amplitude       ", mwf_gauss_ampl
            write (file_id,'(a,es16.8)') "Frequency       ", mwf_gauss_freq
            if (mwf_gauss=='P') write (file_id,'(a,i8)') "Pulse steps       ", mwf_gauss_pulse_time
         else if (mwf_gauss=='S'.or.mwf_gauss=='W') then
            write (file_id,'(a,a)') "MWF freqeuncy broadened Site dependent ", mwf_gauss
            write (file_id,'(a,3es16.8)') "Direction       ", mwf_gauss_dir(1), mwf_gauss_dir(2), mwf_gauss_dir(3)
            write (file_id,'(a,es16.8)') "Amplitude       ", mwf_gauss_ampl
            write (file_id,'(a,es16.8)') "Frequency       ", mwf_gauss_freq
            if (mwf_gauss=='W') write (file_id,'(a,i8)') "Pulse steps       ", mwf_gauss_pulse_time
         endif
         write (file_id,*) " "
      end subroutine print_mwf_info


      subroutine print_initialphase()

         call json_write(file_id,ipmode,key='ip_mode')

         call json_write(file_id,iphfield,3,key='ip_hfield')

         if(mode .eq. 'S' .or. mode .eq. 'R' .or. mode .eq. 'L') then
            call json_write(file_id,(/ipnphase/),1,key='ip_nphase')

            if (ipnphase>0) then
               call json_write(file_id,ipnstep,ipnphase,key='ip_nstep')

               call json_write(file_id,ipdelta_t,ipnphase,key='ip_timestep')

               call json_write(file_id,iplambda1,ipnphase,key='ip_damping')
            end if

         else if (mode .eq. 'M' .or. mode .eq. 'H' .or. mode .eq. 'I') then
            call json_write(file_id,(/ipmcnphase/),1,key='ip_mcanneal')

            if (ipmcnphase>0) then
               call json_write(file_id,ipmcnstep,ipmcnphase,key='ip_mcnstep')

               call json_write(file_id,ipTemp,ipmcnphase,key='ip_temp')
            end if
         end if

      end subroutine print_initialphase

      subroutine print_measurementphase()

         call json_write(file_id,(/mode/),1,key='mode')

         call json_write(file_id,(/temp/),1,key='temp')

         call json_write(file_id,hfield,3,key='hfield')

         if(mode .eq. 'S' .or. mode .eq. 'R' .or. mode .eq. 'L') then
            call json_write(file_id,(/nstep/),1,key='nstep')

            call json_write(file_id,(/delta_t/),1,key='timestep')

            call json_write(file_id,(/mplambda1/),1,key='damping')
         else if (mode .eq. 'M' .or. mode .eq. 'H' .or. mode .eq. 'I') then
            call json_write(file_id,(/mcnstep/),1,key='mcnstep')
         end if

         call json_write(file_id,(/real_time_measure/),1,key='real_time_measure')

         call json_write(file_id,(/do_avrg/),1,key='do_avrg')
         if(do_avrg .eq. 'Y') then
            call json_write(file_id,(/avrg_step/),1,key='avrg_step')
         end if

         call json_write(file_id,(/do_tottraj/),1,key='do_tottraj')
         if(do_tottraj .eq. 'Y') then
            call json_write(file_id,(/tottraj_step/),1,key='tottraj_step')
         end if

         call json_write(file_id,(/do_cumu/),1,key='do_cumu')
         if(do_cumu .eq. 'Y') then
            call json_write(file_id,(/cumu_step/),1,key='cumu_step')
         end if

         call json_write(file_id,(/plotenergy/),1,key='plotenergy')
         !if(plotenergy>0) then
         !   call json_write(file_id,(/ene_step/),1,key='ene_step')
         !end if

         call json_write(file_id,(/do_spintemp/),1,key='do_spintemp')
         if(do_spintemp .eq. 'Y') then
            call json_write(file_id,(/spintemp_step/),1,key='spintemp_step')
         end if
         
         call json_write(file_id,(/ntraj/),1,key='ntraj')
         if(ntraj>0) then
            call json_write(file_id,traj_atom,ntraj,key='traj_atom')
            call json_write(file_id,traj_step,ntraj,key='traj_step')
         end if
      end subroutine print_measurementphase

      subroutine print_correlation_info()

         use Correlation
         use Correlation_utils
         use Omegas

         call json_write(file_id,do_sc,key='do_sc')

         if (do_sc .ne. 'N') then

            if (do_sc .eq. 'Y' .or. do_sc .eq. 'C') then
               call json_write(file_id,(/sc%sc_sep/),1,key='sc_sep')
            end if

            if (do_sc .eq. 'Y' .or. do_sc .eq. 'Q') then
               call json_write(file_id,(/sc%sc_step/),1,key='sc_step')
               call json_write(file_id,(/sc%sc_nstep/),1,key='sc_nstep')
            end if

            call json_write(file_id,sc%do_sc_local_axis,key='do_sc_local_axis')

            call json_write(file_id,(/sc%sc_local_axis_mix/),1,key='sc_local_axis_mix')

            call json_write(file_id,(/sc_window_fun/),1,key='sc_window_fun')

            call json_write(file_id,sc%sc_average,key='sc_average')

            call json_write(file_id,qpoints,key='qpoints')

            call json_write(file_id,qfile,key='qfile')

            call json_write(file_id,do_conv,key='do_conv')

            if (do_conv .ne. 'N') then
               call json_write(file_id,(/sigma_q/),1,key='sigma_q')

               call json_write(file_id,(/sigma_w/),1,key='sigma_w')

               call json_write(file_id,(/LQfactor/),1,key='lorentz_q')

               call json_write(file_id,(/LWfactor/),1,key='lorentz_w')

            end if


         end if

      end subroutine print_correlation_info


   end subroutine prninp

   !> Print yaml summary file for compatibility with GUI
   subroutine print_yaml()
      !
      use InputData
      use prn_averages
      use prn_trajectories
      use Correlation
      use AMS, only : do_ams
      use diamag, only: do_diamag



      !.. Implicit declarations
      implicit none

      !.. Local scalars
      integer :: file_id
      character(len=20) :: filn
      integer, dimension(8) :: times

      call date_and_time(VALUES=times)

      ! Open outputfile and print starting brackets
      write (filn,'(''uppasd.'',a,''.yaml'')') trim(simid)
      file_id=ofileno
      open(file_id, file=filn)
      write(file_id,'(a,a)') "simid: ", simid

      write(file_id,'(a,i4,a,i0.2,a,i0.2)') "date: ", times(1),"-",times(2),"-",times(3)
#if defined(VERSION)
      write (file_id,'(a,a)')  "git_revision: ", VERSION
#endif

      write(file_id,'(a)') "siminfo:"

      call yaml_print("temperature:",temp,file_id)
      call yaml_print("timestep:",delta_t,file_id)
      call yaml_print("damping:",mplambda1,file_id)
      call yaml_print("nstep:",nstep,file_id)

      if (mode .eq. 'S') then
         call yaml_print("mode:",'LLG',file_id)
      else if (mode .eq. 'R') then
         call yaml_print("mode:",'SLD',file_id)
      else if (mode .eq. 'M') then
         call yaml_print("mode:",'M-MC',file_id)
      else if (mode .eq. 'M') then
         call yaml_print("mode:",'H-MC',file_id)
      end if

      call yaml_print("alat:",alat,file_id)

      if (do_sc .eq. 'Y' .or. do_sc .eq. 'Q') then
         call yaml_print("sc_step:",sc%sc_step,file_id)
         call yaml_print("sc_nstep:",sc%sc_nstep,file_id)
      end if
      if (do_sc .eq. 'Y' .or. do_sc .eq. 'C') then
         call yaml_print("sc_sep:",sc%sc_sep,file_id)
      end if

      write(file_id,'(a)') "measurables:"
      ! Averages
      call yes_no_print(do_avrg,"averages: ",'Y',file_id)
      ! Trajectories
      call yes_no_print(ntraj,"trajectories: ",0,file_id)
      ! Moments
      call yes_no_print(do_tottraj,"moments: ",'Y',file_id)
      ! S(q,w)
      call yes_no_print(do_sc,"sqw: ",(/'Y','Q'/),file_id)
      ! AMS
      call yes_no_print(do_ams,"ams: ",'Y',file_id)
      ! nc-AMS
      call yes_no_print(do_diamag,"nc-ams: ",'Y',file_id)
      ! Energies
      call yes_no_print(plotenergy,"totenergy: ",0,file_id)
      ! Cumulants
      call yes_no_print(do_cumu,"cumulants:",'Y',file_id)
      ! Cumulants
      call yes_no_print(do_prnstruct,"coord:",0,file_id)

      close(file_id)

   end subroutine print_yaml

   subroutine yes_no_print_str(variable,label,val,file_id)
      implicit none

      character*1, intent(in) :: variable
      character(*), intent(in) :: label
      character*1, intent(in) :: val
      integer, intent(in) :: file_id

      character*20 :: tlabel
      tlabel=label

      if (variable .eq. val) then
         write(file_id, 10001) adjustl(tlabel), "Yes"
      else
         write(file_id, 10001) adjustl(tlabel), "No "
      end if

      10001 format(4x,a16,a3)
   end subroutine yes_no_print_str

   subroutine yes_no_print_str_arr(variable,label,val,file_id)
      implicit none

      character*1, intent(in) :: variable
      character(*), intent(in) :: label
      character*1,dimension(:), intent(in) :: val
      integer, intent(in) :: file_id

      character*20 :: tlabel
      tlabel=label

      if (any(variable .eq. val)) then
         write(file_id, 10001) adjustl(tlabel), "Yes"
      else
         write(file_id, 10001) adjustl(tlabel), "No "
      end if

      10001 format(4x,a16,a3)
   end subroutine yes_no_print_str_arr

   subroutine yes_no_print_int(variable,label,val,file_id)
      implicit none

      integer, intent(in) :: variable
      character(*), intent(in) :: label
      integer, intent(in) :: val
      integer, intent(in) :: file_id
      
      character*20 :: tlabel
      tlabel=label

      if (variable > val) then
         write(file_id, 10001) adjustl(tlabel), "Yes"
      else
         write(file_id, 10001) adjustl(tlabel), "No "
      end if

      10001 format(4x,a16,a3)
   end subroutine yes_no_print_int

   subroutine yaml_print_real(key,val,file_id)
      implicit none

      character(*), intent(in) :: key
      real(dblprec), intent(in) :: val
      integer, intent(in) :: file_id

      character*16 :: tkey
      character*18 :: tval
      tkey=key
      write(tval,'(g14.8)') val

      write(file_id,'(4x,a,a)') adjustl(tkey), adjustl(tval)

   end subroutine yaml_print_real

   subroutine yaml_print_int(key,val,file_id)
      implicit none

      character(*), intent(in) :: key
      integer, intent(in) :: val
      integer, intent(in) :: file_id

      character*16 :: tkey
      character*8 :: tval
      tkey=key
      write(tval,'(i8)') val

      write(file_id,'(4x,a,a)') adjustl(tkey), adjustl(tval)

   end subroutine yaml_print_int

   subroutine yaml_print_str(key,val,file_id)
      implicit none

      character(*), intent(in) :: key
      character(*), intent(in) :: val
      integer, intent(in) :: file_id

      character*16 :: tkey
      tkey=key

      write(file_id,'(4x,a,a)') adjustl(tkey), val

   end subroutine yaml_print_str
end module PrintInput
