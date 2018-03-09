!> Routine for printing input parameters
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module PrintInput
  use Parameters
  use Profiling
  implicit none
  public


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
    integer :: i, j, k, l, file_id
    character(len=20) :: filn

    !.. Executable statements
    file_id=ofileno

    ! Open outputfile
    write (filn,'(''inp.'',a8,''.out'')') simid
    open(file_id, file=filn)
    call print_header()
    call print_constants()
    call print_input1()
    call print_simulation()
    call print_initmagnetization()
    call print_rotation()
    call print_demag()
    call print_mwf_info()
    call print_initialphase()
    call print_measurementphase()
    close(file_id)


  CONTAINS


    subroutine print_header()
      write (file_id,'(a)') "**************** INPUT DATA                     ****************"
      write (file_id,'(a,a)') "Simulation ID:          ", simid
      write (file_id,*)  " "
    end subroutine print_header

    subroutine print_constants()
      use Constants, only : k_bolt, mub, mry
      write (file_id,'(a)') "**************** Physical constants             ****************"
      write (file_id,'(a,es16.8)') "k_Bolt          ", k_bolt
      write (file_id,'(a,es16.8)') "Bohrmag         ", mub
      write (file_id,'(a,es16.8)') "mRy             ", mry
      write (file_id,'(a)')  " "
    end subroutine print_constants

    subroutine print_input1()
      !
      write (file_id,'(a)') "**************** Input 1                        ****************"
      write (file_id,'(a,3i16)') "N1, N2, N3      ", N1, N2, N3
      write (file_id,'(a,a,a,a,a,a)') "BC                             ", BC1, '               ', &
           BC2, '               ', BC3
      write (file_id,'(a, 3es16.8)') "C1              ",C1(1), C1(2), C1(3)
      write (file_id,'(a, 3es16.8)') "C2              ",C2(1), C2(2), C2(3)
      write (file_id,'(a, 3es16.8)') "C3              ",C3(1), C3(2), C3(3)
      write (file_id,'(a,i16)') "NA              ", NA
      write (file_id,'(a,i16)') "NT              ", NT
      write (file_id,'(a,i16)') "Sym             ", Sym
      write (file_id,'(a,i16)') "Natom           ", Natom
      do i=1,NA
         write (file_id,'(i16,3es16.8)') anumb_inp(i), Bas(1,i),Bas(2,i),Bas(3,i)
      end do
      write (file_id,*) " "
      write (file_id,'(a)') "Exchange"

      ! Sort tensorial exchange (SKKR) input later
      if(do_jtensor/=1) then

         do i=1,NT
            if (do_ralloy==0) then
               write (file_id,'(3i16,2es16.8)') i, atype_inp(i), nn(atype_inp(i)), ammom_inp(i,1,1), Landeg_ch(i,1,1)
               do l=1,NN(atype_inp(i))
                  write (file_id,'(4es16.8)') redcoord(atype_inp(i),l,1), &
                       redcoord(atype_inp(i),l,2), redcoord(atype_inp(i),l,3), jc(atype_inp(i),l,1,1,1)
               end do
            else
               write (file_id,'(4i16)') i, atype_inp(i), nn(i), Nch(i)
               do j=1,Nch(i)
               end do
               do j=1,Nch(i)
                  do k=1,Nch(i)
                     write (file_id,'(4i16,es16.8)') j, k
                     do l=1,NN(i)
                        write (file_id,'(4es16.8)') redcoord(atype_inp(i),l,1), &
                             redcoord(atype_inp(i),l,2), redcoord(atype_inp(i),l,3), jc(i,l,j,k,1)
                     end do
                  end do
               end do
            end if
         end do

      end if

      write (file_id,'(a)')  " "
      write (file_id,'(a)') "Anisotropy"
      do i=1,NT
         write (file_id,'(6es16.8)') anisotropy(i,1,1), anisotropy(i,2,1), anisotropy(i,3,1), anisotropy(i,4,1), &
              anisotropy(i,5,1),anisotropy(i,6,1)
      enddo
      write (file_id,'(a)')  " "
    end subroutine print_input1


    subroutine print_simulation()
      write (file_id,'(a)') "**************** Simulation parameters          ****************"
      write (file_id,'(a,i16)') "LLG             ", llg
      write (file_id,'(a,i16)') "SDE             ", SDEalgh
      write (file_id,'(a,i16)') "Mensemble       ", Mensemble
      write (file_id,'(a)')  " "
    end subroutine print_simulation


    subroutine print_initmagnetization()
      write (file_id,'(a)') "**************** Initial magnetization          ****************"
      write (file_id,'(a, i16)') "initmag        ", initmag
      if(initmag==1) then
         write (file_id,'(a,i16)') "mseed           ", mseed
      endif
      if(initmag==2) then
         write (file_id,'(a,i16)') "mseed           ", mseed
         write (file_id,'(a,es16.8)') "Theta0          ", theta0
         write (file_id,'(a,es16.8)') "Phi0            ", phi0
      endif
      if(initmag==3) then
         do i=1,NA
            write (file_id,'(a,3es16.8)') "Mom             ", aemom_inp(1,i,1,1), aemom_inp(2,i,1,1), aemom_inp(3,i,1,1)
         enddo
      endif
      if(initmag==4) then
         write (file_id,'(a,a)') "File            ", restartfile
      endif
      write (file_id,*) "   "
    end subroutine print_initmagnetization


    subroutine print_rotation()
      write (file_id,'(a)') "**************** Rotation parameters            ****************"
      write (file_id,'(a,i16)') "Rotation        ", roteul
      write (file_id,'(a,3es16.8)') "Rotation Angle  ", rotang(1), rotang(2), rotang(3)
      write (file_id,*) " "
    end subroutine print_rotation


    subroutine print_demag()
      write (file_id,'(a)') "**************** Demagnetization field          ****************"
      write (file_id,'(a,a,a,a,a,a)') "Demag                          ", demag1, '               ', &
           demag2, '               ', demag3
      write (file_id,'(a,a)') "Demag                          ", demag
      write (file_id,'(a,es16.8)') "Volume          ", demagvol
      write (file_id,*) " "
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
      write (file_id,'(a)') "**************** Initial phase                  ****************"
      write (file_id,'(a,a)') "MC, SD, EM                       ", ipmode
      write (file_id,'(a,3es16.8)') "Field           ", iphfield(1), iphfield(2), iphfield(3)
      if(ipmode.eq.'M' .or. ipmode.eq.'H') then
         write (file_id,'(a,i16)') "nphase          ", ipmcnphase
         write (file_id,'(a)') "MC steps, T :"
         do i=1,ipmcnphase
            write (file_id,'(i16,es16.8)') ipmcnstep(i), ipTemp(i)
         enddo
      elseif(ipmode.eq.'S') then
         write (file_id,'(a,i16)') "nphase          ", ipnphase
         write(file_id,'(a)') "Nstep         dt            T        lambda1      lambda2 "
         if (ipnphase>0) then
            do i=1,ipnphase
               write(file_id,'(i10,5es12.5)') ipnstep(i),ipdelta_t(i), iptemp(i), iplambda1(i), iplambda2(i)
            enddo
         endif
      endif
      write (file_id,'(a)') " "
    end subroutine print_initialphase


    subroutine print_measurementphase()
      use Temperature, only : grad

      write (file_id,'(a)') "**************** Measurement phase              ****************"
      write (file_id,'(a,a)') "MC, SD, MEP                         ", mode
      write (file_id,'(a,a)') "Gradient				 ", grad
      if(do_ewald.ne.'N')write (file_id,'(a,a)') "Ewald Summation                ", do_ewald
      if(do_ewald.ne.'N') write(file_id,'(a,es16.8)') "Ewald Parameter                ", Ewald_alpha
      if(grad.eq.'N') then
         write (file_id,'(a,es16.8)') "T               ", Temp
      end if
      write (file_id,'(a,3es16.8)') "Field           ", hfield(1), hfield(2), hfield(3)
      if(mode.eq.'M') then
         write (file_id,'(a,i16)') "Nstep           ", mcnstep
         write (file_id,'(a,2i16)') "Step, Buffer    ", mcavrg_step, mcavrg_buff
      elseif(mode.eq.'S') then
         write (file_id,'(a,2es16.8)') "Lambda          ", mplambda1, mplambda2
         write (file_id,'(a,i16)') "Nstep           ", Nstep
         write (file_id,'(a,es16.8)') "delta_t         ", delta_t
      endif
      write (file_id,'(a)') "   "
      write (file_id,'(a,a,2i16)') "Tottraj                        ", do_tottraj, tottraj_step, tottraj_buff
      write (file_id,'(a,i16)')    "No. trajectories", ntraj
      do i=1,ntraj
         write (file_id,'(a,3i16)')   "Atom Trajectory ", traj_atom(i), traj_step(i), traj_buff(i)
      enddo
   end subroutine print_measurementphase

  end subroutine prninp


end module PrintInput
