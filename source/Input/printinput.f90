!> Routine for printing input parameters
!> @copyright
!> GNU Public License.
module PrintInput
   use Parameters
   use Profiling
   implicit none
   public

   real(dblprec), dimension(3) :: fdum
   integer, dimension(3) :: idum
   character, dimension(3) :: cdum

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

      ! Open outputfile
      write (filn,'(''inp.'',a,''.json'')') trim(simid)
      open(file_id, file=filn)
      write (file_id,'(a)') '{'
      call print_header()
      !!! call print_constants()
      call print_input1()
      call print_simulation()
      call print_initmagnetization()
      call json_key(file_id,'comment')
      write(file_id,'(a)') '"Input data from UppASD simulation"'
      write (file_id,'(a)') "  } "
      close(file_id)
      write (filn,'(''inp.'',a,''.out'')') trim(simid)
      call print_rotation()
      call print_demag()
      call print_mwf_info()
      call print_initialphase()
      call print_measurementphase()
      close(file_id)


   CONTAINS


      subroutine print_header()
         call json_key(file_id,'simid')
         call json_string(file_id,simid)
      end subroutine print_header

      !!! subroutine print_constants()
      !!!    use Constants, only : k_bolt, mub
      !!!    write (file_id,'(a)') "# Physical constants"
      !!!    write (file_id,'(a, /20x ,es16.8)') "Boltzman constant:", k_bolt
      !!!    write (file_id,'(a, /20x ,es16.8)') "Bohr magneton:", mub
      !!!    !write (file_id,'(a,es16.8)') "mRy             ", mry
      !!!    write (file_id,'(a)')  " "
      !!! end subroutine print_constants

      subroutine print_input1()
         !
         ! ncell
         call json_key(file_id,'ncell')
         idum(1)=N1;idum(2)=N2;idum(3)=N3
         call json_int(file_id,idum,3)

         ! bc
         call json_key(file_id,'bc')
         cdum(1)=BC1;cdum(2)=BC2;cdum(3)=BC3
         call json_char(file_id,cdum,3)

         ! cell
         call json_key(file_id,'cell')
         write(file_id,'(a)') '['
         !write (file_id,'(a, 20x,3es15.8/20x ,3es15.8,/20x ,3es15.8)') ' "cell" : ',C1, C2, C3
         call json_float(file_id,C1,3)
         call json_float(file_id,C2,3)
         call json_float(file_id,C3,3,.true.)
         write(file_id,'(25x,a)') '],'

         ! positions
         !!!write (file_id,'(a)') '  "positions" : '
         !!!do i=1,NA
         !!!   write (file_id,'(20x,i0,2x,i0,2x,3es15.8)') anumb_inp(i), atype_inp(i), bas(1:3,i)
         !!!end do

         ! na
         call json_key(file_id,'na')
         call json_int(file_id,(/NA/),1)

         ! nt
         call json_key(file_id,'nt')
         call json_int(file_id,(/NT/),1)

         ! sym
         call json_key(file_id,'sym')
         call json_int(file_id,(/Sym/),1)

         ! natom not used in input
         !write (file_id,'(a, 20x ,i0)') '   "natom" :', Natom

      end subroutine print_input1


      subroutine print_simulation()
         !write (file_id,'(a)') "# Simulation parameters"
         ! llg
         call json_key(file_id,'llg')
         call json_int(file_id,(/llg/),1)

         ! sdealgh
         call json_key(file_id,'sdealgh')
         call json_int(file_id,(/sdealgh/),1)

         ! mensemble
         call json_key(file_id,'mensemble')
         call json_int(file_id,(/mensemble/),1)

      end subroutine print_simulation


      subroutine print_initmagnetization()
         !write (file_id,'(a)') "# Initial magnetization"
         ! initmag
         call json_key(file_id,'initmag')
         call json_int(file_id,(/initmag/),1)
         !!! if(initmag==1) then
         !!!    write (file_id,'(a,i16)') "mseed           ", mseed
         !!! endif
         !!! if(initmag==2) then
         !!!    write (file_id,'(a,i16)') "mseed           ", mseed
         !!!    write (file_id,'(a,es16.8)') "Theta0          ", theta0
         !!!    write (file_id,'(a,es16.8)') "Phi0            ", phi0
         !!! endif
         if(initmag==3) then
            call json_key(file_id,'moments')
            write(file_id,'(a)') '['
            do i=1,NA-1
               call json_float(file_id,aemom_inp(:,i,1,1),3)
            enddo
            call json_float(file_id,aemom_inp(:,i,1,1),3,.true.)
            write(file_id,'(25x,a)') '],'
         endif
         if(initmag==4) then
            call json_key(file_id,'restartfile')
            call json_string(file_id,restartfile)
         endif
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
         if(ham_inp%do_ewald.ne.'N')write (file_id,'(a,a)') "Ewald Summation                ", ham_inp%do_ewald
         if(ham_inp%do_ewald.ne.'N') write(file_id,'(a,es16.8)') "Ewald Parameter                ", ham_inp%Ewald_alpha
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

         subroutine json_key(fileno,key)
            implicit none
            integer, intent(in) :: fileno
            character(*), intent(in) :: key
            !
            write(fileno,'(a20,a)',advance='no') '"'//trim(key)//'"','  :  '
         end subroutine json_key

         subroutine json_string(fileno,val)
            implicit none
            integer, intent(in) :: fileno
            character(*), intent(in) :: val
            !
            write(fileno,'(10x,a,a,a)') ' "',val,'",'
         end subroutine json_string

         subroutine json_char(fileno,val,nel,contflag)
            implicit none
            integer, intent(in) :: fileno
            integer, intent(in) :: nel
            character, dimension(nel), intent(in) :: val
            logical, intent(in), optional :: contflag
            !
            integer :: iel
            logical :: cont = .false.

            !
            if (present(contflag)) cont=contflag
            !
            if(nel==1) then
               write(fileno,'(a,a,a)') ' "',val(1),'",'
            else
               write(fileno,'(10x,a)',advance='no') ' [ '
               do iel=1,nel-1
                  write(fileno,'(a,a,a)',advance='no') '"',val(iel),'" , '
               end do
               if (cont) then 
                  write(fileno,'(a,a,a)') '"',val(iel),'" ] '
               else
                  write(fileno,'(a,a,a)') '"',val(iel),'" ], '
               end if
            end if
         end subroutine json_char

         subroutine json_float(fileno,val,nel,contflag)
            implicit none
            integer, intent(in) :: fileno
            integer, intent(in) :: nel
            real(dblprec), dimension(nel), intent(in) :: val
            logical, intent(in), optional :: contflag
            !
            integer :: iel
            logical :: cont
            !
            if(present(contflag)) then
               cont=contflag
            else
               cont=.false.
            endif
            if(nel==1) then
               write(fileno,'(g14.6,a)') val(1),','
            else
               write(fileno,'(10x,a)',advance='no') ' [ '
               do iel=1,nel-1
                  write(fileno,'(g14.6,a)',advance='no') val(iel),' , '
               end do
               if(cont) then
                  write(fileno,'(g14.6,a)') val(iel),' ]  '
               else
                  write(fileno,'(g14.6,a)') val(iel),' ], '
               end if
            end if
         end subroutine json_float

         subroutine json_int(fileno,val,nel,contflag)
            implicit none
            integer, intent(in) :: fileno
            integer, intent(in) :: nel
            integer, dimension(nel), intent(in) :: val
            logical, intent(in), optional :: contflag
            !
            integer :: iel
            logical :: cont
            if(present(contflag)) then
               cont=contflag
            else
               cont=.false.
            end if
            !
            if(nel==1) then
               write(fileno,'(i8,a)') val(1),','
            else
               write(fileno,'(10x,a)',advance='no') ' [ '
               do iel=1,nel-1
                  write(fileno,'(i8,a)',advance='no') val(iel),' , '
               end do
               if(cont) then
                   write(fileno,'(i8,a)') val(iel),' ]  '
               else
                   write(fileno,'(i8,a)') val(iel),' ], '
                end if
            end if
         end subroutine json_int

!!!          subroutine json_vector(fileno,val,nel)
!!!             implicit none
!!!             integer, intent(in) :: fileno
!!!             integer, intent(in) :: nel
!!!             real(dblprec), dimension(nel), intent(in) :: val
!!!             !
!!!             integer :: iel
!!!             !
!!!             write(fileno,'(a)',advance='no') ' [ '
!!!             do iel=1,nel-1
!!!                write(fileno,'(f12.6,a)',advance='no') val(iel),' , '
!!!             end do
!!!             write(fileno,'(f12.6,a)') val(iel),' ] '
!!!          end subroutine json_vector

end module PrintInput
