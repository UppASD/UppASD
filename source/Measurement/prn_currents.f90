!> This module will calculate the energy and magnetization currents for out of equilibrium situations
!> @author
!> Jonathan Chico
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt

module prn_currents

   use Parameters
   use Profiling

   implicit none

   ! Input parameters
   integer          :: quant_axis            !< Quantization axis for calculating psi
   integer          :: current_step          !< Interval for sampling currents
   integer          :: current_buff          !< Buffer size for currents
   character(len=1) :: do_currents           !< Measure magnon and heat currents

   ! Measurement variables
   integer                                         :: bcount_current   !< Counter of buffer for the currents
   real(dblprec), dimension(:), allocatable        :: indxb_currents   !< Step counter for the currents
   real(dblprec), dimension(:,:,:), allocatable    :: psi_amp_b        !< Buffer for the heat current
   real(dblprec), dimension(:,:,:), allocatable    :: psi_phase_b      !< Buffer for the magnon current
   real(dblprec), dimension(:,:,:), allocatable    :: heat_current_b   !< Buffer for the heat current
   real(dblprec), dimension(:,:,:), allocatable    :: heat_current2_b  !< Buffer for the heat current
   real(dblprec), dimension(:,:,:), allocatable    :: magnon_current_b !< Buffer for the magnon current
   complex(dblprec), dimension(:,:,:), allocatable :: psi_b

   ! Calculation variables
   complex(dblprec), dimension(:,:), allocatable   :: psi_curr           !< Psi variable for the current time step
   complex(dblprec), dimension(:,:), allocatable   :: psi_prev           !< Psi variable for the previous time step
   complex(dblprec), dimension(:,:), allocatable   :: psi_dot            !< Derivative of the psi variable
   complex(dblprec), dimension(:,:), allocatable   :: psi_dot2           !< Derivative of the psi variable
   complex(dblprec), dimension(:,:,:), allocatable :: psi_time
   real(dblprec), dimension(:,:), allocatable      :: heat_current       !< Heat current per site
   real(dblprec), dimension(:,:), allocatable      :: heat_current2      !< Heat current per site
   real(dblprec), dimension(:,:), allocatable      :: magnon_current     !< Magnon current per time
   real(dblprec), dimension(:,:), allocatable      :: ave_heat_current   !< Average heat current per site
   real(dblprec), dimension(:,:), allocatable      :: ave_magnon_current !< Average magnon current per time
   real(dblprec), dimension(:,:), allocatable      :: psi_amp            !< Amplitude of the complex number
   real(dblprec), dimension(:,:), allocatable      :: psi_phase          !< Phase of the complex number

   private

   !private :: bcount_current, indxb_currents, heat_current_b, magnon_current_b,heat_current2_b,psi_amp_b,psi_b,psi_phase_b
   public :: print_currents, flush_currents, allocate_currents, current_init
   public :: do_currents, current_step, current_buff, quant_axis

   contains

   !> Wrapper routine to print all the currents
   subroutine print_currents(Natom,sstep,mstep,Mensemble,emomM,delta_t,real_time_measure,simid,&
              ncoup,nlistsize,nlist,max_no_neigh,coord)

     implicit none

     integer, intent(in) :: Natom         !< Number of atoms in the system
     integer, intent(in) :: sstep         !< Simulation step in logarithmic scale
     integer, intent(in) :: mstep         !< Current simulation step
     integer, intent(in) :: Mensemble     !< Number of ensembles
     integer, intent(in) :: max_no_neigh  !< Maximum number of neighbours for the neighbour lists
     integer, dimension(Natom), intent(in) :: nlistsize          !< Size of neighbour list for Heisenberg exchange couplings
     integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
     real(dblprec), intent(in) :: delta_t !< Time step for real time measurement
     real(dblprec), intent(in), dimension(3,Natom) :: coord !< Coordinates for all the atoms
     real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
     real(dblprec), dimension(max_no_neigh,Natom), intent(in) :: ncoup !< Heisenberg exchange couplings
     character(len=8), intent(in) :: simid             !< Simulation ID
     character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

     if (do_currents=='Y') then

         ! Calculate the psi quantities needed for the currents
         call mag_axes(Natom,sstep,Mensemble,emomM)

         ! Calculate the heat and magnon currents
         call calc_currents(Natom,mstep,Mensemble,max_no_neigh,nlistsize,nlist,ncoup,coord)

         if (mod(sstep-1,current_step)==0) then

            ! Write step to buffer
            call buffer_current(Natom, Mensemble, mstep-1,bcount_current,delta_t,real_time_measure)

            if (bcount_current==current_buff) then

               ! Write buffer to file
               call prn_site_currents(Natom,Mensemble,simid,real_time_measure)
               bcount_current=1
            else
               bcount_current=bcount_current+1
            endif

         endif

     endif

   end subroutine print_currents

   !> This is a subroutine to calculate the heat and magnon currents in the system
   subroutine calc_currents(Natom,Nstep,Mensemble,max_no_neigh,nlistsize,nlist,ncoup,coord)

      implicit none

      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Nstep        !< Current time step
      integer, intent(in) :: Mensemble    !< Number of ensembles
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(Natom), intent(in) :: nlistsize                !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist       !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(3,Natom) :: coord !< Coordinates for all the atoms
      real(dblprec), dimension(max_no_neigh,Natom), intent(in) :: ncoup !< Heisenberg exchange couplings

      !.. Local variables
      integer :: i,j,k

      if (Nstep.eq.1) psi_prev=0.0D0

      ! This is a calculation of the magnon current
      do j=1, Mensemble
         do i=1, Natom
            do k=1, nlistsize(i)
               ! Current flowing in the x-direction sharing the same quantization axis
               if ((coord(1,nlist(k,i)).gt.coord(1,i)).and.abs(coord(2,nlist(k,i))-&
                    coord(2,i))<dbl_tolerance.and.abs(coord(3,nlist(k,i))-coord(3,i))<dbl_tolerance) then
                   magnon_current(1,i)=magnon_current(1,i)+&
                     2*ncoup(k,i)*AIMAG(CONJG(psi_curr(i,j))*(psi_curr(nlist(k,i),j)-psi_curr(i,j)))
               ! Current flowing in the y-direction sharing the same quantization axis
               else if ((coord(2,nlist(k,i)).gt.coord(2,i)).and.abs(coord(1,nlist(k,i))-&
                         coord(1,i))<dbl_tolerance.and.abs(coord(3,nlist(k,i))-coord(3,i))<dbl_tolerance) then
                   magnon_current(2,i)=magnon_current(2,i)+&
                      2*ncoup(k,i)*AIMAG(CONJG(psi_curr(i,j))*(psi_curr(nlist(k,i),j)-psi_curr(i,j)))
               ! Current flowing in the z-direction sharing the same quantization axis
               else if ((coord(3,nlist(k,i)).gt.coord(3,i)).and.abs(coord(1,nlist(k,i))-&
                         coord(1,i))<dbl_tolerance.and.abs(coord(2,nlist(k,i))-coord(2,i))<dbl_tolerance) then
                   magnon_current(3,i)=magnon_current(3,i)+&
                      2*ncoup(k,i)*AIMAG(CONJG(psi_curr(i,j))*(psi_curr(nlist(k,i),j)-psi_curr(i,j)))
               endif
            enddo
         enddo
      enddo

     magnon_current=magnon_current/Mensemble

     ! This is the calculation of the heat currents
     do j=1, Mensemble
        do i=1, Natom
           do k=1, nlistsize(i)
               ! Current flowing in the x-direction sharing the same quantization axis
               if ((coord(1,nlist(k,i)).gt.coord(1,i)).and.abs(coord(2,nlist(k,i))-&
                    coord(2,i))<dbl_tolerance.and.abs(coord(3,nlist(k,i))-coord(3,i))<dbl_tolerance) then
                   heat_current(1,i)=heat_current(1,i)+&
                      2*ncoup(k,i)*REAL(CONJG(psi_dot(i,j))*(psi_curr(nlist(k,i),j)-psi_curr(i,j)))
               ! Current flowing in the y-direction sharing the same quantization axis
               else if ((coord(2,nlist(k,i)).gt.coord(2,i)).and.abs(coord(1,nlist(k,i))-&
                         coord(1,i))<dbl_tolerance.and.abs(coord(3,nlist(k,i))-coord(3,i))<dbl_tolerance) then
                   heat_current(2,i)=heat_current(2,i)+&
                      2*ncoup(k,i)*REAL(CONJG(psi_dot(i,j))*(psi_curr(nlist(k,i),j)-psi_curr(i,j)))
               ! Current flowing in the z-direction sharing the same quantization axis
               else if ((coord(3,nlist(k,i)).gt.coord(3,i)).and.abs(coord(1,nlist(k,i))-&
                         coord(1,i))<dbl_tolerance.and.abs(coord(2,nlist(k,i))-coord(2,i))<dbl_tolerance) then
                   heat_current(3,i)=heat_current(3,i)+&
                      2*ncoup(k,i)*REAL(CONJG(psi_dot(i,j))*(psi_curr(nlist(k,i),j)-psi_curr(i,j)))
               endif

               ! Current flowing in the x-direction sharing the same quantization axis
               if ((coord(1,nlist(k,i)).gt.coord(1,i)).and.abs(coord(2,nlist(k,i))-&
                    coord(2,i))<dbl_tolerance.and.abs(coord(3,nlist(k,i))-coord(3,i))<dbl_tolerance) then
                   heat_current2(1,i)=heat_current2(1,i)+&
                      2*ncoup(k,i)*REAL(CONJG(psi_dot2(i,j))*(psi_prev(nlist(k,i),j)-psi_prev(i,j)))
               ! Current flowing in the y-direction sharing the same quantization axis
               else if ((coord(2,nlist(k,i)).gt.coord(2,i)).and.abs(coord(1,nlist(k,i))-&
                         coord(1,i))<dbl_tolerance.and.abs(coord(3,nlist(k,i))-coord(3,i))<dbl_tolerance) then
                   heat_current2(2,i)=heat_current2(2,i)+&
                      2*ncoup(k,i)*REAL(CONJG(psi_dot2(i,j))*(psi_prev(nlist(k,i),j)-psi_prev(i,j)))
               ! Current flowing in the z-direction sharing the same quantization axis
               else if ((coord(3,nlist(k,i)).gt.coord(3,i)).and.abs(coord(1,nlist(k,i))-&
                         coord(1,i))<dbl_tolerance.and.abs(coord(2,nlist(k,i))-coord(2,i))<dbl_tolerance) then
                   heat_current2(3,i)=heat_current2(3,i)+&
                      2*ncoup(k,i)*REAL(CONJG(psi_dot2(i,j))*(psi_prev(nlist(k,i),j)-psi_prev(i,j)))
               endif

           enddo
        enddo
     enddo

     heat_current=heat_current/Mensemble
     heat_current2=heat_current2/Mensemble

     ! Sum over the currents to be able to do a time average

   end subroutine calc_currents

   !> This routine will sum up the currents over time
   subroutine sum_currents(Natom)

      implicit none

      integer, intent(in) :: Natom  !< Number of atoms in the system

      !.. Local variables
      integer :: i

      do i=1,Natom
         ! Summation over the magnon current for different times
         ave_magnon_current(1,i)=ave_magnon_current(1,i)+magnon_current(1,i)
         ave_magnon_current(2,i)=ave_magnon_current(2,i)+magnon_current(2,i)
         ave_magnon_current(3,i)=ave_magnon_current(3,i)+magnon_current(3,i)
         ! Summation over the heat current for different times
         ave_heat_current(1,i)=ave_heat_current(1,i)+heat_current(1,i)
         ave_heat_current(1,i)=ave_heat_current(1,i)+heat_current(1,i)
         ave_heat_current(1,i)=ave_heat_current(1,i)+heat_current(1,i)

      enddo

   end subroutine sum_currents

   !> Calculation of the psi variables needed for the currents
   subroutine mag_axes(Natom,mstep,Mensemble,emomM)

      implicit none

      integer, intent(in) :: Natom          !< Number of atoms in the system
      integer, intent(in) :: mstep          !< Current simulation step
      integer, intent(in) :: Mensemble      !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM !< Current magnetic moment vector

      ! Local variables
      real(dblprec) :: normM,rtemp,itemp
      complex(dblprec) :: i
      integer :: j, k

      i=(0.0D0,1.0D0)

      if (mstep.eq.1) then
         psi_prev=0.0d0
         psi_time(1,Natom,Mensemble)=0.0d0
      else
         psi_prev=psi_curr
      endif

      ! Calculation of psi with quantization axis along the x-direction
      if (quant_axis==1) then
         do k=1,Natom
            do j=1, Mensemble
               normM=sqrt(emomM(1,k,j)**2+emomM(2,k,j)**2+emomM(3,k,j)**2)
               psi_curr(k,j)=emomM(2,k,j)+i*emomM(3,k,j)
               psi_curr(k,j)=psi_curr(k,j)/sqrt(2*normM*(normM+emomM(1,k,j)))
             !  write(*,*) mstep,k,psi_curr(k,j)
               if ( mstep.eq.2) then
                  psi_time(2,k,j)=psi_curr(k,j)
               else if (mstep.eq.3) then
                  psi_time(3,k,j)=psi_curr(k,j)
               else
                  psi_time(1,k,j)=psi_time(2,k,j)
                  psi_time(2,k,j)=psi_time(3,k,j)
                  psi_time(3,k,j)=psi_curr(k,j)
               endif
            enddo
         enddo
      ! Calculation of psi with quantization axis along the y-direction
      else if (quant_axis==2) then
         do k=1,Natom
            do j=1, Mensemble
               normM=sqrt(emomM(1,k,j)**2+emomM(2,k,j)**2+emomM(3,k,j)**2)
               psi_curr(k,j)=emomM(1,k,j)+i*emomM(3,k,j)
               psi_curr(k,j)=psi_curr(k,j)/sqrt(2*normM*(normM+emomM(2,k,j)))
               if ( mstep.eq.2) then
                  psi_time(2,k,j)=psi_curr(k,j)
               else if (mstep.eq.3) then
                  psi_time(3,k,j)=psi_curr(k,j)
               else
                  psi_time(1,k,j)=psi_time(2,k,j)
                  psi_time(2,k,j)=psi_time(3,k,j)
                  psi_time(3,k,j)=psi_curr(k,j)
               endif
            enddo
         enddo
      ! Calculation of psi with quantization axis along the z-direction
      else if (quant_axis==3) then
          do k=1,Natom
            do j=1, Mensemble
               normM=sqrt(emomM(1,k,j)**2+emomM(2,k,j)**2+emomM(3,k,j)**2)
               psi_curr(k,j)=emomM(1,k,j)+i*emomM(2,k,j)
               psi_curr(k,j)=psi_curr(k,j)/sqrt(2*normM*(normM+emomM(3,k,j)))
            !   write(*,*) mstep,k,psi_curr(k,j)
               if ( mstep.eq.2) then
                  psi_time(2,k,j)=psi_curr(k,j)
               else if (mstep.eq.3) then
                  psi_time(3,k,j)=psi_curr(k,j)
               else
                  psi_time(1,k,j)=psi_time(2,k,j)
                  psi_time(2,k,j)=psi_time(3,k,j)
                  psi_time(3,k,j)=psi_curr(k,j)
               endif
            enddo
         enddo
      endif

      psi_dot=(psi_curr-psi_prev)
      do k=1,Natom
         do j=1,Mensemble
            psi_dot2(k,j)=(psi_time(3,k,j)-psi_time(1,k,j))
            psi_amp(k,j)=abs(psi_curr(k,j))
            rtemp=REAL(psi_curr(k,j))
            itemp=AIMAG(psi_curr(k,j))
            psi_phase(k,j)=atan2(rtemp,itemp)
         enddo
      enddo

   end subroutine mag_axes

   !> Buffer magnon and heat currents
   subroutine buffer_current(Natom, Mensemble, mstep,bcount_current,delta_t,real_time_measure)

     implicit none

     integer, intent(in) :: mstep     !< Current simulation step
     integer, intent(in) :: Natom     !< Number of atoms in system
     integer, intent(in) :: Mensemble !< Number of ensembles
     integer, intent(in) :: bcount_current !< Counter of buffer for currents
     real(dblprec), intent(in) :: delta_t  !< Current time step (used for real time measurements)
     character(len=1), intent(in) :: real_time_measure   !< Real time measurement flag

     !.. Local variables
     integer :: i,k

     do k=1, Mensemble
        do i=1, Natom
           magnon_current_b(1:3,i,bcount_current)=magnon_current(1:3,i)
           if (mstep.eq.1.or.mstep.eq.2) then
              heat_current2_b(1:3,i,bcount_current)=0.0d0
           else
              heat_current2_b(1:3,i,bcount_current)=heat_current2(1:3,i)
           endif
           heat_current_b(1:3,i,bcount_current)=heat_current(1:3,i)
           psi_amp_b(i,k,bcount_current)=psi_amp(i,k)
           psi_phase_b(i,k,bcount_current)=psi_amp(i,k)
           psi_b(i,k,bcount_current)=psi_curr(i,k)
        end do
     end do

     if (real_time_measure=='Y') then
        indxb_currents(bcount_current)=mstep*delta_t
     else
        indxb_currents(bcount_current)=mstep
     endif

   end subroutine buffer_current

  !> Print magnon and heat currents to file
  subroutine prn_site_currents(Natom, Mensemble, simid,real_time_measure)
    !
    !.. Implicit declarations
    implicit none

    integer, intent(in) :: Natom     !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles
    character(len=8), intent(in) :: simid             !< Name of simulation
    character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

    !.. Local variables
    integer :: i, j,k
    character(len=30) :: filn

    !.. Executable statements
    write (filn,'(''magnon_curr.'',a8,''.out'')') simid
    open(ofileno, file=filn, position="append")
    do i=1, bcount_current
       do j=1, Natom
          if (real_time_measure=='Y') then
             write (ofileno,10003) indxb_currents(i), j, magnon_current_b(1,j,i)/(indxb_currents(i)+1), &
                   magnon_current_b(2,j,i)/(indxb_currents(i)+1),&
                   magnon_current_b(3,j,i)/(indxb_currents(i)+1), (magnon_current_b(1,j,i)/(indxb_currents(i)+1))**2+&
                   (magnon_current_b(2,j,i)/(indxb_currents(i)+1))**2+(magnon_current_b(3,j,i)/(indxb_currents(i)+1))**2
          else
             write (ofileno,10002) int(indxb_currents(i)), j, magnon_current_b(1,j,i)/(indxb_currents(i)+1), &
                   magnon_current_b(2,j,i)/(indxb_currents(i)+1),&
                   magnon_current_b(3,j,i)/(indxb_currents(i)+1), (magnon_current_b(1,j,i)/(indxb_currents(i)+1))**2+&
                   (magnon_current_b(2,j,i)/(indxb_currents(i)+1))**2+(magnon_current_b(3,j,i)/(indxb_currents(i)+1))**2
          endif
       end do
    end do
    close(ofileno)

    write (filn,'(''heat_curr.'',a8,''.out'')') simid
    open(ofileno, file=filn, position="append")
    do i=1, bcount_current
       do j=1, Natom
          if (real_time_measure=='Y') then
             write (ofileno,10003) indxb_currents(i), j, heat_current_b(1,j,i)/(indxb_currents(i)+1), &
                   heat_current_b(2,j,i)/(indxb_currents(i)+1),&
                   heat_current_b(3,j,i)/(indxb_currents(i)+1), (heat_current_b(1,j,i)/(indxb_currents(i)+1))**2+&
                  (heat_current_b(2,j,i)/(indxb_currents(i)+1))**2+(heat_current_b(3,j,i)/(indxb_currents(i)+1))**2
          else
             write (ofileno,10002) int(indxb_currents(i)), j, heat_current_b(1,j,i)/(indxb_currents(i)+1), &
                   heat_current_b(2,j,i)/(indxb_currents(i)+1),&
                   heat_current_b(3,j,i)/(indxb_currents(i)+1), (heat_current_b(1,j,i)/(indxb_currents(i)+1))**2+&
                  (heat_current_b(2,j,i)/(indxb_currents(i)+1))**2+(heat_current_b(3,j,i)/(indxb_currents(i)+1))**2
          endif
       end do
    end do
    close(ofileno)

    write (filn,'(''heat_curr2.'',a8,''.out'')') simid
    open(ofileno, file=filn, position="append")
    do i=1, bcount_current
       do j=1, Natom
          if (real_time_measure=='Y') then
             write (ofileno,10003) indxb_currents(i)-2, j, heat_current2_b(1,j,i)/(indxb_currents(i)+1), &
                   heat_current2_b(2,j,i)/(indxb_currents(i)+1),&
                   heat_current2_b(3,j,i)/(indxb_currents(i)+1), (heat_current2_b(1,j,i)/(indxb_currents(i)+1))**2+&
                  (heat_current2_b(2,j,i)/(indxb_currents(i)+1))**2+(heat_current2_b(3,j,i)/(indxb_currents(i)+1))**2
          else
             write (ofileno,10002) int(indxb_currents(i)), j, heat_current2_b(1,j,i)/(indxb_currents(i)+1), &
                   heat_current2_b(2,j,i)/(indxb_currents(i)+1),&
                   heat_current2_b(3,j,i)/(indxb_currents(i)+1), (heat_current2_b(1,j,i)/(indxb_currents(i)+1))**2+&
                  (heat_current2_b(2,j,i)/(indxb_currents(i)+1))**2+(heat_current2_b(3,j,i)/(indxb_currents(i)+1))**2
          endif
       end do
    end do
    close(ofileno)

    write (filn,'(''psi_data.'',a8,''.out'')') simid
    open(ofileno, file=filn, position="append")
    do i=1, bcount_current
       do k=1, Mensemble
          do j=1, Natom
             if (real_time_measure=='Y') then
                write (ofileno,10004) indxb_currents(i), j, k, psi_amp_b(j,k,i),&
                                 psi_phase_b(j,k,i),REAL(psi_b(j,k,i)),AIMAG(psi_b(j,k,i))
             else
                write (ofileno,10001) int(indxb_currents(i)), j, k,psi_amp_b(j,k,i),&
                                 psi_phase_b(j,k,i),REAL(psi_b(j,k,i)),AIMAG(psi_b(j,k,i))
             endif
          end do
       enddo
    end do
    close(ofileno)


    return
    write (*,*) "Error writing the current files"


10002 format (i8,2x,i8,2x,2x, es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3)
10001 format (i8,2x,i8,2x,i8,2x,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3)
10003 format (es16.4,2x,i8,2x,2x,es16.8,es16.8,es16.8,es16.8)
10004 format (es16.4,2x,i8,2x,i8,2x,2x,es16.8,es16.8,es16.8,es16.8)

  end subroutine prn_site_currents

  !> Flush the trajectory measurements, i.e. print to file in the last iteration
  subroutine flush_currents(Natom,Mensemble,real_time_measure,simid)

    implicit none

     integer, intent(in) :: Natom         !< Number of atoms in the system
     integer, intent(in) :: Mensemble     !< Number of ensembles
     character(len=8), intent(in) :: simid             !< Simulation ID
     character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time


    ! All trajectories are printed to file in the last iteration
    if (do_currents=='Y') then
       !Write buffer to file
       bcount_current=bcount_current-1
       call prn_site_currents(Natom, Mensemble, simid,real_time_measure)
    endif

  end subroutine flush_currents

  !> Allocate and initialize the variables needed for printing the currents
  subroutine allocate_currents(Natom,Mensemble,flag)

    implicit none
    !
    integer, intent(in) :: flag      !< Allocate or deallocate (1/-1)
    integer, intent(in) :: Natom     !< Number of atoms in the system
    integer, intent(in) :: Mensemble !< Number of ensembles
    !.. Local variables
    integer :: i_stat, i_all
    !.. Allocate measurement counters
    if(flag>0) then

       if (do_currents=='Y') then

          allocate(heat_current(3,Natom),stat=i_stat)
          call memocc(i_stat,product(shape(heat_current))*kind(heat_current),'heat_current','allocate_currents')
          allocate(heat_current2(3,Natom),stat=i_stat)
          call memocc(i_stat,product(shape(heat_current2))*kind(heat_current2),'heat_current','allocate_currents')

          allocate(magnon_current(3,Natom),stat=i_stat)
          call memocc(i_stat,product(shape(magnon_current))*kind(magnon_current),'magnon_current','allocate_currents')
          allocate(ave_heat_current(3,Natom),stat=i_stat)
          call memocc(i_stat,product(shape(ave_heat_current))*kind(ave_heat_current),'ave_heat_current','allocate_currents')
          allocate(ave_magnon_current(3,Natom),stat=i_stat)
          call memocc(i_stat,product(shape(ave_magnon_current))*kind(ave_magnon_current),'ave_magnon_current','allocate_currents')
          heat_current=0.0D0
          heat_current2=0.0D0
          magnon_current=0.0D0
          ave_heat_current=0.0D0
          ave_magnon_current=0.0D0

          allocate(psi_curr(Natom,Mensemble),stat=i_stat)
          call memocc(i_stat,product(shape(psi_curr))*kind(psi_curr),'psi_curr','allocate_currents')
          allocate(psi_prev(Natom,Mensemble),stat=i_stat)
          call memocc(i_stat,product(shape(psi_prev))*kind(psi_prev),'psi_prev','allocate_currents')
          allocate(psi_dot(Natom,Mensemble),stat=i_stat)
          call memocc(i_stat,product(shape(psi_dot))*kind(psi_dot),'psi_dot','allocate_currents')
          allocate(psi_dot2(Natom,Mensemble),stat=i_stat)
          call memocc(i_stat,product(shape(psi_dot2))*kind(psi_dot2),'psi_dot2','allocate_currents')
          allocate(psi_time(3,Natom,Mensemble),stat=i_stat)
          call memocc(i_stat,product(shape(psi_time))*kind(psi_time),'psi_time','allocate_currents')
          allocate(psi_amp(Natom,Mensemble),stat=i_stat)
          call memocc(i_stat,product(shape(psi_amp))*kind(psi_amp),'psi_amp','allocate_currents')
          allocate(psi_phase(Natom,Mensemble),stat=i_stat)
          call memocc(i_stat,product(shape(psi_phase))*kind(psi_phase),'psi_phase','allocate_currents')
          psi_curr=0.0D0
          psi_prev=0.0D0
          psi_dot=0.0D0
          psi_dot2=0.0D0
          psi_time=0.0d0
          psi_amp=0.0D0
          psi_phase=0.0D0

          allocate(psi_b(Natom,Mensemble,current_buff),stat=i_stat)
          call memocc(i_stat,product(shape(psi_b))*kind(psi_b),'psi_b','allocate_currents')
          allocate(heat_current_b(3,Natom,current_buff),stat=i_stat)
          call memocc(i_stat,product(shape(heat_current_b))*kind(heat_current_b),'heat_current_b','allocate_currents')
          allocate(psi_amp_b(Natom,Mensemble,current_buff),stat=i_stat)
          call memocc(i_stat,product(shape(psi_amp_b))*kind(heat_current_b),'psi_amp_b','allocate_currents')
          allocate(psi_phase_b(Natom,Mensemble,current_buff),stat=i_stat)
          call memocc(i_stat,product(shape(psi_phase_b))*kind(psi_phase_b),'psi_phase_b','allocate_currents')
          allocate(heat_current2_b(3,Natom,current_buff),stat=i_stat)
          call memocc(i_stat,product(shape(heat_current2_b))*kind(heat_current2_b),'heat_current2_b','allocate_currents')
          allocate(magnon_current_b(3,Natom,current_buff),stat=i_stat)
          call memocc(i_stat,product(shape(magnon_current_b))*kind(magnon_current_b),'magnon_current_b','allocate_currents')
          allocate(indxb_currents(current_buff),stat=i_stat)
          call memocc(i_stat,product(shape(indxb_currents))*kind(indxb_currents),'indxb_currents','allocate_currents')
          heat_current_b=0.0D0
          heat_current2_b=0.0D0
          magnon_current_b=0.0D0
          psi_amp_b=0.0d0
          psi_phase_b=0.0d0
          psi_b=0.0D0

       endif

       !.. Initiate current buffer counters
       bcount_current=1
    else

       if (do_currents=='Y') then

          i_all=-product(shape(heat_current))*kind(heat_current)
          deallocate(heat_current,stat=i_stat)
          call memocc(i_stat,i_all,'heat_current','allocate_currents')
          i_all=-product(shape(magnon_current))*kind(magnon_current)
          deallocate(magnon_current,stat=i_stat)
          call memocc(i_stat,i_all,'magnon_current','allocate_currents')
          i_all=-product(shape(ave_heat_current))*kind(ave_heat_current)
          deallocate(ave_heat_current,stat=i_stat)
          call memocc(i_stat,i_all,'ave_heat_current','allocate_currents')
          i_all=-product(shape(ave_magnon_current))*kind(ave_magnon_current)
          deallocate(ave_magnon_current,stat=i_stat)
          call memocc(i_stat,i_all,'ave_magnon_current','allocate_currents')
          i_all=-product(shape(heat_current_b))*kind(heat_current_b)
          deallocate(heat_current_b,stat=i_stat)
          call memocc(i_stat,i_all,'heat_current_b','allocate_currents')
          i_all=-product(shape(magnon_current_b))*kind(magnon_current_b)
          deallocate(magnon_current_b,stat=i_stat)
          call memocc(i_stat,i_all,'magnon_current_b','allocate_currents')
          i_all=-product(shape(psi_phase_b))*kind(psi_phase_b)
          deallocate(psi_phase_b,stat=i_stat)
          call memocc(i_stat,i_all,'psi_phase_b','allocate_currents')
          i_all=-product(shape(psi_amp_b))*kind(psi_amp_b)
          deallocate(psi_amp_b,stat=i_stat)
          call memocc(i_stat,i_all,'psi_amp_b','allocate_currents')
          i_all=-product(shape(indxb_currents))*kind(indxb_currents)
          deallocate(indxb_currents,stat=i_stat)
          call memocc(i_stat,i_all,'indxb_currents','allocate_trajectories')
          i_all=-product(shape(psi_curr))*kind(psi_curr)
          deallocate(psi_curr,stat=i_stat)
          call memocc(i_stat,i_all,'psi_curr','allocate_currents')
          i_all=-product(shape(psi_prev))*kind(psi_prev)
          deallocate(psi_prev,stat=i_stat)
          call memocc(i_stat,i_all,'psi_prev','allocate_currents')
          i_all=-product(shape(psi_dot))*kind(psi_dot)
          deallocate(psi_dot,stat=i_stat)
          call memocc(i_stat,i_all,'psi_dot','allocate_currents')
          i_all=-product(shape(psi_time))*kind(psi_time)
          deallocate(psi_time,stat=i_stat)
          call memocc(i_stat,i_all,'psi_time','allocate_currents')
          i_all=-product(shape(psi_amp))*kind(psi_amp)
          deallocate(psi_amp,stat=i_stat)
          call memocc(i_stat,i_all,'psi_amp','allocate_currents')
          i_all=-product(shape(psi_phase))*kind(psi_phase)
          deallocate(psi_phase,stat=i_stat)
          call memocc(i_stat,i_all,'psi_phase','allocate_currents')
          i_all=-product(shape(psi_b))*kind(psi_b)
          deallocate(psi_b,stat=i_stat)
          call memocc(i_stat,i_all,'psi_b','allocate_currents')

       endif

    end if


  end subroutine allocate_currents

  !> Initialization of variables with default variables for measuring currents
  subroutine current_init()

     implicit none

     quant_axis   = 1
     do_currents   = 'N'
     current_step = 100
     current_buff = 10


  end subroutine current_init


end module prn_currents
