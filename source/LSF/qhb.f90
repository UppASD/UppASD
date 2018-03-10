!> @brief
!> Rescales the temperature used in the Langevin heat bath according to quantum statistics
!> @author
!> Lars Bergqvist, Anders Bergman
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module QHB
   !
   use Parameters
   !
   implicit none


   ! Quantum Heat Bath flag
   character(LEN=1) :: do_qhb                  !< Quantum Heat Bath (N/Q)
   character(LEN=2) :: qhb_mode                !< Temperature control of QHB (TM/TR/MT)
   real(dblprec)    :: tcurie

   public

contains

   !> Rescaling of QHB using quasiharmonic approximation
   subroutine qhb_rescale(temperature,temprescale,domode,mode,mag)
      use Constants, only : k_bolt_ev
      use AMS, only : magdos, tcmfa, tcrpa, msat, magdos_freq

      implicit none

      real(dblprec),intent(in) :: temperature
      real(dblprec),intent(out) :: temprescale
      character(len=1),intent(in) :: domode
      character(len=2),intent(in) :: mode
      real(dblprec),intent(in) :: mag


      real(dblprec) :: emin,emax,chb,qhb,deltae,beta,fac,cfac
      real(dblprec),dimension(magdos_freq) :: mdos,edos,qint
      integer :: k

      mdos=0.d0 ; edos=0.d0
      !      beta=1.d0/3.d0                       !beta critical exponent (3D)
      beta=0.365d0                      !beta critical exponent (3D)
      emin=minval(magdos(:,1))
      emax=maxval(magdos(:,1))
      deltae=(emax-emin)/(magdos_freq-1)
      fac=1.0d0

      if (do_qhb=='Q' .or. do_qhb=='R') then
         if (mode=='TM') then
            if (temperature < tcmfa) then
               fac=(1.d0-temperature/tcmfa)**beta
            else
               fac=-1.d0
            endif
         elseif(mode=='TR') then
            if (temperature < tcrpa) then
               fac=(1.d0-temperature/tcrpa)**beta
            else
               fac=-1.d0
            endif
         elseif(mode=='TC') then
            if (temperature < tcurie) then
               fac=(1.d0-temperature/tcurie)**beta
            else
               fac=-1.d0
            endif
         elseif(mode=='MT') then
            fac=mag/msat
         endif
         if (fac>0) then
            chb=k_bolt_ev*temperature*1000d0
            if(do_qhb=='Q') then
               mdos=magdos(:,2)/fac
               edos=magdos(:,1)*fac
               qint=edos(:)/(exp(edos(:)/chb)-1.d0)
               qint(1)=chb
               qhb=sum(qint(:)*mdos(:))*deltae*fac
               temprescale=qhb/chb
            else
               ! Enforce classical statistics at Tdebye with additional renormalization
               qint=magdos(:,1)/(exp(magdos(:,1)/emax)-1.d0)
               qint(1)=emax
               cfac=(sum(qint(:)*magdos(:,2))*deltae)/emax
               mdos=magdos(:,2)/fac/cfac
               !--------
               edos=magdos(:,1)*fac
               qint=edos(:)/(exp(edos(:)/chb)-1.d0)
               qint(1)=chb
               qhb=sum(qint(:)*mdos(:))*deltae*fac
               temprescale=qhb/chb
               if (temprescale > 1.d0) temprescale=1.d0
            endif
         else
            temprescale=1.d0
         endif
      else    ! use m-DOS from SQW or BLS, direct integration
         chb=k_bolt_ev*temperature*1000
         qint(:)=magdos(:,1)/(exp(magdos(:,1)/chb)-1.d0)
         qint(1)=chb
         qhb=sum(qint*magdos(:,2))*deltae
         temprescale=qhb/chb
      endif

      if (mode=='MT') then
      else
         write(*,3001) 'Temperature rescaling from Quantum Heat Bath:',temprescale, 'T_sim: ', temprescale*temperature
      endif

      3001 format (2x,a,f10.4,2x,a,f7.2)
      3002 format (2x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5)
   end subroutine qhb_rescale

   subroutine init_qhb()
      !
      implicit none
      !
      !Quantum Heat Bath
      do_qhb = 'N'
      qhb_mode= 'TC'
      tcurie=0.01
   end subroutine init_qhb


end module qhb
