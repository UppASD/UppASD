module QHB
   !
   use Parameters
   !
   implicit none


   ! Quantum Heat Bath flag
   character(LEN=1) :: do_qhb                  !< Quantum Heat Bath (N/Q)
   character(LEN=2) :: qhb_mode                !< Temperature control of QHB (TM/TR/MT)
   real(dblprec)    :: tcurie

   ! Mixing scheme flags
   character(LEN=1) :: do_qhb_mix                  !< Do mixing statistics scheme (N/Y)
   character(LEN=2) :: qhb_mix_mode                !< Mixing function (LI/...) only LI/ implemented
   real(dblprec) :: qhb_Tmix = 0                  !< Mix sampling temperature (dE dependent)
   real(dblprec) :: qhb_Tmix_buff = 0
   real(dblprec) :: qhb_Tmix_prn = -1

   public

contains

   !> Rescaling of QHB using quasiharmonic approximation
   subroutine qhb_rescale(temperature,temprescale,temprescalegrad,domode,mode,mag)
      use Constants, only : k_bolt_ev
      use AMS, only : magdos, tcmfa, tcrpa, msat, magdos_freq

      implicit none

      real(dblprec),intent(in) :: temperature
      real(dblprec),intent(out) :: temprescale
      real(dblprec),intent(out) :: temprescalegrad
      character(len=1),intent(in) :: domode
      character(len=2),intent(in) :: mode
      real(dblprec),intent(in) :: mag


      real(dblprec) :: emin,emax,chb,qhb,deltae,beta,fac,cfac,tcrit,ecut,dT,Tg
      real(dblprec),dimension(magdos_freq) :: mdos,edos,qint
      real(dblprec),dimension(-2:2) :: Tresc
      integer :: i

      mdos=0.0_dblprec ; edos=0.0_dblprec ; dT=1.0_dblprec
      beta=0.3650_dblprec         !beta critical exponent (3D)
      emin=minval(magdos(:,1))
      emax=maxval(magdos(:,1))
      deltae=(emax-emin)/(magdos_freq-1)
      fac=1.0_dblprec
      if (mode=='TM') then
         tcrit=tcmfa
      elseif(mode=='TR') then
         tcrit=tcrpa
      elseif(mode=='TC') then
         tcrit=tcurie
      endif

      if (do_qhb=='Q' .or. do_qhb=='R' .or. do_qhb=='P') then
         do i=-2,2
            Tg=temperature+i*dT
            if (Tg <= 0.0_dblprec) Tg=dbl_tolerance
            if (mode=='MT') then
               fac=mag/msat
            else
               if (Tg < tcrit) then
                  fac=(1.0_dblprec-Tg/tcrit)**beta
               else
                  fac=-1.0_dblprec
               endif
            endif
            if (fac>0) then
               chb=k_bolt_ev*Tg*1000
               if(do_qhb=='Q') then
                  mdos=magdos(:,2)/fac
                  edos=magdos(:,1)*fac
                  qint=edos(:)/(exp(edos(:)/chb)-1.0_dblprec)
                  qint(1)=chb
                  qhb=sum(qint(:)*mdos(:))*deltae*fac
                  Tresc(i)=qhb/chb
               elseif(do_qhb=='R') then
                  ! Enforce classical statistics at Tdebye with additional renormalization
                  qint=magdos(:,1)/(exp(magdos(:,1)/emax)-1.0_dblprec)
                  qint(1)=emax
                  cfac=(sum(qint(:)*magdos(:,2))*deltae)/emax
                  mdos=magdos(:,2)/fac/cfac
                  !--------
                  edos=magdos(:,1)*fac*cfac
                  qint=edos(:)/(exp(edos(:)/chb)-1.0_dblprec)
                  qint(1)=chb
                  qhb=sum(qint(:)*mdos(:))*deltae*fac*cfac
                  Tresc(i)=qhb/chb
                  if (Tresc(i) > 1.0_dblprec) Tresc(i)=1.0_dblprec
               else
                  ! Enforce classical statistics at Tdebye with additional renormalization
                  ecut=k_bolt_ev*tcrit*1000.0_dblprec
                  qint=magdos(:,1)/(exp(magdos(:,1)/ecut)-1.0_dblprec)
                  qint(1)=ecut
                  cfac=(sum(qint(:)*magdos(:,2))*deltae)/ecut
!               mdos=magdos(:,2)/fac/cfac
                  mdos=magdos(:,2)/cfac
               !--------
!               edos=magdos(:,1)*fac
                  edos=magdos(:,1)
                  qint=edos(:)/(exp(edos(:)/chb)-1.0_dblprec)
                  qint(1)=chb
                  qhb=sum(qint(:)*mdos(:))*deltae
!               qhb=sum(qint(:)*mdos(:))*deltae*fac
                  Tresc(i)=qhb/chb
                  if (Tresc(i) > 1.0_dblprec) Tresc(i)=1.0_dblprec
               endif
            else
               Tresc(i)=1.0_dblprec
            endif
         enddo
         temprescale=Tresc(0)
         temprescalegrad=max(0.0_dblprec,((1.0_dblprec/12.0_dblprec)*(Tresc(-2)-Tresc(2))+(2.0_dblprec/3.0_dblprec)*(Tresc(1)-Tresc(-1)))/dT)
      else    ! use m-DOS from SQW or BLS, direct integration
         do i=-2,2
            Tg=temperature+i*dT
            if (Tg <= 0.0_dblprec) Tg=dbl_tolerance
            chb=k_bolt_ev*Tg*1000
            qint(:)=magdos(:,1)/(exp(magdos(:,1)/chb)-1.0_dblprec)
            qint(1)=chb
            qhb=sum(qint*magdos(:,2))*deltae
            Tresc(i)=qhb/chb
         enddo
         temprescale=Tresc(0)
         temprescalegrad=max(0.0_dblprec,((1.0_dblprec/12.0_dblprec)*(Tresc(-2)-Tresc(2))+(2.0_dblprec/3.0_dblprec)*(Tresc(1)-Tresc(-1)))/dT)
      endif

      if (mode=='MT') then
      else
         write(*,3001) 'Temperature rescaling from Quantum Heat Bath:',temprescale, 'T_sim: ', temprescale*temperature, 'dx/dt: ',temprescalegrad
      endif
      3001 format (2x,a,f10.4,2x,a,f7.2,2x,a,f8.4)
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

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !> Mixing Quantum-Classic statistics scheme for Metropolis !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   subroutine init_mix()
      !
      implicit none
      do_qhb_mix = 'N'
      qhb_mix_mode= 'LI'
   end subroutine init_mix
   
   subroutine mix_beta(temperature,temprescale,de,beta_new)
    ! Mixing statistics scheme for Metropolis algorithm where the probability
    ! of acceptance of the new state (W) is given by:
    !
    ! W = alpha*exp(de*beta) + (1-alpha)*exp(de*beta/temprescale)
    ! 
    ! So both quantum and classical statistics are considered.
    !
    ! Effectively, a beta_mix is calculated so:
    ! exp(de*beta_mix) = alpha*exp(de*beta)+(1-alpha)*exp(de*beta/temprescale)
  
      use Constants, only : k_bolt ! Ry units
      use AMS, only : tcmfa, tcrpa
  
      real(dblprec),intent(in) :: temperature ! Simulation temperature
      real(dblprec),intent(in) :: temprescale ! Temperature rescaling factor
      real(dblprec),intent(in) :: de ! Energy diference in Ry
      real(dblprec),intent(out) :: beta_new !Beta mix [output]

      real(dblprec) :: tcrit ! Temperature free parameter in mixing
      real(dblprec) :: alpha ! Mixing factor
      real(dblprec) :: beta_qhb, beta_classic
      integer :: i

      beta_qhb=1.0_dblprec/k_bolt/(temprescale*temperature+1.0d-15)
      beta_classic=1.0_dblprec/k_bolt/(temperature+1.0d-15)
      
      ! Mixing functions - control how quantum W and classic W are mixed
      ! 'LI' - Linear mixing
      ! 'SI' - TODO(?) Sigmoide function would imply more free parameters :(
      !
      ! Linear Mixing
      if (qhb_mix_mode=='LI') then
         if (qhb_mode=='TM') then
            tcrit=tcmfa
         elseif(qhb_mode=='TR') then
            tcrit=tcrpa
         elseif(qhb_mode=='TC') then
            tcrit=tcurie
         endif

         if (temperature<=tcrit) then
            alpha=temperature/tcrit
         elseif (temperature>tcrit) then
            alpha=1
         endif

      endif
      ! end of mixing functions block  

      beta_new=log(alpha*(exp(-de*(beta_classic-beta_qhb))-1)+1)
      beta_new=-(beta_new/de)+beta_qhb
      
      ! Calculates and stores Tmix from beta_new
      qhb_Tmix=(1.0_dblprec/k_bolt/beta_new) 

   end subroutine mix_beta

   subroutine qhb_Tmix_cumu(mcmstep,cumu_step)
       integer, intent(in) :: mcmstep, cumu_step

       qhb_Tmix_buff = qhb_Tmix_buff + qhb_Tmix
       
       if (mod(mcmstep, cumu_step) == 0) then
          qhb_Tmix_prn = qhb_Tmix_buff/cumu_step
          qhb_Tmix_buff = 0
       endif

   end subroutine qhb_Tmix_cumu



end module qhb
