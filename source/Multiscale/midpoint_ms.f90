!-------------------------------------------------------------------------------
!  MODULE: Midpoint_ms
!> @brief
!> The semi-implicit midpoint solver for the LLG-equations for multiscale
!> @details Ref: J.H. Mentink et al, J. Phys.: Condens. Matter, 22, 176001 (2010)
!> @authors
!> Johan Mentink
!> Edgar Mendez
!> Nikos Ntallis
!> Manuel Pereiro
!> Jonathan Chico
!> Anders Bergman
!> Nastaran Salehi
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module Midpoint_ms

   use Profiling
   use Parameters
   use MultiscaleDampingBand, only: DampingBandData
   use Multiscale, only : multiscaleBackbuffer, multiscaleBackbufferHead 

   implicit none

   real(dblprec), dimension(:,:,:), allocatable :: btorque_full !< Resulting effective field

   private :: btorque_full

contains


  subroutine applyMultiscaleDapingband(atom,ensemble, dband, demomdt, &
       m1,m2,m3, a1,a2,a3)
    implicit none
    integer, intent(in) :: atom,ensemble
    type(DampingBandData), intent(in) :: dband
    real(dblprec), dimension(3),intent(in) :: demomdt
    real(dblprec), intent(in) :: m1,m2,m3
    real(dblprec), intent(inout) :: a1,a2,a3

    integer :: index
    real(dblprec) :: gamma,sdnorm, ma1,ma2,ma3    
  
    index = dband%interpolation%indices(atom)
    if(index .ne. 0) then
       sdnorm = sum(demomdt**2) ** 0.25_dblprec
       gamma = dband%coefficients(index) * sdnorm
       ma1 = gamma*dband%preinterpolation(1,index,ensemble)
       ma2 = gamma*dband%preinterpolation(2,index,ensemble)
       ma3 = gamma*dband%preinterpolation(3,index,ensemble)

       a1 = a1 - (m2*ma3 - m3*ma2)
       a2 = a2 - (m3*ma1 - m1*ma3)
       a3 = a3 - (m1*ma2 - m2*ma1)
          
    end if
  end subroutine applyMultiscaleDapingband



  subroutine smodeulermpt_ms(Natom, Mensemble, Landeg,bn, lambda1_array, beff, emom, emom2, emomM, mmom, deltat,thermal_field,dband,&
                            STT,do_she,do_sot,btorque,she_btorque,sot_btorque)
           !  nlist,nlistsize,constellationsUnitVec2,unitCellType,OPT_flag,cos_thr)
    use Constants
    use RandomNumbers, only : ranv
    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles
    real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
    real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0d0)
    real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: beff !< Total effective field from application of Hamiltonian
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
    real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
    real(dblprec), intent(in) :: deltat !< Time step
    type(DampingBandData), intent(in) :: dband !< Damping band info (multiscale)
    character(len=1), intent(in) :: STT    !< Treat spin transfer torque
    character(len=1), intent(in) :: do_she !< Treat the SHE spin transfer torque
    character(len=1), intent(in) :: do_sot !< Treat the general SOT model
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: btorque      !< Spin transfer torque
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: she_btorque  !< SHE spin transfer torque
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: sot_btorque  !< Spin orbit torque 
    
    ! ... Local variables ...
    integer :: i, j, ij, ired  
    real(dblprec) :: lldamp

    ! de/dt (damping band)
    real(dblprec),dimension(3) :: dedt
    
    !!$omp threadprivate(e1x,e1y,e1z,etx,ety,etz,s1x,s1y,s1z,f1x,f1y,f1z,a1x,a1y,a1z,b1x,b1y,b1z,Ax,Ay,Az,detAi,a2x,a2y,a2z,dtg,sqrtdtg)

    ! deterministic variables
    real(dblprec) :: a1x, a1y, a1z
    real(dblprec) :: b1x, b1y, b1z
    !
    ! stochastic variables
    real(dblprec) :: s1x, s1y, s1z
    real(dblprec) :: f1x, f1y, f1z
    !
    ! auxilary variables
    real(dblprec) :: Ax, Ay, Az, detAi
    real(dblprec) :: a2x, a2y, a2z
    !
    ! time steps
    real(dblprec) :: dt,sqrtdt
    real(dblprec) :: dtg, sqrtdtg

    !
    ! spins
    real(dblprec)  :: etx, ety, etz
    real(dblprec)  :: e1x, e1y, e1z

!    dt=deltat*bn*gama !dimm. less time
!    sqrtdt=sqrt(dt)

    !!!$omp parallel do default(shared) schedule(guided,128) &
    !!!$omp private(ij,i,j,e1x,e1y,e1z,etx,ety,etz,s1x,s1y,s1z,&
    !!!$omp f1x,f1y,f1z,a1x,a1y,a1z,b1x,b1y,b1z,Ax,Ay,Az,detAi,a2x,a2y,a2z,dtg,sqrtdtg)

    !!$omp parallel do default(shared) schedule(static) private(i,j,ij)
    !!!$omp parallel do default(shared) schedule(guided,128) private(i,j,ij)

    !print*, ranv
    
    !!$omp parallel do default(shared) schedule(guided,128) &
    !!$omp private(ij,i,j,e1x,e1y,e1z,etx,ety,etz,s1x,s1y,s1z,&
    !!$omp f1x,f1y,f1z,a1x,a1y,a1z,b1x,b1y,b1z,Ax,Ay,Az,detAi,a2x,a2y,a2z,dt,dtg,sqrtdtg,lldamp,dedt)
    
    btorque_full=0.0_dblprec
      !
      if(stt/='N') then
         !$omp parallel do default(shared) private(ired,i,j)  schedule(static) collapse(2)
         do j=1,Mensemble
            do i=1,Natom
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               btorque_full(:,i,j)=btorque_full(:,i,j)+btorque(:,i,j)
            end do
         end do
         !$omp end parallel do
      end if

      if(do_she/='N') then
         !$omp parallel do default(shared) private(ired,i,j)  schedule(static) collapse(2)
         do j=1,Mensemble
            do i=1,Natom
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               btorque_full(:,i,j)= btorque_full(:,i,j)+she_btorque(:,i,j)
            end do
         end do
         !$omp end parallel do
      end if
      if(do_sot/='N') then
         !$omp parallel do default(shared) private(ired,i,j)  schedule(static) collapse(2)
         do j=1,Mensemble
            do i=1,Natom
               ! Adding STT, SHE and SOT torques if present (prefactor instead of if-statement)
              btorque_full(:,i,j)= btorque_full(:,i,j)+sot_btorque(:,i,j)
            end do
         end do
         !$omp end parallel do
      end if
    
    do i=1,Natom
       do j=1,Mensemble

       !i=mod(ij-1,Natom)+1
       !j=int((ij-1)/Natom)+1
       lldamp=1.0D0/(1.0D0+lambda1_array(i)*lambda1_array(i))
       dt=deltat*bn*gama*lldamp !dimm. less time
       sqrtdt=sqrt(dt)
       dtg=dt*Landeg(i)
       sqrtdtg=sqrtdt*Landeg(i)
       e1x=emom(1,i,j) !emom is value previous timestep
       e1y=emom(2,i,j) !in step t emom=emom2 due to copym after step f
       e1z=emom(3,i,j)
       b1x=beff(1,i,j) ! effective field approximated with et
       b1y=beff(2,i,j)
       b1z=beff(3,i,j)
       f1x=ranv(1,i,j) !f1=sqrt(2*D)ksi (fluctuations)
       f1y=ranv(2,i,j)
       f1z=ranv(3,i,j)
       
       ! a1 = -b1 - lambda*(e1 cross b1)
       a1x=-b1x-lambda1_array(i)*(e1y*b1z-e1z*b1y)-btorque_full(1,i,j)
       a1y=-b1y-lambda1_array(i)*(e1z*b1x-e1x*b1z)-btorque_full(2,i,j)
       a1z=-b1z-lambda1_array(i)*(e1x*b1y-e1y*b1x)-btorque_full(3,i,j)
       !

           
       if (dband%enable) then
          dedt = predEdt(i,j,dtg)
          call applyMultiscaleDapingband(i,j,dband,dedt, e1x,e1y,e1z, a1x,a1y,a1z)
       endif
       
       ! s1 is stochastic counterpart of a1
       s1x=-f1x-lambda1_array(i)*(e1y*f1z-e1z*f1y)
       s1y=-f1y-lambda1_array(i)*(e1z*f1x-e1x*f1z)
       s1z=-f1z-lambda1_array(i)*(e1x*f1y-e1y*f1x)


       thermal_field(1,i,j)=s1x
       thermal_field(2,i,j)=s1y
       thermal_field(3,i,j)=s1z
      ! print*,  thermal_field(1,i,j),  thermal_field(2,i,j),  thermal_field(3,i,j)
       
       !
       !
       ! semi-implicitness midpoint requires solution of linear system:
       ! A*e2 = At*e1, At=transpose(A) => e2=inv(A)*At*e1
       ! A = I + skew(dt*a1/2 + sqrt(dt)*s1/2)
       ! write A*e2=a2, a2=At*e1 => e2=inv(A)*a2
       ! Ax,Ay,Az off-diagonal components of A
       ! solve with Cramers' rule => define detAi=1/determinant(A)
       !
       Ax=0.5d0*dtg*a1x + 0.5d0*sqrtdtg*s1x 
       Ay=0.5d0*dtg*a1y + 0.5d0*sqrtdtg*s1y 
       Az=0.5d0*dtg*a1z + 0.5d0*sqrtdtg*s1z 

       detAi=1.0d0/(1.0d0+Ax*Ax+Ay*Ay+Az*Az)
       !
       a2x=e1x+e1y*Az-e1z*Ay
       a2y=e1y+e1z*Ax-e1x*Az
       a2z=e1z+e1x*Ay-e1y*Ax
       !
       etx=a2x*(1+Ax*Ax)+a2y*(Ax*Ay+Az)+a2z*(Ax*Az-Ay)
       ety=a2x*(Ay*Ax-Az)+a2y*(1+Ay*Ay)+a2z*(Ay*Az+Ax)
       etz=a2x*(Az*Ax+Ay)+a2y*(Az*Ay-Ax)+a2z*(1+Az*Az)

       ! now use et for writing et'=(e1+et)/2 in emom2
       etx=0.5d0*(e1x+etx*detAi)
       ety=0.5d0*(e1y+ety*detAi)
       etz=0.5d0*(e1z+etz*detAi)
       !
       ! write et'=(e1+et)/2 in emom2
       emom2(1,i,j)=etx
       emom2(2,i,j)=ety
       emom2(3,i,j)=etz

       ! write new emomM
       ! effective field for step f approximated with et
       ! no new variables (e.g emomt) to reduce memory requirements
       emomM(1,i,j)=etx*mmom(i,j)
       emomM(2,i,j)=ety*mmom(i,j)
       emomM(3,i,j)=etz*mmom(i,j)
    enddo
 enddo
    !!$omp end parallel do
    return

  contains
    !! Approximate de/dt for SIB predictor
    function predEdt(atom,ensemble,dt) result (dedt)
      implicit none
      real(dblprec), intent(in) :: dt
      integer, intent(in)       :: atom,ensemble
      real(dblprec),dimension(3) :: dedt
      
      real(dblprec), parameter :: A = 3.0_dblprec / 2.0_dblprec
      real(dblprec), parameter :: B = -2.0_dblprec
      real(dblprec), parameter :: C = 1.0_dblprec / 2.0_dblprec
      
      real(dblprec),dimension(3) :: numerator
      ! current previous and second previous
      real(dblprec),dimension(3) :: emom_0, emom_1, emom_2 
      integer :: current, prev, prev2

      current = multiscaleBackbufferHead
      prev = multiscaleBackbufferHead-1
      prev2 = multiscaleBackbufferHead-2

      if (prev2 < 1) then
         prev2 = prev2 + ubound(multiscaleBackbuffer,4)
         if (prev < 1) then
            prev = prev + ubound(multiscaleBackbuffer,4)
         end if      
      end if
      
      emom_0 = multiscaleBackbuffer(:,atom,ensemble,current)
      emom_1 = multiscaleBackbuffer(:,atom,ensemble,prev)
      emom_2 = multiscaleBackbuffer(:,atom,ensemble,prev2)
      
      numerator = A*emom_0 + B*emom_1 + C*emom_2
      dedt = numerator / dt
      return
    end function predEdt
    
  end subroutine smodeulermpt_ms


  !> Second step of midpoint solver
  subroutine modeulermpf_ms(Natom, Mensemble, Landeg, bn, lambda1_array, beff, emom, emom2, deltat, dband,STT,do_she,do_sot,&
                            btorque,she_btorque,sot_btorque)
    use Constants
    use RandomNumbers, only : ranv
    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles
    real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
    real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0d0)
    real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
    real(dblprec), intent(in) :: deltat !< Time step
    type(DampingBandData), intent(in) :: dband !< Damping band info (multiscale)
    character(len=1), intent(in) :: STT    !< Treat spin transfer torque?
    character(len=1), intent(in) :: do_she !< Treat the SHE spin transfer torque
    character(len=1), intent(in) :: do_sot !< Treat the general SOT model
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: btorque      !< Spin transfer torque
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: she_btorque  !< SHE spin transfer torque
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: sot_btorque  !< Spin orbit torque
    
    ! deterministic variables
    real(dblprec) :: a1x, a1y, a1z
    real(dblprec) :: b1x, b1y, b1z
    !
    ! stochastic variables
    real(dblprec) :: s1x, s1y, s1z
    real(dblprec) :: f1x, f1y, f1z
    !
    ! auxilary variables
    real(dblprec) :: Ax, Ay, Az, detAi
    real(dblprec) :: a2x, a2y, a2z
    !
    ! time steps
    real(dblprec) :: dt,sqrtdt
    real(dblprec) :: dtg, sqrtdtg
    !
    ! spins
    real(dblprec) :: etpx, etpy, etpz
    real(dblprec) :: e1x, e1y, e1z
    !
    ! ... Local variables ...
    integer :: i, j, ij, ired
    real(dblprec) :: lldamp

    ! de/dt (damping band)
    real(dblprec),dimension(3) :: dedt    

    ! scale dt (from evolve.f90)
!    dt=deltat*bn*gama !dimm. less time
!    sqrtdt=sqrt(dt)
    !
    !!!$omp parallel do default(shared) schedule(guided,128) &
    !!!$omp private(ij,i,j,e1x,e1y,e1z,etpx,etpy,etpz,s1x,s1y,&
    !!!$omp s1z,f1x,f1y,f1z,a1x,a1y,a1z,b1x,b1y,b1z,Ax,Ay,Az,detAi,a2x,a2y,a2z,dtg,sqrtdtg)
    !!$omp parallel do default(shared) schedule(guided,128) private(i,j,ij)
    !!$omp parallel do default(shared) schedule(static) private(i,j,ij)

    !!$omp parallel do default(shared) schedule(guided,128) &
    !!$omp private(ij,i,j,e1x,e1y,e1z,etpx,etpy,etpz,s1x,s1y,&
    !!$omp s1z,f1x,f1y,f1z,a1x,a1y,a1z,b1x,b1y,b1z,Ax,Ay,Az,detAi,a2x,a2y,a2z,dt,dtg,sqrtdtg,lldamp,dedt)
    
    btorque_full=0.0_dblprec
      if(stt/='N') then
         !$omp parallel do default(shared) private(ired,i,j)  schedule(static) collapse(2)
         do j=1,Mensemble
            do i=1,Natom
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               btorque_full(:,i,j)=btorque_full(:,i,j)+btorque(:,i,j)
            end do
         end do
         !$omp end parallel do
      end if

      if(do_she/='N') then
         !$omp parallel do default(shared) private(ired,i,j)  schedule(static) collapse(2)
         do j=1,Mensemble
            do i=1,Natom
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               btorque_full(:,i,j)= btorque_full(:,i,j)+she_btorque(:,i,j)
            end do
         end do
         !$omp end parallel do
      end if
      if(do_sot/='N') then
         !$omp parallel do default(shared) private(ired,i,j)  schedule(static) collapse(2)
         do j=1,Mensemble
            do i=1,Natom
               ! Adding STT, SHE and SOT torques if present (prefactor instead of if-statement)
               btorque_full(:,i,j)= btorque_full(:,i,j)+sot_btorque(:,i,j)
            end do
         end do
         !$omp end parallel do
      end if

    do i=1,Natom
      do j=1,Mensemble
       !i=mod(ij-1,Natom)+1
       !j=int((ij-1)/Natom)+1
       lldamp=1.0D0/(1.0D0+lambda1_array(i)*lambda1_array(i))
       dt=deltat*bn*gama*lldamp !dimm. less time
       sqrtdt=sqrt(dt)
       dtg=dt*Landeg(i)
       sqrtdtg=sqrtdt*Landeg(i)
       e1x=emom(1,i,j) ! load e1 back from emom
       e1y=emom(2,i,j) ! emom unchanged by step t
       e1z=emom(3,i,j)
       etpx=emom2(1,i,j) ! load etp=et' back from emom2
       etpy=emom2(2,i,j)
       etpz=emom2(3,i,j)
       b1x=beff(1,i,j) ! effective field approximated with e1
       b1y=beff(2,i,j)
       b1z=beff(3,i,j)
       f1x=ranv(1,i,j) ! f1=sqrt(2*D)*ksi (fluctuations)
       f1y=ranv(2,i,j)
       f1z=ranv(3,i,j)
       
       ! a1 = -b1 - lambda*(et cross b1)  
       a1x=-b1x-lambda1_array(i)*(etpy*b1z-etpz*b1y)-btorque_full(1,i,j)
       a1y=-b1y-lambda1_array(i)*(etpz*b1x-etpx*b1z)-btorque_full(2,i,j)
       a1z=-b1z-lambda1_array(i)*(etpx*b1y-etpy*b1x)-btorque_full(3,i,j)
       !

       ! Multiscale damping band     
       if (dband%enable) then
          !dedt = (emom2(:,i,j) - emom(:,i,j))/dtg
          dedt = corrEdt(i,j,dtg)
          call applyMultiscaleDapingband(i,j,dband,dedt, etpx,etpy,etpz,&
               a1x,a1y,a1z)
       endif
       
       ! s1 is stochastic counterpart of a1
       s1x=-f1x-lambda1_array(i)*(etpy*f1z-etpz*f1y)
       s1y=-f1y-lambda1_array(i)*(etpz*f1x-etpx*f1z)
       s1z=-f1z-lambda1_array(i)*(etpx*f1y-etpy*f1x)

       ! semi-implicitness midpoint requires solution of linear system:
       ! A*e2 = At*e1, At=transpose(A) => e2=inv(A)*At*e1
       ! A = I + skew(dt*a1/2 + sqrt(dt)*s1/2)
       ! write A*e2=a2, a2=At*e1 => e2=inv(A)*a2
       ! Ax,Ay,Az off-diagonal components of A
       ! solve with Cramers' rule => define detAi=1/determinant(A)
       !
       Ax=0.5d0*dtg*a1x+0.5d0*sqrtdtg*s1x
       Ay=0.5d0*dtg*a1y+0.5d0*sqrtdtg*s1y
       Az=0.5d0*dtg*a1z+0.5d0*sqrtdtg*s1z
       detAi=1.0d0/(1.0d0+Ax*Ax+Ay*Ay+Az*Az)
       !
       a2x=e1x+e1y*Az-e1z*Ay
       a2y=e1y+e1z*Ax-e1x*Az
       a2z=e1z+e1x*Ay-e1y*Ax

       ! now use etp simply to write emom2
       etpx=a2x*(1+Ax*Ax)+a2y*(Ax*Ay+Az)+a2z*(Ax*Az-Ay)
       etpy=a2x*(Ay*Ax-Az)+a2y*(1+Ay*Ay)+a2z*(Ay*Az+Ax)
       etpz=a2x*(Az*Ax+Ay)+a2y*(Az*Ay-Ax)+a2z*(1+Az*Az)
       !
       emom2(1,i,j)=etpx*detAi
       emom2(2,i,j)=etpy*detAi
       emom2(3,i,j)=etpz*detAi
       
    end do
 enddo
    !!$omp end parallel do    
    return

  contains
        !! Approximate de/dt for SIB predictor
    function corrEdt(atom,ensemble,dt) result (dedt)
      implicit none
      real(dblprec), intent(in) :: dt
      integer, intent(in)       :: atom,ensemble
      real(dblprec),dimension(3) :: dedt
      
      real(dblprec), parameter :: A = 2.0_dblprec
      real(dblprec), parameter :: B = -3.0_dblprec
      real(dblprec), parameter :: C = 1.0_dblprec 
      
      real(dblprec),dimension(3) :: numerator
      ! current previous and second previous
      real(dblprec),dimension(3) :: emom_0, emom_1, emom_2 
      integer :: current, prev, prev2

      current = multiscaleBackbufferHead
      prev = multiscaleBackbufferHead-1
      prev2 = multiscaleBackbufferHead-2

      if (prev2 < 1) then
         prev2 = prev2 + ubound(multiscaleBackbuffer,4)
         if (prev < 1) then
            prev = prev + ubound(multiscaleBackbuffer,4)
         end if      
      end if
      
      emom_0 = multiscaleBackbuffer(:,atom,ensemble,current)
      emom_1 = multiscaleBackbuffer(:,atom,ensemble,prev)
      emom_2 = multiscaleBackbuffer(:,atom,ensemble,prev2)
      
      numerator = A*emom_0 + B*emom_1 + C*emom_2
      dedt = numerator / dt
      return
    end function corrEdt

    
  end subroutine modeulermpf_ms

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_midpoint_fields
   !> @brief Allocation of auxilary fields for the treatment of STT and SOT based torques
   !----------------------------------------------------------------------------
   subroutine allocate_midpointms_fields(flag,Natom,Mensemble)

      implicit none

      integer, intent(in) :: flag   !< Allocate or deallocate (1/-1)
      integer, intent(in), optional :: Natom !< Number of atoms in system
      integer, intent(in), optional :: Mensemble   !< Number of ensembles

      integer :: i_stat,i_all

      if (flag>0) then
         allocate(btorque_full(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(btorque_full))*kind(btorque_full),'btorque_full','allocate_midpointms_fields')
         btorque_full=0.0_dblprec
      else
         i_all=-product(shape(btorque_full))*kind(btorque_full)
         deallocate(btorque_full,stat=i_stat)
         call memocc(i_stat,i_all,'btorque_full','allocate_midpointms_fields')
      endif

   end subroutine allocate_midpointms_fields

end module midpoint_ms
