!-------------------------------------------------------------------------------
!> MODULE: SphMidpoint
!> @brief
!> The semi-implicit spherical midpoint solver for the LLG-equations
!> Documented in
!> J. Hellsvik, Progress report: Semi-implicit spherical midpoint solver
!> and constructed by combining the solvers in
!> [1] McLachlan et al, Phys. Rev. E 89, 061301(R) (2014) and
!> [2] J. H. Mentink et al, J. Phys.: Condens. Matter, 22, 176001 (2010)
!> @author Johan Hellsvik
!-------------------------------------------------------------------------------
module SphMidpoint
   use Parameters
   use Profiling

   implicit none

   real(dblprec), dimension(:,:,:), allocatable :: emoma  !< Fix point step k unit moment vector
   real(dblprec), dimension(:,:,:), allocatable :: emomb  !< Fix point step k unit moment vector

   real(dblprec), dimension(:,:,:), allocatable :: uveca  !< Fix point step k displacement vector
   real(dblprec), dimension(:,:,:), allocatable :: uvecb  !< Fix point step k displacement vector
   !
   real(dblprec), dimension(:,:,:), allocatable :: vveca  !< Fix point step k velocity vector
   real(dblprec), dimension(:,:,:), allocatable :: vvecb  !< Fix point step k velocity vector
   real(dblprec), dimension(:,:,:), allocatable :: btorque_full !< Resulting effective field

contains

   !----------------------------------------------------------------------------
   !> SUBROUTINE: sibt
   !> @brief
   !> First step of a fixed point iteration variant of the semi-implicit midpoint solver SIB [2]
   !> @author Johan Hellsvik
   !----------------------------------------------------------------------------
   subroutine sibt(Natom,Mensemble,Landeg,bn,lambda1_array,beff,emom,emom2,emomM,   &
      mmom,deltat,thermal_field)

      use Constants
      use RandomNumbers, only : ranv

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0_dblprec)
      real(dblprec), intent(in) :: deltat !< Time step
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      ! .. In/out variables
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: beff !< Total effective field
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field
      !
      ! ... Local variables ...
      integer :: i, j, ij, k
      real(dblprec) :: lldamp
      !
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
      real(dblprec) :: efpdx, efpdy, efpdz
      real(dblprec) :: efpsx, efpsy, efpsz
      !
      ! time steps
      real(dblprec) :: dt,sqrtdt
      real(dblprec) :: dtg, sqrtdtg
      !
      ! spins
      real(dblprec) :: etx, ety, etz
      real(dblprec) :: e1x, e1y, e1z

      real(dblprec) :: efpx, efpy, efpz
      real(dblprec) :: efp0x, efp0y, efp0z
      real(dblprec) :: efp1x, efp1y, efp1z
      !
      ! tolerance for fix point iteration
      real(dblprec) :: epsilon, epsilona
      real(dblprec) :: ediff, ediff2, ediff3, ediff4
      real(dblprec) :: ediffx, ediffy, ediffz
      real(dblprec) :: enorm
      !
      ! Executable statements
      write(*,*) 'entered sibt'
      ediff    = 0.0_dblprec
      ediff2   = 0.0_dblprec
      ediff3   = 0.0_dblprec
      ediff4   = 0.0_dblprec
      ediffx   = 0.0_dblprec
      ediffy   = 0.0_dblprec
      ediffz   = 0.0_dblprec
      epsilon  = 1e-16
      epsilona = 1e-12
      !
      ! Temporarily deactivated omp-threading for easier debugging
      !
      ! Fix point iteration tolerance check is currently performed for the
      ! Cartesian components. A tolerance check that instead
      ! uses the angle between efp0 and efp, calculated from their
      ! scalar product, is possibly a better choice.
      !
      !!!$omp parallel do default(shared) schedule(guided,128) &
      !!!$omp private(ij,i,j,e1x,e1y,e1z,etx,ety,etz,s1x,s1y,s1z, &
      !!!$omp f1x,f1y,f1z,a1x,a1y,a1z,b1x,b1y,b1z,Ax,Ay,Az,detAi,a2x,a2y,a2z,dt,dtg,sqrtdtg,lldamp, &
      !!!$omp ediff, ediff2, ediffx, ediffy, ediffz, k, &
      !!!$omp efp0x, efp0y, efp0z, efpx, efpy, efpz, efpdx, efpdy, efpdz, efpsx, efpsy, efpsz )
      do ij=1,Natom*Mensemble
         i=mod(ij-1,Natom)+1
         j=int((ij-1)/Natom)+1
         lldamp=1.0_dblprec/(1.0_dblprec+lambda1_array(i)*lambda1_array(i))
         dt=deltat*bn*gama*lldamp !dimension less time
         sqrtdt=sqrt(dt)
         dtg=dt*Landeg(i)
         sqrtdtg=sqrtdt*Landeg(i)
         e1x=emom(1,i,j) ! emom is value previous timestep
         e1y=emom(2,i,j) ! in step t emom=emom2 due to copym after step f
         e1z=emom(3,i,j)
         b1x=beff(1,i,j) ! effective field approximated from emom
         b1y=beff(2,i,j)
         b1z=beff(3,i,j)
         f1x=ranv(1,i,j) ! f1=sqrt(2*D)ksi (fluctuations)
         f1y=ranv(2,i,j)
         f1z=ranv(3,i,j)
         !
         ! a1 = -b1 - lambda*(e1 cross b1)
         a1x=-b1x-lambda1_array(i)*(e1y*b1z-e1z*b1y)
         a1y=-b1y-lambda1_array(i)*(e1z*b1x-e1x*b1z)
         a1z=-b1z-lambda1_array(i)*(e1x*b1y-e1y*b1x)
         !
         ! s1 is stochastic counterpart of a1
         s1x=-f1x-lambda1_array(i)*(e1y*f1z-e1z*f1y)
         s1y=-f1y-lambda1_array(i)*(e1z*f1x-e1x*f1z)
         s1z=-f1z-lambda1_array(i)*(e1x*f1y-e1y*f1x)
         !
         thermal_field(1,i,j)=s1x
         thermal_field(2,i,j)=s1y
         thermal_field(3,i,j)=s1z
         !
         ! semi-implicitness midpoint requires solution of linear system:
         ! A*e2 = At*e1, At=transpose(A) => e2=inv(A)*At*e1
         ! A = I + skew(dt*a1/2 + sqrt(dt)*s1/2)
         ! write A*e2=a2, a2=At*e1 => e2=inv(A)*a2
         ! Ax,Ay,Az off-diagonal components of A
         ! solve with Cramers´ rule => define detAi=1/determinant(A)
         !
         Ax=0.5_dblprec*dtg*a1x+0.5_dblprec*sqrtdtg*s1x
         Ay=0.5_dblprec*dtg*a1y+0.5_dblprec*sqrtdtg*s1y
         Az=0.5_dblprec*dtg*a1z+0.5_dblprec*sqrtdtg*s1z
         detAi=1.0_dblprec/(1.0_dblprec+Ax*Ax+Ay*Ay+Az*Az)
         !
         a2x=e1x+e1y*Az-e1z*Ay
         a2y=e1y+e1z*Ax-e1x*Az
         a2z=e1z+e1x*Ay-e1y*Ax
         !
         etx=a2x*(1+Ax*Ax)+a2y*(Ax*Ay+Az)+a2z*(Ax*Az-Ay)
         ety=a2x*(Ay*Ax-Az)+a2y*(1+Ay*Ay)+a2z*(Ay*Az+Ax)
         etz=a2x*(Az*Ax+Ay)+a2y*(Az*Ay-Ax)+a2z*(1+Az*Az)

         ! now use et for writing et´=(e1+et)/2 in emom2
         etx=0.5_dblprec*(e1x+etx*detAi)
         ety=0.5_dblprec*(e1y+ety*detAi)
         etz=0.5_dblprec*(e1z+etz*detAi)
         !
         ! As an alternative to solve the 3 x 3 linear system exactly:
         ! Use of fix-point iteration. This will later also be
         ! used for the semi-implicit spherical midpoint
         ! method for which no analytical solution is
         ! available for the corresponding 3 x 3 NON-linear system.
         !
         efp0x=e1x
         efp0y=e1y
         efp0z=e1z
         write(*,*) 'starts fixed point iteration'
         do k=1,100

            ! Here efp1 is broken out for efficiency. For the spherical
            ! midpoint solver, the corresponding quantity is the u-vector
            ! as defined in [1].
            efp1x = 0.5_dblprec*(e1x+efp0x)
            efp1y = 0.5_dblprec*(e1y+efp0y)
            efp1z = 0.5_dblprec*(e1z+efp0z)

            ! The left factor in the precession cross product
            efpdx = dtg*efp1x
            efpdy = dtg*efp1y
            efpdz = dtg*efp1z

            ! Its stochastic counterpart
            efpsx = sqrtdtg*efp1x
            efpsy = sqrtdtg*efp1y
            efpsz = sqrtdtg*efp1z

            ! Calculation of the next iterate efp from
            ! deterministic: efpd and a1
            ! stochastic: efpd and s1
            ! contributions
            efpx = e1x + efpdy*a1z - efpdz*a1y + efpsy*s1z - efpsz*s1y
            efpy = e1y + efpdz*a1x - efpdx*a1z + efpsz*s1x - efpsx*s1z
            efpz = e1z + efpdx*a1y - efpdy*a1x + efpsx*s1y - efpsy*s1x

            ! Norms used for the fixed point iteration convergence criteria
            ! The ones which are not used should be commented out eventually
            ediff = sqrt( (efp0x-efpx)**2 + (efp0y-efpy)**2 + (efp0z-efpz)**2 )
            enorm = sqrt( efpx**2 + efpy**2 + efpz**2 )
            ediffx = sqrt( (efp0x-efpx)**2 )
            ediffy = sqrt( (efp0y-efpy)**2 )
            ediffz = sqrt( (efp0z-efpz)**2 )
            !write(*,10002) 'k-step ', k, 'ediff  ', ediff, ediffx, ediffy, ediffz

            ! fixed point iteration convergence criteria
            !if (ediff < epsilon) exit
            if (ediffx < epsilona .and. ediffy < epsilona .and. ediffz < epsilona ) exit

            ! Updates efp0 for the next iteration
            efp0x=efpx
            efp0y=efpy
            efp0z=efpz
         end do

         ! After converged fixed point iteration, store (e1 + efp)/2 in efp
         efpx=0.5_dblprec*(e1x+efpx)
         efpy=0.5_dblprec*(e1y+efpy)
         efpz=0.5_dblprec*(e1z+efpz)

         !write(*,*) 'efp step', ki
         !write(*,10001) 'et   :', etx, ety, etz
         !write(*,10001) 'efpt :', efpx, efpy, efpz
         !ediff2 = sqrt( (efpx-etx)**2 + (efpy-ety)**2 + (efpz-etz)**2 )
         !write(*,10002) 'k-step ', k, 'ediff2 ', ediff2
         !enorm2 = sqrt( etx**2 + ety**2 + etz**2 )
         !write(*,10003) 'efpnorm ', enorm, 'etnorm ', enorm2

         ! write efp in emom2.
         ! (* This corresponds to the write et´=(e1+et)/2 in emom2 *)
         emom2(1,i,j)=efpx
         emom2(2,i,j)=efpy
         emom2(3,i,j)=efpz

         ! write new emomM
         ! effective field for step f approximated with efp
         ! no new variables (e.g emomt) to reduce memory requirements
         emomM(1,i,j)=emom2(1,i,j)*mmom(i,j)
         emomM(2,i,j)=emom2(2,i,j)*mmom(i,j)
         emomM(3,i,j)=emom2(3,i,j)*mmom(i,j)

      enddo
      !!!$omp end parallel do

      return

   end subroutine sibt

   !----------------------------------------------------------------------------
   !> SUBROUTINE: sibf
   !> @brief
   !> Second step of a fixed point iteration variant of the semi-implicit midpoint solver SIB [2]
   !> @author Johan Hellsvik
   !----------------------------------------------------------------------------
   subroutine sibf(Natom,Mensemble,Landeg,bn,lambda1_array,beff,emom,emom2,deltat)

      use Constants
      use RandomNumbers, only : ranv

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0_dblprec)
      real(dblprec), intent(in) :: deltat !< Time step
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field
      !from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      !
      ! ... Local variables ...
      integer :: i, j, ij, k
      real(dblprec) :: lldamp
      !
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
      real(dblprec) :: efpdx, efpdy, efpdz
      real(dblprec) :: efpsx, efpsy, efpsz
      !
      ! time steps
      real(dblprec) :: dt,sqrtdt
      real(dblprec) :: dtg, sqrtdtg
      !
      ! spins

      real(dblprec) :: e1x, e1y, e1z
      real(dblprec) :: etpx, etpy, etpz
      real(dblprec) :: efpx, efpy, efpz
      real(dblprec) :: efp0x, efp0y, efp0z
      real(dblprec) :: efp1x, efp1y, efp1z
      !
      ! tolerance for fix point iteration
      real(dblprec) :: epsilon, epsilona
      real(dblprec) :: ediff, ediff2, ediff3, ediff4
      real(dblprec) :: ediffx, ediffy, ediffz
      real(dblprec) :: enorm
      !
      !
      ! Executable statements
      write(*,*) 'entered sibf'
      epsilon = 1e-16
      epsilona = 1e-12
      ediff = 0.0_dblprec
      ediff2 = 0.0_dblprec
      ediff3 = 0.0_dblprec
      ediff4 = 0.0_dblprec
      ediffx = 0.0_dblprec
      ediffy = 0.0_dblprec
      ediffz = 0.0_dblprec
      !
      ! Temporarily deactivated omp-threading for easier debugging
      !
      ! Fix point iteration tolerance check is currently performed for the
      ! Cartesian components. A tolerance check that instead
      ! uses the angle between efp0 and efp, calculated from their
      ! scalar product, is possibly a better choice.
      !
      !!!$omp parallel do default(shared) schedule(guided,128) &
      !!!$omp private(ij,i,j,e1x,e1y,e1z,etx,ety,etz,etpx,etpy,etpz,s1x,s1y,s1z, &
      !!!$omp f1x,f1y,f1z,a1x,a1y,a1z,b1x,b1y,b1z,Ax,Ay,Az,detAi,a2x,a2y,a2z,dt,dtg,sqrtdtg,lldamp, &
      !!!$omp ediff, ediff2, ediffx, ediffy, ediffz, k, &
      !!!$omp efp0x, efp0y, efp0z, efpx, efpy, efpz, efpdx, efpdy, efpdz, efpsx, efpsy, efpsz )
      ! deterministic variables
      do ij=1,Natom*Mensemble
         i=mod(ij-1,Natom)+1
         j=int((ij-1)/Natom)+1
         lldamp=1.0_dblprec/(1.0_dblprec+lambda1_array(i)*lambda1_array(i))
         dt=deltat*bn*gama*lldamp !dimm. less time
         sqrtdt=sqrt(dt)
         dtg=dt*Landeg(i)
         sqrtdtg=sqrtdt*Landeg(i)
         e1x=emom(1,i,j) ! load e1 back from emom
         e1y=emom(2,i,j) ! emom unchanged by step t
         e1z=emom(3,i,j)
         etpx=emom2(1,i,j) ! load etp=et´ back from emom2
         etpy=emom2(2,i,j)
         etpz=emom2(3,i,j)
         b1x=beff(1,i,j) ! effective field approximated from emom2
         b1y=beff(2,i,j)
         b1z=beff(3,i,j)
         f1x=ranv(1,i,j) ! f1=sqrt(2*D)*ksi (fluctuations)
         f1y=ranv(2,i,j)
         f1z=ranv(3,i,j)
         !
         ! a1 = -b1 - lambda*(et cross b1)
         a1x=-b1x-lambda1_array(i)*(etpy*b1z-etpz*b1y)
         a1y=-b1y-lambda1_array(i)*(etpz*b1x-etpx*b1z)
         a1z=-b1z-lambda1_array(i)*(etpx*b1y-etpy*b1x)
         !
         ! s1 is stochastic counterpart of a1
         s1x=-f1x-lambda1_array(i)*(etpy*f1z-etpz*f1y)
         s1y=-f1y-lambda1_array(i)*(etpz*f1x-etpx*f1z)
         s1z=-f1z-lambda1_array(i)*(etpx*f1y-etpy*f1x)
         !
         ! semi-implicitness midpoint requires solution of linear system:
         ! A*e2 = At*e1, At=transpose(A) => e2=inv(A)*At*e1
         ! A = I + skew(dt*a1/2 + sqrt(dt)*s1/2)
         ! write A*e2=a2, a2=At*e1 => e2=inv(A)*a2
         ! Ax,Ay,Az off-diagonal components of A
         ! solve with Cramers´ rule => define detAi=1/determinant(A)
         !
         Ax=0.5_dblprec*dtg*a1x+0.5_dblprec*sqrtdtg*s1x
         Ay=0.5_dblprec*dtg*a1y+0.5_dblprec*sqrtdtg*s1y
         Az=0.5_dblprec*dtg*a1z+0.5_dblprec*sqrtdtg*s1z
         detAi=1.0_dblprec/(1.0_dblprec+Ax*Ax+Ay*Ay+Az*Az)
         !
         a2x=e1x+e1y*Az-e1z*Ay
         a2y=e1y+e1z*Ax-e1x*Az
         a2z=e1z+e1x*Ay-e1y*Ax
         !
         ! now use etp simply to write emom2
         etpx=a2x*(1+Ax*Ax)+a2y*(Ax*Ay+Az)+a2z*(Ax*Az-Ay)
         etpy=a2x*(Ay*Ax-Az)+a2y*(1+Ay*Ay)+a2z*(Ay*Az+Ax)
         etpz=a2x*(Az*Ax+Ay)+a2y*(Az*Ay-Ax)+a2z*(1+Az*Az)
         !
         ! As an alternative to solve the 3 x 3 linear system exactly:
         ! Use of fix-point iteration. This will later also be
         ! used for the semi-implicit spherical midpoint
         ! method for which no analytical solution is
         ! available for the corresponding 3 x 3 NON-linear system.
         !
         efp0x=e1x
         efp0y=e1y
         efp0z=e1z
         write(*,*) 'starts fixed point iteration'
         do k=1,100

            ! Here efp1 is broken out for efficiency. For the spherical
            ! midpoint solver, the corresponding quantity is the u-vector
            ! as defined in [1].
            efp1x = 0.5_dblprec*(e1x+efp0x)
            efp1y = 0.5_dblprec*(e1y+efp0y)
            efp1z = 0.5_dblprec*(e1z+efp0z)

            ! The left factor in the precession cross product
            efpdx = dtg*efp1x
            efpdy = dtg*efp1y
            efpdz = dtg*efp1z

            ! Its stochastic counterpart
            efpsx = sqrtdtg*efp1x
            efpsy = sqrtdtg*efp1y
            efpsz = sqrtdtg*efp1z

            ! Calculation of the next iterate efp from
            ! deterministic: efpd and a1
            ! stochastic: efpd and s1
            ! contributions
            efpx = e1x + efpdy*a1z - efpdz*a1y + efpsy*s1z - efpsz*s1y
            efpy = e1y + efpdz*a1x - efpdx*a1z + efpsz*s1x - efpsx*s1z
            efpz = e1z + efpdx*a1y - efpdy*a1x + efpsx*s1y - efpsy*s1x

            ! Norms used for the fixed point iteration convergence criteria
            ! The ones which are not used should be commented out eventually
            ediff = sqrt( (efp0x-efpx)**2 + (efp0y-efpy)**2 + (efp0z-efpz)**2 )
            enorm = sqrt( efpx**2 + efpy**2 + efpz**2 )
            ediffx = sqrt( (efp0x-efpx)**2 )
            ediffy = sqrt( (efp0y-efpy)**2 )
            ediffz = sqrt( (efp0z-efpz)**2 )
            !write(*,10002) 'k-step ', k, 'ediff  ', ediff, ediffx, ediffy, ediffz

            ! fixed point iteration convergence criteria
            !if (ediff < epsilon) exit
            if (ediffx < epsilona .and. ediffy < epsilona .and. ediffz < epsilona ) exit

            ! Updates efp0 for the next iteration
            efp0x=efpx
            efp0y=efpy
            efp0z=efpz
         end do

         !write(*,*) 'efp step', ki
         !write(*,10001) 'etp  :', etpx*detAi, etpy*detAi, etpz*detAi
         !write(*,10001) 'efpt :', efpx, efpy, efpz
         !ediff2 = sqrt( (efpx-etpx*detAi)**2 + (efpy-etpy*detAi)**2 + (efpz-etpz*detAi)**2 )
         !write(*,10002) 'k-step ', k, 'ediff2 ', ediff2
         !enorm2 = sqrt( (etpx*detAi)**2 + (etpy*detAi)**2 + (etpz*detAi)**2 )
         !write(*,10003) 'efpnorm ', enorm, 'etnorm ', enorm2

         ! write efp in emom2.
         emom2(1,i,j)=efpx
         emom2(2,i,j)=efpy
         emom2(3,i,j)=efpz

      end do
      !!!$omp end parallel do

      return


   end subroutine sibf

   !----------------------------------------------------------------------------
   !> SUBROUTINE: sispht
   !> @brief
   !> First step of a semi-implicit variant of the spherical midpoint solver [1]
   !> @author Johan Hellsvik
   !----------------------------------------------------------------------------
   subroutine sispht(Natom,Mensemble,Landeg,bn,lambda1_array,beff,emom,emom2,emomM, &
      mmom,deltat,thermal_field)

      use Constants
      use RandomNumbers, only : ranv

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0_dblprec)
      real(dblprec), intent(in) :: deltat !< Time step
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      ! .. In/out variables
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: beff !< Total effective field
      !from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field
      !
      ! ... Local variables ...
      integer :: i, j, ij, k
      real(dblprec) :: lldamp
      !
      ! deterministic variables
      real(dblprec) :: a1x, a1y, a1z
      real(dblprec) :: b1x, b1y, b1z
      !
      ! stochastic variables
      real(dblprec) :: s1x, s1y, s1z
      real(dblprec) :: f1x, f1y, f1z
      !
      ! auxilary variables
      real(dblprec) :: efpdx, efpdy, efpdz
      real(dblprec) :: efpsx, efpsy, efpsz
      !
      ! time steps
      real(dblprec) :: dt,sqrtdt
      real(dblprec) :: dtg, sqrtdtg
      !
      ! spins
      real(dblprec) :: e1x, e1y, e1z
      real(dblprec) :: efpx, efpy, efpz
      real(dblprec) :: efp0x, efp0y, efp0z
      real(dblprec) :: efpux, efpuy, efpuz, unorm
      !
      ! tolerance for fix point iteration
      real(dblprec) :: epsilon, epsilona
      real(dblprec) :: ediff, ediff2, ediff3, ediff4
      real(dblprec) :: ediffx, ediffy, ediffz
      real(dblprec) :: enorm
      !
      !
      ! Executable statements
      !write(*,*) 'entered sispht'
      ediff    = 0.0_dblprec
      ediff2   = 0.0_dblprec
      ediff3   = 0.0_dblprec
      ediff4   = 0.0_dblprec
      ediffx   = 0.0_dblprec
      ediffy   = 0.0_dblprec
      ediffz   = 0.0_dblprec
      epsilon  = 1e-16
      epsilona = 1e-12
      !
      ! Temporarily deactivated omp-threading for easier debugging
      !
      ! Fix point iteration tolerance check is currently performed for the
      ! Cartesian components. A tolerance check that instead
      ! uses the angle between efp0 and efp, calculated from their
      ! scalar product, is possibly a better choice.
      !
      !!!$omp parallel do default(shared) schedule(guided,128) &
      !!!$omp private(ij,i,j,e1x,e1y,e1z,etx,ety,etz,s1x,s1y,s1z, &
      !!!$omp f1x,f1y,f1z,a1x,a1y,a1z,b1x,b1y,b1z,Ax,Ay,Az,detAi,a2x,a2y,a2z,dt,dtg,sqrtdtg,lldamp, &
      !!!$omp ediff, ediff2, ediffx, ediffy, ediffz, k, &
      !!!$omp efp0x, efp0y, efp0z, efpx, efpy, efpz, efpdx, efpdy, efpdz, efpsx, efpsy, efpsz )
      do ij=1,Natom*Mensemble
         i=mod(ij-1,Natom)+1
         j=int((ij-1)/Natom)+1
         lldamp=1.0_dblprec/(1.0_dblprec+lambda1_array(i)*lambda1_array(i))
         dt=deltat*bn*gama*lldamp !dimension less time
         sqrtdt=sqrt(dt)
         dtg=dt*Landeg(i)
         sqrtdtg=sqrtdt*Landeg(i)
         e1x=emom(1,i,j) ! emom is value previous timestep
         e1y=emom(2,i,j) ! in step t emom=emom2 due to copym after step f
         e1z=emom(3,i,j)
         b1x=beff(1,i,j) ! effective field approximated from emom
         b1y=beff(2,i,j)
         b1z=beff(3,i,j)
         f1x=ranv(1,i,j) ! f1=sqrt(2*D)ksi (fluctuations)
         f1y=ranv(2,i,j)
         f1z=ranv(3,i,j)
         !
         ! a1 = -b1 - lambda*(e1 cross b1)
         a1x=-b1x-lambda1_array(i)*(e1y*b1z-e1z*b1y)
         a1y=-b1y-lambda1_array(i)*(e1z*b1x-e1x*b1z)
         a1z=-b1z-lambda1_array(i)*(e1x*b1y-e1y*b1x)
         !
         ! s1 is stochastic counterpart of a1
         s1x=-f1x-lambda1_array(i)*(e1y*f1z-e1z*f1y)
         s1y=-f1y-lambda1_array(i)*(e1z*f1x-e1x*f1z)
         s1z=-f1z-lambda1_array(i)*(e1x*f1y-e1y*f1x)
         !
         thermal_field(1,i,j)=s1x
         thermal_field(2,i,j)=s1y
         thermal_field(3,i,j)=s1z
         !
         ! semi-implicitness midpoint requires solution of linear system:
         ! A*e2 = At*e1, At=transpose(A) => e2=inv(A)*At*e1
         ! A = I + skew(dt*a1/2 + sqrt(dt)*s1/2)
         ! write A*e2=a2, a2=At*e1 => e2=inv(A)*a2
         ! Ax,Ay,Az off-diagonal components of A
         ! solve with Cramers´ rule => define detAi=1/determinant(A)
         !
         !Ax=0.5_dblprec*dtg*a1x+0.5_dblprec*sqrtdtg*s1x
         !Ay=0.5_dblprec*dtg*a1y+0.5_dblprec*sqrtdtg*s1y
         !Az=0.5_dblprec*dtg*a1z+0.5_dblprec*sqrtdtg*s1z
         !detAi=1.0_dblprec/(1.0_dblprec+Ax*Ax+Ay*Ay+Az*Az)
         !
         !a2x=e1x+e1y*Az-e1z*Ay
         !a2y=e1y+e1z*Ax-e1x*Az
         !a2z=e1z+e1x*Ay-e1y*Ax
         !
         !etx=a2x*(1+Ax*Ax)+a2y*(Ax*Ay+Az)+a2z*(Ax*Az-Ay)
         !ety=a2x*(Ay*Ax-Az)+a2y*(1+Ay*Ay)+a2z*(Ay*Az+Ax)
         !etz=a2x*(Az*Ax+Ay)+a2y*(Az*Ay-Ax)+a2z*(1+Az*Az)

         ! now use et for writing et´=(e1+et)/2 in emom2
         !etx=0.5_dblprec*(e1x+etx*detAi)
         !ety=0.5_dblprec*(e1y+ety*detAi)
         !etz=0.5_dblprec*(e1z+etz*detAi)
         !
         ! As an alternative to solve the 3 x 3 linear system exactly:
         ! Use of fix-point iteration. This will later also be
         ! used for the semi-implicit spherical midpoint
         ! method for which no analytical solution is
         ! available for the corresponding 3 x 3 NON-linear system.
         !
         efp0x=e1x
         efp0y=e1y
         efp0z=e1z
         !write(*,*) 'starts fixed point iteration'
         do k=1,100

            ! The u-vector as defined in [1]
            unorm =  sqrt( (e1x+efp0x)**2 + (e1y+efp0y)**2 + (e1z+efp0z)**2 )
            efpux = (e1x+efp0x) / unorm
            efpuy = (e1y+efp0y) / unorm
            efpuz = (e1z+efp0z) / unorm

            ! The left factor in the precession cross product
            efpdx = dtg*efpux
            efpdy = dtg*efpuy
            efpdz = dtg*efpuz

            ! Its stochastic counterpart
            efpsx = sqrtdtg*efpux
            efpsy = sqrtdtg*efpuy
            efpsz = sqrtdtg*efpuz

            ! Calculation of the next iterate efp from
            ! deterministic: efpd and a1
            ! stochastic: efpd and s1
            ! contributions
            efpx = e1x + efpdy*a1z - efpdz*a1y + efpsy*s1z - efpsz*s1y
            efpy = e1y + efpdz*a1x - efpdx*a1z + efpsz*s1x - efpsx*s1z
            efpz = e1z + efpdx*a1y - efpdy*a1x + efpsx*s1y - efpsy*s1x

            ! Norms used for the fixed point iteration convergence criteria
            ! The ones which are not used should be commented out eventually
            ediff = sqrt( (efp0x-efpx)**2 + (efp0y-efpy)**2 + (efp0z-efpz)**2 )
            enorm = sqrt( efpx**2 + efpy**2 + efpz**2 )
            ediffx = sqrt( (efp0x-efpx)**2 )
            ediffy = sqrt( (efp0y-efpy)**2 )
            ediffz = sqrt( (efp0z-efpz)**2 )
            !write(*,10002) 'k-step ', k, 'ediff  ', ediff, ediffx, ediffy, ediffz

            ! fixed point iteration convergence criteria
            !if (ediff < epsilon) exit
            if (ediffx < epsilona .and. ediffy < epsilona .and. ediffz < epsilona ) exit

            ! Updates efp0 for the next iteration
            efp0x=efpx
            efp0y=efpy
            efp0z=efpz
         end do

         ! After converged fixed point iteration, store (e1 + efp) / |e1 + efp| in efpu
         unorm =  sqrt( (e1x+efp0x)**2 + (e1y+efp0y)**2 + (e1z+efp0z)**2 )
         efpux = (e1x+efp0x) / unorm
         efpuy = (e1y+efp0y) / unorm
         efpuz = (e1z+efp0z) / unorm
         ! After converged fixed point iteration, store (e1 + efp)/2 in efp
         !efpx=0.5_dblprec*(e1x+efpx)
         !efpy=0.5_dblprec*(e1y+efpy)
         !efpz=0.5_dblprec*(e1z+efpz)

         !write(*,*) 'efp step', ki
         !write(*,10001) 'et   :', etx, ety, etz
         !write(*,10001) 'efpt :', efpx, efpy, efpz
         !ediff2 = sqrt( (efpx-etx)**2 + (efpy-ety)**2 + (efpz-etz)**2 )
         !write(*,10002) 'k-step ', k, 'ediff2 ', ediff2
         !enorm2 = sqrt( etx**2 + ety**2 + etz**2 )
         !write(*,10003) 'efpnorm ', enorm, 'etnorm ', enorm2

         ! write efp in emom2.
         ! (* This corresponds to the write et´=(e1+et)/2 in emom2 *)
         emom2(1,i,j)=efpux
         emom2(2,i,j)=efpuy
         emom2(3,i,j)=efpuz

         ! write new emomM
         ! effective field for step f approximated with efp
         ! no new variables (e.g emomt) to reduce memory requirements
         emomM(1,i,j)=emom2(1,i,j)*mmom(i,j)
         emomM(2,i,j)=emom2(2,i,j)*mmom(i,j)
         emomM(3,i,j)=emom2(3,i,j)*mmom(i,j)

      enddo
      !!!$omp end parallel do

      return

   end subroutine sispht

   !-----------------------------------------------------------------------------
   !> SUBROUTINE: sisphf
   !> @brief
   !> Second step of a semi-implicit variant of the spherical midpoint solver [1]
   !> @author Johan Hellsvik
   !-----------------------------------------------------------------------------
   subroutine sisphf(Natom,Mensemble,Landeg,bn,lambda1_array,beff,emom,emom2,deltat)

      use Constants
      use RandomNumbers, only : ranv

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0_dblprec)
      real(dblprec), intent(in) :: deltat !< Time step
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field
      !from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      !
      ! ... Local variables ...
      integer :: i, j, ij, k
      real(dblprec) :: lldamp
      !
      ! deterministic variables
      real(dblprec) :: a1x, a1y, a1z
      real(dblprec) :: b1x, b1y, b1z
      !
      ! stochastic variables
      real(dblprec) :: s1x, s1y, s1z
      real(dblprec) :: f1x, f1y, f1z
      !
      ! auxilary variables
      real(dblprec) :: efpdx, efpdy, efpdz
      real(dblprec) :: efpsx, efpsy, efpsz
      !
      ! time steps
      real(dblprec) :: dt,sqrtdt
      real(dblprec) :: dtg, sqrtdtg
      !
      ! spins
      real(dblprec) :: e1x, e1y, e1z
      real(dblprec) :: etpx, etpy, etpz
      real(dblprec) :: efpx, efpy, efpz
      real(dblprec) :: efp0x, efp0y, efp0z
      real(dblprec) :: efpux, efpuy, efpuz, unorm
      !
      ! tolerance for fix point iteration
      real(dblprec) :: epsilon, epsilona
      real(dblprec) :: ediff, ediff2, ediff3, ediff4
      real(dblprec) :: ediffx, ediffy, ediffz
      real(dblprec) :: enorm
      !
      !
      ! Executable statements
      !write(*,*) 'entered sisphf'
      epsilon = 1e-16
      epsilona = 1e-12
      ediff = 0.0_dblprec
      ediff2 = 0.0_dblprec
      ediff3 = 0.0_dblprec
      ediff4 = 0.0_dblprec
      ediffx = 0.0_dblprec
      ediffy = 0.0_dblprec
      ediffz = 0.0_dblprec
      !
      ! Temporarily deactivated omp-threading for easier debugging
      !
      ! Fix point iteration tolerance check is currently performed for the
      ! Cartesian components. A tolerance check that instead
      ! uses the angle between efp0 and efp, calculated from their
      ! scalar product, is possibly a better choice.
      !
      !!!$omp parallel do default(shared) schedule(guided,128) &
      !!!$omp private(ij,i,j,e1x,e1y,e1z,etx,ety,etz,etpx,etpy,etpz,s1x,s1y,s1z, &
      !!!$omp f1x,f1y,f1z,a1x,a1y,a1z,b1x,b1y,b1z,Ax,Ay,Az,detAi,a2x,a2y,a2z,dt,dtg,sqrtdtg,lldamp, &
      !!!$omp ediff, ediff2, ediffx, ediffy, ediffz, k, &
      !!!$omp efp0x, efp0y, efp0z, efpx, efpy, efpz, efpdx, efpdy, efpdz, efpsx, efpsy, efpsz )
      ! deterministic variables
      do ij=1,Natom*Mensemble
         i=mod(ij-1,Natom)+1
         j=int((ij-1)/Natom)+1
         lldamp=1.0_dblprec/(1.0_dblprec+lambda1_array(i)*lambda1_array(i))
         dt=deltat*bn*gama*lldamp !dimm. less time
         sqrtdt=sqrt(dt)
         dtg=dt*Landeg(i)
         sqrtdtg=sqrtdt*Landeg(i)
         e1x=emom(1,i,j) ! load e1 back from emom
         e1y=emom(2,i,j) ! emom unchanged by step t
         e1z=emom(3,i,j)
         etpx=emom2(1,i,j) ! load etp=et´ back from emom2
         etpy=emom2(2,i,j)
         etpz=emom2(3,i,j)
         b1x=beff(1,i,j) ! effective field approximated from emom2
         b1y=beff(2,i,j)
         b1z=beff(3,i,j)
         f1x=ranv(1,i,j) ! f1=sqrt(2*D)*ksi (fluctuations)
         f1y=ranv(2,i,j)
         f1z=ranv(3,i,j)
         !
         ! a1 = -b1 - lambda*(et cross b1)
         a1x=-b1x-lambda1_array(i)*(etpy*b1z-etpz*b1y)
         a1y=-b1y-lambda1_array(i)*(etpz*b1x-etpx*b1z)
         a1z=-b1z-lambda1_array(i)*(etpx*b1y-etpy*b1x)
         !
         ! s1 is stochastic counterpart of a1
         s1x=-f1x-lambda1_array(i)*(etpy*f1z-etpz*f1y)
         s1y=-f1y-lambda1_array(i)*(etpz*f1x-etpx*f1z)
         s1z=-f1z-lambda1_array(i)*(etpx*f1y-etpy*f1x)
         !
         ! semi-implicitness midpoint requires solution of linear system:
         ! A*e2 = At*e1, At=transpose(A) => e2=inv(A)*At*e1
         ! A = I + skew(dt*a1/2 + sqrt(dt)*s1/2)
         ! write A*e2=a2, a2=At*e1 => e2=inv(A)*a2
         ! Ax,Ay,Az off-diagonal components of A
         ! solve with Cramers´ rule => define detAi=1/determinant(A)
         !
         !Ax=0.5_dblprec*dtg*a1x+0.5_dblprec*sqrtdtg*s1x
         !Ay=0.5_dblprec*dtg*a1y+0.5_dblprec*sqrtdtg*s1y
         !Az=0.5_dblprec*dtg*a1z+0.5_dblprec*sqrtdtg*s1z
         !detAi=1.0_dblprec/(1.0_dblprec+Ax*Ax+Ay*Ay+Az*Az)
         !
         !a2x=e1x+e1y*Az-e1z*Ay
         !a2y=e1y+e1z*Ax-e1x*Az
         !a2z=e1z+e1x*Ay-e1y*Ax
         !
         ! now use etp simply to write emom2
         !etpx=a2x*(1+Ax*Ax)+a2y*(Ax*Ay+Az)+a2z*(Ax*Az-Ay)
         !etpy=a2x*(Ay*Ax-Az)+a2y*(1+Ay*Ay)+a2z*(Ay*Az+Ax)
         !etpz=a2x*(Az*Ax+Ay)+a2y*(Az*Ay-Ax)+a2z*(1+Az*Az)
         !
         ! As an alternative to solve the 3 x 3 linear system exactly:
         ! Use of fix-point iteration. This will later also be
         ! used for the semi-implicit spherical midpoint
         ! method for which no analytical solution is
         ! available for the corresponding 3 x 3 NON-linear system.
         !
         efp0x=e1x
         efp0y=e1y
         efp0z=e1z
         !write(*,*) 'starts fixed point iteration'
         do k=1,100

            ! The u-vector as defined in [1]
            unorm =  sqrt( (e1x+efp0x)**2 + (e1y+efp0y)**2 + (e1z+efp0z)**2 )
            efpux = (e1x+efp0x) / unorm
            efpuy = (e1y+efp0y) / unorm
            efpuz = (e1z+efp0z) / unorm

            ! The left factor in the precession cross product
            efpdx = dtg*efpux
            efpdy = dtg*efpuy
            efpdz = dtg*efpuz

            ! Its stochastic counterpart
            efpsx = sqrtdtg*efpux
            efpsy = sqrtdtg*efpuy
            efpsz = sqrtdtg*efpuz

            ! Calculation of the next iterate efp from
            ! deterministic: efpd and a1
            ! stochastic: efpd and s1
            ! contributions
            efpx = e1x + efpdy*a1z - efpdz*a1y + efpsy*s1z - efpsz*s1y
            efpy = e1y + efpdz*a1x - efpdx*a1z + efpsz*s1x - efpsx*s1z
            efpz = e1z + efpdx*a1y - efpdy*a1x + efpsx*s1y - efpsy*s1x

            ! Norms used for the fixed point iteration convergence criteria
            ! The ones which are not used should be commented out eventually
            ediff = sqrt( (efp0x-efpx)**2 + (efp0y-efpy)**2 + (efp0z-efpz)**2 )
            enorm = sqrt( efpx**2 + efpy**2 + efpz**2 )
            ediffx = sqrt( (efp0x-efpx)**2 )
            ediffy = sqrt( (efp0y-efpy)**2 )
            ediffz = sqrt( (efp0z-efpz)**2 )
            !write(*,10002) 'k-step ', k, 'ediff  ', ediff, ediffx, ediffy, ediffz

            ! fixed point iteration convergence criteria
            !if (ediff < epsilon) exit
            if (ediffx < epsilona .and. ediffy < epsilona .and. ediffz < epsilona ) exit

            ! Updates efp0 for the next iteration
            efp0x=efpx
            efp0y=efpy
            efp0z=efpz
         end do

         !write(*,*) 'efp step', ki
         !write(*,10001) 'etp  :', etpx*detAi, etpy*detAi, etpz*detAi
         !write(*,10001) 'efpt :', efpx, efpy, efpz
         !ediff2 = sqrt( (efpx-etpx*detAi)**2 + (efpy-etpy*detAi)**2 + (efpz-etpz*detAi)**2 )
         !write(*,10002) 'k-step ', k, 'ediff2 ', ediff2
         !enorm2 = sqrt( (etpx*detAi)**2 + (etpy*detAi)**2 + (etpz*detAi)**2 )
         !write(*,10003) 'efpnorm ', enorm, 'etnorm ', enorm2

         ! write efp in emom2.
         emom2(1,i,j)=efpx
         emom2(2,i,j)=efpy
         emom2(3,i,j)=efpz

      end do
      !!!$omp end parallel do

      return

   end subroutine sisphf

   !----------------------------------------------------------------------------
   !> SUBROUTINE: imp_fp
   !> @brief
   !> Fix-point iteration step of implicit midpoint solver
   !> @author Johan Hellsvik
   !> @note Jonathan Chico: Added terms dealing with STT, SHE and SOT fields
   !----------------------------------------------------------------------------
   subroutine imp_fp(Natom,Mensemble,Landeg,bn,lambda1_array,beff,emom,emom2,emomM, &
      mmom,deltat,SDEalgh,epsilon,converged,interrupt,mioninv,eeff,   &
      uvec,uvec2,vvec,vvec2)

      use Constants
      use SpinTorques, only : STT,do_sot,do_she,btorque,she_btorque,sot_btorque
      use RandomNumbers, only : ranv

      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: SDEalgh   !< Solver for equations of motion
      integer, intent(in) :: Mensemble !< Number of ensembles
      logical, intent(in) :: interrupt !< Force interrupt
      real(dblprec), intent(in) :: bn        !< Scaling factor for LLG equation (parameter=1.0_dblprec)
      real(dblprec), intent(in) :: deltat    !< Time step
      real(dblprec), intent(in) :: epsilon   !< Fix-point converged criteria
      real(dblprec), dimension(Natom), intent(in) :: Landeg          !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array   !< Damping parameter
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Magnitude of magnetic moments
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mioninv  !< Inverse ionic mass
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff   !< Total effective force from application of Hamiltonian
      ! .. In/out variables
      logical, intent(inout) :: converged !< Fix-point iteration converged or not?
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: beff   !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: vvec   !< Current ionic velocity
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: vvec2  !< Final (or temporary) ionic velocity
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: uvec2  !< Final (or temporary) ionic displacement
      !
      ! ... Local variables ...
      integer :: i, j, ij
      real(dblprec) :: lldamp
      !
      ! deterministic variables
      real(dblprec) :: a1x, a1y, a1z
      real(dblprec) :: b1x, b1y, b1z
      real(dblprec) :: bsttx, bstty, bsttz
      !
      ! stochastic variables
      real(dblprec) :: s1x, s1y, s1z
      real(dblprec) :: f1x, f1y, f1z
      !
      ! auxilary variables
      real(dblprec) :: efpdx, efpdy, efpdz
      real(dblprec) :: efpsx, efpsy, efpsz
      !
      ! time steps
      real(dblprec) :: dt,sqrtdt
      real(dblprec) :: dtg, sqrtdtg
      !
      ! spins
      real(dblprec) :: e1x, e1y, e1z
      real(dblprec) :: etpx, etpy, etpz
      real(dblprec) :: efpx, efpy, efpz
      real(dblprec) :: unorm
      !
      ! tolerance for fix point iteration
      real(dblprec) :: epsilonv
      real(dblprec) :: ediffm, ediffu, ediffv
      real(dblprec) :: enormm, enormu, enormv
      !
      ! lattice variables
      real(dblprec) amuinv
      real(dblprec) angstrominv
      !
      ! Executable statements
      !write(*,*) 'entered imp_fp'
      epsilonv = epsilon * 1e10
      ediffm = 0.0_dblprec
      ediffu = 0.0_dblprec
      ediffv = 0.0_dblprec

      btorque_full=0.0_dblprec
      if(stt/='N') then
         !$omp parallel do default(shared) private(i,j)  schedule(static) collapse(2)
         do j=1,Mensemble
            do i=1,Natom
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               btorque_full(:,i,j)=btorque_full(:,i,j)+btorque(:,i,j)
            end do
         end do
         !$omp end parallel do
      end if

      if(do_she/='N') then
         !$omp parallel do default(shared) private(i,j)  schedule(static) collapse(2)
         do j=1,Mensemble
            do i=1,Natom
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               btorque_full(:,i,j)= btorque_full(:,i,j)+she_btorque(:,i,j)
            end do
         end do
         !$omp end parallel do
      end if
      if(do_sot/='N') then
         !$omp parallel do default(shared) private(i,j)  schedule(static) collapse(2)
         do j=1,Mensemble
            do i=1,Natom
               ! Adding STT, SHE and SOT torques if present (prefactor instead of if-statement)
               btorque_full(:,i,j)= btorque_full(:,i,j)+sot_btorque(:,i,j)
            end do
         end do
         !$omp end parallel do
      end if
      !dt = deltat
      amuinv = 1.0_dblprec/amu
      angstrominv = 1.0_dblprec/angstrom
      !
      ! Temporarily deactivated omp-threading for easier debugging
      !
      !$omp parallel do default(shared) schedule(guided,128) &
      !$omp private(ij, i, j, e1x, e1y, e1z, etpx, etpy, etpz, s1x, s1y, s1z, &
      !$omp f1x, f1y, f1z, a1x, a1y, a1z, b1x, b1y, b1z, dt, dtg, sqrtdtg, lldamp, &
      !$omp efpx, efpy, efpz, efpdx, efpdy, efpdz, efpsx, efpsy, efpsz,bsttx,bstty,bsttz )
      do ij=1,Natom*Mensemble
         i=mod(ij-1,Natom)+1
         j=int((ij-1)/Natom)+1

         lldamp=1.0_dblprec/(1.0_dblprec+lambda1_array(i)*lambda1_array(i))
         dt=deltat*bn*gama*lldamp !dimm. less time
         sqrtdt=sqrt(dt)
         dtg=dt*Landeg(i)
         sqrtdtg=sqrtdt*Landeg(i)

         e1x=emom(1,i,j) ! load e1 back from emom
         e1y=emom(2,i,j) !
         e1z=emom(3,i,j)
         etpx=emom2(1,i,j) ! load etpx back from emom2
         etpy=emom2(2,i,j)
         etpz=emom2(3,i,j)

         b1x=beff(1,i,j) ! effective field approximated from emomM ( i.e. from emomM )
         b1y=beff(2,i,j)
         b1z=beff(3,i,j)
         bsttx=btorque_full(1,i,j)
         bstty=btorque_full(2,i,j)
         bsttz=btorque_full(3,i,j)
         f1x=ranv(1,i,j) ! f1=sqrt(2*D)*ksi (fluctuations)
         f1y=ranv(2,i,j)
         f1z=ranv(3,i,j)

         ! a1 = -b1 - lambda*(et cross b1)
         a1x=-bsttx-b1x-lambda1_array(i)*(etpy*b1z-etpz*b1y)
         a1y=-bstty-b1y-lambda1_array(i)*(etpz*b1x-etpx*b1z)
         a1z=-bsttz-b1z-lambda1_array(i)*(etpx*b1y-etpy*b1x)
         !
         ! s1 is stochastic counterpart of a1
         s1x=-f1x-lambda1_array(i)*(etpy*f1z-etpz*f1y)
         s1y=-f1y-lambda1_array(i)*(etpz*f1x-etpx*f1z)
         s1z=-f1z-lambda1_array(i)*(etpx*f1y-etpy*f1x)
         !
         ! The left factor in the precession cross product
         efpdx = dtg*etpx
         efpdy = dtg*etpy
         efpdz = dtg*etpz

         ! Its stochastic counterpart
         efpsx = sqrtdtg*etpx
         efpsy = sqrtdtg*etpy
         efpsz = sqrtdtg*etpz

         ! Calculation of the next iterate efp from
         ! deterministic: efpd and a1
         ! stochastic: efps and s1
         ! contributions
         efpx = e1x + efpdy*a1z - efpdz*a1y + efpsy*s1z - efpsz*s1y
         efpy = e1y + efpdz*a1x - efpdx*a1z + efpsz*s1x - efpsx*s1z
         efpz = e1z + efpdx*a1y - efpdy*a1x + efpsx*s1y - efpsy*s1x

         !write(*,10001) 'efp  :', efpx, efpy, efpz

         emoma(1,i,j)=efpx
         emoma(2,i,j)=efpy
         emoma(3,i,j)=efpz

         !uvec(1,i,j) = uvec(1,i,j) + vvec(1,i,j) * dt + mry * eeff(1,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt2
         !uvec(2,i,j) = uvec(2,i,j) + vvec(2,i,j) * dt + mry * eeff(2,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt2
         !uvec(3,i,j) = uvec(3,i,j) + vvec(3,i,j) * dt + mry * eeff(3,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt2

         uveca(1,i,j) = uvec(1,i,j) + deltat * vvec2(1,i,j)
         uveca(2,i,j) = uvec(2,i,j) + deltat * vvec2(2,i,j)
         uveca(3,i,j) = uvec(3,i,j) + deltat * vvec2(3,i,j)

         !vvec(1,i,j) = vvec(1,i,j) + mry * ( eeff(1,i,j) + e2eff(1,i,j) ) * angstrominv**2 * amuinv * mioninv(i,j) * dth
         !vvec(2,i,j) = vvec(2,i,j) + mry * ( eeff(2,i,j) + e2eff(2,i,j) ) * angstrominv**2 * amuinv * mioninv(i,j) * dth
         !vvec(3,i,j) = vvec(3,i,j) + mry * ( eeff(3,i,j) + e2eff(3,i,j) ) * angstrominv**2 * amuinv * mioninv(i,j) * dth

         vveca(1,i,j) = vvec(1,i,j) + deltat * mry * eeff(1,i,j) * angstrominv**2 * amuinv * mioninv(i,j)
         vveca(2,i,j) = vvec(2,i,j) + deltat * mry * eeff(2,i,j) * angstrominv**2 * amuinv * mioninv(i,j)
         vveca(3,i,j) = vvec(3,i,j) + deltat * mry * eeff(3,i,j) * angstrominv**2 * amuinv * mioninv(i,j)

      end do
      !$omp end parallel do

      do ij=1,Natom*Mensemble
         i=mod(ij-1,Natom)+1
         j=int((ij-1)/Natom)+1

         ! Norms used for the fixed point iteration convergence criteria
         ! Pythagorean 2-norm not useful for combination of components
         ! with large difference in magnitude
         ediffm = ediffm + ( emoma(1,i,j) - emomb(1,i,j) )**2 &
         + ( emoma(2,i,j) - emomb(2,i,j) )**2 &
         + ( emoma(3,i,j) - emomb(3,i,j) )**2
         ediffu = ediffu + ( uveca(1,i,j) - uvecb(1,i,j) )**2 &
         + ( uveca(2,i,j) - uvecb(2,i,j) )**2 &
         + ( uveca(3,i,j) - uvecb(3,i,j) )**2
         ediffv = ediffv + ( vveca(1,i,j) - vvecb(1,i,j) )**2 &
         + ( vveca(2,i,j) - vvecb(2,i,j) )**2 &
         + ( vveca(3,i,j) - vvecb(3,i,j) )**2

      end do

      enormm = sqrt(ediffm)
      enormu = sqrt(ediffu)
      enormv = sqrt(ediffv)

      ! fixed point iteration convergence criteria
      if ( enormm < epsilon .and. enormu < epsilon .and. enormv < epsilonv) converged = .true.
      if (converged .or. interrupt) then
         !if ( enormm < epsilon .and. enormu < epsilon .and. enormv < epsilon) then

         ! Converged
         !write(*,10001) 'Converged'
         !write(*,10005) 'enormm ', enormm, '  enormu ', enormu, '  enormv ', enormv
         !converged = .true.

         !$omp parallel do default(shared) schedule(guided,128) &
         !$omp private(ij,i,j)
         do ij=1,Natom*Mensemble
            i=mod(ij-1,Natom)+1
            j=int((ij-1)/Natom)+1

            emomb(1,i,j) = emoma(1,i,j)
            emomb(2,i,j) = emoma(2,i,j)
            emomb(3,i,j) = emoma(3,i,j)
            emom2(1,i,j) = emoma(1,i,j)
            emom2(2,i,j) = emoma(2,i,j)
            emom2(3,i,j) = emoma(3,i,j)
            emom(1,i,j) = emoma(1,i,j)
            emom(2,i,j) = emoma(2,i,j)
            emom(3,i,j) = emoma(3,i,j)
            emomM(1,i,j) = emom2(1,i,j)*mmom(i,j)
            emomM(2,i,j) = emom2(2,i,j)*mmom(i,j)
            emomM(3,i,j) = emom2(3,i,j)*mmom(i,j)

            uvecb(1,i,j) = uveca(1,i,j)
            uvecb(2,i,j) = uveca(2,i,j)
            uvecb(3,i,j) = uveca(3,i,j)
            uvec2(1,i,j) = uveca(1,i,j)
            uvec2(2,i,j) = uveca(2,i,j)
            uvec2(3,i,j) = uveca(3,i,j)
            uvec(1,i,j) = uveca(1,i,j)
            uvec(2,i,j) = uveca(2,i,j)
            uvec(3,i,j) = uveca(3,i,j)

            vvecb(1,i,j) = vveca(1,i,j)
            vvecb(2,i,j) = vveca(2,i,j)
            vvecb(3,i,j) = vveca(3,i,j)
            vvec2(1,i,j) = vveca(1,i,j)
            vvec2(2,i,j) = vveca(2,i,j)
            vvec2(3,i,j) = vveca(3,i,j)
            vvec(1,i,j) = vveca(1,i,j)
            vvec(2,i,j) = vveca(2,i,j)
            vvec(3,i,j) = vveca(3,i,j)

         end do
         !$omp end parallel do

      else

         ! Not converged
         !write(*,10001) 'Not converged'
         !write(*,10005) 'enormm ', enormm, '  enormu ', enormu, '  enormv ', enormv
         converged = .false.

         !$omp parallel do default(shared) schedule(guided,128) &
         !$omp private(ij,i,j)
         do ij=1,Natom*Mensemble
            i=mod(ij-1,Natom)+1
            j=int((ij-1)/Natom)+1

            emomb(1,i,j) = emoma(1,i,j)
            emomb(2,i,j) = emoma(2,i,j)
            emomb(3,i,j) = emoma(3,i,j)

            ! This if statement should eventually be moved, by e.g. duplicating
            ! subroutines to one imp and one sph_imp variant.
            if ( SDEalgh== 21) then
               !write(*,*) 'Calculates midpoint epfx'
               ! Calculate the midpoint efpx = ( e1+ efp ) / 2
               emom2(1,i,j) = (emom(1,i,j)+emoma(1,i,j))/2
               emom2(2,i,j) = (emom(2,i,j)+emoma(2,i,j))/2
               emom2(3,i,j) = (emom(3,i,j)+emoma(3,i,j))/2
            elseif (SDEalgh == 22) then
               ! Calculate the u-vector as defined in [1]
               !write(*,*) 'Calculates spherical midpoint epf'
               unorm =  sqrt( (emom(1,i,j)+emoma(1,i,j))**2 &
               + (emom(2,i,j)+emoma(2,i,j))**2 &
               + (emom(3,i,j)+emoma(3,i,j))**2 )
               emom2(1,i,j) = (emom(1,i,j)+emoma(1,i,j)) / unorm
               emom2(2,i,j) = (emom(2,i,j)+emoma(2,i,j)) / unorm
               emom2(3,i,j) = (emom(3,i,j)+emoma(3,i,j)) / unorm
            end if

            emomM(1,i,j) = emom2(1,i,j)*mmom(i,j)
            emomM(2,i,j) = emom2(2,i,j)*mmom(i,j)
            emomM(3,i,j) = emom2(3,i,j)*mmom(i,j)

            uvecb(1,i,j) = uveca(1,i,j)
            uvecb(2,i,j) = uveca(2,i,j)
            uvecb(3,i,j) = uveca(3,i,j)

            uvec2(1,i,j) = (uvec(1,i,j)+uveca(1,i,j))/2
            uvec2(2,i,j) = (uvec(2,i,j)+uveca(2,i,j))/2
            uvec2(3,i,j) = (uvec(3,i,j)+uveca(3,i,j))/2

            vvecb(1,i,j) = vveca(1,i,j)
            vvecb(2,i,j) = vveca(2,i,j)
            vvecb(3,i,j) = vveca(3,i,j)

            vvec2(1,i,j) = (vvec(1,i,j)+vveca(1,i,j))/2
            vvec2(2,i,j) = (vvec(2,i,j)+vveca(2,i,j))/2
            vvec2(3,i,j) = (vvec(3,i,j)+vveca(3,i,j))/2

         end do
         !$omp end parallel do

      end if

   end subroutine imp_fp

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_emomfp
   !> @brief Allocate fixpoint iteration directions of magnetic moments,
   !> displacements, and velocity arrays
   !----------------------------------------------------------------------------
   subroutine allocate_emomfp(Natom, Mensemble, flag)

      implicit none

      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(emoma(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(emoma))*kind(emoma),'emoma','allocate_emomfp')
         allocate(emomb(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(emomb))*kind(emomb),'emomb','allocate_emomfp')
         allocate(uveca(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(uveca))*kind(uveca),'uveca','allocate_emomfp')
         allocate(uvecb(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(uvecb))*kind(uvecb),'uvecb','allocate_emomfp')
         allocate(vveca(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(vveca))*kind(vveca),'vveca','allocate_emomfp')
         allocate(vvecb(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(vvecb))*kind(vvecb),'vvecb','allocate_emomfp')
         allocate(btorque_full(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(btorque_full))*kind(btorque_full),'btorque_full','allocate_emomfp')
         btorque_full=0.0_dblprec
      else
         i_all=-product(shape(emoma))*kind(emoma)
         deallocate(emoma,stat=i_stat)
         call memocc(i_stat,i_all,'emoma','allocate_emomfp')
         i_all=-product(shape(emomb))*kind(emomb)
         deallocate(emomb,stat=i_stat)
         call memocc(i_stat,i_all,'emomb','allocate_emomfp')
         i_all=-product(shape(uveca))*kind(uveca)
         deallocate(uveca,stat=i_stat)
         call memocc(i_stat,i_all,'uveca','allocate_emomfp')
         i_all=-product(shape(uvecb))*kind(uvecb)
         deallocate(uvecb,stat=i_stat)
         call memocc(i_stat,i_all,'uvecb','allocate_emomfp')
         i_all=-product(shape(vveca))*kind(vveca)
         deallocate(vveca,stat=i_stat)
         call memocc(i_stat,i_all,'vveca','allocate_emomfp')
         i_all=-product(shape(vvecb))*kind(vvecb)
         deallocate(vvecb,stat=i_stat)
         call memocc(i_stat,i_all,'vvecb','allocate_emomfp')
         i_all=-product(shape(btorque_full))*kind(btorque_full)
         deallocate(btorque_full,stat=i_stat)
         call memocc(i_stat,i_all,'btorque_full','allocate_emomfp')
      end if

   end subroutine allocate_emomfp

end module sphmidpoint
