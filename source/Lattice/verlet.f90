!> Velocity-Verlet solver for lattice degrees of freedom and
!! the Grønbech-Jensen Farago (G-JF) velocity Verlet solver for Langevin LD
!! N. Grønbech-Jensen and O. Farago, Molecular Physics 111, 983 (2013)
!! N. Grønbech-Jensen, N. R. Hayre, and O. Farago, Comp. Phys. Comm 185, 524 (2014)
!! 

module Verlet
   use Parameters
   use Profiling
   use RandomDrivers, only : lattranv

   implicit none


contains


   !> Update of displacements
   subroutine u_vverlet(Natom, Mensemble, mioninv, eeff, uvec, vvec, deltat)

      use Constants

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mioninv   !< Inverse ionic mass
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: vvec   !< Current ionic velocity
      real(dblprec), intent(in) :: deltat !< Time step

      ! ... Local variables ...
      integer :: i, j, ij
      real(dblprec) dt, dt2, dt0
      real(dblprec) amuinv
      real(dblprec) angstrominv

      ! Units for Lattice Dynamics

      ! --- Units for input ---
      ! harmonic force constants    ll_inptens           mRyd/angstrom**2 (standard)
      ! harmonic force constants    ll_inptens_phonopy   ev/angstrom**2 (phonopy)
      ! ternary force constants     lll_inptens          mRyd/angstrom**3
      ! quartic force constants     llll_inptens         mRyd/angstrom**4
      ! bilinear ML coupling        ml_inptens           mRyd/angstrom
      ! ternary MML coupling        mml_inptens          mRyd/angstrom
      ! quartic MMLL coupling       mmll_inptens         mRyd/angstrom**2
      ! forces: mRyd/angstrom       (at present external electric field is not available)
      ! magnetic field: B
      ! atomic masses: amu
      ! time: s
      ! displacements: angstrom
      ! velocities: angstrom/s

      ! --- Units for output
      ! harmonic force constants:   mRyd/angstrom**2
      ! ternary force constants:    mRyd/angstrom**3
      ! quartic force constants:    mRyd/angstrom**4
      ! forces: mRyd/angstrom
      ! magnetic field: B
      ! atomic masses: amu
      ! time: s
      ! displacements: angstrom
      ! velocities: angstrom/s
      ! energies: mRyd to be consistent with SD output

      ! --- Internal units ---
      ! harmonic force constants    ll_tens              mRyd/angstrom**2
      ! ternary force constants     lll_tens             mRyd/angstrom**3
      ! quartic force constants     llll_tens            mRyd/angstrom**4
      ! bilinear ML coupling        ml_tens              mRyd/angstrom
      ! ternary MML coupling        mml_tens             mRyd/angstrom
      ! quartic MMLL coupling       mmll_tens            mRyd/angstrom**2
      ! forces: mRyd/angstrom
      ! magnetic field: B
      ! atomic masses: amu
      ! time: s
      ! displacements: angstrom
      ! velocities: angstrom/s
      ! energies: mRyd

      dt0 = 1.0_dblprec
      dt = deltat*dt0
      dt2 = 0.5_dblprec*dt**2
      amuinv = 1.0_dblprec/amu
      angstrominv = 1.0_dblprec/angstrom

      ! Note outdated after change to have mRyd as internal energy unit throughout
      ! -- AB note --
      ! eeff is in eV/angstrom**2 thus it is rescaled by k_force = eV/angstrom**2 to get J/m
      ! then the mass in mioninv is rescaled from amu by k_mass = 1/amu 
      ! what about vvec? 
      !
      !print *,'D',mioninv(1,1),ev * angstrominv**2 * amu , eeff(1,1,1)

      !$omp parallel do default(shared) schedule(static) private(i, j, ij)
      do ij=1,Natom*Mensemble

         i=mod(ij-1,Natom)+1
         j=int((ij-1)/Natom)+1

         uvec(1,i,j) = uvec(1,i,j) + vvec(1,i,j) * dt + mry * eeff(1,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt2
         uvec(2,i,j) = uvec(2,i,j) + vvec(2,i,j) * dt + mry * eeff(2,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt2
         uvec(3,i,j) = uvec(3,i,j) + vvec(3,i,j) * dt + mry * eeff(3,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt2
         !uvec(1,i,j) = uvec(1,i,j) + vvec(1,i,j) * dt + ev * eeff(1,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt2
         !uvec(2,i,j) = uvec(2,i,j) + vvec(2,i,j) * dt + ev * eeff(2,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt2
         !uvec(3,i,j) = uvec(3,i,j) + vvec(3,i,j) * dt + ev * eeff(3,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt2

      end do
      !$omp end parallel do


      10004 format (i8,30es16.8)

   end subroutine u_vverlet


   !> Update of velocities
   subroutine v_vverlet(Natom, Mensemble, mioninv, eeff, e2eff, vvec, deltat)

      use Constants

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mioninv   !< Inverse ionic mass
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: e2eff !< Next step total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: vvec   !< Current ionic velocity
      real(dblprec), intent(in) :: deltat !< Time step

      ! ... Local variables ...
      integer :: i, j, ij
      real(dblprec) dt, dth, dt0
      real(dblprec) amuinv
      real(dblprec) angstrominv

      dt0 = 1.0_dblprec
      dt = deltat*dt0
      dth = 0.5_dblprec*dt
      amuinv = 1.0_dblprec/amu
      angstrominv = 1.0_dblprec/angstrom

      !$omp parallel do default(shared) schedule(static) private(i, j, ij)
      do ij=1,Natom*Mensemble

         i=mod(ij-1,Natom)+1
         j=int((ij-1)/Natom)+1

         vvec(1,i,j) = vvec(1,i,j) + mry * ( eeff(1,i,j) + e2eff(1,i,j) ) * angstrominv**2 * amuinv * mioninv(i,j) * dth
         vvec(2,i,j) = vvec(2,i,j) + mry * ( eeff(2,i,j) + e2eff(2,i,j) ) * angstrominv**2 * amuinv * mioninv(i,j) * dth
         vvec(3,i,j) = vvec(3,i,j) + mry * ( eeff(3,i,j) + e2eff(3,i,j) ) * angstrominv**2 * amuinv * mioninv(i,j) * dth
         !vvec(1,i,j) = vvec(1,i,j) + ev * ( eeff(1,i,j) + e2eff(1,i,j) ) * angstrominv**2 * amuinv * mioninv(i,j) * dth
         !vvec(2,i,j) = vvec(2,i,j) + ev * ( eeff(2,i,j) + e2eff(2,i,j) ) * angstrominv**2 * amuinv * mioninv(i,j) * dth
         !vvec(3,i,j) = vvec(3,i,j) + ev * ( eeff(3,i,j) + e2eff(3,i,j) ) * angstrominv**2 * amuinv * mioninv(i,j) * dth

      end do
      !$omp end parallel do

      10004 format (i8,30es16.8)

   end subroutine v_vverlet


   !> Gronbech-Jensen and Farago: Update of displacements
   subroutine u_gjfverlet(Natom, Mensemble, mioninv, eeff, uvec, vvec, deltat, lattdamp)

      use Constants

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mioninv   !< Inverse ionic mass
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: vvec   !< Current ionic velocity
      real(dblprec), intent(in) :: deltat !< Time step
      real(dblprec), intent(in) :: lattdamp !< Ionic damping constant

      ! ... Local variables ...
      integer :: i, j, ij
      real(dblprec) dt, dth, dt2, dt0
      real(dblprec) amuinv
      real(dblprec) angstrominv
      real(dblprec) bfac

      dt0 = 1.0_dblprec
      dt = deltat*dt0
      dt2 = 0.5_dblprec*dt**2
      dth = 0.5_dblprec*dt
      amuinv = 1.0_dblprec/amu
      angstrominv = 1.0_dblprec/angstrom

      !$omp parallel do default(shared) schedule(static) private(i, j, ij)
      do ij=1,Natom*Mensemble

         i=mod(ij-1,Natom)+1
         j=int((ij-1)/Natom)+1

         ! bfac should be dimensionless. lattdamp has units of kg / s.
         bfac = 1.0_dblprec / ( 1.0_dblprec + lattdamp * dth * mioninv(i,j) * amuinv ) 
         !write(*,'(a,f16.10,es16.8)') 'bfac  ', bfac, bfac

         uvec(1,i,j) = uvec(1,i,j) + bfac *  vvec(1,i,j) * dt &
            + bfac * mry * eeff(1,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt2 &
            !+ bfac * ev * eeff(1,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt2 &
         + bfac * angstrominv * amuinv * mioninv(i,j) * dth * lattranv(1,i,j)

         uvec(2,i,j) = uvec(2,i,j) + bfac *  vvec(2,i,j) * dt &
            + bfac * mry * eeff(2,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt2 &
            !+ bfac * ev * eeff(2,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt2 &
         + bfac * angstrominv * amuinv * mioninv(i,j) * dth * lattranv(2,i,j)

         uvec(3,i,j) = uvec(3,i,j) + bfac * vvec(3,i,j) * dt &
            + bfac * mry * eeff(3,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt2 &
            !+ bfac * ev * eeff(3,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt2 &
         + bfac * angstrominv * amuinv * mioninv(i,j) * dth * lattranv(3,i,j)

      end do
      !$omp end parallel do


      10004 format (i8,30es16.8)

   end subroutine u_gjfverlet


   !> Gronbech-Jensen and Farago: Update of velocities
   subroutine v_gjfverlet(Natom, Mensemble, mioninv, eeff, e2eff, uvec, uvec2, vvec, deltat, lattdamp)

      use Constants

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mioninv   !< Inverse ionic mass
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: e2eff !< Next step total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: uvec2  !< Final (or temporary) ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: vvec   !< Current ionic velocity
      real(dblprec), intent(in) :: deltat !< Time step
      real(dblprec), intent(in) :: lattdamp !< Ionic damping constant

      ! ... Local variables ...
      integer :: i, j, ij
      real(dblprec) dt, dth, dt0
      real(dblprec) amuinv
      real(dblprec) angstrominv
      real(dblprec) a1, a2
      real(dblprec) afac, bfac

      dt0 = 1.0_dblprec
      dt = deltat*dt0
      dth = 0.5_dblprec*dt
      amuinv = 1.0_dblprec/amu
      angstrominv = 1.0_dblprec/angstrom

      !$omp parallel do default(shared) schedule(static) private(i, j, ij)
      do ij=1,Natom*Mensemble

         i=mod(ij-1,Natom)+1
         j=int((ij-1)/Natom)+1

         ! afac and bfac should be dimensionless. lattdamp has units of kg / s.
         ! The bfac is calculated also in u_gjfverlet, consider to use a module array for this variable
         bfac = 1.0_dblprec / ( 1.0_dblprec + lattdamp * dth * amuinv * mioninv(i,j) )
         a1 = 1.0_dblprec - lattdamp * dth * amuinv * mioninv(i,j)
         a2 = 1.0_dblprec + lattdamp * dth * amuinv * mioninv(i,j)
         afac = a1 / a2

         !write(*,'(a,f16.10,a,f16.10)') 'lattdamp ', lattdamp
         !write(*,'(a,f16.10,a,f16.10)') 'afac  ', afac, ' bfac  ', bfac
         !write(*,'(a,3f16.10)') 'lattranv ', lattranv(1:3,i,j)

         !write(*,'(a,4f16.10)') 'vvecx', afac * vvec(1,i,j), &
         !      amuinv * mioninv(i,j) * ( afac * eeff(1,i,j) + e2eff(1,i,j) ) * dth, &
         !      amuinv * mioninv(i,j) * bfac * lattranv(1,i,j)

         vvec(1,i,j) = afac * vvec(1,i,j) & 
            + angstrominv**2 * amuinv * mioninv(i,j) * mry * ( afac * eeff(1,i,j) + e2eff(1,i,j) ) * dth &
            !+ angstrominv**2 * amuinv * mioninv(i,j) * ev * ( afac * eeff(1,i,j) + e2eff(1,i,j) ) * dth &
         + angstrominv * amuinv * mioninv(i,j) * bfac * lattranv(1,i,j)

         vvec(2,i,j) = afac * vvec(2,i,j) & 
            + angstrominv**2 * amuinv * mioninv(i,j) * mry * ( afac * eeff(2,i,j) + e2eff(2,i,j) ) * dth &
            !+ angstrominv**2 * amuinv * mioninv(i,j) * ev * ( afac * eeff(2,i,j) + e2eff(2,i,j) ) * dth &
         + angstrominv * amuinv * mioninv(i,j) * bfac * lattranv(2,i,j)

         vvec(3,i,j) = afac * vvec(3,i,j) & 
            + angstrominv**2 * amuinv * mioninv(i,j) * mry * ( afac * eeff(3,i,j) + e2eff(3,i,j) ) * dth &
            !+ angstrominv**2 * amuinv * mioninv(i,j) * ev * ( afac * eeff(3,i,j) + e2eff(3,i,j) ) * dth &
         + angstrominv * amuinv * mioninv(i,j) * bfac * lattranv(3,i,j)

      end do
      !$omp end parallel do

      10004 format (i8,30es16.8)

   end subroutine v_gjfverlet


   !> Copies uvec2 to uvec, and vvec2 to vvec. A temporary construction,
   !> which is not necessary and should be removed eventually.
   subroutine copy_vverlet(Natom, Mensemble, uvec, uvec2, vvec, vvec2)
      use Constants

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: uvec2  !< Final (or temporary) ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: vvec   !< Current ionic velocity
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: vvec2   !< Final (or temporary) ionic velocity

      ! ... Local variables ...
      integer :: i, j, ij

      do ij=1,Natom*Mensemble

         i=mod(ij-1,Natom)+1
         j=int((ij-1)/Natom)+1

         uvec(1,i,j) = uvec2(1,i,j)
         uvec(2,i,j) = uvec2(2,i,j)
         uvec(3,i,j) = uvec2(3,i,j)
         vvec(1,i,j) = vvec2(1,i,j)
         vvec(2,i,j) = vvec2(2,i,j)
         vvec(3,i,j) = vvec2(3,i,j)

      end do

   end subroutine copy_vverlet


   !> Set average displacement to zero
   !SLDTODO Change to set position of mass centrum to zero
   subroutine set_avrgu0(Natom, Mensemble, uvec, vvec)
      use Constants

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: vvec   !< Current ionic velocity

      ! ... Local variables ...
      integer :: i, j

      real(dblprec) :: avrgux, avrguy, avrguz

      do j=1,Mensemble

         avrgux = 0.0_dblprec
         avrguy = 0.0_dblprec
         avrguz = 0.0_dblprec

         do i=1,Natom

            avrgux = avrgux + uvec(1,i,j)
            avrguy = avrguy + uvec(2,i,j)
            avrguz = avrguz + uvec(3,i,j)

         end do

         avrgux = avrgux / Natom
         avrguy = avrguy / Natom
         avrguz = avrguz / Natom

         do i=1,Natom

            uvec(1,i,j) = uvec(1,i,j) - avrgux
            uvec(2,i,j) = uvec(2,i,j) - avrguy
            uvec(3,i,j) = uvec(3,i,j) - avrguz

         end do

      end do

   end subroutine set_avrgu0


   !> Set average linear momenta to zero
   !SLDTODO right now works with velocities, not momenta!
   subroutine set_avrgp0(Natom, Mensemble, uvec, vvec)
      use Constants

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: vvec   !< Current ionic velocity

      ! ... Local variables ...
      integer :: i, j

      real(dblprec) :: avrgpx, avrgpy, avrgpz

      do j=1,Mensemble

         avrgpx = 0.0_dblprec
         avrgpy = 0.0_dblprec
         avrgpz = 0.0_dblprec

         do i=1,Natom

            avrgpx = avrgpx + vvec(1,i,j)
            avrgpy = avrgpy + vvec(2,i,j)
            avrgpz = avrgpz + vvec(3,i,j)

         end do

         avrgpx = avrgpx / Natom
         avrgpy = avrgpy / Natom
         avrgpz = avrgpz / Natom

         do i=1,Natom

            vvec(1,i,j) = vvec(1,i,j) - avrgpx
            vvec(2,i,j) = vvec(2,i,j) - avrgpy
            vvec(3,i,j) = vvec(3,i,j) - avrgpz

         end do

      end do

   end subroutine set_avrgp0


   !> Rescaling of velocities
   subroutine rescale_vel(Natom, Mensemble, uvec, vvec, alpha)
      use Constants

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: vvec   !< Current ionic velocity
      real(dblprec), intent(in) :: alpha   !< Rescale factor

      ! ... Local variables ...
      integer :: i, j


      do j=1,Mensemble

         do i=1,Natom

            vvec(1,i,j) = vvec(1,i,j) * alpha
            vvec(2,i,j) = vvec(2,i,j) * alpha
            vvec(3,i,j) = vvec(3,i,j) * alpha

         end do

      end do

   end subroutine rescale_vel


  !> Update of displacements
  subroutine u1_vverlet(Natom, Mensemble, mioninv, eeff, uvec, vvec, deltat)

    use Constants

    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles
    real(dblprec), dimension(Natom,Mensemble), intent(in) :: mioninv   !< Inverse ionic mass
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff !< Total effective field from application of Hamiltonian
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: uvec   !< Current ionic displacement
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: vvec   !< Current ionic velocity
    real(dblprec), intent(in) :: deltat !< Time step

    ! ... Local variables ...
    integer :: i, j, ij
    real(dblprec) dt, dt2, dt0
    real(dblprec) amuinv
    real(dblprec) angstrominv

    dt0 = 1.0_dblprec
    dt = deltat*dt0
    dt2 = 0.5_dblprec*dt**2
    amuinv = 1.0_dblprec/amu
    angstrominv = 1.0_dblprec/angstrom

    !$omp parallel do default(shared) schedule(static) private(i, j, ij)
    do ij=1,Natom*Mensemble

       i=mod(ij-1,Natom)+1
       j=int((ij-1)/Natom)+1

       uvec(1,i,j) = uvec(1,i,j) + vvec(1,i,j) * dt
       uvec(2,i,j) = uvec(2,i,j) + vvec(2,i,j) * dt
       uvec(3,i,j) = uvec(3,i,j) + vvec(3,i,j) * dt

    end do
    !$omp end parallel do


10004 format (i8,30es16.8)

  end subroutine u1_vverlet


  !> Update of velocities
  subroutine v1_vverlet(Natom, Mensemble, mioninv, eeff, e2eff, vvec, deltat)

    use Constants

    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles
    real(dblprec), dimension(Natom,Mensemble), intent(in) :: mioninv   !< Inverse ionic mass
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff !< Total effective field from application of Hamiltonian
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: e2eff !< Next step total effective field from application of Hamiltonian
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: vvec   !< Current ionic velocity
    real(dblprec), intent(in) :: deltat !< Time step

    ! ... Local variables ...
    integer :: i, j, ij
    real(dblprec) dt, dth, dt0
    real(dblprec) amuinv
    real(dblprec) angstrominv

    dt0 = 1.0_dblprec
    dt = deltat*dt0
    dth = 0.5_dblprec*dt
    amuinv = 1.0_dblprec/amu
    angstrominv = 1.0_dblprec/angstrom

    !$omp parallel do default(shared) schedule(static) private(i, j, ij)
    do ij=1,Natom*Mensemble

       i=mod(ij-1,Natom)+1
       j=int((ij-1)/Natom)+1

       vvec(1,i,j) = vvec(1,i,j) + mry * eeff(1,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt
       vvec(2,i,j) = vvec(2,i,j) + mry * eeff(2,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt
       vvec(3,i,j) = vvec(3,i,j) + mry * eeff(3,i,j) * angstrominv**2 * amuinv * mioninv(i,j) * dt

    end do
    !$omp end parallel do

10004 format (i8,30es16.8)

  end subroutine v1_vverlet


end module Verlet

