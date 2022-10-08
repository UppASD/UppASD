!!!program main
!!!   implicit none
!!!
!!!   call threetemp
!!!
!!!end program main
module temperature_3tm

   use Parameters
   use Profiling
   !integer, parameter :: dblprec = selected_real_kind(15, 307)  !< define precision for double reals

   real(dblprec) :: ps
   real(dblprec) :: gamma_C_elec
   real(dblprec) :: C_elec
   real(dblprec) :: C_latt
   real(dblprec) :: C_spin
   real(dblprec) :: Gel
   real(dblprec) :: Ges
   real(dblprec) :: Gsl
   real(dblprec) :: P_pulse
   real(dblprec) :: t0_pulse
   real(dblprec) :: sigma_pulse
   real(dblprec) :: G_cool
   !
   real(dblprec) :: Temp_final
   real(dblprec) :: Temp_latt_init  !< Initial lattice temperature for 3TM
   real(dblprec) :: Temp_spin_init  !< Initial spin temperature for 3TM
   real(dblprec) :: Temp_elec_init  !< Initial electron temperature for 3TM
   !
   character :: do_cs_temp          !< Use variable Cv for spins (Read from file)
   character :: do_cl_temp          !< Use variable Cv for ions  (Read from file)
   character :: do_ce_temp          !< Use variable Cv for electrons  (Read from file)
   character :: do_ct_temp          !< Use variable Cv for everything  (Read from file)

   integer   :: cs_ntemp            !< Number of read temperatures for Cs(t)
   integer   :: cl_ntemp            !< Number of read temperatures for Cl(t)
   integer   :: ce_ntemp            !< Number of read temperatures for Ce(t)

   integer   :: print_3tm_step      !< Interval for printing temperatures to file

   character(len=35) :: csfile  !< Specific heat for spins
   character(len=35) :: clfile  !< Specific heat for vibrations
   character(len=35) :: ctfile  !< Specific heat for all subsystems

   real(dblprec), dimension(:,:), allocatable   :: cs_array      !< spin Cv as an array
   real(dblprec), dimension(:,:), allocatable   :: cl_array      !< lattice Cv as an array
   real(dblprec), dimension(:,:), allocatable   :: ce_array      !< electron Cv as an array

   real(dblprec), dimension(3) :: xp

   private

   public :: set_temperature_3tm_defaults, read_parameters_3tm, threetemp_single
   public :: threetemp_print, set_initial_temp_3tm, init_3tm_cv, threetemp_elec
   public :: threetemp_elec_print, threetemp_elec_init, unify_3tm_params


contains 

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: set_temperature_3tm_defaults
   !> @brief Initialize default values for 3TM
   !
   !> @author Anastasiia Pervishko, Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine set_temperature_3tm_defaults()

      implicit none

      ! Original values:
      !   ps=1.0E-12, C_elec=6.0E+3, C_latt=2.2E+6, C= 0.7E+6, Gel=8.0E+17, Ges=6.0E+17, Gsl=0.3E+17, pi=3.14
      !ps=1.0e-12_dblprec
      ps=1.0_dblprec
      C_elec=6.0e3_dblprec
      C_latt=2.2e6_dblprec
      C_spin= 0.7e6_dblprec
      Gel=8.0e17_dblprec
      Ges=6.0e17_dblprec
      Gsl=0.3e17_dblprec
      P_pulse=1700.0_dblprec
      t0_pulse=1.0e-12_dblprec
      !t0_pulse=2.0e-12_dblprec
      !sigma_pulse=0.01_dblprec
      sigma_pulse=0.02e-12_dblprec
      G_cool=5.0e10_dblprec
      !
      Temp_final=300.0_dblprec
      Temp_spin_init=-1.0_dblprec
      Temp_latt_init=-1.0_dblprec
      Temp_elec_init=-1.0_dblprec
      !
      do_ct_temp = 'N'
      do_cs_temp = 'N'
      do_cl_temp = 'N'
      do_ce_temp = 'N'

      print_3tm_step = 100

   end subroutine set_temperature_3tm_defaults

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: threetemp_demo
   !> @brief Routine for timestep evolution of temperatures according to 3TM
   !
   !> @author Anastasiia Pervishko
   !-----------------------------------------------------------------------------
   subroutine threetemp_demo()

      implicit none
      integer, parameter :: n=3
      real(dblprec) :: ti, tf, dt, tmax
      real(dblprec), dimension(3) :: xi, xf
      integer :: i

      open (unit=16, file='temperature.dat')


      ps=1.0e-12
      ti =    0.0               ! initial time (ps)
      xi(1) = 300.0             ! initial Te
      xi(2) = 300.0             ! initial TL
      xi(3) = 300.0             ! initial Ts

      dt   = 0.001               ! timestep (ps)
      tmax = 10.0                ! final time (ps)

      write (16,102) ti, xi(1), xi(2), xi(3)

      do while (ti <= tmax)
         tf = ti + dt

         call rk4n(f_3tm_var,ti, tf, xi, xf, n)

         write(16,102) tf, xf(1), xf(2), xf(3)

         ti = tf
         do i = 1,n
            xi(i) = xf(i)
         end do
      end do
      ps=1.0

      102 format(4(1pe12.3))
   end subroutine threetemp_demo

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: threetemp_demo
   !> @brief Routine for timestep evolution of temperatures according to 3TM
   !
   !> @author Anders Bergman, Anastasiia Pervishko
   !-----------------------------------------------------------------------------
   subroutine threetemp_print(rstep,mstep,delta_t,simid)

      implicit none

      integer, intent(in) :: rstep !< Starting simulation step
      integer, intent(in) :: mstep !< Number of simulation steps
      real(dblprec), intent(in) :: delta_t !< Timestep
      character(len=8), intent(in) :: simid !< simulation name

      integer, parameter :: n=3
      real(dblprec) :: ti, tf
      real(dblprec), dimension(3) :: xi, xf
      integer :: i,istep
      character(len=30) :: filn

      write (filn,'(''temperature_3tm.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position='append')


      ti =    rstep*delta_t    ! initial time (ps)
      xi(1) = Temp_elec_init      ! initial Te
      xi(2) = Temp_latt_init      ! initial TL
      xi(3) = Temp_spin_init      ! initial Ts


      write (ofileno,102) ti, xi(1), xi(2), xi(3)

      do istep=rstep,rstep+mstep
         tf = ti + delta_t

         !print '(i8,g20.8,3f16.4)',istep,ti,xi
         call rk4n(f_3tm_var,ti, tf, xi, xf, n)

         write(ofileno,102) tf, xf(1), xf(2), xf(3)

         ti = tf
         do i = 1,n
            xi(i) = xf(i)
         end do
      end do

      close(ofileno)

      102 format(4(1pe16.6))
   end subroutine threetemp_print

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: threetemp_demo
   !> @brief Routine for timestep evolution of temperatures according to 3TM
   !
   !> @author Anders Bergman, Anastasiia Pervishko
   !-----------------------------------------------------------------------------
   subroutine threetemp_print_elec(rstep,mstep,delta_t,simid)

      implicit none

      integer, intent(in) :: rstep !< Starting simulation step
      integer, intent(in) :: mstep !< Number of simulation steps
      real(dblprec), intent(in) :: delta_t !< Timestep
      character(len=8), intent(in) :: simid !< simulation name

      integer, parameter :: n=3
      real(dblprec) :: ti, tf
      real(dblprec), dimension(3) :: xi, xf
      integer :: i,istep
      character(len=30) :: filn

      write (filn,'(''temperature_3tm.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position='append')


      ti =    rstep*delta_t    ! initial time (ps)
      xi(1) = Temp_elec_init      ! initial Te
      xi(2) = Temp_latt_init      ! initial TL
      xi(3) = Temp_spin_init      ! initial Ts


      write (ofileno,102) ti, xi(1), xi(2), xi(3)

      do istep=rstep,rstep+mstep
         tf = ti + delta_t
         xf=xi
         call threetemp_elec(tf,delta_t,xf(3), xf(2), xf(1))

         !call eh_elec(f_3tm_var,ti, tf, xi, xf, n)

         write(ofileno,102) tf, xf(1), xf(2), xf(3)

         ti = tf
         do i = 1,n
            xi(i) = xf(i)
         end do
      end do

      close(ofileno)

      102 format(4(1pe16.6))
   end subroutine threetemp_print_elec

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: threetemp_single
   !> @brief Routine for single timestep evolution of temperatures according to 3TM
   !
   !> @author Anastasiia Pervishko, Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine threetemp_single(t_in,dt,Temp_s,Temp_l,Temp_e)

      implicit none

      real(dblprec), intent(in) :: t_in        !< Current time
      real(dblprec), intent(in) :: dt          !< Time step
      real(dblprec), intent(inout) :: Temp_s   !< Spin temperature
      real(dblprec), intent(inout) :: Temp_l   !< Lattice temperature
      real(dblprec), intent(inout) :: Temp_e   !< Electronic temperature

      integer, parameter :: n=3
      real(dblprec) ::  t_fin
      real(dblprec), dimension(3) :: xi, xf

      xi(1) = Temp_e            ! initial Te
      xi(2) = Temp_l            ! initial TL
      xi(3) = Temp_s            ! initial Ts


      t_fin = t_in + dt

      call rk4n(f_3tm_var,t_in, t_fin, xi, xf, n)

      Temp_e =  xf(1)           ! initial Te
      Temp_l =  xf(2)           ! initial TL
      Temp_s =  xf(3)           ! initial Ts

   end subroutine threetemp_single

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: threetemp_elec
   !> @brief Routine for single timestep evolution of temperatures according to 3TM
   !
   !> @author Anastasiia Pervishko, Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine threetemp_elec(t_in,dt,Temp_s,Temp_l,Temp_e)

      implicit none

      real(dblprec), intent(in) :: t_in        !< Current time
      real(dblprec), intent(in) :: dt          !< Time step
      real(dblprec), intent(inout) :: Temp_s   !< Spin temperature
      real(dblprec), intent(inout) :: Temp_l   !< Lattice temperature
      real(dblprec), intent(inout) :: Temp_e   !< Electronic temperature

      integer, parameter :: n=3
      real(dblprec) ::  t_fin, A, B
      real(dblprec), dimension(3) :: xi, xf

      xi(1) = Temp_e            ! initial Te
      xi(2) = Temp_l            ! initial Tl
      xi(3) = Temp_s            ! initial Ts


!     write(2222,'(2x,g20.6,3f20.8)') t_in,xi

      t_fin = t_in + dt

      A=f_C_latt(xi(2))/f_C_elec(xi(1))
      B=f_C_spin(xi(3))/f_C_elec(xi(1))


      call eh_elec(f_3tm_ele, t_in, t_fin, xi, xp, xf, n)

      !write(*,'(a,2g20.6,3f20.8)') 'Te',Temp_e, xf(1)

      !- A*(x(2)-xp(2)) - B*(x(3) - xp(3))
      !write(*,'(a,2g20.6,3f20.8)') 'AB',A*(Temp_l-xp(2)),B*(Temp_s-xp(3))

      !Temp_e =  xf(1)           ! initial Te
      Temp_e =  xf(1) - A*(Temp_l-xp(2)) - B*(Temp_s-xp(3))
      !!! print '(3f12.6)',xi
      !!! print '(3f12.6)',xp
      !!! print '(3f12.6)',xf
      !!! print *,'------------------------------------------------------------------------------'
      !!   P_pulse/(1.0e-12/ps)*exp(-(t - t0_pulse/ps)**(2)/(2.0_dblprec*(sigma_pulse/ps)**2)) &
      !!- ps*G_cool*(x(1)-Temp_final)

      xp = xi

   end subroutine threetemp_elec

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: threetemp_elec_print
   !> @brief Routine for printing sub-system temperatures
   !
   !> @author Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine threetemp_elec_print(mstep,simid,Temp_s,Temp_l,Temp_e)

      implicit none

      integer, intent(in) :: mstep !< Number of simulation steps
      character(len=8), intent(in) :: simid !< simulation name
      real(dblprec), intent(inout) :: Temp_s   !< Spin temperature
      real(dblprec), intent(inout) :: Temp_l   !< Lattice temperature
      real(dblprec), intent(inout) :: Temp_e   !< Electronic temperature

      character(len=30) :: filn

      if (mod(mstep,print_3tm_step)==0) then
         write (filn,'(''temperature_3tm.'',a,''.out'')') trim(simid)
         open(ofileno,file=filn, position='append')

         write (ofileno,102) mstep, Temp_s,Temp_l,Temp_e

         close(ofileno)

      end if

      return

      102 format(1x,i8,3f12.6)
   end subroutine threetemp_elec_print

!!!    !-----------------------------------------------------------------------------
!!!    ! SUBROUTINE: f_3tm
!!!    !> @brief Temperature functions for the three-temperature model with constant
!!!    !         heat capacities
!!!    !> @author Anastasiia Pervishko
!!!    !-----------------------------------------------------------------------------
!!!    subroutine f_3tm(t, x, dx,n)
!!!       implicit none
!!!       real(dblprec) ::  t
!!!       integer :: n
!!!       real(dblprec), dimension(n) :: x
!!!       real(dblprec), dimension(n) :: dx
!!!       real(dblprec) :: A, B, D, F, G, H
!!! 
!!!  !   real, parameter:: A=ps*Gel/C_elec, B=ps*Ges/C_elec, D=ps*Gel/C_latt,F=ps*Gsl/C_latt, G=ps*Ges/C_spin, H=ps*Gsl/C_spin
!!!       Gel_ce=Gel/f_C_elec(x(1)) !/8e3    ! Electron-Lattice coupling
!!!       Ges_ce=Ges/f_C_elec(x(1)) !/8e3    ! Electron-Spin coupling
!!!       Gel_cl=Gel/f_C_latt(x(2)) !/8e4    ! Lattice-Electron coupling
!!!       Gsl_cl=Gsl/f_C_latt(x(2)) !/8e4    ! Lattice-Spin coupling
!!!       Ges_cs=Ges/f_C_spin(x(3)) !/8e4    ! Spin-Electron coupling
!!!       Gsl_cs=Gsl/f_C_spin(x(3)) !/8e4    ! Spin-lattice coupling
!!!       Gep_ce=P_pulse/f_C_elec(x(1))      ! Electron-Pulse coupling
!!!       Gec_ce=G_cool/f_C_elec(x(1))       ! Electron-Ambience coupling
!!!       
!!!       dx(1) = -Gel_ce*(x(1)-x(2)) -Ges_ce*(x(1)-x(3)) &
!!!       +  Gep_ce*exp(-(t - t0_pulse)**(2)/(2.0_dblprec*sigma_pulse**2)) &
!!!       -  Gec_ce*(x(1)-Temp_final)
!!!       dx(2) = Gel_cl*(x(1)-x(2)) - Gsl_cl*(x(2)-x(3))
!!!       dx(3) = Ges_cs*(x(1)-x(3)) + Gsl_cs*(x(2)-x(3))
!!! 
!!!       A=ps*Gel/(gamma_C_elec*x(1))
!!!       B=ps*Ges/(gamma_C_elec*x(1))
!!!       D=ps*Gel/C_latt
!!!       F=ps*Gsl/C_latt
!!!       G=ps*Ges/C_spin
!!!       H=ps*Gsl/C_spin
!!! 
!!!       dx(1) = -A*(x(1)-x(2)) - B*(x(1)-x(3)) + &
!!!          P_pulse/(1.0e-12/ps)*exp(-(t - t0_pulse/ps)**(2)/(2.0_dblprec*(sigma_pulse/ps)**2)) &
!!!       - ps*G_cool*(x(1)-Temp_final)
!!!       dx(2) = D*(x(1)-x(2)) - F*(x(2)-x(3))
!!!       dx(3) = G*(x(1)-x(3)) + H*(x(2)-x(3))
!!! 
!!!    end subroutine f_3tm

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: f_3tm_var
   !> @brief Temperature functions for the three-temperature model for variable
   !         heat capacities
   !> @author Anastasiia Pervishko, Anders Bergman
   !> @note: Update Nov 2021: 
   !> 1) Removed prefactors for coupling constants
   !> 2) Removed division with T_e in dT_e/dt expression 
   !> (Here we should not assume C_e = K * T_e 
   !> 3) Removed time factor ps 
   !-----------------------------------------------------------------------------
   subroutine f_3tm_var(t, x, dx,n)
      implicit none
      real(dblprec) ::  t
      integer :: n
      real(dblprec), dimension(n) :: x
      real(dblprec), dimension(n) :: dx
      real(dblprec) :: Gel_ce, Ges_ce, Gel_cl, Gsl_cl, Ges_cs, Gsl_cs, Gep_ce, Gec_ce

      Gel_ce=Gel/f_C_elec(x(1)) !/8e3    ! Electron-Lattice coupling
      Ges_ce=Ges/f_C_elec(x(1)) !/8e3    ! Electron-Spin coupling
      Gel_cl=Gel/f_C_latt(x(2)) !/8e4    ! Lattice-Electron coupling
      Gsl_cl=Gsl/f_C_latt(x(2)) !/8e4    ! Lattice-Spin coupling
      Ges_cs=Ges/f_C_spin(x(3)) !/8e4    ! Spin-Electron coupling
      Gsl_cs=Gsl/f_C_spin(x(3)) !/8e4    ! Spin-lattice coupling
      Gep_ce=P_pulse/f_C_elec(x(1))      ! Electron-Pulse coupling
      Gec_ce=G_cool/f_C_elec(x(1))       ! Electron-Ambience coupling
      
      !!! print '(a,6g14.6)', 'T_ele ', dx(1)
      !!! print '(a,6g14.6)', 'Couplings:',Gel, Ges, Gsl, P_pulse, G_cool
      !!! print '(a,6g14.6)', 'Capacities:',f_C_elec(x(1)), f_C_latt(x(2)), f_C_spin(x(3))
      !!! print '(a,6g14.6)', 'Effective couplings:',Gel_ce, Ges_ce, Gep_ce, Gec_ce
      dx(1) = -Gel_ce*(x(1)-x(2)) -Ges_ce*(x(1)-x(3)) &
      +  Gep_ce*exp(-(t - t0_pulse)**(2)/(2.0_dblprec*sigma_pulse**2)) &
      -  Gec_ce*(x(1)-Temp_final)
      dx(2) = Gel_cl*(x(1)-x(2)) - Gsl_cl*(x(2)-x(3))
      dx(3) = Ges_cs*(x(1)-x(3)) + Gsl_cs*(x(2)-x(3))

   end subroutine f_3tm_var

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: f_3tm_ele
   !> @brief Temperature functions for the three-temperature model for variable
   !         heat capacities
   !> @author Anastasiia Pervishko, Anders Bergman
   !> @note Edit Nov 21: Removed ps time scaling
   !-----------------------------------------------------------------------------
   subroutine f_3tm_ele(t, x, xp, dx,n)
      implicit none
      real(dblprec), intent(in) ::  t
      real(dblprec), dimension(n) :: x
      real(dblprec), dimension(n) :: xp
      real(dblprec), dimension(n) :: dx
      integer :: n

      real(dblprec) :: Gep_ce, Gec_ce

      Gep_ce=P_pulse/f_C_elec(x(1))      ! Electron-Pulse coupling
      Gec_ce=G_cool/f_C_elec(x(1))       ! Electron-Ambience coupling

      dx(1) = Gep_ce*exp(-(t - t0_pulse)**(2)/(2.0_dblprec*sigma_pulse**2)) - Gec_ce*(x(1)-Temp_final)  

   end subroutine f_3tm_ele


   !-----------------------------------------------------------------------------
   ! SUBROUTINE: rk4n
   !> @brief Runge-Kutta solver for function fcn
   !
   !> @author Anastasiia Pervishko
   !-----------------------------------------------------------------------------
   subroutine rk4n(fcn,ti, tf, xi, xf, n)
      implicit none
      real(dblprec) :: ti, tf
      integer :: n 
      real(dblprec), dimension(n) :: xi
      real(dblprec), dimension(n) :: xf
      external :: fcn

      integer  :: j
      real(dblprec) h, t
      real(dblprec), dimension(n) ::  x
      real(dblprec), dimension(n) ::  dx
      real(dblprec), dimension(n) ::  k1,k2,k3,k4

      h = tf-ti
      t = ti

      call fcn(t, xi, dx, n)
      do j=1,n
         k1(j) = h*dx(j)
         x(j)  = xi(j) + k1(j)/2.0_dblprec
      end do

      call fcn(t+h/2.0_dblprec, x, dx, n)
      do j=1,n
         k2(j) = h*dx(j)
         x(j)  = xi(j) + k2(j)/2.0_dblprec
      end do

      call fcn(t+h/2.0_dblprec, x, dx, n)
      do j=1,n
         k3(j) = h*dx(j)
         x(j)  = xi(j) + k3(j)
      end do

      call fcn(t+h, x, dx, n)
      do j=1,n
         k4(j) = h*dx(j)
         xf(j) = xi(j) + k1(j)/6.0_dblprec+k2(j)/3.0_dblprec+k3(j)/3.0_dblprec+k4(j)/6.0_dblprec
      end do

   end subroutine rk4n

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: eh_elec
   !> @brief Euler-Heun  solver for function fcn. Custom made to only evolve
   !>        the first of all variables (application: T_elec in 3TM)
   !
   !> @author Anders Bergman (based of rk4n from Anastasiia Pervishko)
   !-----------------------------------------------------------------------------
   subroutine eh_elec(fcn,ti, tf, xi, xp, xf, n)
      implicit none
      real(dblprec) :: ti, tf
      integer :: n 
      real(dblprec), dimension(n) :: xi
      real(dblprec), dimension(n) :: xp
      real(dblprec), dimension(n) :: xf
      external :: fcn

      integer  :: j
      real(dblprec) h, t
      real(dblprec), dimension(n) ::  x
      real(dblprec), dimension(n) ::  dx_pred
      real(dblprec), dimension(n) ::  dx_corr
      real(dblprec), dimension(n) ::  k1,k2,k3,k4

      h = tf-ti
      t = ti

      ! Calculate predictor derivatives
      call fcn(t, xi, xp, dx_pred, n)

      ! Create predictor
      x(1) = xi(1) + h*dx_pred(1) 
      x(2:n) = xi(2:n)

      ! Call calculate corrector derivatives
      call fcn(t+h/2.0_dblprec, x, xp, dx_corr, n)

      xf(1) = xi(1) + h/2.0_dblprec * (dx_pred(1) + dx_corr(1))
      xf(2:n) = x(2:n)

   end subroutine eh_elec

   function f_C_latt(temp) result(C_l)
      !
      use Math_Functions, only : f_interp_1d
      !
      implicit none
      ! 
      real(dblprec), intent(in) :: temp
      real(dblprec) :: C_l
      !
      if (do_cl_temp/='Y') then
         C_l=C_latt
      else
         C_l=f_interp_1d(temp,cl_array(:,1),cl_array(:,2),cl_ntemp)
      end if

      return
      !
   end function f_C_latt

   function f_C_spin(temp) result(C_s)
      !
      use Math_Functions, only : f_interp_1d
      !
      implicit none
      ! 
      real(dblprec), intent(in) :: temp
      real(dblprec) :: C_s
      !
      if (do_cs_temp/='Y') then
         C_s=C_spin
      else
         C_s=f_interp_1d(temp,cs_array(:,1),cs_array(:,2),cs_ntemp)
      end if

      return
      !
   end function f_C_spin

   function f_C_elec(temp) result(C_e)
      !
      use Math_Functions, only : f_interp_1d
      !
      implicit none
      ! 
      real(dblprec), intent(in) :: temp
      real(dblprec) :: C_e
      !
      if (do_ce_temp/='Y') then
         C_e=gamma_C_elec*temp
      else
         C_e=f_interp_1d(temp,ce_array(:,1),ce_array(:,2),ce_ntemp)
      end if

      return
      !
   end function f_C_elec


   !!!    !-----------------------------------------------------------------------------
   !!!    ! SUBROUTINE: allocate_temp
   !!!    !> @brief Subroutine to allocate the arrays that contain the temperature for
   !!!    !> the three-temperature model
   !!!    !> @author Anders Bergman
   !!!    !-----------------------------------------------------------------------------
   !!!    subroutine allocate_3tm(natom, ip_nphase,flag)
   !!! 
   !!!       use Profiling 
   !!! 
   !!!       implicit none
   !!! 
   !!!       integer, intent(in) :: natom
   !!!       integer, intent(in) :: ip_nphase
   !!!       integer, intent(in) :: flag
   !!! 
   !!!       integer :: i_stat,i_all
   !!! 
   !!!       if(flag>=0) then
   !!!          allocate(temp_s_array(natom),stat=i_stat)
   !!!          call memocc(i_stat,product(shape(temp_s_array))*kind(temp_s_array),'temp_s_array','allocate_3tm')
   !!!          allocate(temp_l_array(natom),stat=i_stat)
   !!!          call memocc(i_stat,product(shape(temp_l_array))*kind(temp_l_array),'temp_l_array','allocate_3tm')
   !!!          allocate(temp_e_array(natom),stat=i_stat)
   !!!          call memocc(i_stat,product(shape(temp_e_array))*kind(temp_e_array),'temp_e_array','allocate_3tm')
   !!!          if(ip_nphase>0) then
   !!!             allocate(iptemp_s_array(natom,ip_nphase),stat=i_stat)
   !!!             call memocc(i_stat,product(shape(iptemp_s_array))*kind(iptemp_s_array),'iptemp_s_array','allocate_3tm')
   !!!             allocate(iptemp_l_array(natom,ip_nphase),stat=i_stat)
   !!!             call memocc(i_stat,product(shape(iptemp_l_array))*kind(iptemp_l_array),'iptemp_l_array','allocate_3tm')
   !!!             allocate(iptemp_e_array(natom,ip_nphase),stat=i_stat)
   !!!             call memocc(i_stat,product(shape(iptemp_e_array))*kind(iptemp_e_array),'iptemp_e_array','allocate_3tm')
   !!!          end if
   !!!       else
   !!!          i_all=-product(shape(temp_s_array))*kind(temp_s_array)
   !!!          deallocate(temp_s_array,stat=i_stat)
   !!!          call memocc(i_stat,i_all,'temp_s_array','allocate_3tm')
   !!!          i_all=-product(shape(temp_l_array))*kind(temp_l_array)
   !!!          deallocate(temp_l_array,stat=i_stat)
   !!!          call memocc(i_stat,i_all,'temp_l_array','allocate_3tm')
   !!!          i_all=-product(shape(temp_e_array))*kind(temp_e_array)
   !!!          deallocate(temp_e_array,stat=i_stat)
   !!!          call memocc(i_stat,i_all,'temp_e_array','allocate_3tm')
   !!!          if(ip_nphase>0) then
   !!!             i_all=-product(shape(iptemp_s_array))*kind(iptemp_s_array)
   !!!             deallocate(iptemp_s_array,stat=i_stat)
   !!!             call memocc(i_stat,i_all,'iptemp_s_array','allocate_3tm')
   !!!             i_all=-product(shape(iptemp_l_array))*kind(iptemp_l_array)
   !!!             deallocate(iptemp_l_array,stat=i_stat)
   !!!             call memocc(i_stat,i_all,'iptemp_l_array','allocate_3tm')
   !!!             i_all=-product(shape(iptemp_e_array))*kind(iptemp_e_array)
   !!!             deallocate(iptemp_e_array,stat=i_stat)
   !!!             call memocc(i_stat,i_all,'iptemp_e_array','allocate_3tm')
   !!!          end if
   !!!       end if
   !!! 
   !!!    end subroutine allocate_3tm

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: set_initial_temp_3tm
   !> @brief Set initial temperatures for 3TM 
   !
   !> @author Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine set_initial_temp_3tm(Temp,Temp_s,Temp_l,Temp_e)

      implicit none

      real(dblprec), intent(in) :: Temp   !< Global temperature for non-3TM simulations
      real(dblprec), intent(out) :: Temp_s !< Global temperature for non-3TM simulations
      real(dblprec), intent(out) :: Temp_l !< Global temperature for non-3TM simulations
      real(dblprec), intent(out) :: Temp_e !< Global temperature for non-3TM simulations

      ! If initial temperatures are not specified then 
      ! set all initial temperaturees to global Temp
      if (Temp_latt_init<0.0_dblprec) Temp_latt_init=Temp
      if (Temp_spin_init<0.0_dblprec) Temp_spin_init=Temp
      if (Temp_elec_init<0.0_dblprec) Temp_elec_init=Temp

      Temp_l=Temp_latt_init
      Temp_s=Temp_spin_init
      Temp_e=Temp_elec_init


      xp(1) = Temp_e
      xp(2) = Temp_l
      xp(3) = Temp_s

   end subroutine set_initial_temp_3tm

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: set_initial_temp_3tm
   !> @brief Set initial temperatures for 3TM 
   !
   !> @author Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine threetemp_elec_init(simid,Temp,Temp_s,Temp_l,Temp_e)

      implicit none

      character(len=8), intent(in) :: simid !< simulation name
      real(dblprec), intent(in) :: Temp   !< Global temperature for non-3TM simulations
      real(dblprec), intent(in) :: Temp_s !< Global temperature for non-3TM simulations
      real(dblprec), intent(in) :: Temp_l !< Global temperature for non-3TM simulations
      real(dblprec), intent(inout) :: Temp_e !< Global temperature for non-3TM simulations
      character(len=30) :: filn

      if (Temp_elec_init<0.0_dblprec) Temp_elec_init=Temp

      Temp_e = Temp_elec_init
      xp(1) = Temp_e
      xp(2) = Temp_l
      xp(3) = Temp_s

      write (filn,'(''temperature_3tm.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position='append')

      write(ofileno,'(a)') "#  Timestep    T_spin      T_latt      T_elec"

      write (ofileno,102) 0, Temp_s,Temp_l,Temp_e

      102 format(1x,i8,3f12.6)
   end subroutine threetemp_elec_init
   !---------------------------------------------------------------------------
   ! SUBROUTINE: read_parameters_3tm
   !> @brief
   !> Read input parameters for the three-temperature model
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_parameters_3tm(ifile)

      use FileParser
      use ErrorHandling

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword, cache
      integer :: rd_len, i_err, i_errb
      logical :: comment

      do
         10      continue
         ! Read file character for character until first whitespace
         keyword=""
         call bytereader(keyword,rd_len,ifile,i_errb)

         ! converting Capital letters
         call caps2small(keyword)

         ! check for comment markers (currently % and #)
         comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&
            (scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1.or.&
            (scan(trim(keyword),'!')==1))

         if (comment) then
            read(ifile,*)
         else
            ! Parse keyword
            keyword=trim(keyword)
            select case(keyword)
            ! This is the flags for the ElkGeometry module

         case('gamma_ce')
            read(ifile,*,iostat=i_err) gamma_C_elec
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('ce')
            read(ifile,*,iostat=i_err) C_elec
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('cl')
            read(ifile,*,iostat=i_err) C_latt
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('cs')
            read(ifile,*,iostat=i_err) C_spin
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('gel')
            read(ifile,*,iostat=i_err) Gel
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('gsl')
            read(ifile,*,iostat=i_err) Gsl
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('ges')
            read(ifile,*,iostat=i_err) Ges
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('p_pulse')
            read(ifile,*,iostat=i_err) P_pulse
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('t0_pulse')
            read(ifile,*,iostat=i_err) t0_pulse
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('sigma_pulse')
            read(ifile,*,iostat=i_err) sigma_pulse
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('temp_spin_init')
            read(ifile,*,iostat=i_err) Temp_spin_init
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('temp_latt_init')
            read(ifile,*,iostat=i_err) Temp_latt_init
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('temp_elec_init')
            read(ifile,*,iostat=i_err) Temp_elec_init
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('temp_final')
            read(ifile,*,iostat=i_err) Temp_final
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('g_cool')
            read(ifile,*,iostat=i_err) G_cool
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_ct_temp')
            read(ifile,*,iostat=i_err) do_ct_temp
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_cs_temp')
            read(ifile,*,iostat=i_err) do_cs_temp
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_cl_temp')
            read(ifile,*,iostat=i_err) do_cl_temp
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('cv_file')
            do_ct_temp='Y'
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            ctfile=trim(adjustl(cache))
            call ErrorHandling_check_file_exists(ctfile, &
               'Please specify cv_file <file> where <file> is a valid heat capacity file')

         case('cv_spinfile')
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            csfile=trim(adjustl(cache))
            call ErrorHandling_check_file_exists(csfile, &
               'Please specify cv_spinfile <file> where <file> is a valid heat capacity file')

         case('cv_lattfile')
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            clfile=trim(adjustl(cache))
            call ErrorHandling_check_file_exists(clfile, &
               'Please specify cv_lattfile <file> where <file> is a valid heat capacity file')

         case('print_3tm_step')
            read(ifile,*,iostat=i_err) print_3tm_step
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         end select
      end if

      ! End of file
      if (i_errb==20) goto 20
      ! End of row
      if (i_errb==10) goto 10
   end do

   20   continue

   rewind(ifile)
   return
end subroutine read_parameters_3tm

subroutine init_3tm_cv(flag)
   !
   implicit none
   !
   integer, intent(in) :: flag
   !
   integer :: itemp, i_stat, i_all
   real(dblprec) :: tdum,cdum
   !

   if (flag>0) then
      if (do_ct_temp=='Y') then
            open(ifileno,file=ctfile)

            ! Find the number of spin temperatures 
            cs_ntemp=0
            do 
               read(ifileno,*,end=20) tdum, cdum
               cs_ntemp=cs_ntemp+1
            end do
            20         continue
            cl_ntemp=cs_ntemp
            ce_ntemp=cs_ntemp

            rewind(ifileno)

            ! Allocate array for spin Cv
            allocate(cs_array(cs_ntemp,2),stat=i_stat)
            call memocc(i_stat,product(shape(cs_array))*kind(cs_array),'cs_array','init_3tm_cv')
            allocate(cl_array(cs_ntemp,2),stat=i_stat)
            call memocc(i_stat,product(shape(cl_array))*kind(cl_array),'cl_array','init_3tm_cv')
            allocate(ce_array(cs_ntemp,2),stat=i_stat)
            call memocc(i_stat,product(shape(ce_array))*kind(ce_array),'ce_array','init_3tm_cv')

            ! Read all input spin temperatures and heat capacities
            ! array(:,1) contains temperatures, array(:,2) contains heat capacities
            do itemp=1,cs_ntemp
               read(ifileno,*) cs_array(itemp,1), cl_array(itemp,2) , cs_array(itemp,2), ce_array(itemp,2)
               cl_array(itemp,1)=cs_array(itemp,1)
               ce_array(itemp,1)=ce_array(itemp,1)
            end do
            close(ifileno)
            do_cs_temp='Y'
            do_cl_temp='Y'
            do_ce_temp='Y'
            print *,' Total CV file read'
      else
         if (do_cs_temp=='Y') then
            open(ifileno,file=csfile)

            ! Find the number of spin temperatures 
            cs_ntemp=0
            do 
               read(ifileno,*,end=30) tdum, cdum
               cs_ntemp=cs_ntemp+1
            end do
            30         continue
            rewind(ifileno)

            ! Allocate array for spin Cv
            allocate(cs_array(cs_ntemp,2),stat=i_stat)
            call memocc(i_stat,product(shape(cs_array))*kind(cs_array),'cs_array','init_3tm_cv')

            ! Read all input spin temperatures and heat capacities
            do itemp=1,cs_ntemp
               read(ifileno,*) cs_array(itemp,1), cs_array(itemp,2)
            end do
            close(ifileno)
            print *,' Spin CV file read'
         end if

         if (do_cl_temp=='Y') then
            open(ifileno,file=clfile)
            ! Find the number of lattice temperatures 
            cl_ntemp=0
            do 
               read(ifileno,*,end=40) tdum, cdum
               cl_ntemp=cl_ntemp+1
            end do
            40         continue

            ! Allocate array for lattice Cv
            allocate(cl_array(cl_ntemp,2),stat=i_stat)
            call memocc(i_stat,product(shape(cl_array))*kind(cl_array),'cl_array','init_3tm_cv')

            ! Read all input spin temperatures and heat capacities
            do itemp=1,cs_ntemp
               read(ifileno,*) cl_array(itemp,1), cl_array(itemp,2)
            end do
            print *,' Lattice CV file read'
         end if
      end if

   else

      if (do_cs_temp=='Y') then
         ! Deallocate arrays
         i_all=-product(shape(cs_array))*kind(cs_array)
         deallocate(cs_array,stat=i_stat)
         call memocc(i_stat,i_all,'cs_array','init_3tm_cv')
      end if

      if (do_cl_temp=='Y') then
         i_all=-product(shape(cl_array))*kind(cl_array)
         deallocate(cl_array,stat=i_stat)
         call memocc(i_stat,i_all,'cl_array','init_3tm_cv')
      end if
   end if

end subroutine init_3tm_cv


subroutine unify_3tm_params(C1,C2,C3,alat,NA)
   !
   use Constants 
   !
   implicit none
   !
   real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
   real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
   real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
   real(dblprec), intent(in) :: alat   !< Lattice parameter
   integer, intent(in) :: NA  !< Number of atoms in unit cell

   real(dblprec) :: cell_vol, density, cell_area, p_dens


   !!!cell_vol=(C1(1)*C2(2)*C3(3)-C1(1)*C2(3)*C3(2)+&
   !!!   C1(2)*C2(3)*C3(1)-C1(2)*C2(1)*C3(3)+&
   !!!   C1(3)*C2(1)*C3(2)-C1(3)*C2(2)*C3(1))*alat**3

   !!!cell_area = (C1(2)*C2(3)-C1(3)*C2(2))**2 &
   !!!          + (C1(3)*C2(1)-C1(1)*C2(3))**2 &
   !!!          + (C1(1)*C2(2)-C1(2)*C2(1))**2 
   !!!cell_area = sqrt(cell_area)*alat**2

   !!!density = NA/cell_vol
   density = 4.0d0/alat**3 ! hard check for fcc
   p_dens = density * 22.0d-9   ! hard fix for 22nm thickness

   if(allocated(ce_array)) ce_array(:,2) = ce_array(:,2) * k_bolt
   if(allocated(cl_array)) cl_array(:,2) = cl_array(:,2) * k_bolt
   if(allocated(cs_array)) cs_array(:,2) = cs_array(:,2) * k_bolt

   !print *,'Rescaling 3TM couplings', density, p_dens

   ! Pulse intensity given per area
   P_pulse = P_pulse  / p_dens / (sigma_pulse*sqrt(2.0d0*pi))
   ! Remaing parameters given per volume
   G_cool = G_cool / density
   C_latt = C_latt / density
   C_spin = C_spin / density
   C_elec = C_elec / density
   gamma_C_elec = gamma_C_elec / density
   Gel = Gel / density
   Ges = Ges / density
   Gsl = Gsl / density

!!!! P_pulse     70            #Power/heat of the pulse  (J/m^2)
!!!! t0_pulse    1.0e-12       #Center of the pulse  (s)
!!!! sigma_pulse 0.02e-12      #Width of the pulse (s)
!!!! G_cool      1.0e1       #Cooling rate (only on electron bath) (J/s)
!!!! Temp_final  300           #End temperature  (K)
!!!! gamma_Ce    6.0e3         #Heat capacity prefactor for electrons  (J/m^3/K^2)(because C_e = gamma * Te)
!!!! Cl          2.2e6         #Heat capacity for lattice  (J/m^3/K) (constant)
!!!! Cs          0.7e6         #Heat capacity for spins  (J/m^3/K) (constant)
!!!! Gel         8.0e17        #Electron-lattice transfer rate (W/m^3/K)
!!!! Ges         6.0e17        #Electron-spin transfer rate (W/m^3/K)
!!!! Gsl         0.3e17        #Spin-lattice transfer rate (W/m^3/K)



   return

end subroutine unify_3tm_params

end module temperature_3tm
