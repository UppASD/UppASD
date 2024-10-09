module chroma

   use Parameters
   use Profiling
   use Constants, only : pi
   use RandomNumbers

   !real(dblprec) ::
   !
   character :: do_quant_colour = 'N'

   integer   :: nterms

   character(len=35) :: qcfile  !< Specific heat for spins
   character(len=35) :: qccomment 

   ! Lorentzian factors
   real(dblprec), dimension(:), allocatable   :: c_i
   real(dblprec), dimension(:), allocatable   :: gamma_i
   real(dblprec), dimension(:), allocatable   :: omega_i

   ! Data arrays
   real(dblprec), dimension(:,:), allocatable   :: xsi 
   real(dblprec), dimension(:,:), allocatable   :: xsi_dot
   real(dblprec), dimension(:,:), allocatable   :: xsi_0
   real(dblprec), dimension(:,:), allocatable   :: xsi_dot0
   real(dblprec), dimension(:,:), allocatable   :: eta
   real(dblprec), dimension(:,:), allocatable   :: eta_0

   ! Work arrays
   real(dblprec), dimension(:), allocatable :: xsi_pred
   real(dblprec), dimension(:), allocatable :: xsidot_pred
   real(dblprec), dimension(:), allocatable :: xsidot_fp
   real(dblprec), dimension(:), allocatable :: xsidot_fc

   ! Input parameters
   real(dblprec) :: qc_gamma  = 0.5_dblprec
   real(dblprec) :: qc_omega  = 50.0_dblprec
   real(dblprec) :: qc_fac    = 1.0_dblprec
   real(dblprec) :: qc_dt     = 0.0001_dblprec
   integer :: qc_nruns = 1000000

   private

   public :: do_quant_colour
   public :: read_parameters_qc
   public :: initialize_qc
   public :: qc_generator, qc_single_generator

   !! subroutine set_qc_defaults()
   !! end subroutine set_temperature_3tm_defaults
   !! subroutine allocate_qc_data(Natom, Mensemble, flag)
   !! end subroutine allocate_qc_data
   !! subroutine read_parameters_qc(ifile)


contains 

   subroutine qc_single_generator(rand_vec,ndim,dt)

      use Constants, only : pi, gama, hbar, k_bolt, mub
      use Temperature, only : temp_array
      use Damping, only : lambda1_array
      use InputData, only : delta_t
      use Momentdata, only : mmomi

      implicit none

      integer, intent(in) :: ndim
      real(dblprec), dimension(ndim), intent(inout) :: rand_vec
      real(dblprec), intent(in) :: dt
      
      real(dblprec) :: dtau

      integer :: iterm

      ! Parameters in THz
      dtau = dt / 1.0e-12_dblprec
      !print *,'D',dt,dtau
      call f_rk4_deg2_omp(ndim,dtau,xsi,xsi_dot,xsi_0,xsi_dot0)
      !call f_heun_deg2_omp(ndim,dtau,xsi,xsi_dot,xsi_0,xsi_dot0)


      rand_vec = 0.0d0
      do iterm = 1,nterms
         rand_vec = rand_vec +  xsi(:,iterm)*c_i(iterm)
      end do
      !rand_vec = rand_vec *k_bolt*Temp_array(1)*sqrt(2.0_dblprec*lambda1_array(1)/mmomi(1,1))/mub
      !print *,k_bolt*Temp_array(1)*sqrt(2.0_dblprec*lambda1_array(1)/mmomi(1,1))/mub

   end subroutine qc_single_generator
      
   subroutine qc_generator(rand_vec,ndim,dt)

      use Constants, only : pi, gama, hbar, k_bolt, mub
      use Temperature, only : temp_array
      use Damping, only : lambda1_array
      use InputData, only : delta_t
      use Momentdata, only : mmomi

      implicit none

      integer, intent(in) :: ndim
      real(dblprec), dimension(ndim), intent(inout) :: rand_vec
      real(dblprec), intent(in) :: dt
      
      real(dblprec) :: dtau

      integer :: iterm

      ! Dimensionless time
      dtau = dt * Temp_array(1) * k_bolt / hbar
      !print *,'D',dt,dtau
      call f_rk4_deg2_omp(ndim,dtau,xsi,xsi_dot,xsi_0,xsi_dot0)
      !call f_heun_deg2_omp(ndim,dtau,xsi,xsi_dot,xsi_0,xsi_dot0)


      rand_vec = 0.0d0
      do iterm = 1,nterms
         rand_vec = rand_vec +  xsi(:,iterm)*c_i(iterm) !*2.00_dblprec
      end do
      !rand_vec = rand_vec *k_bolt*Temp_array(1)*sqrt(2.0_dblprec*lambda1_array(1)/mmomi(1,1))/mub
      !print *,k_bolt*Temp_array(1)*sqrt(2.0_dblprec*lambda1_array(1)/mmomi(1,1))/mub

   end subroutine qc_generator
      
   subroutine qc_warmup(Natom,Mensemble,nruns,dt,rand_vec)
      use Constants, only : pi, gama, hbar, k_bolt, mub
      use Temperature, only : temp_array
      use Damping, only : lambda1_array

      implicit none

      integer, intent(in) :: Natom
      integer, intent(in) :: Mensemble
      integer, intent(in) :: nruns
      real(dblprec), intent(in) :: dt
      real(dblprec), dimension(3*Natom*Mensemble), intent(inout) :: rand_vec

      real(dblprec) :: dtau
      integer :: irun, iterm
      character(len=30) :: filn

      xsi = 0.0d0
      xsi_dot = 0.0d0
      xsi_0 = 0.0d0
      xsi_dot0 = 0.0d0
      eta_0 = 0.0d0


      filn='colournoise.out'
      open(ofileno,file=filn, position='append')


      dtau = dt * Temp_array(1) * k_bolt / hbar
      do irun=1,nruns/2
         !call f_heun_deg2_omp(3*Natom*Mensemble,dtau,xsi,xsi_dot,xsi_0,xsi_dot0)
         call f_rk4_deg2_omp(3*Natom*Mensemble,dtau,xsi,xsi_dot,xsi_0,xsi_dot0)
         rand_vec = 0.0d0
         do iterm = 1,nterms
            rand_vec = rand_vec +  xsi(:,iterm)*c_i(iterm)
         end do
      end do
      do irun=irun,nruns
         !call f_heun_deg2_omp(3*Natom*Mensemble,dtau,xsi,xsi_dot,xsi_0,xsi_dot0)
         call f_rk4_deg2_omp(3*Natom*Mensemble,dtau,xsi,xsi_dot,xsi_0,xsi_dot0)
         rand_vec = 0.0d0
         do iterm = 1,nterms
            rand_vec = rand_vec +  xsi(:,iterm)*c_i(iterm)
         end do
         !write(ofileno,'(3f12.6)') rand_vec(:,1)
         write(ofileno,'(9g20.10)') rand_vec(1),eta(1,1:2),(c_i(iterm),xsi(1,iterm),iterm=1,2)
      end do

      close(ofileno)

   end subroutine qc_warmup
      
   !-----------------------------------------------------------------------------
   ! SUBROUTINE: f_heun_deg2_omp
   !> @brief Heun solver for coloured noise generation
   !
   !> @author Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine f_heun_deg2_omp(ndim,dt,xsi_fun,xsidot_fun,xsi_old,xsidot_old)

      implicit none

      integer, intent(in) :: ndim
      real(dblprec), intent(in) :: dt

      real(dblprec), dimension(ndim,nterms), intent(inout) :: xsi_fun
      real(dblprec), dimension(ndim,nterms), intent(inout) :: xsidot_fun
      real(dblprec), dimension(ndim,nterms), intent(inout) :: xsi_old
      real(dblprec), dimension(ndim,nterms), intent(inout) :: xsidot_old

      integer :: iterm, i
      real(dblprec) :: dthalf, sqgamma, omega2

      !outrand = 0.0d0
      dthalf = dt/2.0_dblprec
      do iterm = 1,nterms


         call rng_gaussian(eta,nterms*ndim,1.0_dblprec)
         !eta(:,iterm) = sqrt(2.0*gamma_i(iterm)) * eta((:,iterm)
         sqgamma = sqrt(2.0_dblprec*gamma_i(iterm))
         !print *,'dt',dt
         !sqgamma = sqrt(2.0_dblprec*gamma_i(iterm)/dt)
         omega2 = omega_i(iterm)**2

         !$omp parallel do default(shared) private(i) schedule(static)
         do i = 1,ndim
            xsi_pred(i) = xsi_old(i,iterm) + xsidot_old(i,iterm) * dt

            xsidot_fp(i) = eta_0(i,iterm)*sqgamma - omega2*xsi_old(i,iterm) - gamma_i(iterm)*xsidot_old(i,iterm)
            xsidot_pred(i) = xsidot_old(i,iterm) + dt*xsidot_fp(i)

            xsidot_fc(i) = eta(i,iterm)*sqgamma - omega2*xsi_pred(i) - gamma_i(iterm)*xsidot_pred(i)

            xsi_fun(i,iterm) = xsi_old(i,iterm) + dthalf*xsidot_old(i,iterm) + dthalf*xsidot_pred(i)

            xsidot_fun(i,iterm) = xsidot_old(i,iterm) + xsidot_fp(i)*dthalf + xsidot_fc(i)*dthalf

            xsi_old(i,iterm) = xsi_fun(i,iterm)
            xsidot_old(i,iterm) = xsidot_fun(i,iterm)

           ! outrand(i) = outrand(i) + c_i(iterm)*xsi_fun(i,iterm)

            eta_0(i,iterm) = eta(i,iterm)
         end do

      end do
      !write(ofileno,'(9g20.10)') outrand(1),eta(1,:),(c_i(iterm)*xsi_fun(i,iterm),iterm=1,nterms)

   end subroutine f_heun_deg2_omp

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: f_heun_deg2_omp
   !> @brief Heun solver for coloured noise generation
   !
   !> @author Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine f_rk4_deg2_omp(ndim,dt,xsi_fun,xsidot_fun,xsi_old,xsidot_old)

      implicit none

      integer, intent(in) :: ndim
      real(dblprec), intent(in) :: dt

      real(dblprec), dimension(ndim,nterms), intent(inout) :: xsi_fun
      real(dblprec), dimension(ndim,nterms), intent(inout) :: xsidot_fun
      real(dblprec), dimension(ndim,nterms), intent(inout) :: xsi_old
      real(dblprec), dimension(ndim,nterms), intent(inout) :: xsidot_old

      integer :: iterm, i
      real(dblprec) :: dthalf, sqgamma, omega2

      real(dblprec) :: k1, k1_dot, k2, k2_dot, k3, k3_dot, k4, k4_dot, eta_half

      !outrand = 0.0d0
      dthalf = dt/2.0_dblprec
      call rng_gaussian(eta,nterms*ndim,1.0_dblprec)
      do iterm = 1,nterms


         !call rng_gaussian(eta,ndim,1.0_dblprec)
         !eta(:,iterm) = sqrt(2.0*gamma_i(iterm)) * eta((:,iterm)
         sqgamma = sqrt(2.0_dblprec*gamma_i(iterm))
         !print *,'dt',dt
         !sqgamma = sqrt(2.0_dblprec*gamma_i(iterm)/dt)
         omega2 = omega_i(iterm)**2

         !$omp parallel do default(shared) private(i,eta_half,k1,k1_dot,k2,k2_dot,k3,k3_dot,k4,k4_dot) schedule(static)
         do i = 1,ndim
            !print *,'--------------------------------'
            !print '(a,6g20.10)','Xsi, Xsi_dot ',xsi_old(i,iterm), xsidot_old(i,iterm)

            eta_half = 0.5_dblprec*eta_0(i,iterm) + 0.5_dblprec*eta(i,iterm)
            !print '(a,6g20.10)','Sigma,Gamma ',sqgamma,omega2,gamma_i(iterm),eta_0(i,iterm),eta_half,eta(i,iterm)

            k1     =  xsidot_old(i,iterm)
            k1_dot =  eta_0(i,iterm)*sqgamma - omega2*xsi_old(i,iterm) - gamma_i(iterm)*xsidot_old(i,iterm)
            !k1_dot =  eta_0(i,iterm)*sqgamma - omega2*xsi_old(i,iterm) - gamma_i(iterm)*xsidot_old(i,iterm)

            k2     =  xsidot_old(i,iterm) + dthalf * k1_dot
            k2_dot =  eta_0(i,iterm)*sqgamma - omega2*(xsi_old(i,iterm)+dthalf*k1) - gamma_i(iterm)*(xsidot_old(i,iterm)+dthalf*k1_dot)
            !k2_dot =  eta_half*sqgamma - omega2*(xsi_old(i,iterm)+dthalf*k1) - gamma_i(iterm)*(xsidot_old(i,iterm)+dthalf*k1_dot)

            k3     =  xsidot_old(i,iterm) + dthalf * k2_dot
            k3_dot =  eta_0(i,iterm)*sqgamma - omega2*(xsi_old(i,iterm)+dthalf*k2) - gamma_i(iterm)*(xsidot_old(i,iterm)+dthalf*k2_dot)
            !k3_dot =  eta_half*sqgamma - omega2*(xsi_old(i,iterm)+dthalf*k2) - gamma_i(iterm)*(xsidot_old(i,iterm)+dthalf*k2_dot)

            k4     =  xsidot_old(i,iterm) + dt * k3_dot
            k4_dot =  eta_0(i,iterm)*sqgamma - omega2*(xsi_old(i,iterm)+dt*k3) - gamma_i(iterm)*(xsidot_old(i,iterm)+dt*k3_dot)
            !k4_dot =  eta(i,iterm)*sqgamma - omega2*(xsi_old(i,iterm)+dt*k3) - gamma_i(iterm)*(xsidot_old(i,iterm)+dt*k3_dot)

            xsi_fun(i,iterm) = xsi_old(i,iterm) + dt/6.0_dblprec*(k1+2.0_dblprec*k2+2.0_dblprec*k3+k4)
            xsidot_fun(i,iterm) = xsidot_old(i,iterm) + dt/6.0_dblprec*(k1_dot+2.0_dblprec*k2_dot+2.0_dblprec*k3_dot+k4_dot)
            !print '(2i)',i,iterm
            !print '(a,6g20.10)','k     ',k1,k2,k3,k4,xsi_fun(i,iterm),xsi_old(i,iterm)
            !print '(a,6g20.10)','k_dot ',k1_dot,k2_dot,k3_dot,k4_dot,xsidot_fun(i,iterm),xsidot_old(i,iterm)

            xsi_old(i,iterm) = xsi_fun(i,iterm)
            xsidot_old(i,iterm) = xsidot_fun(i,iterm)
            !print '(a,6g20.10)','xsi, xsi_dot ',xsi_old(i,iterm), xsidot_old(i,iterm)

           ! outrand(i) = outrand(i) + c_i(iterm)*xsi_fun(i,iterm)

            eta_0(i,iterm) = eta(i,iterm)
         end do
         !$omp end parallel do

      end do
      !print *,'--------------------------',nterms, ndim
      !print *,gamma_i
      !print *,omega_i
      !print *,xsi_old(5,1), xsi_old(5,2)
      !print *,xsidot_old(5,1), xsidot_old(5,2)
      !write(ofileno,'(9g20.10)') outrand(1),eta(1,:),(c_i(iterm)*xsi_fun(i,iterm),iterm=1,nterms)

   end subroutine f_rk4_deg2_omp

   subroutine initialize_qc(Natom,Mensemble,flag)

      implicit none

      integer, intent(in) :: Natom
      integer, intent(in) :: Mensemble
      integer, intent(in) :: flag

      real(dblprec), dimension(3*Natom*Mensemble) :: rand_vec


      if(flag>=0) then  
         if (do_quant_colour=='Y') then
            call set_qc_defaults()
         else if (do_quant_colour=='S') then
            call set_qc_pulse()
         else if (do_quant_colour=='F') then
            call read_qcfile()
         end if
      end if

      !print *, 'Allocating QC data', nterms
      call allocate_qc_data(Natom, Mensemble, flag)

      write(*,'(1x,A)') '--------------------------------------'
      write(*,'(1x,*(A))') ' Initializing ', &
             char(27)//'[31m'//'c'//char(27)//'[32m'//'o'//char(27)//'[33m'//'l'//&
             char(27)//'[34m'//'o'//char(27)//'[35m'//'u'//char(27)//'[36m'//'r'//&
             char(27)//'[37m'//'e'//char(27)//'[31m'//'d'//char(27)//'[0m', &
             ' noise'
      write(*,'(1x,A,*(F10.4))') '  C_i:     ', c_i
      write(*,'(1x,A,*(F10.4))') '  Gamma_i: ', gamma_i
      write(*,'(1x,A,*(F10.4))') '  Omega_i: ', omega_i
      write(*,'(1x,A,I10)') '  Nruns  : ', qc_nruns
      write(*,'(1x,A,4x,E10.4)') '  qc_dt  : ', qc_dt
      write(*,'(1x,A)') '--------------------------------------'

      write(*,'(1x,A,I0,A)', advance='no') 'QC initialized, now warmup for ', qc_nruns, ' steps'
      call qc_warmup(Natom,Mensemble,qc_nruns,qc_dt,rand_vec)
      write(*,'(A)') ' done'
      !call qc_warmup(Natom,Mensemble,qc_nruns,qc_dt,rand_vec)

   end subroutine initialize_qc

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: read_qcfile
   !> @brief Read parameters for quantum coloured noise
   !
   !> @author Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine read_qcfile()

      use Profiling 

      implicit none

      integer :: iterm, i_stat

      print *,'Opening qcfile: ',qcfile,'.'
      open(ifileno, file=adjustl(qcfile))

      read(ifileno,'(a)') qccomment
      read(ifileno,*) nterms

      allocate(c_i(nterms),stat=i_stat)
      call memocc(i_stat,product(shape(c_i))*kind(c_i),'c_i','read_qcfile')
      allocate(omega_i(nterms),stat=i_stat)
      call memocc(i_stat,product(shape(omega_i))*kind(omega_i),'omega_i','read_qcfile')
      allocate(gamma_i(nterms),stat=i_stat)
      call memocc(i_stat,product(shape(gamma_i))*kind(gamma_i),'gamma_i','read_qcfile')

      do iterm = 1, nterms
         read(ifileno,*) c_i(iterm), omega_i(iterm), gamma_i(iterm)
         print '(i4,3f12.6)',iterm,c_i(iterm), omega_i(iterm), gamma_i(iterm)
      end do

   end subroutine read_qcfile

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: set_qc_defaults
   !> @brief Initialize default values for quantum noise 
   !
   !> @author Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine set_qc_defaults()

      use Profiling 


      implicit none

      integer :: i_stat

      nterms = 2

      allocate(c_i(nterms),stat=i_stat)
      call memocc(i_stat,product(shape(c_i))*kind(c_i),'c_i','set_qc_defaults')
      allocate(omega_i(nterms),stat=i_stat)
      call memocc(i_stat,product(shape(omega_i))*kind(omega_i),'omega_i','set_qc_defaults')
      allocate(gamma_i(nterms),stat=i_stat)
      call memocc(i_stat,product(shape(gamma_i))*kind(gamma_i),'gamma_i','set_qc_defaults')

      c_i(1) = 1.8315d0
      c_i(2) = 0.3429d0
      omega_i(1) = 2.7189d0
      omega_i(2) = 1.2223d0
      gamma_i(1) = 5.0142d0
      gamma_i(2) = 3.2974d0


   end subroutine set_qc_defaults

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: set_qc_pulse
   !> @brief Initialize default values for pumped coloured noise 
   !
   !> @author Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine set_qc_pulse()

      use Profiling 


      implicit none

      integer :: i_stat

      nterms = 1

      allocate(c_i(nterms),stat=i_stat)
      call memocc(i_stat,product(shape(c_i))*kind(c_i),'c_i','set_qc_defaults')
      allocate(omega_i(nterms),stat=i_stat)
      call memocc(i_stat,product(shape(omega_i))*kind(omega_i),'omega_i','set_qc_defaults')
      allocate(gamma_i(nterms),stat=i_stat)
      call memocc(i_stat,product(shape(gamma_i))*kind(gamma_i),'gamma_i','set_qc_defaults')

      omega_i(1) = qc_omega
      gamma_i(1) = qc_gamma
      c_i(1) = sqrt(qc_fac*qc_gamma**2/pi)

   end subroutine set_qc_pulse

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_qc_data
   !> @brief Subroutine to allocate the arrays needed for quantum noise
   !> @author Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine allocate_qc_data(Natom, Mensemble, flag)

      use Profiling 

      implicit none

      integer, intent(in) :: Natom
      integer, intent(in) :: Mensemble
      integer, intent(in) :: flag

      integer :: i_stat,i_all

      if(flag>=0) then
         allocate(xsi(3*Natom*Mensemble,nterms),stat=i_stat)
         call memocc(i_stat,product(shape(xsi))*kind(xsi),'xsi','allocate_qc_data')
         allocate(xsi_dot(3*Natom*Mensemble,nterms),stat=i_stat)
         call memocc(i_stat,product(shape(xsi_dot))*kind(xsi_dot),'xsi_dot','allocate_qc_data')
         allocate(xsi_dot0(3*Natom*Mensemble,nterms),stat=i_stat)
         call memocc(i_stat,product(shape(xsi_dot0))*kind(xsi_dot0),'xsi_dot0','allocate_qc_data')
         allocate(xsi_0(3*Natom*Mensemble,nterms),stat=i_stat)
         call memocc(i_stat,product(shape(xsi_0))*kind(xsi_0),'xsi_0','allocate_qc_data')
         allocate(eta(3*Natom*Mensemble,nterms),stat=i_stat)
         call memocc(i_stat,product(shape(eta))*kind(eta),'eta','allocate_qc_data')
         allocate(eta_0(3*Natom*Mensemble,nterms),stat=i_stat)
         call memocc(i_stat,product(shape(eta_0))*kind(eta_0),'eta_0','allocate_qc_data')
         !
         allocate(xsi_pred(3*Natom*Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(xsi_pred))*kind(xsi_pred),'xsi_pred','allocate_qc_data')
         allocate(xsidot_pred(3*Natom*Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(xsidot_pred))*kind(xsidot_pred),'xsidot_pred','allocate_qc_data')
         allocate(xsidot_fp(3*Natom*Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(xsidot_fp))*kind(xsidot_fp),'xsidot_fp','allocate_qc_data')
         allocate(xsidot_fc(3*Natom*Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(xsidot_fc))*kind(xsidot_fc),'xsidot_fc','allocate_qc_data')
      else
         i_all=-product(shape(xsi))*kind(xsi)
         deallocate(xsi,stat=i_stat)
         call memocc(i_stat,i_all,'xsi','allocate_qc_data')
         i_all=-product(shape(xsi_dot))*kind(xsi_dot)
         deallocate(xsi_dot,stat=i_stat)
         call memocc(i_stat,i_all,'xsi_dot','allocate_qc_data')
         i_all=-product(shape(xsi_dot0))*kind(xsi_dot0)
         deallocate(xsi_dot0,stat=i_stat)
         call memocc(i_stat,i_all,'xsi_dot0','allocate_qc_data')
         i_all=-product(shape(xsi_0))*kind(xsi_0)
         deallocate(xsi_0,stat=i_stat)
         call memocc(i_stat,i_all,'xsi_0','allocate_qc_data')

         i_all=-product(shape(eta))*kind(eta)
         deallocate(eta,stat=i_stat)
         call memocc(i_stat,i_all,'eta','allocate_qc_data')
         i_all=-product(shape(eta_0))*kind(eta_0)
         deallocate(eta_0,stat=i_stat)
         call memocc(i_stat,i_all,'eta_0','allocate_qc_data')

         i_all=-product(shape(xsi_pred))*kind(xsi_pred)
         deallocate(xsi_pred,stat=i_stat)
         call memocc(i_stat,i_all,'xsi_pred','allocate_qc_data')
         i_all=-product(shape(xsidot_pred))*kind(xsidot_pred)
         deallocate(xsidot_pred,stat=i_stat)
         call memocc(i_stat,i_all,'xsidot_pred','allocate_qc_data')
         i_all=-product(shape(xsidot_fp))*kind(xsidot_fp)
         deallocate(xsidot_fp,stat=i_stat)
         call memocc(i_stat,i_all,'xsidot_fp','allocate_qc_data')
         i_all=-product(shape(xsidot_fc))*kind(xsidot_fc)
         deallocate(xsidot_fc,stat=i_stat)
         call memocc(i_stat,i_all,'xsidot_fc','allocate_qc_data')


         i_all=-product(shape(c_i))*kind(c_i)
         deallocate(c_i,stat=i_stat)
         call memocc(i_stat,i_all,'c_i','allocate_qc_data')
         i_all=-product(shape(omega_i))*kind(omega_i)
         deallocate(omega_i,stat=i_stat)
         call memocc(i_stat,i_all,'omega_i','allocate_qc_data')
         i_all=-product(shape(gamma_i))*kind(gamma_i)
         deallocate(gamma_i,stat=i_stat)
         call memocc(i_stat,i_all,'gamma_i','allocate_qc_data')
      end if

   end subroutine allocate_qc_data

   !---------------------------------------------------------------------------
   ! SUBROUTINE: read_parameters_qc
   !> @brief
   !> Read input parameters for quantum colored noise
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_parameters_qc(ifile)

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

            case('do_quant_colour')
               read(ifile,*,iostat=i_err) do_quant_colour
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('qcfile')
               read(ifile,*,iostat=i_err) qcfile
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('qc_omega')
               read(ifile,*,iostat=i_err) qc_omega
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('qc_gamma')
               read(ifile,*,iostat=i_err) qc_gamma
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('qc_fac')
               read(ifile,*,iostat=i_err) qc_fac
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('qc_dt')
               read(ifile,*,iostat=i_err) qc_dt
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('qc_nruns')
               read(ifile,*,iostat=i_err) qc_nruns
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
   end subroutine read_parameters_qc

end module chroma

!!!    !-----------------------------------------------------------------------------
!!!    ! SUBROUTINE: f_euler_deg2
!!!    !> @brief Generic forward Euler solver for ODEs of second order
!!!    !> @note g_fun = f_fun' and g_old/f_old are data for previous iteration
!!!    !
!!!    !> @author Anders Bergman
!!!    !-----------------------------------------------------------------------------
!!!    subroutine f_euler_deg2(ndim,dt,f_fun,g_fun,f_old,g_old)
!!! 
!!!       implicit none
!!! 
!!!       integer, intent(in) :: ndim
!!!       real(dblprec), intent(in) :: dt
!!! 
!!!       real(dblprec), dimension(ndim,nterms), intent(inout) :: f_fun
!!!       real(dblprec), dimension(ndim,nterms), intent(inout) :: g_fun
!!!       real(dblprec), dimension(ndim,nterms), intent(inout) :: f_old
!!!       real(dblprec), dimension(ndim,nterms), intent(inout) :: g_old
!!! 
!!!       !real(dblprec), dimension(ndim) :: outrand
!!! 
!!!       integer :: iterm
!!! 
!!!       !print *,'Eta check: ',allocated(eta)
!!!       do iterm = 1, nterms
!!!          call rng_gaussian(eta(:,iterm),ndim,2.0_dblprec*gamma_i(iterm))
!!!          !call fill_rngarray(eta,ndim)
!!! 
!!!          g_fun(:,iterm) = g_old(:,iterm) &
!!!             + eta(:,iterm)*dt - omega_i(iterm)**2*f_old(:,iterm)*dt - gamma_i(iterm)*g_old(:,iterm)*dt
!!!          f_fun(:,iterm) = f_old(:,iterm) + g_fun(:,iterm) * dt
!!! 
!!!          f_old = f_fun
!!!          g_old = g_fun
!!! 
!!!          !outrand = outrand + c_i(iterm)*f_fun(:,iterm)
!!!       end do
!!! 
!!!       !write(ofileno,'(9g20.10)') outrand(1),eta(1,:),(c_i(iterm)*f_fun(:,iterm),iterm=1,nterms)
!!!       !write(ofileno,'(5g20.10)') c_i(1)*f_fun(1,1),eta(1)
!!! 
!!!    end subroutine f_euler_deg2
!!! 
!!!    !-----------------------------------------------------------------------------
!!!    ! SUBROUTINE: f_heun_deg2
!!!    !> @brief Generic Heun solver for ODEs of second order
!!!    !> @note g_fun = f_fun' and g_old/f_old are data for previous iteration
!!!    !
!!!    !> @author Anders Bergman
!!!    !-----------------------------------------------------------------------------
!!!    subroutine f_heun_deg2(ndim,dt,xsi_fun,xsidot_fun,xsi_old,xsidot_old)
!!! 
!!!       implicit none
!!! 
!!!       integer, intent(in) :: ndim
!!!       real(dblprec), intent(in) :: dt
!!! 
!!!       real(dblprec), dimension(ndim,nterms), intent(inout) :: xsi_fun
!!!       real(dblprec), dimension(ndim,nterms), intent(inout) :: xsidot_fun
!!!       real(dblprec), dimension(ndim,nterms), intent(inout) :: xsi_old
!!!       real(dblprec), dimension(ndim,nterms), intent(inout) :: xsidot_old
!!! 
!!!       !real(dblprec), dimension(ndim) :: xsi_pred
!!!       !real(dblprec), dimension(ndim) :: xsidot_pred
!!!       !real(dblprec), dimension(ndim) :: xsidot_fp
!!!       !real(dblprec), dimension(ndim) :: xsidot_fc
!!! 
!!!       !real(dblprec), dimension(ndim) :: outrand
!!! 
!!!       integer :: iterm
!!! 
!!!       !outrand = 0.0d0
!!!       do iterm = 1,nterms
!!! 
!!! 
!!!          call rng_gaussian(eta,ndim,1.0_dblprec)
!!!          !eta(:,iterm) = sqrt(2.0*gamma_i(iterm)) * eta((:,iterm)
!!! 
!!!          xsi_pred = xsi_old(:,iterm) + xsidot_old(:,iterm) * dt
!!! 
!!!          xsidot_fp = eta_0(:,iterm)*sqrt(2.0*gamma_i(iterm)) - omega_i(iterm)**2*xsi_old(:,iterm) - gamma_i(iterm)*xsidot_old(:,iterm)
!!!          xsidot_pred = xsidot_old(:,iterm) + dt*xsidot_fp
!!! 
!!!          !call rng_gaussian(eta(:,iterm),ndim,gamma_i(iterm))
!!! 
!!!          xsidot_fc = eta(:,iterm)*sqrt(2.0*gamma_i(iterm)) - omega_i(iterm)**2*xsi_pred - gamma_i(iterm)*xsidot_pred
!!! 
!!!          xsi_fun(:,iterm) = xsi_old(:,iterm) + 0.5d0*dt*xsidot_old(:,iterm) + 0.5d0*dt*xsidot_pred
!!! 
!!!          xsidot_fun(:,iterm) = xsidot_old(:,iterm) + 0.5d0*xsidot_fp*dt + 0.5d0*xsidot_fc*dt
!!! 
!!!          xsi_old(:,iterm) = xsi_fun(:,iterm)
!!!          xsidot_old(:,iterm) = xsidot_fun(:,iterm)
!!! 
!!!          !outrand = outrand + c_i(iterm)*xsi_fun(:,iterm)
!!! 
!!!          eta_0(:,iterm) = eta(:,iterm)
!!! 
!!!       end do
!!!       !write(ofileno,'(9g20.10)') outrand(1),eta(1,:),(c_i(iterm)*xsi_fun(1,iterm),iterm=1,nterms)
!!! 
!!!    end subroutine f_heun_deg2
!!!    !-----------------------------------------------------------------------------
!!!    ! SUBROUTINE: f_euler_deg2_omp
!!!    !> @brief Euler solver for coloured noise generation
!!!    !
!!!    !> @author Anders Bergman
!!!    !-----------------------------------------------------------------------------
!!!    subroutine f_euler_deg2_omp(ndim,dt,xsi_fun,xsidot_fun,xsi_old,xsidot_old)
!!! 
!!!       implicit none
!!! 
!!!       integer, intent(in) :: ndim
!!!       real(dblprec), intent(in) :: dt
!!! 
!!!       real(dblprec), dimension(ndim,nterms), intent(inout) :: xsi_fun
!!!       real(dblprec), dimension(ndim,nterms), intent(inout) :: xsidot_fun
!!!       real(dblprec), dimension(ndim,nterms), intent(inout) :: xsi_old
!!!       real(dblprec), dimension(ndim,nterms), intent(inout) :: xsidot_old
!!! 
!!!       integer :: iterm, i
!!!       real(dblprec) :: sqgamma, omega2, xsidot_fps
!!! 
!!!       do iterm = 1,nterms
!!! 
!!! 
!!!          call rng_gaussian(eta,ndim,1.0_dblprec)
!!!          sqgamma = sqrt(2.0_dblprec*gamma_i(iterm)/dt)
!!!          omega2 = omega_i(iterm)**2
!!! 
!!!          !$omp parallel do default(shared) private(i,xsidot_fps) schedule(static)
!!!          do i = 1,ndim
!!! 
!!!             xsidot_fps = eta_0(i,iterm)*sqgamma - omega2*xsi_old(i,iterm) - gamma_i(iterm)*xsidot_old(i,iterm)
!!!             xsidot_fun(i,iterm) = xsidot_old(i,iterm) + xsidot_fps*dt + xsidot_fc(i)*dt
!!!             xsi_fun(i,iterm) = xsi_old(i,iterm) + xsidot_fun(i,iterm) * dt
!!! 
!!!             xsi_old(i,iterm) = xsi_fun(i,iterm)
!!!             xsidot_old(i,iterm) = xsidot_fun(i,iterm)
!!! 
!!!             eta_0(i,iterm) = eta(i,iterm)
!!!          end do
!!! 
!!!       end do
!!! 
!!!    end subroutine f_euler_deg2_omp
!!! 
