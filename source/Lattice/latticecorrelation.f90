!> Data and routines for calculate connected ionic displacement correlation function U(r,t) and Fourier transforms U(q,w)
!> @author
!> Johan Hellsvik
module LatticeCorrelation
   use Parameters
   use Profiling
   !
   implicit none
   !

   !Displacement-displacement correlation parameters from input file
   character(len=1) :: do_uc        !< Perform spin-correlation sampling (Y/N/C)
   character(len=1) :: do_ur        !< Perform spin-correlation sampling in real space directly (Y/N)
   !character(len=2) :: do_conv      !< Flag for convolution of Magnon DOS Lorentzian(LY/LG/LW/N) Gaussian(GY/GW/GQ/N)
   !character(len=1) :: qpoints      !< Flag for q-point generation (F=file,A=automatic,C=full cell)
   character(len=1) :: do_uc_proj   !< Measure sublattice projection of U(q,w) (Y/N/C)
   !character(len=1) :: do_uc_local_axis   !< Perform SQW along local quantization axis (Y/N
   character(len=1) :: do_uqt_traj   !< Measure time trajectory of U(q,t) (Y/N)
   !character(len=1) :: do_connected !< Perform the connected part U(q,w)
   !character(LEN=35) :: qfile       !< File name for q-points
   integer :: uc_sep                !< Separation between averages
   integer :: uc_mode               !< Spin correlation mode (0-3)
   integer :: uc_step               !< Separation between sampling steps
   integer :: uc_nstep              !< Number of steps to sample
   !real(dblprec) :: sigma_q         !< Sigma parameter in Q for the Gaussian convolution
   !real(dblprec) :: sigma_w         !< Sigma parameter in W for the Gaussian convolution
   !real(dblprec) :: LWfactor        !< Gamma parameter in W for Lorentzian convolution
   !real(dblprec) :: LQfactor        !< Gamma parameter in Q for Lorentzian convolution

   ! Working variables to perform the printing of the correlation
   integer :: lattnq  !< Number of q-points to sample
   integer :: lattnw  !< Number of lattice frequencies to sample
   integer :: uc_samp_done !< Flag to keep track of if S(q) sampling is done (for do_uc='C')
   integer :: uc_samp_done_proj !< Flag to keep track of if S(q) sampling is done (for do_uc_proj='C')
   real(dblprec) :: qfac
   real(dblprec) :: wfac
   real(dblprec) :: r_max       ! Max length of distance vectors
   real(dblprec) :: nrinv
   real(dblprec) :: na2inv
   real(dblprec) :: nscinv
   real(dblprec) :: nrscinv
   real(dblprec) :: scnstepinv
   real(dblprec), dimension(:), allocatable :: lattw                 !< Lattice Frequencies
   real(dblprec), dimension(:), allocatable :: r_norm                ! Length of distance vectors
   real(dblprec), dimension(:), allocatable :: q_weight              !< Weights of the q points for DOS calculations
   real(dblprec), dimension(:,:), allocatable :: lattq               !< q-points
   real(dblprec), dimension(:,:), allocatable :: corr_k              ! Correlation in q G(k)
   real(dblprec), dimension(:,:), allocatable :: corr_s              ! Correlation in r G(r)
   real(dblprec), dimension(:,:), allocatable :: corr_sr             ! Correlation in r G(r) calculated directly
   real(dblprec), dimension(:,:), allocatable :: corr_kt0            ! Temporary for g(q,t0)
   real(dblprec), dimension(:,:,:), allocatable :: corr_s_proj       ! Correlation in q or r (G(k),G(r)) (sublattice)
   real(dblprec), dimension(:,:,:,:), allocatable :: corr_k_proj     ! Correlation in q G(k) (sublattice)
   real(dblprec), dimension(:,:,:,:), allocatable :: corr_ss_proj    ! Correlation in r G(r) (sublattice)
   complex(dblprec), dimension(:,:), allocatable :: corr_st          ! Temporary for g(r)
   complex(dblprec), dimension(:,:), allocatable :: corr_vst         ! Temporary for g(r)
   complex(dblprec), dimension(:,:,:), allocatable :: sqw            !< U(q,w)
   complex(dblprec), dimension(:,:,:), allocatable :: corr_kt        ! Correlation for G(k,t)
   complex(dblprec), dimension(:,:,:), allocatable :: corr_kw        ! Correlation for G(k,w)
   complex(dblprec), dimension(:,:,:), allocatable :: corr_vkt       ! Correlation for G(k,t)
   complex(dblprec), dimension(:,:,:), allocatable :: corr_vkw       ! Correlation for G(k,w)
   complex(dblprec), dimension(:,:,:), allocatable :: corr_st_proj   ! Temporary for g(r)
   complex(dblprec), dimension(:,:,:,:), allocatable :: corr_kt_proj ! Correlation for G(k,t)
   complex(dblprec), dimension(:,:,:,:), allocatable :: corr_kw_proj ! Correlation for G(k,w)

   real(dblprec), dimension(:,:), allocatable :: phonon_dos !< Phonon density of states
   integer,dimension(2) :: qmin  !< Index of smallest wave vector and gamma point
   real(dblprec), dimension(:), allocatable :: udeltat_corr ! Array storing delta_t for each sample
   real(dblprec), dimension(:), allocatable :: ucstep_arr  ! Array storing uc_step for each sample

   integer :: sc_window_fun  !< Choice of FFT window function (1=box, 2=Hann, 3=Hamming, 4=Blackman-Harris)
   integer :: uc_tidx

   public

   private :: qfac, wfac, r_max, nrinv, na2inv, nscinv, nrscinv, scnstepinv, r_norm, q_weight, corr_k, corr_s, &
      corr_sr, corr_kt0, corr_s_proj, corr_k_proj, corr_ss_proj, corr_st, sqw, corr_kt, &
      corr_kt_proj, corr_kw_proj, phonon_dos, sc_window_fun


contains


   !--------------------------------------------------------------------------
   !
   ! DESCRIPTION
   !> @brief
   !> Driver for correlation function calculations
   !---------------------------------------------------------------------------------
   subroutine lattcorrelation_wrapper(Natom, Mensemble, coord, simid, uvec, vvec, &
         mstep, delta_t, NT, atype, Nchmax, achtype, flag, flag_p )
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: uvec  !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: vvec  !< Current ionic velocity
      integer, intent(in) :: mstep  !< Current simulation step
      real(dblprec), intent(in) :: delta_t               !< Time step
      integer, intent(in) :: NT  !< Number of types of atoms
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      integer, intent(in) :: Nchmax  !< Number of chemical types 
      integer, dimension(Natom), intent(in) :: achtype !< Type of atom
      integer, intent(inout) :: flag  !< Setup, sample, or print 
      integer, intent(inout) :: flag_p  !< Setup, sample, or print 
      !
      if(flag==0) uc_tidx=0

      !Lattice correlation
      ! Sample displacements for correlation functions
      ! Sample S(r) through S(q)
      if(do_uc=='C') then

         if(mod(mstep-1,uc_sep)==0 .or. flag==2) call lattcalc_gk(Natom, Mensemble, coord, simid, &
            uvec, do_uqt_traj, mstep, flag)
         if (do_uc_proj=='C') then
            if(mod(mstep-1,uc_sep)==0) call lattcalc_gk_proj(Natom, Mensemble, coord, &
               simid, uvec, flag_p, NT, atype)
         endif

         ! S(k,w)       
      else if(do_uc=='Q') then

         if(mod(mstep-1,uc_step)==1.and.uc_tidx<=uc_nstep) then
            !!!AB restruc reinstate later
            !!!if(adaptive_time_flag.and.deltat_correction_flag) then
            !!!   adapt_step = mstep+adapt_to_sc_ratio*(sc_step-1) ! Potential adjustment of delta t 
            !!! occuring next time at iteration adapt_ste p
            !!!   deltat_correction_flag = .false.
            !!!end if
            uc_tidx=uc_tidx+1
            udeltat_corr(uc_tidx) = delta_t ! Save current time step
            ucstep_arr(uc_tidx) = uc_step  ! Save current sampling period
            call lattcalc_gkt(Natom, Mensemble, coord, simid, uvec, vvec, flag)
            !if (do_uc_proj=='Q') then
            !   call lattcalc_gkt_proj(Natom, Mensemble, coord, simid, emomM, flag_p, atype, nt)
            !end if
            !if (do_uc_projch=='Q') then
            !   call lattcalc_gkt_projch(Natom, Mensemble, coord, simid, emomM, flag_p, achtype, nchmax)
            !endif

         else if (flag==2) then

            call lattcalc_gkt(Natom, Mensemble, coord, simid, uvec, vvec, flag)
            !if (do_uc_proj=='Q') then
            !   call lattcalc_gkt_proj(Natom, Mensemble, coord, simid, emomM, flag_p, atype, nt)
            !end if
            !if (do_uc_projch=='Q') then
            !   call lattcalc_gkt_projch(Natom, Mensemble, coord, simid, emomM, flag_p, achtype, nchmax)
            !endif

         end if

      end if

      ! Sample U(r) directly for non periodic systems
      !if (do_ur=='Y') then
      !   if (mod(mstep-1,uc_sep)==0) call lattcalc_ur(Natom, Mensemble, coord, simid, uvec, flag)
      !endif

   end subroutine lattcorrelation_wrapper


   !> @brief
   !> Initialization of the variables needed for printing the correlation
   subroutine lattcorrelation_init()

      implicit none

      !Spin-spin correlation
      uc_sep       = 100
      uc_mode      = 2
      uc_step      = 1
      uc_nstep     = 1
      !sigma_q      = 0.0_dblprec
      !sigma_w      = 0.0_dblprec
      !LWfactor     = 0.0_dblprec
      !LQfactor     = 0.0_dblprec
      do_uc        = 'N'
      do_ur        = 'N'
      !qfile        = 'qfile'
      !qpoints      = 'F'
      !do_conv      = 'N'
      !do_uc_proj   = 'N'
      !do_uc_local_axis   = 'N'
      do_uqt_traj   = 'N'
      !do_connected = 'N'

   end subroutine lattcorrelation_init


   !> Calculate U(q) for a cell filling mesh for transform to U(r)
   subroutine lattcalc_gk(Natom, Mensemble, coord, simid, uvec, do_uqt_traj, mstep, cgk_flag)

      use Constants, only : pi

      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: uvec  !< Current ionic displacement
      character(len=1), intent(in) :: do_uqt_traj !< Measure time trajectory of U(q,t) (Y/N)
      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(inout) :: cgk_flag  !< Allocate or deallocate (1/-1)
      !
      integer :: iq,r,l,i_stat,i_all
      character(len=30) :: filn
      real(dblprec) :: epowqr!,s0,sp
      real(dblprec) :: qdr,nainv,k_min,qfac,wfac
      real(dblprec), dimension(3) :: s0,sp,cl

      complex(dblprec), dimension(3) :: cl_step

      if(cgk_flag==0) then
         ! First call, allocate and clear arrays
         allocate(corr_s(3,lattnq),stat=i_stat)
         call memocc(i_stat,product(shape(corr_s))*kind(corr_s),'corr_s','lattcalc_gk')
         allocate(corr_k(3,lattnq),stat=i_stat)
         call memocc(i_stat,product(shape(corr_k))*kind(corr_k),'corr_k','lattcalc_gk')
         allocate(corr_kt0(3,lattnq),stat=i_stat)
         call memocc(i_stat,product(shape(corr_kt0))*kind(corr_kt0),'corr_kt0','lattcalc_gk')
         corr_k=0.0_dblprec
         cgk_flag=1
         uc_samp_done=0
      end if

      qfac=2._dblprec*pi
      wfac=1._dblprec
      nainv=1.0_dblprec/Natom

      if (cgk_flag==1) then
         ! Calculate s(k) for the current iteration and add to average of G(k)
         corr_s=0.0_dblprec
         corr_kt0=0.0_dblprec

         !if (do_connected=='Y') call calc_mavrg_vec(Natom,Mensemble,emomM,mavrg_vec,mavg_axis)
         !SLDTODO insert OMP statements
         !!!$omp parallel do default(shared) private(r,iq,l,qdr,epowqr) schedule(static)
         do iq=1,lattnq
            do l=1,Mensemble
               corr_s(:,iq)=0.0_dblprec
               do r=1,Natom
                  qdr=lattq(1,iq)*coord(1,r)+lattq(2,iq)*coord(2,r)+lattq(3,iq)*coord(3,r)
                  !SLDTODO Consider to use exponential Fourier transform instead
                  epowqr=cos(qfac*qdr)*nainv
                  !write(*,*) 'uvec(1,r,l)', uvec(1,r,l)
                  !write(*,*) 'uvec(2,r,l)', uvec(2,r,l)
                  !write(*,*) 'uvec(3,r,l)', uvec(3,r,l)
                  corr_s(1,iq)=corr_s(1,iq)+epowqr*uvec(1,r,l)!-epowqr*mavrg_vec(1)
                  corr_s(2,iq)=corr_s(2,iq)+epowqr*uvec(2,r,l)!-epowqr*mavrg_vec(2)
                  corr_s(3,iq)=corr_s(3,iq)+epowqr*uvec(3,r,l)!-epowqr*mavrg_vec(3)
               end do
               !write(*,*) 'corr_s(1:3,iq)', iq, corr_s(1:3,iq)
               corr_k(1,iq)=corr_k(1,iq) + corr_s(1,iq)
               corr_k(2,iq)=corr_k(2,iq) + corr_s(2,iq)
               corr_k(3,iq)=corr_k(3,iq) + corr_s(3,iq)
               corr_kt0(1,iq)=corr_kt0(1,iq) + corr_s(1,iq)
               corr_kt0(2,iq)=corr_kt0(2,iq) + corr_s(2,iq)
               corr_kt0(3,iq)=corr_kt0(3,iq) + corr_s(3,iq)
            end do
         end do
         !!!$omp end parallel do

         if(do_uqt_traj == 'Y') then
            ! Write U(q,t0)
            corr_kt0=corr_kt0/Mensemble
            write (filn,'(''uqt0.'',a,''.out'')') trim(simid)
            open(200,file=filn, position='append')
            do iq=1,lattnq
               write(200,'(i10, i10,3f10.4,4es18.8)') mstep, iq,(lattq(l,iq),l=1,3), &
                  !write(200,'(i10, i10,3f10.4,4f18.8)') mstep, iq,(lattq(l,iq),l=1,3), &
               (((corr_kt0(l,iq))),l=1,3), &
                  sqrt(corr_kt0(1,iq)**2+corr_kt0(2,iq)**2+corr_kt0(3,iq)**2)
            end do
            close(200)
         end if

         uc_samp_done=uc_samp_done+1
      end if

      if (cgk_flag==2) then
         i_all=-product(shape(corr_kt0))*kind(corr_kt0)
         deallocate(corr_kt0,stat=i_stat)
         call memocc(i_stat,i_all,'corr_kt0','lattcalc_gk')
         ! Finish sampling and write S(q)
         corr_k=corr_k/(uc_samp_done*Mensemble)

         !Write G(k)
         write (filn,'(''uq.'',a,''.out'')') trim(simid)
         open(200,file=filn,status='replace')
         do iq=1,lattnq
            write(200,'(i10,3f10.4,4es18.8)') iq,(lattq(l,iq),l=1,3),(((corr_k(l,iq))),l=1,3), &
               sqrt(corr_k(1,iq)**2+corr_k(2,iq)**2+corr_k(3,iq)**2)
            !write(200,'(i10,3f10.4,4f18.8)') iq,(lattq(l,iq),l=1,3),(((corr_k(l,iq))),l=1,3), &
            !     sqrt(corr_k(1,iq)**2+corr_k(2,iq)**2+corr_k(3,iq)**2)
         end do
         close(200)

         ! Calculate the correlation length following the Katzgraber recipe
         !SLDTODO For the displacement correlation, no substraction of the mean displacement
         !SLDTODO has been done, i.e. we do not have the connected correlation function.
         !SLDTODO how to define/calculate the displacement correlation?
         s0=corr_k(:,qmin(1))
         sp=corr_k(:,qmin(2))
         k_min=sqrt(lattq(1,qmin(2))**2+lattq(2,qmin(2))**2+lattq(3,qmin(2))**2)
         cl_step = (s0/sp-1.0_dblprec)/2.0_dblprec/sin(qfac*k_min/2.0_dblprec)
         cl=real(sqrt(cl_step))
         write(*,'(2x,a20,2x,f11.5,2x,f11.5,2x,f11.5)') &
            'Displacement Correlation lengths:',cl(1),cl(2),cl(3)

         !write(*,*) s0/sp-1.0_dblprec,sin(qfac*k_min/2.0_dblprec);stop

         ! Transform G(k) to G(r)
         i_all=-product(shape(corr_s))*kind(corr_s)
         deallocate(corr_s,stat=i_stat)
         call memocc(i_stat,i_all,'corr_s','lattcalc_gk')
         allocate(corr_s(3,natom),stat=i_stat)
         call memocc(i_stat,product(shape(corr_s))*kind(corr_s),'corr_s','lattcalc_gk')
         corr_s=0.0_dblprec

         !$omp parallel do default(shared) private(r,iq,qdr,epowqr) schedule(static)
         do r=1,Natom
            do iq=1,lattnq
               qdr=lattq(1,iq)*coord(1,r)+lattq(2,iq)*coord(2,r)+lattq(3,iq)*coord(3,r)
               epowqr=cos(-qfac*qdr)
               corr_s(1,r)=corr_s(1,r)+epowqr*corr_k(1,iq)
               corr_s(2,r)=corr_s(2,r)+epowqr*corr_k(2,iq)
               corr_s(3,r)=corr_s(3,r)+epowqr*corr_k(3,iq)
            end do
         end do
         !$omp end parallel do

         ! Write G(r)
         write (filn,'(''ur.'',a,''.out'')') trim(simid)
         open(200,file=filn,status='replace')
         do r=1,Natom
            write(200,'(i10,3f10.4,4f18.8)') r,(coord(l,r),l=1,3),(((corr_s(l,r))),l=1,3),&
               sqrt(corr_s(1,r)**2+corr_s(2,r)**2+corr_s(3,r)**2)
         end do
         close(200)

         ! Write G(|r|)
         write (filn,'(''ura.'',a,''.out'')') trim(simid)
         open(200,file=filn,status='replace')
         do r=1,Natom
            write(200,'(7f18.8)') sqrt((coord(1,r)-coord(1,1))**2+(coord(2,r)-coord(2,1))**2+(coord(3,r)-coord(3,1))**2),&
               (((corr_s(l,r))),l=1,3),&
               sqrt(corr_s(1,r)**2+corr_s(2,r)**2+corr_s(3,r)**2)
         end do
         close(200)

         ! Deallocate arrays
         i_all=-product(shape(corr_k))*kind(corr_k)
         deallocate(corr_k,stat=i_stat)
         call memocc(i_stat,i_all,'corr_k','lattcalc_gk')

         i_all=-product(shape(corr_s))*kind(corr_s)
         deallocate(corr_s,stat=i_stat)
         call memocc(i_stat,i_all,'corr_s','lattcalc_gk')
      end if
      return
   end subroutine lattcalc_gk


   subroutine copy_q(nq, q, flag)

      integer :: nq  !< Number of q-points to sample
      real(dblprec), dimension(3,nq):: q  !< q-points
      integer :: flag

      integer :: i_stat,i_all

      if(flag==1) then
         lattnq = nq
         allocate(lattq(3,lattnq),stat=i_stat)
         call memocc(i_stat,product(shape(lattq))*kind(lattq),'lattq','copy_q')
         lattq = q
      else if (flag==-1) then
         i_all=-product(shape(lattq))*kind(lattq)
         deallocate(lattq,stat=i_stat)
         call memocc(i_stat,i_all,'lattq','copy_q')
      end if

   end subroutine copy_q


   !> Calculate suitable values of frequencies for U(q,t) -> U(q,w) transform
   subroutine lattset_w(delta_t, uc_step, uc_nstep)
      !
      use Constants, only : pi, hbar_mev
      !
      implicit none
      !
      real(dblprec), intent(in) :: delta_t !< Time step
      integer, intent(in) :: uc_step !< Separation between sampling steps
      integer, intent(in) :: uc_nstep !< Number of steps to sample
      !
      integer :: i_stat,j
      real(dblprec) :: dt !< Time step
      real(dblprec) :: dww
      real(dblprec) :: emin,emax
      !
      dt=uc_step*delta_t
      lattnw = uc_nstep
      allocate(lattw(lattnw),stat=i_stat)
      call memocc(i_stat,product(shape(lattw))*kind(lattw),'lattw','lattset_w')
      dww = 2*pi/(lattnw*dt)

      ! Previous convention was j*ddw
      do j=1,lattnw
         lattw(j)=(j-1)*dww
      end do

      emin=hbar_mev*(lattw(2)-lattw(1))
      emax=hbar_mev*lattw(lattnw)
      write(*,'(1x,a,f6.3,a,f6.1,a)') 'Phonon sampling between ',emin,' meV and ',emax,' meV.'
   end subroutine lattset_w


   !> Calculate correction to the sampling frequencies, w
   !> Perform correction to all frequencies >= uc_tidx
   subroutine lattcalc_corr_w(deltat_upd, uc_step, uc_nstep, uc_tidx)

      implicit none

      real(dblprec), intent(in) :: deltat_upd                       !< Updated time step
      integer, intent(in) :: uc_step                                !< Separation between sampling steps
      integer, intent(in) :: uc_nstep                               !< Number of steps to sample
      integer, intent(in) :: uc_tidx                                !< Sampling time interval indicator

      integer :: j
      real(dblprec) :: dt, pi
      real(dblprec) :: dww

      pi = 4.*atan(1.0_dblprec)

      dt = uc_step*deltat_upd ! Time interval between new samples
      lattnw = uc_nstep
      dww = 2*pi/(lattnw*dt)

      ! Perform correction of the current/future sampling frequencies
      do j=uc_tidx+1,lattnw ! The addition of 1 is due to function calling at the end of the previous sampling period
         lattw(j)=j*dww
      end do

   end subroutine lattcalc_corr_w


   !> Calculate U(q) for a cell filling mesh for transform to U(r) (sublattice projection)
   !! @todo Add non-diagonal components
   subroutine lattcalc_gk_proj(Natom, Mensemble, coord, simid, uvec, flag, NT, atype)
      !
      use Constants, only : pi
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: uvec  !< Current ionic displacement
      integer, intent(inout) :: flag  !< Allocate or deallocate (1/-1)
      integer, intent(in) :: NT  !< Number of types of atoms
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      !
      integer :: iq,r,l,i_stat,i_all,k,m,ic
      integer,dimension(nt) :: ncounter
      integer,dimension(natom,nt) :: indxcoord
      character(len=40) :: filn
      real(dblprec) :: epowqr!,s0,sp
      real(dblprec) :: qdr
      real(dblprec), dimension(nt) :: nainv

      if(flag==0) then

         ! First call, allocate and clear arrays
         allocate(corr_s_proj(3,lattnq,nt),stat=i_stat)
         call memocc(i_stat,product(shape(corr_s_proj))*kind(corr_s_proj),'corr_s_proj','calc_gk_proj')
         allocate(corr_k_proj(3,lattnq,nt,nt),stat=i_stat)
         call memocc(i_stat,product(shape(corr_k_proj))*kind(corr_k_proj),'corr_k_proj','calc_gk_proj')
         corr_k_proj=0.0_dblprec
         corr_s_proj=0.0_dblprec
         flag=1
         uc_samp_done_proj=0
      end if

      qfac=2._dblprec*pi
      wfac=1._dblprec
      ncounter=0
      indxcoord=0
      do k=1,nt
         do r=1,Natom
            if (k==atype(r)) then
               ncounter(k)=ncounter(k)+1
               indxcoord(ncounter(k),k)=r
            end if
         enddo
      enddo

      nainv(:)=1.0_dblprec/ncounter(:)
      if (flag==1) then

         ! Calculate s(k) for the current iteration and add to average of G(k)
         corr_s_proj=0.0_dblprec

         !$omp parallel do default(shared) private(r,iq,l,qdr,epowqr,k) schedule(static)
         do iq=1,lattnq
            do l=1,Mensemble
               corr_s_proj(:,iq,:)=0.0_dblprec
               do r=1,Natom
                  k=atype(r)
                  qdr=lattq(1,iq)*coord(1,r)+lattq(2,iq)*coord(2,r)+lattq(3,iq)*coord(3,r)
                  !SLDTODO Consider to use exponential Fourier transform instead
                  epowqr=cos(qfac*qdr)*nainv(k)
                  corr_s_proj(1,iq,k)=corr_s_proj(1,iq,k)+epowqr!*emomM(1,r,l)
                  corr_s_proj(2,iq,k)=corr_s_proj(2,iq,k)+epowqr!*emomM(2,r,l)
                  corr_s_proj(3,iq,k)=corr_s_proj(3,iq,k)+epowqr!*emomM(3,r,l)
               end do
               do m=1,nt
                  do k=1,nt
                     corr_k_proj(1,iq,k,m)=corr_k_proj(1,iq,k,m) + (corr_s_proj(1,iq,k)*corr_s_proj(1,iq,m))
                     corr_k_proj(2,iq,k,m)=corr_k_proj(2,iq,k,m) + (corr_s_proj(2,iq,k)*corr_s_proj(2,iq,m))
                     corr_k_proj(3,iq,k,m)=corr_k_proj(3,iq,k,m) + (corr_s_proj(3,iq,k)*corr_s_proj(3,iq,m))
                  enddo
               enddo
            end do
         end do
         !$omp end parallel do

         uc_samp_done_proj=uc_samp_done_proj+1
      end if

      if (flag==2) then

         ! Finish sampling and write S(q)
         corr_k_proj=corr_k_proj/(uc_samp_done_proj*Mensemble)
         !Write G(k)
         do k=1,nt
            do m=k,nt
               if(k>10) then
                  write (filn,'(''projuq.'',i2,''.'',i2,''.'',a,''.out'')') k,m,trim(simid)
               else
                  write (filn,'(''projuq.'',i1,''.'',i1,''.'',a,''.out'')') k,m,trim(simid)
               end if
               open(299, file=filn)
               do iq=1,lattnq
                  write(299,'(i5,3f10.4,4f18.8)') iq,(lattq(l,iq),l=1,3),(((corr_k_proj(l,iq,k,m))),l=1,3), &
                     sqrt(corr_k_proj(1,iq,k,m)**2+corr_k_proj(2,iq,k,m)**2+corr_k_proj(3,iq,k,m)**2)
               enddo
               close(299)
            end do
         end do

         ! Transform G(k) to G(r)
         i_all=-product(shape(corr_s_proj))*kind(corr_s_proj)
         deallocate(corr_s_proj,stat=i_stat)
         call memocc(i_stat,i_all,'corr_s_proj','calc_gk_proj')
         allocate(corr_ss_proj(3,natom,nt,nt),stat=i_stat)
         call memocc(i_stat,product(shape(corr_ss_proj))*kind(corr_ss_proj),'corr_ss_proj','calc_gk_proj')
         corr_ss_proj=0.0_dblprec

         !$omp parallel do default(shared) private(r,iq,qdr,epowqr,k,ic,m) schedule(static)
         do r=1,Natom
            m=atype(r)
            do k=1,nt
               do iq=1,lattnq
                  ic=indxcoord(1,k)
                  qdr=lattq(1,iq)*(coord(1,r)-coord(1,ic))+lattq(2,iq)*(coord(2,r)-coord(2,ic))+ &
                     lattq(3,iq)*(coord(3,r)-coord(3,ic))
                  epowqr=cos(-qfac*qdr)
                  corr_ss_proj(1,r,k,m)=corr_ss_proj(1,r,k,m)+epowqr*corr_k_proj(1,iq,k,m)
                  corr_ss_proj(2,r,k,m)=corr_ss_proj(2,r,k,m)+epowqr*corr_k_proj(2,iq,k,m)
                  corr_ss_proj(3,r,k,m)=corr_ss_proj(3,r,k,m)+epowqr*corr_k_proj(3,iq,k,m)
               end do
            end do
         end do
         !$omp end parallel do

         ! Write G(r)
         do k=1,nt
            ic=indxcoord(1,k)
            do m=k,nt
               if(k>10) then
                  write (filn,'(''projur.'',i2,''.'',i2,''.'',a,''.out'')') k,m,trim(simid)
               else
                  write (filn,'(''projur.'',i1,''.'',i1,''.'',a,''.out'')') k,m,trim(simid)
               end if
               open(298, file=filn)
               do r=1,Natom
                  if (m==atype(r)) then
                     write(298,'(i5,3f10.4,4f18.8)') r,(coord(1,r)-coord(1,ic)),(coord(2,r)-coord(2,ic)), &
                        (coord(3,r)-coord(3,ic)),(((corr_ss_proj(l,r,k,m))),l=1,3),&
                        sqrt(corr_ss_proj(1,r,k,m)**2+corr_ss_proj(2,r,k,m)**2+corr_ss_proj(3,r,k,m)**2)
                  end if
               enddo
               close(298)
            end do
         end do

         ! Write G(|r|)
         do k=1,nt
            ic=indxcoord(1,k)
            do m=k,nt
               if(k>10) then
                  write (filn,'(''projura.'',i2,''.'',i2,''.'',a,''.out'')') k,m,trim(simid)
               else
                  write (filn,'(''projura.'',i1,''.'',i1,''.'',a,''.out'')') k,m,trim(simid)
               end if
               open(297, file=filn)
               do r=1,Natom
                  if (m==atype(r)) then
                     write(297,'(7f18.8)') sqrt( (coord(1,r)-coord(1,ic))**2+(coord(2,r)-coord(2,ic))**2+ &
                        (coord(3,r)-coord(3,ic))**2),&
                        (((corr_ss_proj(l,r,k,m))),l=1,3),&
                        sqrt(corr_ss_proj(1,r,k,m)**2+corr_ss_proj(2,r,k,m)**2+corr_ss_proj(3,r,k,m)**2)
                  end if
               enddo
               close(297)
            end do
         end do

         ! Deallocate arrays
         i_all=-product(shape(corr_k_proj))*kind(corr_k_proj)
         deallocate(corr_k_proj,stat=i_stat)
         call memocc(i_stat,i_all,'corr_k_proj','calc_gk_proj')

         i_all=-product(shape(corr_ss_proj))*kind(corr_ss_proj)
         deallocate(corr_ss_proj,stat=i_stat)
         call memocc(i_stat,i_all,'corr_ss_proj','calc_gk_proj')
      end if
      return
   end subroutine lattcalc_gk_proj


   !> Calculate U(q,t) for obtaining U(q,w) after FT
   subroutine lattcalc_gkt(Natom, Mensemble, coord, simid, uvec, vvec, flag)
      !subroutine lattcalc_gkt(Natom, Mensemble, coord, simid, uvec, vvec, uc_tidx, uc_nstep, &
      !           flag,q_flag,do_connected,do_conv,sigma_q,sigma_w,LQfactor,LWfactor)
      !
      use Constants
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: uvec  !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: vvec  !< Current ionic velocity
      !integer, intent(inout) :: uc_tidx ! < Sampling time index
      !integer, intent(in) :: uc_nstep !< Number of steps to sample
      integer, intent(inout) :: flag  !< Setup, sample, or print
      !character, intent(in) :: q_flag
      !character(len=2), intent(in) :: do_conv
      !real(dblprec), intent(in) :: sigma_q, sigma_w
      !real(dblprec), intent(in) :: LQfactor,LWfactor
      !real(dblprec), dimension(uc_nstep+1), intent(in) :: udeltat_corr
      !character(len=1), intent(in) :: do_connected
      !
      !
      integer :: iq,iw,step,r,l,i_stat,i_all,j
      character(len=30) :: filn, filn2
      complex(dblprec) :: epowqr, i, iqfac,tt, epowwt
      real(dblprec) :: qdr, nainv

      real(dblprec), dimension(uc_nstep+1) :: dt
      !
      i=(0.0_dblprec,1.0_dblprec)

      if(flag==0) then
         ! First call, allocate and clear arrays
         allocate(corr_st(3,lattnq),stat=i_stat)
         call memocc(i_stat,product(shape(corr_st))*kind(corr_st),'corr_st','lattcalc_gkt')
         !
         allocate(corr_kt(3,lattnq,uc_nstep+1),stat=i_stat)
         call memocc(i_stat,product(shape(corr_kt))*kind(corr_kt),'corr_kt','lattcalc_gkt')
         corr_kt=0.0_dblprec
         !
         allocate(corr_vst(3,lattnq),stat=i_stat)
         call memocc(i_stat,product(shape(corr_vst))*kind(corr_vst),'corr_vst','lattcalc_gkt')
         !
         allocate(corr_vkt(3,lattnq,uc_nstep+1),stat=i_stat)
         call memocc(i_stat,product(shape(corr_vkt))*kind(corr_vkt),'corr_vkt','lattcalc_gkt')
         corr_vkt=0.0_dblprec
         !
         allocate(r_norm(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(r_norm))*kind(r_norm),'r_norm','lattcalc_gkt')

         r_max=0.0_dblprec
         do r=1,Natom
            r_norm(r)=coord(1,r)**2+coord(2,r)**2+coord(3,r)**2
            if(r_max<r_norm(r)) r_max=r_norm(r)
         end do

         flag=1
         uc_samp_done=0
         qfac=2.0_dblprec*pi
      end if

      nainv=1.0_dblprec/Natom

      ! Calculate g(k) for the current iteration and add to G(k,t)
      if (flag==1) then
         !SLDTODO Are the variables corr_st and corr_vst needed?
         corr_st=0.0_dblprec
         corr_vst=0.0_dblprec
         iqfac=i*qfac
         !$omp parallel do default(shared) private(r,iq,l,qdr,epowqr) schedule(static)
         do iq=1,lattnq
            do l=1,Mensemble
               do r=1,Natom
                  qdr=lattq(1,iq)*coord(1,r)+lattq(2,iq)*coord(2,r)+lattq(3,iq)*coord(3,r)
                  epowqr=exp(iqfac*qdr)*nainv
                  !write(*,*) 'qdr ', qdr, ' epowqr ', epowqr
                  corr_kt(1,iq,uc_tidx)=corr_kt(1,iq,uc_tidx)+epowqr*uvec(1,r,l)
                  corr_kt(2,iq,uc_tidx)=corr_kt(2,iq,uc_tidx)+epowqr*uvec(2,r,l)
                  corr_kt(3,iq,uc_tidx)=corr_kt(3,iq,uc_tidx)+epowqr*uvec(3,r,l)
                  !write(*,*) 'corr_kt(1,iq,uc_tidx) ', corr_kt(1,iq,uc_tidx)
                  corr_vkt(1,iq,uc_tidx)=corr_vkt(1,iq,uc_tidx)+epowqr*vvec(1,r,l)
                  corr_vkt(2,iq,uc_tidx)=corr_vkt(2,iq,uc_tidx)+epowqr*vvec(2,r,l)
                  corr_vkt(3,iq,uc_tidx)=corr_vkt(3,iq,uc_tidx)+epowqr*vvec(3,r,l)
                  !write(*,*) 'corr_vkt(1,iq,uc_tidx) ', corr_vkt(1,iq,uc_tidx)
               end do
            end do
         end do
         !$omp end parallel do
      end if

      ! Print G(k,t)
      if (flag==4) then
         write (filn,'(''uqt.'',a,''.out'')') trim(simid)
         open(200,file=filn, position='append')
         write (filn2,'(''vqt.'',a,''.out'')') trim(simid)
         open(201,file=filn2, position='append')
         do j=1,uc_nstep
            do iq=1,lattnq

               write(200,'(i10, i10,3f10.4,4es18.8)') j, iq,(lattq(l,iq),l=1,3), &
                  !write(200,'(i10, i10,3f10.4,4f18.8)') mstep, iq,(lattq(l,iq),l=1,3), &
               !(((corr_kt(l,iq,j))),l=1,3), &
               !sqrt(corr_kt(1,iq,j)**2+corr_kt(2,iq,j)**2+corr_kt(3,iq,j)**2)
               abs(corr_kt(1,iq,j)), abs(corr_kt(2,iq,j)), abs(corr_kt(3,iq,j)), &
                  sqrt(corr_kt(1,iq,j)**2+corr_kt(2,iq,j)**2+corr_kt(3,iq,j)**2)
               !corr_kt(1,iq,uc_tidx)=corr_kt(1,iq,uc_tidx)+epowqr*uvec(1,r,l)
               !corr_kt(2,iq,uc_tidx)=corr_kt(2,iq,uc_tidx)+epowqr*uvec(2,r,l)
               !corr_kt(3,iq,uc_tidx)=corr_kt(3,iq,uc_tidx)+epowqr*uvec(3,r,l)

               write(201,'(i10, i10,3f10.4,4es18.8)') j, iq,(lattq(l,iq),l=1,3), &
                  abs(corr_vkt(1,iq,j)), abs(corr_vkt(2,iq,j)), abs(corr_vkt(3,iq,j)), &
                  sqrt(corr_vkt(1,iq,j)**2+corr_vkt(2,iq,j)**2+corr_vkt(3,iq,j)**2)

            end do
         end do
         close(200)
         close(201)
      end if

      ! Final operations, transform and print
      if (flag==2) then

         ! Allocate arrays
         i_all=-product(shape(corr_st))*kind(corr_st)
         deallocate(corr_st,stat=i_stat)
         call memocc(i_stat,i_all,'corr_st','lattcalc_gkt')
         !
         allocate(corr_kw(3,lattnq,lattnw),stat=i_stat)
         call memocc(i_stat,product(shape(corr_kw))*kind(corr_kw),'corr_kw','lattcalc_gkt')
         corr_kw=0.0_dblprec

         i_all=-product(shape(corr_vst))*kind(corr_vst)
         deallocate(corr_vst,stat=i_stat)
         call memocc(i_stat,i_all,'corr_vst','lattcalc_gkt')
         !
         allocate(corr_vkw(3,lattnq,lattnw),stat=i_stat)
         call memocc(i_stat,product(shape(corr_vkw))*kind(corr_vkw),'corr_vkw','lattcalc_gkt')
         corr_vkw=0.0_dblprec

         ! Finish sampling and transform (k,t)->(k,w)
         !corr_kt=corr_kt/(sct_samp_done*Mensemble)
         !corr_kt=abs(corr_kt*corr_kt)

         if(uc_tidx.GT.uc_nstep) then
            uc_tidx = uc_nstep
         end if

         j=1
         !write(*,*) 'allocated? ucstep_arr ', allocated(ucstep_arr)
         !write(*,*) 'allocated? udeltat_corr ', allocated(udeltat_corr)
         !write(*,*) 'allocated? dt ', allocated(dt)
         do while(j.LE.uc_tidx)
            !write(*,*) 'j ', j, ' uc_tidx ', uc_tidx
            !write(*,*) 'shape(ucstep_arr) ', shape(ucstep_arr)
            !write(*,*) 'shape(udeltat_corr) ', shape(udeltat_corr)
            !write(*,*) 'shape(dt) ', shape(dt)
            dt(j) = ucstep_arr(j)*udeltat_corr(j)
            j = j+1
         end do
         !SLDTODO Fix udeltat_corr
         write(*,*) 'in latticecorrelation (latt)'
         write(*,*) 'dt(1:10) ', dt(1:10)
         !      dt=1.0e-1 !!!! SLDTODO HARD CODED !!!!
         write(*,*) 'dt(1:10) ', dt(1:10)

         wfac=1.0_dblprec
         corr_kw=0.0_dblprec
         corr_vkw=0.0_dblprec

         !!!$omp parallel do default(shared) private(iw,iq,step,tt,epowwt) schedule(static)
         write(*,*) 'lattnq ', lattnq, ' lattnw ', lattnw
         do iw=1,lattnw
            do step=1,uc_tidx
               tt=i*wfac*(step-1)*dt(step)
               epowwt=exp(lattw(iw)*tt)!*(0.54_dblprec-0.46*cos(2.0_dblprec*pi*(step-1)/(uc_tidx-1)))
               !epowwt=exp(w(iw)*tt)
               do iq=1,lattnq

                  corr_kw(1,iq,iw)=corr_kw(1,iq,iw)+epowwt*corr_kt(1,iq,step)
                  corr_kw(2,iq,iw)=corr_kw(2,iq,iw)+epowwt*corr_kt(2,iq,step)
                  corr_kw(3,iq,iw)=corr_kw(3,iq,iw)+epowwt*corr_kt(3,iq,step)

                  corr_vkw(1,iq,iw)=corr_vkw(1,iq,iw)+epowwt*corr_vkt(1,iq,step)
                  corr_vkw(2,iq,iw)=corr_vkw(2,iq,iw)+epowwt*corr_vkt(2,iq,step)
                  corr_vkw(3,iq,iw)=corr_vkw(3,iq,iw)+epowwt*corr_vkt(3,iq,step)

               enddo
               !write(*,*) '---------'
               !write(*,*) 'iq,iw ', iq,iw
               !write(*,*) 'corr_kw(1,iq,iw) ', corr_kw(1,iq,iw)
            enddo
            !corr_kw(:,iq,iw)=abs(corr_kw(:,iq,iw))!**2
         enddo
         !!!$omp end parallel do


         ! Write U(q,w)
         write (filn,'(''uqw.'',a,''.out'')') trim(simid)
         write (filn2,'(''vqw.'',a,''.out'')') trim(simid)
         open(9, file=filn)
         open(10, file=filn2)
         do iq=1,Lattnq
            do iw=1,Lattnw

               write (9,10005) iq,lattq(1,iq), lattq(2,iq),lattq(3,iq),iw, &
                  abs(corr_kw(1,iq,iw)),abs(corr_kw(2,iq,iw)),abs(corr_kw(3,iq,iw)), &
                  abs(corr_kw(1,iq,iw)**2+corr_kw(2,iq,iw)**2+corr_kw(3,iq,iw)**2)**0.5_dblprec

               write (10,10005) iq,lattq(1,iq), lattq(2,iq),lattq(3,iq),iw, &
                  abs(corr_vkw(1,iq,iw)),abs(corr_vkw(2,iq,iw)),abs(corr_vkw(3,iq,iw)), &
                  abs(corr_vkw(1,iq,iw)**2+corr_vkw(2,iq,iw)**2+corr_vkw(3,iq,iw)**2)**0.5_dblprec

            end do
         end do
         close(9)
         close(10)

         ! Write U(q,w) decomposed in real and imaginary parts
         write (filn,'(''cuqw.'',a,''.out'')') trim(simid)
         write (filn2,'(''cvqw.'',a,''.out'')') trim(simid)
         open(29, file=filn)
         open(30, file=filn2)
         do iq=1,Lattnq
            do iw=1,Lattnw

               write (29,10005) iq, lattq(1,iq), lattq(2,iq), lattq(3,iq),iw, &
                  abs(corr_kw(1,iq,iw)),atan2(aimag(corr_kw(1,iq,iw)),real(corr_kw(1,iq,iw))), &
                  abs(corr_kw(2,iq,iw)),atan2(aimag(corr_kw(2,iq,iw)),real(corr_kw(2,iq,iw))), &
                  abs(corr_kw(3,iq,iw)),atan2(aimag(corr_kw(3,iq,iw)),real(corr_kw(3,iq,iw)))

               write (30,10005) iq, lattq(1,iq), lattq(2,iq), lattq(3,iq),iw, &
                  abs(corr_vkw(1,iq,iw)),atan2(aimag(corr_vkw(1,iq,iw)),real(corr_vkw(1,iq,iw))), &
                  abs(corr_vkw(2,iq,iw)),atan2(aimag(corr_vkw(2,iq,iw)),real(corr_vkw(2,iq,iw))), &
                  abs(corr_vkw(3,iq,iw)),atan2(aimag(corr_vkw(3,iq,iw)),real(corr_vkw(3,iq,iw)))

            end do
         end do
         close(29)
         close(30)

         !if (do_conv.eq.'LQ' .or. do_conv.eq.'LW' .or. do_conv.eq.'LW' .or. &
         !    do_conv.eq.'GQ' .or. do_conv.eq.'GW' .or. do_conv.eq.'GY') then
         !   ! Calculate the convolution for the phonon DOS
         !   call lattcalc_dos_conv(do_conv,sigma_q,sigma_w,LQfactor,LWfactor,corr_kw)
         !endif

         ! Calculate and write the phonon DOS
         allocate(phonon_dos(3,lattnw/2),stat=i_stat)
         call memocc(i_stat,product(shape(phonon_dos))*kind(phonon_dos),'phonon_dos','lattcalc_gkt')
         !
         write (filn,'(''phondos.'',a,''.out'')') trim(simid)
         open(9, file=filn)
         !       write (9,'(a)') "#    E(meV)            DOS_x           DOS_y           DOS_z          DOS_tot"
         phonon_dos=0.0_dblprec
         !if (q_flag=='I'.or.q_flag=='B') then
         !    do iw=1,Lattnw/2
         !       do iq=1,Lattnq
         !          phonon_dos(1,iw)=phonon_dos(1,iw)+abs(corr_kw(1,iq,iw))*q_weight(iq) ! Put a weight for the IR BZ points
         !          phonon_dos(2,iw)=phonon_dos(2,iw)+abs(corr_kw(2,iq,iw))*q_weight(iq)
         !          phonon_dos(3,iw)=phonon_dos(3,iw)+abs(corr_kw(3,iq,iw))*q_weight(iq)
         !       end do
         !       write (9,10006) hbar_mev*lattw(iw),phonon_dos(1,iw),phonon_dos(2,iw), phonon_dos(3,iw), &
         !            sqrt(phonon_dos(1,iw)**2+phonon_dos(2,iw)**2+phonon_dos(3,iw)**2)
         !    end do
         !else
         do iw=1,Lattnw/2
            do iq=1,Lattnq
               phonon_dos(1,iw)=phonon_dos(1,iw)+abs(corr_kw(1,iq,iw)) ! Put a weight for the IR BZ points
               phonon_dos(2,iw)=phonon_dos(2,iw)+abs(corr_kw(2,iq,iw))
               phonon_dos(3,iw)=phonon_dos(3,iw)+abs(corr_kw(3,iq,iw))
            end do
            write (9,10006) hbar_mev*lattw(iw),phonon_dos(1,iw),phonon_dos(2,iw), phonon_dos(3,iw), &
               sqrt(phonon_dos(1,iw)**2+phonon_dos(2,iw)**2+phonon_dos(3,iw)**2)
         end do
         !endif
         close(9)

         ! Deallocate arrays
         i_all=-product(shape(corr_kt))*kind(corr_kt)
         deallocate(corr_kt,stat=i_stat)
         call memocc(i_stat,i_all,'corr_kt','lattcalc_gkt')
         !
         i_all=-product(shape(corr_kw))*kind(corr_kw)
         deallocate(corr_kw,stat=i_stat)
         call memocc(i_stat,i_all,'corr_kw','lattcalc_gkt')
         !
         i_all=-product(shape(corr_vkt))*kind(corr_vkt)
         deallocate(corr_vkt,stat=i_stat)
         call memocc(i_stat,i_all,'corr_vkt','lattcalc_gkt')
         !
         i_all=-product(shape(corr_vkw))*kind(corr_vkw)
         deallocate(corr_vkw,stat=i_stat)
         call memocc(i_stat,i_all,'corr_vkw','lattcalc_gkt')
         !
         i_all=-product(shape(phonon_dos))*kind(phonon_dos)
         deallocate(phonon_dos,stat=i_stat)
         call memocc(i_stat,i_all,'phonon_dos','lattcalc_gkt')
         !
         i_all=-product(shape(r_norm))*kind(r_norm)
         deallocate(r_norm,stat=i_stat)
         call memocc(i_stat,i_all,'r_norm','lattcalc_gkt')
         !
         !i_all=-product(shape(m_proj))*kind(m_proj)
         !deallocate(m_proj,stat=i_stat)
         !call memocc(i_stat,i_all,'m_proj','lattcalc_gkt')
         !
      end if
      return
      !
      10005 format (i7,3f10.6,2x,i7,2x,9es16.8)
      10006 format (5G16.8)
      !
   end subroutine lattcalc_gkt


   !SLDTODO Adopt to LD dynamics
   !< Perform only spatial correlation in real space to be able to deal with non periodic systems
   subroutine lattcalc_ur(Natom, Mensemble, coord, simid, emomM, cr_flag)
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      integer, intent(inout) :: cr_flag  !< Allocate or deallocate (1/-1)
      !
      integer :: iatom,r,l,i_stat,i_all
      character(len=30) :: filn
      real(dblprec), dimension(3,Natom) :: connected
      real(dblprec) :: nainv

      if(cr_flag==0) then
         ! First call, allocate and clear arrays
         allocate(corr_sr(3,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(corr_sr))*kind(corr_s),'corr_sr','calc_sr')
         corr_sr=0.0_dblprec
         cr_flag=1
         uc_samp_done=0
         !return
      end if

      nainv=1.0_dblprec/Natom

      if (cr_flag==1) then
         ! Calculate s(k) for the current iteration and add to average of G(k)

         !$omp parallel do default(shared) private(iatom,l) schedule(static)
         do l=1,Mensemble
            do iatom=1,Natom
               connected(1,iatom)=connected(1,iatom)+emomM(1,iatom,l)
               connected(2,iatom)=connected(2,iatom)+emomM(2,iatom,l)
               connected(3,iatom)=connected(3,iatom)+emomM(3,iatom,l)
            enddo
         enddo
         !$omp end parallel do

         !$omp parallel do default(shared) private(r,iatom,l) schedule(static)
         do l=1,Mensemble
            do iatom=1,Natom             
               do r=1,Natom
                  corr_sr(1,iatom)=corr_sr(1,iatom)+emomM(1,iatom,l)*emomM(1,r,l)*nainv**2-connected(1,iatom)*connected(1,r)/Mensemble**2
                  corr_sr(2,iatom)=corr_sr(2,iatom)+emomM(2,iatom,l)*emomM(1,r,l)*nainv**2-connected(2,iatom)*connected(2,r)/Mensemble**2
                  corr_sr(3,iatom)=corr_sr(3,iatom)+emomM(3,iatom,l)*emomM(1,r,l)*nainv**2-connected(3,iatom)*connected(3,r)/Mensemble**2
               end do
            enddo
         end do
         !$omp end parallel do

         uc_samp_done=uc_samp_done+1
      end if

      if (cr_flag==2) then
         ! Finish sampling and write S(q)
         corr_sr=corr_sr/(uc_samp_done*Mensemble)

         ! Write G(r)
         write (filn,'(''dir_sr.'',a,''.out'')') trim(simid)
         open(200,file=filn,status='replace')
         do r=1,Natom
            write(200,'(i10,3f10.4,4f18.8)') r,(coord(l,r),l=1,3),(((corr_sr(l,r))),l=1,3),&
               sqrt(corr_sr(1,r)**2+corr_sr(2,r)**2+corr_sr(3,r)**2)
         end do
         close(200)

         ! Write G(|r|)
         write (filn,'(''dir_sra.'',a,''.out'')') trim(simid)
         open(200,file=filn,status='replace')
         do r=1,Natom
            write(200,'(7f18.8)') sqrt((coord(1,r)-coord(1,1))**2+(coord(2,r)-coord(2,1))**2+(coord(3,r)-coord(3,1))**2),&
               (((corr_sr(l,r))),l=1,3),&
               sqrt(corr_sr(1,r)**2+corr_sr(2,r)**2+corr_sr(3,r)**2)
         end do
         close(200)

         ! Deallocate arrays
         i_all=-product(shape(corr_sr))*kind(corr_s)
         deallocate(corr_sr,stat=i_stat)
         call memocc(i_stat,i_all,'corr_sr','calc_sr')
      end if
      return

   end subroutine lattcalc_ur


   !< Perform convolutions of the magnon DOS
   subroutine lattcalc_dos_conv(do_conv,sigma_q,sigma_w,LQfactor,LWfactor,corr_kw)

      use Constants

      implicit none

      character(len=2), intent(in) :: do_conv
      real(dblprec), intent(in) :: sigma_q, sigma_w
      real(dblprec), intent(in) :: LQfactor,LWfactor

      complex(dblprec), dimension(3,lattnq,lattnw), intent(inout) :: corr_kw

      integer :: iq,iw,j,t, conv_range_w, conv_range_q
      integer :: u, qq
      real(dblprec) :: conv_cutoff_w, facw1, facw2, sfacww ! Variables for the convolution in frequency
      real(dblprec) :: conv_cutoff_q, facq1, facq2, sfacqq ! Variables for the convolution in qpoints

      ! Variables for the convolutions in frequency and qpoints domain, this are optional for the magnon DOS calculations
      ! This are the variables for qpoints convolutions
      if (do_conv=='GQ') then
         ! Variables for the reciprocal space for Gaussian function
         facq1 = -1.0_dblprec/(2*sigma_q**2)
         facq2 = 1.0_dblprec/(sqrt(2*pi)*sigma_q*lattnq)
         conv_cutoff_q=0.01_dblprec
         conv_range_q=int(sqrt(log(conv_cutoff_q)/facq1)+0.5_dblprec)

      else if (do_conv=='LQ') then
         ! Variables for the reciprocal space for Lorentz distribution
         facq1=-1.0_dblprec/(LQfactor**2)
         facq2=LQfactor/(Lattnq*pi)
         conv_cutoff_q=0.01_dblprec
         conv_range_q=int(sqrt(log(conv_cutoff_q)/facq1)+0.5_dblprec)

      else if (do_conv=='GW') then
         ! Variables for the frequency convolutions for Gaussian function
         facw1 = -1.0_dblprec/(2*sigma_w**2)
         facw2 = 1.0_dblprec/(sqrt(2*pi)*sigma_w*lattnw)
         conv_cutoff_w=0.01_dblprec
         conv_range_w=int(sqrt(log(conv_cutoff_w)/facw1)+0.5_dblprec)

      else if (do_conv=='LW') then
         ! Variables for the frequency convolution with Lorentz distribution
         facw1=-1.0_dblprec/(LWfactor**2)
         facw2=LWfactor/(Lattnw*pi)
         conv_cutoff_w=0.01_dblprec
         conv_range_w=int(sqrt(log(conv_cutoff_w)/facw1)+0.5_dblprec)

      else if (do_conv=='GY') then
         ! Variables for both qpoints and frequencies convolutions for Gaussian function
         facq1 = -1.0_dblprec/(2*sigma_q**2)
         facw1 = -1.0_dblprec/(2*sigma_w**2)
         facq2 = 1.0_dblprec/(sqrt(2*pi)*sigma_q*lattnq)
         facw2 = 1.0_dblprec/(sqrt(2*pi)*sigma_w*lattnw)
         conv_cutoff_q=0.01_dblprec
         conv_cutoff_w=0.01_dblprec
         conv_range_q=int(sqrt(log(conv_cutoff_q)/facq1)+0.5_dblprec)
         conv_range_w=int(sqrt(log(conv_cutoff_q)/facw1)+0.5_dblprec)

      else if (do_conv=='LY') then
         ! Variables for both qpoints and frequencies convolutions for Lorentz distribution
         facq1 = -1.0_dblprec/(LQfactor**2)
         facw1 = -1.0_dblprec/(LWfactor**2)
         facq2 = LQfactor/(Lattnq*pi)
         facw2 = LWfactor/(Lattnw*pi)
         conv_cutoff_q=0.01_dblprec
         conv_cutoff_w=0.01_dblprec
         conv_range_q=int(sqrt(log(conv_cutoff_q)/facq1)+0.5_dblprec)
         conv_range_w=int(sqrt(log(conv_cutoff_q)/facw1)+0.5_dblprec)

      endif

      ! Do the convolution of the corr_kw to smooth out the magnon DOS
      if (do_conv=='GQ') then
         do iq=1,lattnq
            !$omp parallel do default(shared) private(iw,u,qq,sfacqq)
            do iw=1,lattnw/2
               if (iw==1.and.iq==1) write(*,*) '-------------------------> CONV_GQ'
               ! Convolution with Gaussian resolution function
               do u = max(1,iq-conv_range_q), min(lattnq-1,iq+conv_range_q)
                  qq=(iq-u)
                  sfacqq=exp(facq1*qq**2)*facq2 ! Factor which controls the convolution
                  ! Convolution over the reciprocal space
                  corr_kw(1,iq,iw)=corr_kw(1,iq,iw)+corr_kw(1,iq,iw)*sfacqq
                  corr_kw(2,iq,iw)=corr_kw(2,iq,iw)+corr_kw(2,iq,iw)*sfacqq
                  corr_kw(3,iq,iw)=corr_kw(3,iq,iw)+corr_kw(3,iq,iw)*sfacqq
               enddo
            enddo
            !$omp end parallel do
         enddo
      else if (do_conv=='GW') then
         do iq=1,Lattnq
            !$omp parallel do default(shared) private(iw,j,t,sfacww)
            do iw=1, lattnw/2
               if (iw==1.and.iq==1) write(*,*) '-------------------------> CONV_GW'
               ! Convolution with Gaussian resolution function
               do j = max(1,iw-conv_range_w), min(lattnw-1,iw+conv_range_w)
                  t=(iw-j)
                  sfacww=exp(facw1*t**2)*facw2 ! This is the parameter that controls the convolution (Should the sum be over the absolute values?)
                  ! Convolution over frequency space with a gaussian resolution function
                  corr_kw(1,iq,iw)=corr_kw(1,iq,iw)+corr_kw(1,iq,iw)*sfacww
                  corr_kw(2,iq,iw)=corr_kw(2,iq,iw)+corr_kw(2,iq,iw)*sfacww
                  corr_kw(3,iq,iw)=corr_kw(3,iq,iw)+corr_kw(3,iq,iw)*sfacww
               enddo
            enddo
            !$omp end parallel do
         enddo
      else if (do_conv=='GY') then
         do iq=1,Lattnq
            !$omp parallel do default(shared) private(iw,u,qq,sfacqq)
            do iw=1, lattnw/2
               if (iw==1.and.iq==1) write(*,*) '-------------------------> CONV_GY'
               ! Convolution with Gaussian resolution function
               do u = max(1,iq-conv_range_q), min(lattnq-1,iq+conv_range_q)
                  qq=(iq-u)
                  sfacqq=exp(facq1*qq**2)*facq2 ! This parameter controls the convolution
                  ! Convolution over the qpoints
                  corr_kw(1,iq,iw)=corr_kw(1,iq,iw)+corr_kw(1,iq,iw)*sfacqq
                  corr_kw(2,iq,iw)=corr_kw(2,iq,iw)+corr_kw(2,iq,iw)*sfacqq
                  corr_kw(3,iq,iw)=corr_kw(3,iq,iw)+corr_kw(3,iq,iw)*sfacqq
               enddo
            enddo
            !$omp end parallel do
         enddo
         do iq=1,Lattnq
            !$omp parallel do default(shared) private(iw,j,t,sfacww)
            do iw=1, lattnw/2
               if (iw==1.and.iq==1) write(*,*) '-------------------------> CONV_GY'
               ! Convolution with Gaussian resolution function
               do j = max(1,iw-conv_range_w), min(lattnw-1,iw+conv_range_w)
                  t=(iw-j)
                  sfacww=exp(facw1*t**2)*facw2 ! This is the parameter that controls the convolution (Should the sum be over the absolute values?)
                  ! Convolution over frequency space with a gaussian resolution function
                  corr_kw(1,iq,iw)=corr_kw(1,iq,iw)+corr_kw(1,iq,iw)*sfacww
                  corr_kw(2,iq,iw)=corr_kw(2,iq,iw)+corr_kw(2,iq,iw)*sfacww
                  corr_kw(3,iq,iw)=corr_kw(3,iq,iw)+corr_kw(3,iq,iw)*sfacww
               enddo
               !             write(*,*) "CORR", iq, iw, corr_kw_dos(1:3,iq,iw)
            enddo
            !$omp end parallel do
         enddo
      else if (do_conv=='LY') then
         do iq=1,Lattnq
            !$omp parallel do default(shared) private(iw,u,qq,sfacqq)
            do iw=1, lattnw/2
               if (iw==1.and.iq==1) write(*,*) '-------------------------> CONV_LY'
               ! Convolution with Gaussian resolution function
               do u = max(1,iq-conv_range_q), min(lattnq-1,iq+conv_range_q)
                  qq=(iq-u)
                  sfacqq=facq1/(LQfactor**2 +qq**2) ! This parameter controls the convolution
                  ! Convolution over the qpoints
                  corr_kw(1,iq,iw)=corr_kw(1,iq,iw)+corr_kw(1,iq,iw)*sfacqq
                  corr_kw(2,iq,iw)=corr_kw(2,iq,iw)+corr_kw(2,iq,iw)*sfacqq
                  corr_kw(3,iq,iw)=corr_kw(3,iq,iw)+corr_kw(3,iq,iw)*sfacqq
               enddo
            enddo
            !$omp end parallel do
         enddo
         do iq=1,Lattnq
            !$omp parallel do default(shared) private(iw,j,t,sfacww)
            do iw=1, Lattnw/2
               if (iw==1.and.iq==1) write(*,*) '--------------------------> CONV_LY'
               do j = max(1,iw-conv_range_w), min(lattnw-1,iw+conv_range_w)
                  t=(iw-j)
                  sfacww=facw1/(LWfactor**2+t**2)
                  corr_kw(1,iq,iw)=corr_kw(1,iq,iw)+corr_kw(1,iq,iw)*sfacww
                  corr_kw(2,iq,iw)=corr_kw(2,iq,iw)+corr_kw(2,iq,iw)*sfacww
                  corr_kw(3,iq,iw)=corr_kw(3,iq,iw)+corr_kw(3,iq,iw)*sfacww
               enddo
            enddo
            !$omp end parallel do
         enddo
      else if (do_conv=='LW') then
         do iq=1,Lattnq
            !$omp parallel do default(shared) private(iw,j,t,sfacww)
            do iw=1, Lattnw/2
               if (iw==1.and.iq==1) write(*,*) '--------------------------> CONV_LW'
               do j = max(1,iw-conv_range_w), min(lattnw-1,iw+conv_range_w)
                  t=(iw-j)
                  sfacww=facw1/(LWfactor**2+t**2)
                  corr_kw(1,iq,iw)=corr_kw(1,iq,iw)+corr_kw(1,iq,iw)*sfacww
                  corr_kw(2,iq,iw)=corr_kw(2,iq,iw)+corr_kw(2,iq,iw)*sfacww
                  corr_kw(3,iq,iw)=corr_kw(3,iq,iw)+corr_kw(3,iq,iw)*sfacww
               enddo
            enddo
            !$omp end parallel do
         enddo
      else if (do_conv=='LQ') then
         do iq=1,Lattnq
            !$omp parallel do default(shared) private(iw,u,qq,sfacqq)
            do iw=1, lattnw/2
               if (iw==1.and.iq==1) write(*,*) '-------------------------> CONV_LQ'
               ! Convolution with Gaussian resolution function
               do u = max(1,iq-conv_range_q), min(lattnq-1,iq+conv_range_q)
                  qq=(iq-u)
                  sfacqq=facq1/(LQfactor**2 +qq**2) ! This parameter controls the convolution
                  ! Convolution over the qpoints
                  corr_kw(1,iq,iw)=corr_kw(1,iq,iw)+corr_kw(1,iq,iw)*sfacqq
                  corr_kw(2,iq,iw)=corr_kw(2,iq,iw)+corr_kw(2,iq,iw)*sfacqq
                  corr_kw(3,iq,iw)=corr_kw(3,iq,iw)+corr_kw(3,iq,iw)*sfacqq
               enddo
            enddo
            !$omp end parallel do
         enddo
      endif

      return

   end subroutine lattcalc_dos_conv


   !> Allocation of adapt timestep
   subroutine lattallocate_deltatcorr(uc_nstep, allocate_flag)

      implicit none

      integer, intent(in) :: uc_nstep                                 !< Number of steps to sample
      logical, intent(in) :: allocate_flag                            !< Allocate/deallocate

      integer :: i_stat, i_all

      if(allocate_flag) then
         if(.not.allocated(udeltat_corr)) then
            allocate(udeltat_corr(uc_nstep+1),stat=i_stat)
            call memocc(i_stat,product(shape(udeltat_corr))*kind(udeltat_corr),'udeltat_corr','lattallocate_deltatcorr')
         end if
         if(.not.allocated(ucstep_arr)) then
            allocate(ucstep_arr(uc_nstep+1),stat=i_stat)
            call memocc(i_stat,product(shape(ucstep_arr))*kind(ucstep_arr),'ucstep_arr','lattallocate_deltatcorr')
         end if
      else
         if(allocated(udeltat_corr)) then
            i_all=-product(shape(udeltat_corr))*kind(udeltat_corr)
            deallocate(udeltat_corr,stat=i_stat)
            call memocc(i_stat,i_all,'udeltat_corr','lattallocate_deltatcorr')
         end if
         if(allocated(ucstep_arr)) then
            i_all=-product(shape(ucstep_arr))*kind(ucstep_arr)
            deallocate(ucstep_arr,stat=i_stat)
            call memocc(i_stat,i_all,'ucstep_arr','lattallocate_deltatcorr')
         end if
      end if

   end subroutine lattallocate_deltatcorr


end module LatticeCorrelation
