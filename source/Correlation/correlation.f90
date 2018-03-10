!-------------------------------------------------------------------------------
! MODULE: CORRELATION
!> @brief Data and routines for calculate connected spin correlation function \f$ \mathbf{S}\left(\mathbf{r},t\right)\f$
!> and Fourier transforms \f$\mathbf{S}\left(\mathbf{q},\omega\right)\f$
!> @details
!> In this routine magnon disperion relations, magnon density of states, etc. are calculated from the
!> time and spatially displaced correlation function \f$ C^k \left(\mathbf{r}-\mathbf{r'},t\right)\f$ which is defined as
!> \f$ C^k (\mathbf{r}-\mathbf{r'},t) = \langle m^k_{\mathbf{r}}(t) m^k_{\mathbf{r'}}(0) \rangle - \langle m^k_{\mathbf{r}}(t) \rangle \langle m^k_{\mathbf{r'}}(0) \rangle\f$.
!> Using this one can then calculate the dynamical structure factor \f$ \mathbf{S}\left(\mathbf{q},\omega\right)\f$, via a Fourier transforms
!> of the correlation function
!> \f$ S^k(\mathbf{q},\omega) = \frac{1}{\sqrt{2\pi}N} \sum_{\mathbf{r},\mathbf{r'}} e^{i\mathbf{q}\cdot(\mathbf{r}-\mathbf{r'})} \int_{-\infty}^{\infty} e^{i\omega t} C^k (\mathbf{r}-\mathbf{r'},t) dt\f$
!> @author
!> A. Bergman, L. Bergqvist, J. Hellsvik, J. Chico
!! @todo Automatic generation of q-points for do_sc="Y"
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License.
!! See http://www.gnu.org/copyleft/gpl.txt
!-------------------------------------------------------------------------------
module Correlation
   use Parameters
   use Profiling
   !
   implicit none
   !

   !Spin-spin correlation parameters from input file
   character(len=1) :: do_sc        !< Perform spin-correlation sampling (Y/N/C)
   character(len=1) :: do_sc_bimag  !< Perform spin-correlation sampling of bi-magnons (Y/N)
   character(len=1) :: do_sr        !< Perform spin-correlation sampling in real space direclty (Y/N)
   character(len=2) :: do_conv      !< Flag for convolution of Magnon DOS Lorentzian(LY/LG/LW/N) Gaussian(GY/GW/GQ/N)
   character(len=1) :: qpoints      !< Flag for q-point generation (F=file,A=automatic,C=full cell)
   character(len=1) :: do_sc_proj   !< Measure sublattice projection of S(q,w) (Y/N/C)
   character(len=1) :: do_sc_projch   !< Measure chemical sublattice projection of S(q,w) (Y/N/C)
   character(len=1) :: do_sc_local_axis   !< Perform SQW along local quantization axis (Y/N)
   character(len=1) :: do_sc_proj_axis   !< Perform projected SQW along local quantization axis (Y/N)
   character(len=1) :: do_sc_dosonly !< Do not print s(q,w), only Magnon DOS (Y/N)
   character(len=1) :: do_sc_complex !< Print the complex values s(q,w) (Y/N)
   character(len=1) :: do_qt_traj   !< Measure time trajectory of S(q,t) (Y/N)
   character(len=1) :: do_connected !< Perform the connected part S(q,w)
   character(LEN=35) :: qfile       !< File name for q-points
   integer :: sc_sep                !< Separation between averages
   integer :: sc_mode               !< Spin correlation mode (0-3)
   integer :: sc_step               !< Separation between sampling steps
   integer :: sc_nstep              !< Number of steps to sample
   real(dblprec) :: sigma_q         !< Sigma parameter in Q for the Gaussian convolution
   real(dblprec) :: sigma_w         !< Sigma parameter in W for the Gaussian convolution
   real(dblprec) :: LWfactor        !< Gamma parameter in W for Lorentzian convolution
   real(dblprec) :: LQfactor        !< Gamma parameter in Q for Lorentzian convolution

   ! Working variables to perform the printing of the correlation
   integer :: nq  !< Number of q-points to sample
   integer :: nw  !< Number of frequencies to sample
   integer :: sc_samp_done !< Flag to keep track of if S(q) sampling is done (for do_sc='C')
   integer :: sc_samp_done_sr !< Flag to keep track of if S(r) sampling is done (for do_sr='Y')
   integer :: sc_samp_done_proj !< Flag to keep track of if S(q) sampling is done (for do_sc_proj='C')
   real(dblprec) :: qfac
   real(dblprec) :: wfac
   real(dblprec) :: r_max       ! Max length of distance vectors
   real(dblprec) :: nrinv
   real(dblprec) :: na2inv
   real(dblprec) :: nscinv
   real(dblprec) :: nrscinv
   real(dblprec) :: scnstepinv
   real(dblprec), dimension(:), allocatable :: w                     !< Frequencies
   real(dblprec), dimension(:), allocatable :: r_norm                ! Length of distance vectors
   real(dblprec), dimension(:), allocatable :: q_weight              !< Weights of the q points for DOS calculations
   real(dblprec), dimension(:,:), allocatable :: q                   !< q-points
   real(dblprec), dimension(:,:), allocatable :: corr_k              ! Correlation in q G(k)
   real(dblprec), dimension(:,:), allocatable :: corr_s              ! Correlation in r G(r)
   real(dblprec), dimension(:,:), allocatable :: corr_sr             ! Correlation in r G(r) calculated directly
   real(dblprec), dimension(:,:), allocatable :: corr_kt0            ! Temporary for g(q,t0)
   real(dblprec), dimension(:,:,:), allocatable :: corr_s_proj       ! Correlation in q or r (G(k),G(r)) (sublattice)
   real(dblprec), dimension(:,:,:,:), allocatable :: corr_k_proj     ! Correlation in q G(k) (sublattice)
   real(dblprec), dimension(:,:,:,:), allocatable :: corr_ss_proj    ! Correlation in r G(r) (sublattice)
   complex(dblprec), dimension(:,:,:), allocatable :: sqw            !< S(q,w)
   complex(dblprec), dimension(:,:,:), allocatable :: corr_kt        ! Correlation for G(k,t)
   complex(dblprec), dimension(:,:,:), allocatable :: corr_kw        ! Correlation for G(k,w)
   complex(dblprec), dimension(:,:), allocatable :: corr_bkt        ! Correlation for B(k,t)
   complex(dblprec), dimension(:,:), allocatable :: corr_bkw        ! Correlation for B(k,w)
   complex(dblprec), dimension(:,:,:), allocatable :: corr_st_proj   ! Temporary for g(r)
   complex(dblprec), dimension(:,:,:,:), allocatable :: corr_kt_proj ! Correlation for G(k,t)
   complex(dblprec), dimension(:,:,:,:), allocatable :: corr_kw_proj ! Correlation for G(k,w)

   real(dblprec), dimension(:,:), allocatable :: magnon_dos !< Magnon density of states
   real(dblprec), dimension(:), allocatable :: bimagnon_dos !< Bi-magnon density of states
   integer,dimension(2) :: qmin  !< Index of smallest wave vector and gamma point

   real(dblprec), dimension(:), allocatable :: deltat_corr ! Array storing delta_t for each sample
   real(dblprec), dimension(:), allocatable :: scstep_arr  ! Array storing sc_step for each sample

   real(dblprec), dimension(:,:), allocatable :: mavg_axis      !< Internal array for finding the global direction of magnetization
   real(dblprec), dimension(:,:,:), allocatable ::  bm_proj      !< Array containing bi-moments projected to parallel/perpendicular directions

   real(dblprec), dimension(:,:,:), allocatable :: mavg_local_axis!< Internal array for finding the local average direction of magnetization
   real(dblprec), dimension(:,:,:,:), allocatable :: mavg_local_rotmat !< Internal array for finding rotatons to project the local average direction of magnetization
   real(dblprec), dimension(:,:,:), allocatable :: mort_axis    !< Orthogonal (and parallel ) components of magnetization
   real(dblprec), dimension(:,:,:), allocatable ::  m_proj      !< Array containing moments projected to parallel/perpendicular directions
   !
   real(dblprec), dimension(:,:,:), allocatable :: mavg_proj_axis!< Internal array for finding the local average direction of magnetization
   real(dblprec), dimension(:,:,:,:), allocatable :: mavg_proj_rotmat !< Internal array for finding rotatons to project the local average direction of magnetization
   real(dblprec), dimension(:,:,:), allocatable :: mort_axis_p    !< Orthogonal (and parallel ) components of magnetization
   real(dblprec), dimension(:,:,:), allocatable ::  m_proj_p      !< Array containing moments projected to parallel/perpendicular directions
   !
   integer :: sc_window_fun  !< Choice of FFT window function (1=box, 2=Hann, 3=Hamming, 4=Blackman-Harris)
   real(dblprec) :: sc_local_axis_mix  !< Mixing parameter for updating the local quantization axis. Rec. value: 0.1-0.5

   !
   real(dblprec), dimension(3) :: r_mid !< Local variable containing the center coordinate for the system. Used for FT in calc_gk

   integer :: sc_tidx

   private

   public :: do_sc, q, nq, calc_corr_w, w, sc_nstep, do_conv,sigma_q,sigma_w,LQfactor,LWfactor, do_sc_proj
   public :: setup_qcoord, read_q, set_w, qpoints
   public :: correlation_init, read_parameters_correlation,  allocate_deltatcorr, deallocate_q
   public :: correlation_wrapper

contains

   !---------------------------------------------------------------------------
   ! SUBROUTINE: correlation_wrapper
   !> @brief
   !> Driver for correlation function calculations
   !---------------------------------------------------------------------------
   subroutine correlation_wrapper(Natom, Mensemble, coord, simid, emomM, mstep, delta_t, NT, atype, Nchmax, achtype, flag, flag_p )
      !
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      integer, intent(in) :: mstep  !< Current simulation step
      real(dblprec), intent(in) :: delta_t               !< Time step
      integer, intent(in) :: NT  !< Number of types of atoms
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      integer, intent(in) :: Nchmax  !< Number of chemical types
      integer, dimension(Natom), intent(in) :: achtype !< Type of atom
      integer, intent(inout) :: flag  !< Setup, sample, or print
      integer, intent(inout) :: flag_p  !< Setup, sample, or print
      !
      if(flag==0) sc_tidx=0

      ! Sample S(r) through S(q)
      if(do_sc=='C') then

         if(mod(mstep-1,sc_sep)==0.or.flag==2) call calc_gk(Natom, Mensemble, coord, simid, emomM, mstep, flag)
         if (do_sc_proj=='C') then
            if(mod(mstep-1,sc_sep)==0.or.flag_p==2) call calc_gk_proj(Natom, Mensemble, coord, simid, emomM, flag_p, NT, atype)
         endif
         if (do_sc_projch=='C') then
            if(mod(mstep-1,sc_sep)==0.or.flag_p==2) call calc_gk_projch(Natom, Mensemble,coord, simid, emomM, flag_p,Nchmax,achtype)
         endif

      else if(do_sc=='Q') then
         ! S(k,w)
         if(mod(mstep-1,sc_step)==1.and.sc_tidx<=sc_nstep) then
            sc_tidx=sc_tidx+1
            deltat_corr(sc_tidx) = delta_t ! Save current time step
            scstep_arr(sc_tidx) = sc_step  ! Save current sampling period
            call calc_gkt(Natom, Mensemble, coord, simid, emomM, flag)
            if (do_sc_proj=='Q') then
               call calc_gkt_proj(Natom, Mensemble, coord, simid, emomM, flag_p, atype, nt)
            end if
            if (do_sc_projch=='Q') then
               call calc_gkt_projch(Natom, Mensemble, coord, simid, emomM, flag_p, achtype, nchmax)
            endif
         else if (flag==2) then
            call calc_gkt(Natom, Mensemble, coord, simid, emomM, flag)
            if (do_sc_proj=='Q') then
               call calc_gkt_proj(Natom, Mensemble, coord, simid, emomM, flag_p, atype, nt)
            end if
            if (do_sc_projch=='Q') then
               call calc_gkt_projch(Natom, Mensemble, coord, simid, emomM, flag_p, achtype, nchmax)
            endif
         end if
      end if

      ! Sample S(r) directly for non periodic systems
      if (do_sr=='Y') then
         if (mod(mstep-1,sc_sep)==0) call calc_sr(Natom, Mensemble, coord, simid, emomM, flag)
      endif

      !
   end subroutine correlation_wrapper

   !--------------------------------------------------------------------------
   ! SUBROUTINE: correlation_init
   !> @brief
   !> Driver for correlation function calculations
   !---------------------------------------------------------------------------------
   subroutine correlation_init()

      implicit none

      !Spin-spin correlation
      sc_sep       = 100
      sc_mode      = 2
      sc_step      = 1
      sc_nstep     = 1
      sigma_q      = 0.0_dblprec
      sigma_w      = 0.0_dblprec
      LWfactor     = 0.0_dblprec
      LQfactor     = 0.0_dblprec
      do_sc        = 'N'
      do_sc_bimag  = 'N'
      do_sr        = 'N'
      qfile        = 'qfile'
      qpoints      = 'F'
      do_conv      = 'N'
      do_sc_proj   = 'N'
      do_sc_projch   = 'N'
      do_sc_local_axis   = 'N'
      sc_local_axis_mix = 0.0_dblprec
      do_sc_proj_axis   = 'N'
      do_sc_dosonly   = 'N'
      do_sc_complex   = 'N'
      do_qt_traj   = 'N'
      do_connected = 'N'
      sc_window_fun = 1

   end subroutine correlation_init

   !-----------------------------------------------------------------------------
   !> @brief Calculate \f$ \mathbf{S}\left(\mathbf{q}\right)\f$ for a cell filling
   !> mesh for transform to \f$ \mathbf{S}\left(\mathbf{r}\right) \f$
   !-----------------------------------------------------------------------------
   subroutine calc_gk(Natom, Mensemble, coord, simid, emomM, mstep, cgk_flag)

      use Constants, only : pi

      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(inout) :: cgk_flag  !< Allocate or deallocate (1/-1)
      !
      integer :: iq,r,l,i_stat,i_all
      character(len=30) :: filn
      real(dblprec) :: epowqr
      real(dblprec) :: qdr,nainv,k_min,qfac,wfac
      real(dblprec), dimension(3) :: s0,sp,cl
      real(dblprec), dimension(3) :: mavrg_vec
      complex(dblprec), dimension(3) :: cl_step

      if(.not.(do_sc=='C')) return

      if(cgk_flag==0) then
         ! First call, allocate and clear arrays
         allocate(corr_s(3,nq),stat=i_stat)
         call memocc(i_stat,product(shape(corr_s))*kind(corr_s),'corr_s','calc_gk')
         allocate(corr_k(3,nq),stat=i_stat)
         call memocc(i_stat,product(shape(corr_k))*kind(corr_k),'corr_k','calc_gk')
         allocate(corr_kt0(3,nq),stat=i_stat)
         call memocc(i_stat,product(shape(corr_kt0))*kind(corr_kt0),'corr_kt0','calc_gk')
         allocate(mavg_axis(3,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(mavg_axis))*kind(mavg_axis),'mavg_axis','calc_gk')
         mavg_axis=0.0d0
         corr_k=0.0d0
         cgk_flag=1
         sc_samp_done=0
         call find_rmid(r_mid,coord,Natom)
         return
      end if

      qfac=2.d0*pi
      wfac=1.d0
      nainv=1.0d0/Natom

      if(do_connected/='Y') then ! Keep the average magnetization for subtraction below if do_connected
         ! otherwise put it to zero (in order to remove a number of if-statements in-loop)
         mavrg_vec=0.0d0
      end if


      if (cgk_flag==1) then

         if(mod(mstep-1,sc_sep)==0) then
            ! Calculate s(k) for the current iteration and add to average of G(k)
            corr_s=0.0d0
            corr_kt0=0.0d0

            if (do_connected=='Y') call calc_mavrg_vec(Natom,Mensemble,emomM,mavrg_vec,mavg_axis)
            !$omp parallel do default(shared) private(r,iq,l,qdr,epowqr) schedule(static)
            do iq=1,nq
               do l=1,Mensemble
                  corr_s(:,iq)=0.0d0
                  do r=1,Natom
                     qdr=q(1,iq)*(coord(1,r)-r_mid(1))+q(2,iq)*(coord(2,r)-r_mid(2))+q(3,iq)*(coord(3,r)-r_mid(3))
                     epowqr=cos(qfac*qdr)*nainv*sc_window_fac(sc_window_fun,iq,nq)
                     corr_s(1,iq)=corr_s(1,iq)+epowqr*emomM(1,r,l)-epowqr*mavrg_vec(1)
                     corr_s(2,iq)=corr_s(2,iq)+epowqr*emomM(2,r,l)-epowqr*mavrg_vec(2)
                     corr_s(3,iq)=corr_s(3,iq)+epowqr*emomM(3,r,l)-epowqr*mavrg_vec(3)
                  end do
                  corr_k(1,iq)=corr_k(1,iq) + abs(corr_s(1,iq))**2
                  corr_k(2,iq)=corr_k(2,iq) + abs(corr_s(2,iq))**2
                  corr_k(3,iq)=corr_k(3,iq) + abs(corr_s(3,iq))**2
                  corr_kt0(1,iq)=corr_kt0(1,iq) + corr_s(1,iq)
                  corr_kt0(2,iq)=corr_kt0(2,iq) + corr_s(2,iq)
                  corr_kt0(3,iq)=corr_kt0(3,iq) + corr_s(3,iq)
               end do
            end do
            !$omp end parallel do

            if(do_qt_traj == 'Y') then
               ! Write S(q,t0)
               corr_kt0=corr_kt0/Mensemble
               write (filn,'(''sqt0.'',a8,''.out'')') simid
               open(ofileno,file=filn, position='append')
               do iq=1,nq
                  write(ofileno,'(i10, i10,3f10.4,4f18.8)') mstep, iq,(q(l,iq),l=1,3), &
                     (((corr_kt0(l,iq))),l=1,3), &
                     sqrt(corr_kt0(1,iq)**2+corr_kt0(2,iq)**2+corr_kt0(3,iq)**2)
               end do
               close(ofileno)
            end if

            sc_samp_done=sc_samp_done+1
         end if
      end if

      if (cgk_flag==2) then
         i_all=-product(shape(corr_kt0))*kind(corr_kt0)
         deallocate(corr_kt0,stat=i_stat)
         call memocc(i_stat,i_all,'corr_kt0','calc_gk')
         ! Finish sampling and write S(q)
         corr_k=corr_k/(sc_samp_done*Mensemble)

         !Write G(k)
         write (filn,'(''sq.'',a8,''.out'')') simid
         open(ofileno,file=filn,status='replace')
         do iq=1,nq
            write(ofileno,'(i10,3f10.4,5f18.8)') iq,(q(l,iq),l=1,3),(((corr_k(l,iq))),l=1,3), &
               sqrt(corr_k(1,iq)**2+corr_k(2,iq)**2+corr_k(3,iq)**2),corr_k(1,iq)+corr_k(2,iq)+corr_k(3,iq)
         end do
         close(ofileno)


         ! Calculate the correlation length following the Katzgraber recipe
         s0=corr_k(:,qmin(1))
         sp=corr_k(:,qmin(2))
         k_min=sqrt(q(1,qmin(2))**2+q(2,qmin(2))**2+q(3,qmin(2))**2)
         cl_step = (s0/sp-1.0d0)/4.0d0/sin((qfac*k_min/2.0d0))**2
         cl=real(sqrt(cl_step))
         write(*,'(2x,a20,2x,f11.5,2x,f11.5,2x,f11.5)') &
            'Correlation lengths:',cl(1),cl(2),cl(3)

         ! Transform G(k) to G(r)
         i_all=-product(shape(corr_s))*kind(corr_s)
         deallocate(corr_s,stat=i_stat)
         call memocc(i_stat,i_all,'corr_s','calc_gk')
         allocate(corr_s(3,natom),stat=i_stat)
         call memocc(i_stat,product(shape(corr_s))*kind(corr_s),'corr_s','calc_gk')
         corr_s=0.0d0

         !$omp parallel do default(shared) private(r,iq,qdr,epowqr) schedule(static)
         do r=1,Natom
            do iq=1,nq
               qdr=q(1,iq)*(coord(1,r)-r_mid(1))+q(2,iq)*(coord(2,r)-r_mid(2))+q(3,iq)*(coord(3,r)-r_mid(3))
               epowqr=cos(-qfac*qdr)*sc_window_fac(sc_window_fun,iq,nq)
               corr_s(1,r)=corr_s(1,r)+epowqr*corr_k(1,iq)
               corr_s(2,r)=corr_s(2,r)+epowqr*corr_k(2,iq)
               corr_s(3,r)=corr_s(3,r)+epowqr*corr_k(3,iq)
            end do
         end do
         !$omp end parallel do

         ! Write G(r)
         write (filn,'(''sr.'',a8,''.out'')') simid
         open(ofileno,file=filn,status='replace')
         do r=1,Natom
            write(ofileno,'(i10,3f10.4,5f18.8)') r,(coord(l,r)-r_mid(l),l=1,3),(((corr_s(l,r))),l=1,3),&
               sqrt(corr_s(1,r)**2+corr_s(2,r)**2+corr_s(3,r)**2),corr_s(1,r)+corr_s(2,r)+corr_s(3,r)
         end do
         close(ofileno)

         ! Write G(|r|)
         write (filn,'(''sra.'',a8,''.out'')') simid
         open(ofileno,file=filn,status='replace')
         do r=1,Natom
            write(ofileno,'(7f18.8)') sqrt( (coord(1,r)-r_mid(1))**2+(coord(2,r)-r_mid(2))**2+(coord(3,r)-r_mid(3))**2),&
               (((corr_s(l,r))),l=1,3),&
               sqrt(corr_s(1,r)**2+corr_s(2,r)**2+corr_s(3,r)**2),corr_s(1,r)+corr_s(2,r)+corr_s(3,r)
         end do
         close(ofileno)

         ! Deallocate arrays
         i_all=-product(shape(corr_k))*kind(corr_k)
         deallocate(corr_k,stat=i_stat)
         call memocc(i_stat,i_all,'corr_k','calc_gk')

         i_all=-product(shape(corr_s))*kind(corr_s)
         deallocate(corr_s,stat=i_stat)
         call memocc(i_stat,i_all,'corr_s','calc_gk')
      end if
      return
   end subroutine calc_gk


   !-----------------------------------------------------------------------------
   ! SUBROUTINE: deallocate_q
   !> @brief Deallocate arrays for the q-points
   !-----------------------------------------------------------------------------
   subroutine deallocate_q()
      !
      implicit none
      !
      !
      integer :: i_all,i_stat
      !
      if(allocated(q)) then
         i_all=-product(shape(q))*kind(q)
         deallocate(q,stat=i_stat)
         call memocc(i_stat,i_all,'q','deallocate_q')
      end if
      if(do_sc=='C'.or.do_sc=='D') then

         if (allocated(q_weight)) then
            i_all=-product(shape(q_weight))*kind(q_weight)
            deallocate(q_weight,stat=i_stat)
            call memocc(i_stat,i_all,'q_weight','deallocate_q')
         endif

         if (do_sc=='D') then
            i_all=-product(shape(w))*kind(w)
            deallocate(w,stat=i_stat)
            call memocc(i_stat,i_all,'w','deallocate_q')
         end if
      end if
   end subroutine deallocate_q

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: read_q
   !> @brief Read q-points from file for \f$ \mathbf{S}\left(\mathbf{r},t\right) \rightarrow \mathbf{S}\left(\mathbf{q},t\right)\f$ transform
   !-----------------------------------------------------------------------------
   subroutine read_q(C1,C2,C3)
      !
      use Sorting, only : qsort
      implicit none
      !
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      integer :: iq
      integer :: i_stat,i_all
      integer, dimension(:), allocatable :: ia
      real(dblprec), dimension(3) :: b1,r1
      real(dblprec), dimension(3) :: b2,r2
      real(dblprec), dimension(3) :: b3,r3
      real(dblprec), dimension(3) :: q_tmp
      real(dblprec), dimension(:), allocatable :: dq
      real(dblprec) :: c1r1, c2r2, c3r3

      open(ifileno, file=adjustl(qfile))
      read (ifileno,*) nq
      allocate(q(3,nq),stat=i_stat)
      call memocc(i_stat,product(shape(q))*kind(q),'q','read_q')
      allocate(dq(nq),stat=i_stat)
      call memocc(i_stat,product(shape(dq))*kind(dq),'dq','read_q')
      allocate(ia(nq),stat=i_stat)
      call memocc(i_stat,product(shape(ia))*kind(ia),'ia','read_q')
      if (qpoints=='F') then
         do iq=1,nq
            read (ifileno,*) q(1,iq), q(2,iq), q(3,iq)
         enddo
      else if (qpoints=='I') then
         allocate(q_weight(nq),stat=i_stat)
         call memocc(i_stat,product(shape(q_weight))*kind(q_weight),'q_weight','read_q')
         do iq=1,nq
            read(ifileno,*) q(1,iq), q(2,iq), q(3,iq), q_weight(iq)
         enddo
      else  ! qpoints==D
         ! Calculate reciprocal lattice vectors
         ! r1 = C2xC3
         r1(1)=C2(2)*C3(3)-C2(3)*C3(2)
         r1(2)=C2(3)*C3(1)-C2(1)*C3(3)
         r1(3)=C2(1)*C3(2)-C2(2)*C3(1)
         ! r2 = C3xC1
         r2(1)=C3(2)*C1(3)-C3(3)*C1(2)
         r2(2)=C3(3)*C1(1)-C3(1)*C1(3)
         r2(3)=C3(1)*C1(2)-C3(2)*C1(1)
         ! r3 = C1xC2
         r3(1)=C1(2)*C2(3)-C1(3)*C2(2)
         r3(2)=C1(3)*C2(1)-C1(1)*C2(3)
         r3(3)=C1(1)*C2(2)-C1(2)*C2(1)
         ! cell volume C1*(C2xC3)
         c1r1=C1(1)*r1(1)+C1(2)*r1(2)+C1(3)*r1(3)
         c2r2=C2(1)*r2(1)+C2(2)*r2(2)+C2(3)*r2(3)
         c3r3=C3(1)*r3(1)+C3(2)*r3(2)+C3(3)*r3(3)
         ! b1=(2pi)*r1/(C1*r1)
         b1(1)=r1(1)/c1r1
         b1(2)=r1(2)/c1r1
         b1(3)=r1(3)/c1r1
         ! b2=(2pi)*r2/(C1*r1)
         b2(1)=r2(1)/c2r2
         b2(2)=r2(2)/c2r2
         b2(3)=r2(3)/c2r2
         ! b3=(2pi)*r3/(C1*r1)
         b3(1)=r3(1)/c3r3
         b3(2)=r3(2)/c3r3
         b3(3)=r3(3)/c3r3

         if (qpoints=='B') then
            allocate(q_weight(nq),stat=i_stat)
            call memocc(i_stat,product(shape(q_weight))*kind(q_weight),'q_weight','read_q')
            do iq=1,nq
               read (ifileno,*) q_tmp(1), q_tmp(2) , q_tmp(3), q_weight(iq)
               q(1,iq)=q_tmp(1)*b1(1)+q_tmp(2)*b2(1)+q_tmp(3)*b3(1)
               q(2,iq)=q_tmp(1)*b1(2)+q_tmp(2)*b2(2)+q_tmp(3)*b3(2)
               q(3,iq)=q_tmp(1)*b1(3)+q_tmp(2)*b2(3)+q_tmp(3)*b3(3)
            enddo

         else
            do iq=1,nq
               read (ifileno,*) q_tmp(1), q_tmp(2) , q_tmp(3)
               q(1,iq)=q_tmp(1)*b1(1)+q_tmp(2)*b2(1)+q_tmp(3)*b3(1)
               q(2,iq)=q_tmp(1)*b1(2)+q_tmp(2)*b2(2)+q_tmp(3)*b3(2)
               q(3,iq)=q_tmp(1)*b1(3)+q_tmp(2)*b2(3)+q_tmp(3)*b3(3)
            enddo
         endif
      endif

      do iq=1,nq
         dq(iq)=(q(1,iq)**2+q(2,iq)**2+q(3,iq)**2)
      enddo
      qmin(1)=minloc(dq,1)
      qmin(2)=minloc(dq,1,abs(dq-dq(qmin(1)))> 0.d0)

      close(ifileno)
      open(ofileno,file='qpoints.out',status='replace')
      do iq=1,nq
         write(ofileno,'(i6,3f14.6)') iq,q(1,iq),q(2,iq),q(3,iq)
      end do
      close(ofileno)

      i_all=-product(shape(dq))*kind(dq)
      deallocate(dq,stat=i_stat)
      call memocc(i_stat,i_all,'dq','read_q')
      i_all=-product(shape(ia))*kind(ia)
      deallocate(ia,stat=i_stat)
      call memocc(i_stat,i_all,'ia','read_q')

   end subroutine read_q


   !-----------------------------------------------------------------------------
   ! SUBROUTINE: set_w
   !> @brief Calculate suitable values of frequencies for \f$ \mathbf{S}\left(\mathbf{q},t\right) \rightarrow \mathbf{S}\left(\mathbf{q},\omega\right)\f$ transform
   !> @todo Change setup to smarter algorithm with respect to the exchange strength of the
   !> system
   !-----------------------------------------------------------------------------
   subroutine set_w(delta_t)
      !
      use Constants, only : pi, hbar_mev
      !
      implicit none
      !
      real(dblprec), intent(in) :: delta_t !< Time step
      !
      integer :: i_stat,j
      real(dblprec) :: dt !< Time step
      real(dblprec) :: dww
      real(dblprec) :: emin,emax
      !
      dt=sc_step*delta_t
      nw = sc_nstep
      allocate(w(nw),stat=i_stat)
      call memocc(i_stat,product(shape(w))*kind(w),'w','set_w')
      dww = 2*pi/(nw*dt)

      ! Previous convention was j*ddw
      do j=1,nw
         w(j)=(j-1)*dww
      end do

      emin=hbar_mev*(w(2)-w(1))
      emax=hbar_mev*w(nw)
      write(*,'(1x,a,f6.3,a,f6.1,a)') 'Spin wave sampling between ',emin,' meV and ',emax,' meV.'
   end subroutine set_w

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_corr_w
   !> @brief Calculate correction to the sampling frequencies, w
   !> Perform correction to all frequencies >= sc_tidx
   !-----------------------------------------------------------------------------
   subroutine calc_corr_w(deltat_upd)

      implicit none

      real(dblprec), intent(in) :: deltat_upd                       !< Updated time step

      integer :: j
      real(dblprec) :: dt, pi
      real(dblprec) :: dww

      pi = 4.*atan(1.0d0)

      dt = sc_step*deltat_upd ! Time interval between new samples
      nw = sc_nstep
      dww = 2*pi/(nw*dt)

      ! Perform correction of the current/future sampling frequencies
      do j=sc_tidx+1,nw ! The addition of 1 is due to function calling at the end of the previous sampling period
         w(j)=j*dww
      end do

   end subroutine calc_corr_w

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: setup_qcoord
   !> @brief Automatic setup of q-point mesh for calculating \f$\mathbf{S}\left(\mathbf{q}\right)\rightarrow\mathbf{S}\left(\mathbf{r}\right)\f$
   !-----------------------------------------------------------------------------
   subroutine setup_qcoord(N1,N2,N3,C1,C2,C3)

      use Sorting, only : qsort
      !
      !
      implicit none
      !
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      !or (P)lane 2D extended cell
      !
      integer :: iq,xq,yq,zq
      integer :: i_stat, i_all
      integer,dimension(:), allocatable :: ia
      real(dblprec), dimension(3) :: b1,r1
      real(dblprec), dimension(3) :: b2,r2
      real(dblprec), dimension(3) :: b3,r3
      real(dblprec), dimension(:), allocatable :: dq
      real(dblprec) :: c1r1, c2r2, c3r3

      ! Calculate reciprocal lattice vectors
      ! r1 = C2xC3
      r1(1)=C2(2)*C3(3)-C2(3)*C3(2)
      r1(2)=C2(3)*C3(1)-C2(1)*C3(3)
      r1(3)=C2(1)*C3(2)-C2(2)*C3(1)
      ! r2 = C3xC1
      r2(1)=C3(2)*C1(3)-C3(3)*C1(2)
      r2(2)=C3(3)*C1(1)-C3(1)*C1(3)
      r2(3)=C3(1)*C1(2)-C3(2)*C1(1)
      ! r3 = C1xC2
      r3(1)=C1(2)*C2(3)-C1(3)*C2(2)
      r3(2)=C1(3)*C2(1)-C1(1)*C2(3)
      r3(3)=C1(1)*C2(2)-C1(2)*C2(1)
      ! cell volume C1*(C2xC3)
      c1r1=C1(1)*r1(1)+C1(2)*r1(2)+C1(3)*r1(3)
      c2r2=C2(1)*r2(1)+C2(2)*r2(2)+C2(3)*r2(3)
      c3r3=C3(1)*r3(1)+C3(2)*r3(2)+C3(3)*r3(3)
      ! b1=(2pi)*r1/(C1*r1)
      b1(1)=r1(1)/c1r1
      b1(2)=r1(2)/c1r1
      b1(3)=r1(3)/c1r1
      ! b2=(2pi)*r2/(C1*r1)
      b2(1)=r2(1)/c2r2
      b2(2)=r2(2)/c2r2
      b2(3)=r2(3)/c2r2
      ! b3=(2pi)*r3/(C1*r1)
      b3(1)=r3(1)/c3r3
      b3(2)=r3(2)/c3r3
      b3(3)=r3(3)/c3r3
      !
      if(allocated(q)) then
         i_all=-product(shape(q))*kind(q)
         deallocate(q,stat=i_stat)
         call memocc(i_stat,i_all,'q','setup_qcoord')
      end if

      if(qpoints=='C') then
         nq=0
         do zq=-(N3-1)/2,(N3)/2
            do yq=-(N2-1)/2,(N2)/2
               do xq=-(N1-1)/2,(N1)/2
                  nq=nq+1
               end do
            end do
         end do
         allocate(ia(nq),stat=i_stat)
         call memocc(i_stat,product(shape(ia))*kind(ia),'ia','setup_qcoord')
         allocate(q(3,nq),stat=i_stat)
         call memocc(i_stat,product(shape(q))*kind(q),'q','setup_qcoord')
         allocate(dq(nq),stat=i_stat)
         call memocc(i_stat,product(shape(dq))*kind(dq),'dq','setup_qcoord')
         iq=0

         do zq=-(N3-1)/2,(N3)/2
            do yq=-(N2-1)/2,(N2)/2
               do xq=-(N1-1)/2,(N1)/2
                  iq=iq+1
                  q(:,iq)=xq/(1.0d0*N1)*b1+yq/(1.0d0*N2)*b2+zq/(1.0d0*N3)*b3
               end do
            end do
         end do
         do iq=1,nq
            dq(iq)=q(1,iq)**2+q(2,iq)**2+q(3,iq)**2
         enddo
         qmin(1)=minloc(dq,1)
         qmin(2)=minloc(dq,1,abs(dq-dq(qmin(1)))> 0.d0)

      else if(qpoints=='X') then    ! Like 'C' but for double size of q-grid
         nq=0
         do zq=-(N3),(N3)
            do yq=-(N2),(N2)
               do xq=-(N1),(N1)
                  nq=nq+1
               end do
            end do
         end do

         allocate(ia(nq),stat=i_stat)
         call memocc(i_stat,product(shape(ia))*kind(ia),'ia','setup_qcoord')
         allocate(q(3,nq),stat=i_stat)
         call memocc(i_stat,product(shape(q))*kind(q),'q','setup_qcoord')
         allocate(dq(nq),stat=i_stat)
         call memocc(i_stat,product(shape(dq))*kind(dq),'dq','setup_qcoord')
         iq=0

         do zq=-(N3),(N3)
            do yq=-(N2),(N2)
               do xq=-(N1),(N1)
                  iq=iq+1
                  q(:,iq)=xq/(1.0d0*N1)*b1+yq/(1.0d0*N2)*b2+zq/(1.0d0*N3)*b3
                  dq(iq)=q(1,iq)**2+q(2,iq)**2+q(3,iq)**2
               end do
            end do
         end do
         qmin(1)=minloc(dq,1)
         qmin(2)=minloc(dq,1,abs(dq-dq(qmin(1)))> 0.d0)

      else if(qpoints=='H') then    ! Like 'C' but for double size of q-grid
         nq=0
         do zq=-(N3),(N3)
            do yq=-(N2),(N2)
               do xq=-(N1),(N1)
                  nq=nq+1
               end do
            end do
         end do

         allocate(ia(nq),stat=i_stat)
         call memocc(i_stat,product(shape(ia))*kind(ia),'ia','setup_qcoord')
         allocate(q(3,nq),stat=i_stat)
         call memocc(i_stat,product(shape(q))*kind(q),'q','setup_qcoord')
         allocate(dq(nq),stat=i_stat)
         call memocc(i_stat,product(shape(dq))*kind(dq),'dq','setup_qcoord')
         iq=0

         do zq=-(N3),(N3)
            do yq=-(N2),(N2)
               do xq=-(N1),(N1)
                  iq=iq+1
                  q(:,iq)=xq/(2.0d0*N1)*b1+yq/(2.0d0*N2)*b2+zq/(2.0d0*N3)*b3
                  dq(iq)=q(1,iq)**2+q(2,iq)**2+q(3,iq)**2
               end do
            end do
         end do
         qmin(1)=minloc(dq,1)
         qmin(2)=minloc(dq,1,abs(dq-dq(qmin(1)))> 0.d0)

      else if(qpoints=='P') then
         nq=(4*N1+1)*(4*N3+1)
         allocate(ia(nq),stat=i_stat)
         call memocc(i_stat,product(shape(ia))*kind(ia),'ia','setup_qcoord')
         allocate(q(3,nq),stat=i_stat)
         call memocc(i_stat,product(shape(q))*kind(q),'q','setup_qcoord')
         allocate(dq(nq),stat=i_stat)
         call memocc(i_stat,product(shape(dq))*kind(dq),'dq','setup_qcoord')
         iq=0
         do zq=-2*N3,2*N3
            yq=0
            do xq=-2*N3,2*N3
               iq=iq+1
               q(:,iq)=xq/(1.0d0*N1)*b1+yq/(1.0d0*N2)*b2+zq/(1.0d0*N3)*b3
               dq(iq)=q(1,iq)**2+q(2,iq)**2+q(3,iq)**2
            end do
         end do
         qmin(1)=minloc(dq,1)
         qmin(2)=minloc(dq,1,abs(dq-dq(qmin(1)))> 0.d0)

      else if(qpoints=='A') then
         ! Currently G-x-y-G-z-y
         ! where x,y,z are the directions of the reciprocal lattice vectors
         nq=N1+(N2-1)+(N2-1)+(N3-1)+(N2-1)+(N3-1)
         allocate(q(3,nq),stat=i_stat)
         call memocc(i_stat,product(shape(q))*kind(q),'q','setup_qcoord')
         allocate(dq(nq),stat=i_stat)
         call memocc(i_stat,product(shape(dq))*kind(dq),'dq','setup_qcoord')
         q=0.0d0
         iq=0
         xq=0
         yq=0
         zq=0
         qmin(1)=1
         ! G->x
         do xq=0,N1-1
            iq=iq+1
            q(:,iq)=xq/(1.0d0*N1)*b1+yq/(1.0d0*N2)*b2+zq/(1.0d0*N3)*b3
         end do
         ! x-> y
         do yq=1,N2-1
            !       xq=(N1-1)-yq*(N1-1)/(N2-1)
            iq=iq+1
            q(:,iq)=xq/(1.0d0*N1)*b1+yq/(1.0d0*N2)*b2+zq/(1.0d0*N3)*b3
         end do
         ! xy->G
         do yq=N2-2,0,-1
            xq=yq/(N2-1)
            iq=iq+1
            q(:,iq)=xq/(1.0d0*N1)*b1+yq/(1.0d0*N2)*b2+zq/(1.0d0*N3)*b3
         end do
         ! G->y
         do yq=1,N2-1
            iq=iq+1
            q(:,iq)=xq/(1.0d0*N1)*b1+yq/(1.0d0*N2)*b2+zq/(1.0d0*N3)*b3
         end do
         ! y->z
         do zq=0,N3-1
            yq=N2-1-zq*(N2-1)/(N3-1)
            iq=iq+1
            q(:,iq)=xq/(1.0d0*N1)*b1+yq/(1.0d0*N2)*b2+zq/(1.0d0*N3)*b3
         end do
         ! z-G
         do zq=N3-2,1,-1
            iq=iq+1
            q(:,iq)=xq/(1.0d0*N1)*b1+yq/(1.0d0*N2)*b2+zq/(1.0d0*N3)*b3
         end do
      end if
      open(ofileno,file='qpoints.out',status='replace')
      do iq=1,nq
         write(ofileno,'(i6,3f14.6)') iq,q(1,iq),q(2,iq),q(3,iq)
      end do
      close(ofileno)
      if (qpoints=='A') then
         do iq=1,nq
            dq(iq)=q(1,iq)**2+q(2,iq)**2+q(3,iq)**2
         enddo
         qmin(1)=minloc(dq,1)
         qmin(2)=minloc(dq,1,abs(dq-dq(qmin(1)))> 0.d0)
      endif

      if (qpoints=='P' .or. qpoints=='C') then
         i_all=-product(shape(ia))*kind(ia)
         deallocate(ia,stat=i_stat)
         call memocc(i_stat,i_all,'ia','setup_qcoord')

         i_all=-product(shape(dq))*kind(dq)
         deallocate(dq,stat=i_stat)
         call memocc(i_stat,i_all,'dq','setup_qcoord')
      end if
   end subroutine setup_qcoord


   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_gk_proj
   !> Calculate \f$ \mathbf{S}\left(\mathbf{q}\right)\f$ for a cell filling mesh for transform to \f$\mathbf{S}\left(\mathbf{r}\right)\f$ (sublattice projection)
   !> @todo Add non-diagonal components
   !-----------------------------------------------------------------------------
   subroutine calc_gk_proj(Natom, Mensemble, coord, simid, emomM, flag, NT, atype)
      !
      use Constants, only : pi
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
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
         allocate(corr_s_proj(3,nq,nt),stat=i_stat)
         call memocc(i_stat,product(shape(corr_s_proj))*kind(corr_s_proj),'corr_s_proj','calc_gk_proj')
         allocate(corr_k_proj(3,nq,nt,nt),stat=i_stat)
         call memocc(i_stat,product(shape(corr_k_proj))*kind(corr_k_proj),'corr_k_proj','calc_gk_proj')
         corr_k_proj=0.0d0
         corr_s_proj=0.0d0
         flag=1
         sc_samp_done_proj=0
      end if

      qfac=2.d0*pi
      wfac=1.d0
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

      nainv(:)=1.0d0/ncounter(:)
      if (flag==1) then

         ! Calculate s(k) for the current iteration and add to average of G(k)
         corr_s_proj=0.0d0

         !$omp parallel do default(shared) private(r,iq,l,qdr,epowqr,k) schedule(static)
         do iq=1,nq
            do l=1,Mensemble
               corr_s_proj(:,iq,:)=0.0d0
               do r=1,Natom
                  k=atype(r)
                  qdr=q(1,iq)*coord(1,r)+q(2,iq)*coord(2,r)+q(3,iq)*coord(3,r)
                  epowqr=cos(qfac*qdr)*nainv(k)
                  corr_s_proj(1,iq,k)=corr_s_proj(1,iq,k)+epowqr*emomM(1,r,l)
                  corr_s_proj(2,iq,k)=corr_s_proj(2,iq,k)+epowqr*emomM(2,r,l)
                  corr_s_proj(3,iq,k)=corr_s_proj(3,iq,k)+epowqr*emomM(3,r,l)
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

         sc_samp_done_proj=sc_samp_done_proj+1
      end if

      if (flag==2) then

         ! Finish sampling and write S(q)
         corr_k_proj=corr_k_proj/(sc_samp_done_proj*Mensemble)
         !Write G(k)
         do k=1,nt
            do m=k,nt
               if(k>10) then
                  write (filn,'(''projsq.'',i2,''.'',i2,''.'',a8,''.out'')') k,m,simid
               else
                  write (filn,'(''projsq.'',i1,''.'',i1,''.'',a8,''.out'')') k,m,simid
               end if
               open(ofileno, file=filn)
               do iq=1,nq
                  write(ofileno,'(i5,3f10.4,4f18.8)') iq,(q(l,iq),l=1,3),(((corr_k_proj(l,iq,k,m))),l=1,3), &
                     sqrt(corr_k_proj(1,iq,k,m)**2+corr_k_proj(2,iq,k,m)**2+corr_k_proj(3,iq,k,m)**2)
               enddo
               close(ofileno)
            end do
         end do

         ! Transform G(k) to G(r)
         i_all=-product(shape(corr_s_proj))*kind(corr_s_proj)
         deallocate(corr_s_proj,stat=i_stat)
         call memocc(i_stat,i_all,'corr_s_proj','calc_gk_proj')
         allocate(corr_ss_proj(3,natom,nt,nt),stat=i_stat)
         call memocc(i_stat,product(shape(corr_ss_proj))*kind(corr_ss_proj),'corr_ss_proj','calc_gk_proj')
         corr_ss_proj=0.0d0

         !$omp parallel do default(shared) private(r,iq,qdr,epowqr,k,ic,m) schedule(static)
         do r=1,Natom
            m=atype(r)
            do k=1,nt
               do iq=1,nq
                  ic=indxcoord(1,k)
                  qdr=q(1,iq)*(coord(1,r)-coord(1,ic))+q(2,iq)*(coord(2,r)-coord(2,ic))+ &
                     q(3,iq)*(coord(3,r)-coord(3,ic))
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
                  write (filn,'(''projsr.'',i2,''.'',i2,''.'',a8,''.out'')') k,m,simid
               else
                  write (filn,'(''projsr.'',i1,''.'',i1,''.'',a8,''.out'')') k,m,simid
               end if
               open(ofileno, file=filn)
               do r=1,Natom
                  if (m==atype(r)) then
                     write(ofileno,'(i5,3f10.4,4f18.8)') r,(coord(1,r)-coord(1,ic)),(coord(2,r)-coord(2,ic)), &
                        (coord(3,r)-coord(3,ic)),(((corr_ss_proj(l,r,k,m))),l=1,3),&
                        sqrt(corr_ss_proj(1,r,k,m)**2+corr_ss_proj(2,r,k,m)**2+corr_ss_proj(3,r,k,m)**2)
                  end if
               enddo
               close(ofileno)
            end do
         end do

         ! Write G(|r|)
         do k=1,nt
            ic=indxcoord(1,k)
            do m=k,nt
               if(k>10) then
                  write (filn,'(''projsra.'',i2,''.'',i2,''.'',a8,''.out'')') k,m,simid
               else
                  write (filn,'(''projsra.'',i1,''.'',i1,''.'',a8,''.out'')') k,m,simid
               end if
               open(ofileno, file=filn)
               do r=1,Natom
                  if (m==atype(r)) then
                     write(ofileno,'(7f18.8)') sqrt( (coord(1,r)-coord(1,ic))**2+(coord(2,r)-coord(2,ic))**2+ &
                        (coord(3,r)-coord(3,ic))**2),&
                        (((corr_ss_proj(l,r,k,m))),l=1,3),&
                        sqrt(corr_ss_proj(1,r,k,m)**2+corr_ss_proj(2,r,k,m)**2+corr_ss_proj(3,r,k,m)**2)
                  end if
               enddo
               close(ofileno)
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
   end subroutine calc_gk_proj

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_gk_projch
   !> @brief Calculate \f$ \mathbf{S}\left(\mathbf{q}\right)\f$ for a cell filling mesh for transform to \f$\mathbf{S}\left(\mathbf{r}\right)\f$ (sublattice projection)
   !-----------------------------------------------------------------------------
   subroutine calc_gk_projch(Natom, Mensemble, coord, simid, emomM, flag, Nchmax, achtype)
      !
      use Constants, only : pi
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      integer, intent(inout) :: flag  !< Allocate or deallocate (1/-1)
      integer, intent(in) :: Nchmax  !< Number of chemical types
      integer, dimension(Natom), intent(in) :: achtype !< Type of atom
      !
      integer :: iq,r,l,i_stat,i_all,k,m,ic
      integer,dimension(nchmax) :: ncounter
      integer,dimension(natom,nchmax) :: indxcoord
      character(len=40) :: filn
      real(dblprec) :: epowqr!,s0,sp
      real(dblprec) :: qdr
      real(dblprec), dimension(nchmax) :: nainv

      if(flag==0) then
         ! First call, allocate and clear arrays
         allocate(corr_s_proj(3,nq,nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(corr_s_proj))*kind(corr_s_proj),'corr_s_proj','calc_gk_projch')
         allocate(corr_k_proj(3,nq,nchmax,nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(corr_k_proj))*kind(corr_k_proj),'corr_k_proj','calc_gk_projch')
         corr_k_proj=0.0d0
         corr_s_proj=0.0d0
         flag=1
         sc_samp_done_proj=0
      end if

      qfac=2.d0*pi
      wfac=1.d0
      ncounter=0
      indxcoord=0
      do k=1,nchmax
         do r=1,Natom
            if (k==achtype(r)) then
               ncounter(k)=ncounter(k)+1
               indxcoord(ncounter(k),k)=r
            end if
         enddo
      enddo

      nainv(:)=1.0d0/ncounter(:)
      if (flag==1) then

         ! Calculate s(k) for the current iteration and add to average of G(k)
         corr_s_proj=0.0d0

         !$omp parallel do default(shared) private(r,iq,l,qdr,epowqr,k) schedule(static)
         do iq=1,nq
            do l=1,Mensemble
               corr_s_proj(:,iq,:)=0.0d0
               do r=1,Natom
                  k=achtype(r)
                  qdr=q(1,iq)*(coord(1,r)-coord(1,qmin(1)))+q(2,iq)*(coord(2,r)-coord(2,qmin(1)))+q(3,iq)*(coord(3,r)-coord(3,qmin(1)))
                  epowqr=cos(qfac*qdr)*nainv(k)
                  corr_s_proj(1,iq,k)=corr_s_proj(1,iq,k)+epowqr*emomM(1,r,l)
                  corr_s_proj(2,iq,k)=corr_s_proj(2,iq,k)+epowqr*emomM(2,r,l)
                  corr_s_proj(3,iq,k)=corr_s_proj(3,iq,k)+epowqr*emomM(3,r,l)
               end do
               do m=1,nchmax
                  do k=1,nchmax
                     corr_k_proj(1,iq,k,m)=corr_k_proj(1,iq,k,m) + (corr_s_proj(1,iq,k)*corr_s_proj(1,iq,m))
                     corr_k_proj(2,iq,k,m)=corr_k_proj(2,iq,k,m) + (corr_s_proj(2,iq,k)*corr_s_proj(2,iq,m))
                     corr_k_proj(3,iq,k,m)=corr_k_proj(3,iq,k,m) + (corr_s_proj(3,iq,k)*corr_s_proj(3,iq,m))
                  enddo
               enddo
            end do
         end do
         !$omp end parallel do

         sc_samp_done_proj=sc_samp_done_proj+1
      end if

      if (flag==2) then

         ! Finish sampling and write S(q)
         corr_k_proj=corr_k_proj/(sc_samp_done_proj*Mensemble)
         !Write G(k)
         do k=1,nchmax
            do m=k,nchmax
               if(k>10) then
                  write (filn,'(''projchsq.'',i2,''.'',i2,''.'',a8,''.out'')') k,m,simid
               else
                  write (filn,'(''projchsq.'',i1,''.'',i1,''.'',a8,''.out'')') k,m,simid
               end if
               open(ofileno, file=filn)
               do iq=1,nq
                  write(ofileno,'(i5,3f10.4,5f18.8)') iq,(q(l,iq),l=1,3),(((corr_k_proj(l,iq,k,m))),l=1,3), &
                     sqrt(corr_k_proj(1,iq,k,m)**2+corr_k_proj(2,iq,k,m)**2+corr_k_proj(3,iq,k,m)**2), &
                     corr_k_proj(1,iq,k,m)+corr_k_proj(2,iq,k,m)+corr_k_proj(3,iq,k,m)
               enddo
               close(ofileno)
            end do
         end do

         ! Transform G(k) to G(r)
         i_all=-product(shape(corr_s_proj))*kind(corr_s_proj)
         deallocate(corr_s_proj,stat=i_stat)
         call memocc(i_stat,i_all,'corr_s_proj','calc_gk_proj')
         allocate(corr_ss_proj(3,natom,nchmax,nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(corr_ss_proj))*kind(corr_ss_proj),'corr_ss_proj','calc_gk_projch')
         corr_ss_proj=0.0d0

         !$omp parallel do default(shared) private(r,iq,qdr,epowqr,k,ic,m) schedule(static)
         do r=1,Natom
            m=achtype(r)
            do k=1,nchmax
               do iq=1,nq
                  ic=indxcoord(1,k)
                  qdr= q(1,iq)*( coord(1,r)-coord(1,ic)-coord(1,qmin(1)) )+ &
                     q(2,iq)*( coord(2,r)-coord(2,ic)-coord(2,qmin(1)) )+ &
                     q(3,iq)*( coord(3,r)-coord(3,ic)-coord(3,qmin(1)) )
                  epowqr=cos(-qfac*qdr)
                  corr_ss_proj(1,r,k,m)=corr_ss_proj(1,r,k,m)+epowqr*corr_k_proj(1,iq,k,m)
                  corr_ss_proj(2,r,k,m)=corr_ss_proj(2,r,k,m)+epowqr*corr_k_proj(2,iq,k,m)
                  corr_ss_proj(3,r,k,m)=corr_ss_proj(3,r,k,m)+epowqr*corr_k_proj(3,iq,k,m)
               end do
            end do
         end do
         !$omp end parallel do

         ! Write G(r)
         do k=1,nchmax
            ic=indxcoord(1,k)
            do m=k,nchmax
               if(k>10) then
                  write (filn,'(''projchsr.'',i2,''.'',i2,''.'',a8,''.out'')') k,m,simid
               else
                  write (filn,'(''projchsr.'',i1,''.'',i1,''.'',a8,''.out'')') k,m,simid
               end if
               open(ofileno, file=filn)
               do r=1,Natom
                  if (m==achtype(r)) then
                     write(ofileno,'(i5,3f10.4,5f18.8)') r,(coord(1,r)-coord(1,ic)-coord(1,qmin(1))), &
                        (coord(2,r)-coord(2,ic)-coord(2,qmin(1))), &
                        (coord(3,r)-coord(3,ic)-coord(3,qmin(1))),(((corr_ss_proj(l,r,k,m))),l=1,3),&
                        sqrt(corr_ss_proj(1,r,k,m)**2+corr_ss_proj(2,r,k,m)**2+corr_ss_proj(3,r,k,m)**2), &
                        corr_ss_proj(1,r,k,m)+corr_ss_proj(2,r,k,m)+corr_ss_proj(3,r,k,m)
                  end if
               enddo
               close(ofileno)
            end do
         end do

         ! Write G(|r|)
         do k=1,nchmax
            ic=indxcoord(1,k)
            do m=k,nchmax
               if(k>10) then
                  write (filn,'(''projchsra.'',i2,''.'',i2,''.'',a8,''.out'')') k,m,simid
               else
                  write (filn,'(''projchsra.'',i1,''.'',i1,''.'',a8,''.out'')') k,m,simid
               end if
               open(ofileno, file=filn)
               do r=1,Natom
                  if (m==achtype(r)) then
                     write(ofileno,'(8f18.8)') sqrt( (coord(1,r)-coord(1,ic)-coord(1,qmin(1)))**2+ &
                        (coord(2,r)-coord(2,ic)-coord(2,qmin(1)))**2+ &
                        (coord(3,r)-coord(3,ic)-coord(3,qmin(1)))**2),&
                        (((corr_ss_proj(l,r,k,m))),l=1,3),&
                        sqrt(corr_ss_proj(1,r,k,m)**2+corr_ss_proj(2,r,k,m)**2+corr_ss_proj(3,r,k,m)**2),&
                        corr_ss_proj(1,r,k,m)+corr_ss_proj(2,r,k,m)+corr_ss_proj(3,r,k,m)
                  end if
               enddo
               close(ofileno)
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
   end subroutine calc_gk_projch

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_gkt
   !> @brief Calculate \f$\mathbf{S}\left(\mathbf{q},t\right)\f$ for obtaining \f$\mathbf{S}\left(\mathbf{q},\omega\right)\f$ after FT
   !-----------------------------------------------------------------------------
   subroutine calc_gkt(Natom, Mensemble, coord, simid, emomM, flag)
      !
      use Constants
      use BLS, only : gramms
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      integer, intent(inout) :: flag  !< Setup, sample, or print
      !
      integer :: iq,iw,step,r,l,i_stat,i_all,j
      character(len=30) :: filn
      complex(dblprec) :: epowqr, i, iqfac,tt, epowwt
      real(dblprec) :: qdr,nainv,mavg_norm, mm, win_fac
      real(dblprec), dimension(3) :: mavrg_vec
      real(dblprec), dimension(sc_nstep+1) :: dt
      real(dblprec), dimension(3,3) :: local_rotmat
      real(dblprec), dimension(3) :: k_rod,v_rod,v_para,v_perp
      !
      i=(0.0d0,1.0d0)


      if(flag==0) then
         ! First call, allocate and clear arrays

         allocate(corr_kt(3,nq,sc_nstep+1),stat=i_stat)
         call memocc(i_stat,product(shape(corr_kt))*kind(corr_kt),'corr_kt','calc_gkt')
         corr_kt=0.0d0
         !
         allocate(r_norm(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(r_norm))*kind(r_norm),'r_norm','calc_gkt')
         r_max=0.0d0
         do r=1,Natom
            r_norm(r)=coord(1,r)**2+coord(2,r)**2+coord(3,r)**2
            if(r_max<r_norm(r)) r_max=r_norm(r)
         end do
         !
         allocate(mavg_axis(3,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(mavg_axis))*kind(mavg_axis),'mavg_axis','calc_gkt')
         mavg_axis=0.0d0
         !
         allocate(mort_axis(3,3,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(mort_axis))*kind(mort_axis),'mort_axis','calc_gkt')
         mort_axis=0.0d0
         !
         if (do_sc_local_axis=='Y') then
            allocate(mavg_local_axis(3,Natom,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(mavg_local_axis))*kind(mavg_local_axis),'mavg_local_axis','calc_gkt')
            call calc_mavrg_vec(Natom,Mensemble,emomM,mavrg_vec,mavg_axis)
            mavg_local_axis=emomM
            allocate(mavg_local_rotmat(3,3,Natom,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(mavg_local_rotmat))*kind(mavg_local_rotmat),'mavg_local_rotmat','calc_gkt')
            call find_local_rotmat(Natom*Mensemble,emomM,mavg_local_rotmat)
         else
            call calc_mavrg_vec(Natom,Mensemble,emomM,mavrg_vec,mavg_axis)
            do l=1,Mensemble
               mavg_axis(:,l)=mavg_axis(:,l)/Natom
               mavg_norm=sum(mavg_axis(:,l)*mavg_axis(:,l))**0.5d0
               if(mavg_norm>1.0d-2) then
                  mavg_axis(:,l)=mavg_axis(:,l)/mavg_norm
               else
                  mavg_axis(:,l)=(/0.0d0,0.0d0,1.0d0/)
               end if
            end do
            call gramms(mavg_axis,mort_axis,Mensemble)
         end if
         !
         !
         allocate(m_proj(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(m_proj))*kind(m_proj),'m_proj','calc_gk')
         !

         flag=1
         sc_samp_done=0
         qfac=2.0d0*pi
      end if

      nainv=1.0d0/Natom

      ! Calculate g(k) for the current iteration and add to G(k,t)
      if (flag==1) then

         !corr_st=0.0d0
         iqfac=i*qfac

         call calc_mavrg_vec(Natom,Mensemble,emomM,mavrg_vec,mavg_axis)
         do l=1,Mensemble
            mavg_axis(:,l)=mavg_axis(:,l)/Natom
            mavg_norm=sum(mavg_axis(:,l)*mavg_axis(:,l))**0.5d0
            if(mavg_norm>1.0d-2) then
               mavg_axis(:,l)=mavg_axis(:,l)/mavg_norm
            else
               mavg_axis(:,l)=(/0.0d0,0.0d0,1.0d0/)
            end if
         end do

         if(do_connected/='Y'.and.do_sc_local_axis/='Y') then
            ! Keep the average magnetization for subtraction below if do_connected
            ! otherwise put it to zero (in order to remove a number of if-statements in-loop)
            mavrg_vec=0.0d0
         end if


         if(do_sc_local_axis/='Y') then
            !$omp parallel do default(shared) private(l,r) schedule(static) collapse(2)
            do l=1,Mensemble
               do r=1,Natom
                  m_proj(:,r,l)=emomM(:,r,l)-mavrg_vec(:)
               end do
            end do
            !$omp end parallel do
         else
            !! Sample current moment directions and weight in for average
            !$omp parallel do schedule(static) default(shared) private(l,r,local_rotmat) collapse(2)
            do l=1,Mensemble
               do r=1,Natom
                  mavg_local_axis(:,r,l)=mavg_local_axis(:,r,l)*(1.0_dblprec-sc_local_axis_mix)+emomM(:,r,l)*sc_local_axis_mix
               end do
            end do
            !$omp end parallel do
            ! Find the proper rotation axes
            call find_local_rotmat(Natom*Mensemble,mavg_local_axis,mavg_local_rotmat)
            !$omp parallel do default(shared) private(l,r,local_rotmat) schedule(static) collapse(2)
            do l=1,Mensemble
               do r=1,Natom
                  local_rotmat=mavg_local_rotmat(:,:,r,l)
                  m_proj(1,r,l)=local_rotmat(1,1)*emomM(1,r,l)+local_rotmat(2,1)*emomM(2,r,l)+local_rotmat(3,1)*emomM(3,r,l)
                  m_proj(2,r,l)=local_rotmat(1,2)*emomM(1,r,l)+local_rotmat(2,2)*emomM(2,r,l)+local_rotmat(3,2)*emomM(3,r,l)
                  m_proj(3,r,l)=local_rotmat(1,3)*emomM(1,r,l)+local_rotmat(2,3)*emomM(2,r,l)+local_rotmat(3,3)*emomM(3,r,l)
               end do
            end do
            !$omp end parallel do
         end if

         win_fac=1.0d0*nainv
         !$omp parallel do default(shared) private(r,iq,l,qdr,epowqr) schedule(static)
         do iq=1,nq
            do l=1,Mensemble
#if _OPENMP >= 201307
               !$omp simd   private(qdr,epowqr)
#endif
               do r=1,Natom
                  ! No need to window the q-transform
                  qdr=q(1,iq)*coord(1,r)+q(2,iq)*coord(2,r)+q(3,iq)*coord(3,r)
                  epowqr=exp(iqfac*qdr)*win_fac
                  !DIR$ vector always aligned
                  corr_kt(:,iq,sc_tidx)=corr_kt(:,iq,sc_tidx)+epowqr*m_proj(:,r,l)
               end do
            end do
         end do
         !$omp end parallel do
      end if

      ! Final operations, transform and print
      if (flag==2) then

         ! Allocate arrays
         allocate(corr_kw(3,nq,nw),stat=i_stat)
         call memocc(i_stat,product(shape(corr_kw))*kind(corr_kw),'corr_kw','calc_gkt')
         corr_kw=0.0d0

         ! Finish sampling and transform (k,t)->(k,w)
         if(sc_tidx.GT.sc_nstep) then
            sc_tidx = sc_nstep
         end if

         j=1
         do while(j.LE.sc_tidx)
            dt(j) = scstep_arr(j)*deltat_corr(j)
            j = j+1
         end do

         wfac=1.0d0
         corr_kw=0.0d0

         !$omp parallel do default(shared) private(iw,iq,step,tt,epowwt) schedule(static)
         do iw=1,nw
            do step=1,sc_tidx
               tt=i*wfac*(step-1)*dt(step)
               ! Apply window function as determined by 'sc_window_fun'
               epowwt=exp(w(iw)*tt)*sc_window_fac(sc_window_fun,step,sc_tidx)
               !
#if _OPENMP >= 201307
               !$omp simd
#endif
               do iq=1,nq
                  corr_kw(:,iq,iw)=corr_kw(:,iq,iw)+epowwt*corr_kt(:,iq,step)
               enddo
            enddo
         enddo
         !$omp end parallel do

         if (do_conv.eq.'LQ' .or. do_conv.eq.'LW' .or. do_conv.eq.'LW' .or. &
            do_conv.eq.'GQ' .or. do_conv.eq.'GW' .or. do_conv.eq.'GY') then
            ! Calculate the convolution for the magnon DOS
            call calc_dos_conv(corr_kw)
         endif

         if(do_sc_dosonly/="Y") then
            ! Write S(q,w)
            write (filn,'(''sqw.'',a8,''.out'')') simid
            open(ofileno, file=filn)
            do iq=1,Nq
               do iw=1,Nw/2
                  write (ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, &
                     abs(corr_kw(1,iq,iw)),abs(corr_kw(2,iq,iw)),abs(corr_kw(3,iq,iw)), &
                     abs(corr_kw(1,iq,iw)**2+corr_kw(2,iq,iw)**2+corr_kw(3,iq,iw)**2)**0.5d0
               end do
            end do
            close(ofileno)

            if(do_sc_complex=='Y') then
               ! Write S(q,w) decomposed in real and imaginary parts
               write (filn,'(''csqw.'',a8,''.out'')') simid
               open(ofileno, file=filn)
               do iq=1,Nq
                  do iw=1,Nw/2
                     write (ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, &
                        abs(corr_kw(1,iq,iw)),atan2(aimag(corr_kw(1,iq,iw)),real(corr_kw(1,iq,iw))), &
                        abs(corr_kw(2,iq,iw)),atan2(aimag(corr_kw(2,iq,iw)),real(corr_kw(2,iq,iw))), &
                        abs(corr_kw(3,iq,iw)),atan2(aimag(corr_kw(3,iq,iw)),real(corr_kw(3,iq,iw)))
                  end do
               end do
               close(ofileno)
            end if

         end if


         ! Calculate and write the magnon DOS
         allocate(magnon_dos(3,nw/2),stat=i_stat)
         call memocc(i_stat,product(shape(magnon_dos))*kind(magnon_dos),'magnon_dos','calc_gkt')
         !
         write (filn,'(''swdos.'',a8,''.out'')') simid
         open(ofileno, file=filn)
         magnon_dos=0.0d0
         if (qpoints=='I'.or.qpoints=='B') then
            do iw=1,Nw/2
               do iq=1,Nq
                  magnon_dos(1,iw)=magnon_dos(1,iw)+abs(corr_kw(1,iq,iw))*q_weight(iq) ! Put a weight for the IR BZ points
                  magnon_dos(2,iw)=magnon_dos(2,iw)+abs(corr_kw(2,iq,iw))*q_weight(iq)
                  magnon_dos(3,iw)=magnon_dos(3,iw)+abs(corr_kw(3,iq,iw))*q_weight(iq)
               end do
               write (ofileno,10006) hbar_mev*w(iw),magnon_dos(1,iw),magnon_dos(2,iw), magnon_dos(3,iw), &
                  sqrt(magnon_dos(1,iw)**2+magnon_dos(2,iw)**2+magnon_dos(3,iw)**2)
            end do
         else
            do iw=1,Nw/2
               do iq=1,Nq
                  magnon_dos(1,iw)=magnon_dos(1,iw)+abs(corr_kw(1,iq,iw)) ! Put a weight for the IR BZ points
                  magnon_dos(2,iw)=magnon_dos(2,iw)+abs(corr_kw(2,iq,iw))
                  magnon_dos(3,iw)=magnon_dos(3,iw)+abs(corr_kw(3,iq,iw))
               end do
               write (ofileno,10006) hbar_mev*w(iw),magnon_dos(1,iw),magnon_dos(2,iw), magnon_dos(3,iw), &
                  sqrt(magnon_dos(1,iw)**2+magnon_dos(2,iw)**2+magnon_dos(3,iw)**2)
            end do
         endif
         close(ofileno)
         write (filn,'(''quasiblsdos.'',a8,''.out'')') simid
         open(ofileno, file=filn)
         magnon_dos=0.0d0
         if (qpoints=='I'.or.qpoints=='B') then
            do iw=1,Nw/2
               do iq=1,Nq
                  magnon_dos(1,iw)=magnon_dos(1,iw)+real(corr_kw(1,iq,iw)*conjg(corr_kw(1,iq,iw)))*q_weight(iq) ! Put a weight for the IR BZ points
                  magnon_dos(2,iw)=magnon_dos(2,iw)+real(corr_kw(2,iq,iw)*conjg(corr_kw(2,iq,iw)))*q_weight(iq)
                  magnon_dos(3,iw)=magnon_dos(3,iw)+real(corr_kw(3,iq,iw)*conjg(corr_kw(3,iq,iw)))*q_weight(iq)
               end do
               write (ofileno,10006) hbar_mev*w(iw),magnon_dos(1,iw),magnon_dos(2,iw), magnon_dos(3,iw), &
                  sqrt(magnon_dos(1,iw)**2+magnon_dos(2,iw)**2+magnon_dos(3,iw)**2)
            end do
         else
            do iw=1,Nw/2
               do iq=1,Nq
                  magnon_dos(1,iw)=magnon_dos(1,iw)+real(corr_kw(1,iq,iw)*conjg(corr_kw(1,iq,iw)))
                  magnon_dos(2,iw)=magnon_dos(2,iw)+real(corr_kw(2,iq,iw)*conjg(corr_kw(2,iq,iw)))
                  magnon_dos(3,iw)=magnon_dos(3,iw)+real(corr_kw(3,iq,iw)*conjg(corr_kw(3,iq,iw)))
               end do
               write (ofileno,10006) hbar_mev*w(iw),magnon_dos(1,iw),magnon_dos(2,iw), magnon_dos(3,iw), &
                  sqrt(magnon_dos(1,iw)**2+magnon_dos(2,iw)**2+magnon_dos(3,iw)**2)
            end do
         endif
         close(ofileno)

         ! Deallocate arrays
         i_all=-product(shape(corr_kt))*kind(corr_kt)
         deallocate(corr_kt,stat=i_stat)
         call memocc(i_stat,i_all,'corr_kt','calc_gkt')
         !
         i_all=-product(shape(corr_kw))*kind(corr_kw)
         deallocate(corr_kw,stat=i_stat)
         call memocc(i_stat,i_all,'corr_kw','calc_gkt')
         !
         i_all=-product(shape(magnon_dos))*kind(magnon_dos)
         deallocate(magnon_dos,stat=i_stat)
         call memocc(i_stat,i_all,'magnon_dos','calc_gkt')
         !
         i_all=-product(shape(r_norm))*kind(r_norm)
         deallocate(r_norm,stat=i_stat)
         call memocc(i_stat,i_all,'r_norm','calc_gkt')
         !
         if (do_sc_local_axis=='Y') then
            i_all=-product(shape(mavg_local_axis))*kind(mavg_local_axis)
            deallocate(mavg_local_axis,stat=i_stat)
            call memocc(i_stat,i_all,'mavg_local_axis','calc_gkt')
            i_all=-product(shape(mavg_local_rotmat))*kind(mavg_local_rotmat)
            deallocate(mavg_local_rotmat,stat=i_stat)
            call memocc(i_stat,i_all,'mavg_local_rotmat','calc_gkt')
         end if
         !
         i_all=-product(shape(mort_axis))*kind(mort_axis)
         deallocate(mort_axis,stat=i_stat)
         call memocc(i_stat,i_all,'mort_axis','calc_gkt')
         !
         i_all=-product(shape(m_proj))*kind(m_proj)
         deallocate(m_proj,stat=i_stat)
         call memocc(i_stat,i_all,'m_proj','calc_gkt')
         !
      end if
      return
      !
      10005 format (i7,3f10.6,2x,i7,2x,9es16.8)
      10006 format (5G16.8)
      !
   end subroutine calc_gkt

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_gkt_projch
   !> @brief Calculate \f$\mathbf{S}\left(\mathbf{q},t\right)\f$ for obtaining \f$\mathbf{S}\left(\mathbf{q},\omega\right)\f$ after FT
   !-----------------------------------------------------------------------------
   subroutine calc_gkt_projch(Natom, Mensemble, coord, simid, emomM, flag, achtype, nchmax)
      !
      use Constants
      use BLS, only : gramms
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      integer, intent(inout) :: flag  !< Setup, sample, or print
      integer, dimension(Natom), intent(in) :: achtype !< Chemical Type of atom
      integer, intent(in) :: nchmax  !< Number of chemical types
      !
      integer :: iq,iw,step,r,j,l,i_stat,i_all,k
      integer,dimension(nchmax) :: ncounter
      character(len=30) :: filn
      real(dblprec), dimension(sc_nstep) :: dt
      real(dblprec) :: qdr
      real(dblprec),dimension(nchmax) :: nainv
      real(dblprec),dimension(3,3) :: local_rotmat
      complex(dblprec) :: i, iqfac, iwfac, tt, epowqr,epowwt
      !
      i=(0.0d0,1.0d0)
      !
      if(flag==0) then
         ! First call, allocate and clear arrays
         allocate(corr_st_proj(3,nq,nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(corr_st_proj))*kind(corr_st_proj),'corr_st_proj','calc_gkt_projch')
         !
         allocate(corr_kt_proj(3,nq,sc_nstep+1,nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(corr_kt_proj))*kind(corr_kt_proj),'corr_kt_proj','calc_gkt_projch')
         corr_kt_proj=0.0d0
         !
         if (do_sc_proj_axis=='Y') then
            allocate(mavg_proj_axis(3,Nchmax,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(mavg_proj_axis))*kind(mavg_proj_axis),'mavg_proj_axis','calc_gkt_projch')
            allocate(mavg_proj_rotmat(3,3,Nchmax,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(mavg_proj_rotmat))*kind(mavg_proj_rotmat),'mavg_proj_rotmat','calc_gkt_projch')
            ! Sample current moment directions and weight in for average
            mavg_proj_axis=0.0d0
            do l=1,Mensemble
               do r=1,Natom
                  k=achtype(r)
                  mavg_proj_axis(:,k,l)=mavg_proj_axis(:,k,l)+emomM(:,r,l)
               end do
            end do
            call find_local_rotmat(nchmax*Mensemble,mavg_proj_axis,mavg_proj_rotmat)
         end if
         !
         allocate(mort_axis_p(3,3,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(mort_axis_p))*kind(mort_axis_p),'mort_axis_p','calc_gkt_projch')
         mort_axis_p=0.0d0
         !
         allocate(m_proj_p(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(m_proj_p))*kind(m_proj_p),'m_proj_p','calc_gkt_projch')
         !
         !
         flag=1
         sc_samp_done=0
         qfac=2.0d0*pi
      end if

      ncounter=0
      do k=1,nchmax
         do r=1,Natom
            if (k==achtype(r)) then
               ncounter(k)=ncounter(k)+1
            end if
         enddo
      enddo
      nainv(:)=1.0d0/ncounter(:)

      ! Calculate g(k) for the current iteration and add to G(k,t)
      if (flag==1) then
         iqfac=i*qfac
         corr_st_proj=0.0d0
         if(do_sc_proj_axis/='Y') then
            do l=1,Mensemble
               do r=1,Natom
                  m_proj_p(1,r,l)=emomM(1,r,l)
                  m_proj_p(2,r,l)=emomM(2,r,l)
                  m_proj_p(3,r,l)=emomM(3,r,l)
               end do
            end do
         else
            !$omp parallel do default(shared) private(l,r,k,local_rotmat) schedule(static) collapse(2)
            do l=1,Mensemble
               do r=1,Natom
                  k=achtype(r)
                  local_rotmat=mavg_proj_rotmat(:,:,1,l)
                  m_proj_p(1,r,l)=local_rotmat(1,1)*emomM(1,r,l)+local_rotmat(2,1)*emomM(2,r,l)+local_rotmat(3,1)*emomM(3,r,l)
                  m_proj_p(2,r,l)=local_rotmat(1,2)*emomM(1,r,l)+local_rotmat(2,2)*emomM(2,r,l)+local_rotmat(3,2)*emomM(3,r,l)
                  m_proj_p(3,r,l)=local_rotmat(1,3)*emomM(1,r,l)+local_rotmat(2,3)*emomM(2,r,l)+local_rotmat(3,3)*emomM(3,r,l)
               end do
            end do
            !$omp end parallel do
         end if

         !$omp parallel do default(shared) private(r,k,iq,l,qdr,epowqr) schedule(static)
         do iq=1,nq
            do l=1,Mensemble
               corr_st_proj(:,iq,:)=(0.0d0,0.0d0)
               do r=1,Natom
                  k=achtype(r)
                  qdr=q(1,iq)*coord(1,r)+q(2,iq)*coord(2,r)+q(3,iq)*coord(3,r)
                  epowqr=exp(iqfac*qdr)*nainv(k)
                  corr_st_proj(1,iq,k)=corr_st_proj(1,iq,k)+epowqr*m_proj_p(1,r,l)
                  corr_st_proj(2,iq,k)=corr_st_proj(2,iq,k)+epowqr*m_proj_p(2,r,l)
                  corr_st_proj(3,iq,k)=corr_st_proj(3,iq,k)+epowqr*m_proj_p(3,r,l)
               end do
               do k=1,nchmax
                  corr_kt_proj(1,iq,sc_tidx,k)=corr_kt_proj(1,iq,sc_tidx,k) + corr_st_proj(1,iq,k)
                  corr_kt_proj(2,iq,sc_tidx,k)=corr_kt_proj(2,iq,sc_tidx,k) + corr_st_proj(2,iq,k)
                  corr_kt_proj(3,iq,sc_tidx,k)=corr_kt_proj(3,iq,sc_tidx,k) + corr_st_proj(3,iq,k)
               end do
            end do
         end do
         !$omp end parallel do
      end if

      ! Final operations, transform and print
      if (flag==2) then

         ! Allocate arrays
         i_all=-product(shape(corr_st_proj))*kind(corr_st_proj)
         deallocate(corr_st_proj,stat=i_stat)
         call memocc(i_stat,i_all,'corr_st_proj','calc_gkt_projch')
         !
         allocate(corr_kw_proj(3,nq,nw,nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(corr_kw_proj))*kind(corr_kw_proj),'corr_kw_proj','calc_gkt_projch')
         corr_kw_proj=0.0d0

         ! Finish sampling and transform (k,t)->(k,w)

         if(sc_tidx.GT.sc_nstep) then
            sc_tidx = sc_nstep
         end if

         j=1
         do while(j.LE.sc_tidx)
            dt(j) = scstep_arr(j)*deltat_corr(j)
            j = j+1
         end do

         wfac=1.0d0
         iwfac=i*wfac
         !$omp parallel do default(shared) private(iw,iq,k,step,tt,epowwt) schedule(static)
         do iw=1,nw
            do k=1,nchmax
               do step=1,sc_tidx
                  tt=iwfac*(step-1)*dt(step)
                  epowwt=exp(w(iw)*tt)
                  do iq=1,nq
                     corr_kw_proj(1,iq,iw,k)=corr_kw_proj(1,iq,iw,k)+epowwt*corr_kt_proj(1,iq,step,k)
                     corr_kw_proj(2,iq,iw,k)=corr_kw_proj(2,iq,iw,k)+epowwt*corr_kt_proj(2,iq,step,k)
                     corr_kw_proj(3,iq,iw,k)=corr_kw_proj(3,iq,iw,k)+epowwt*corr_kt_proj(3,iq,step,k)
                  enddo
               enddo
            end do
         enddo
         !$omp end parallel do

         ! Write S(q,w)
         do k=1,nchmax
            if(k>10) then
               write (filn,'(''projchsqw.'',i2,''.'',a8,''.out'')') k,simid
            else
               write (filn,'(''projchsqw.'',i1,''.'',a8,''.out'')') k,simid
            end if
            open(ofileno, file=filn)
            do iq=1,Nq
               do iw=1,Nw
                  write (ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, &
                     abs(corr_kw_proj(1,iq,iw,k)),abs(corr_kw_proj(2,iq,iw,k)),abs(corr_kw_proj(3,iq,iw,k)), &
                     abs((corr_kw_proj(1,iq,iw,k)**2+corr_kw_proj(2,iq,iw,k)**2+corr_kw_proj(3,iq,iw,k)**2)**0.5d0)
               end do
            end do
            close(ofileno)
         end do

         ! Deallocate arrays
         i_all=-product(shape(corr_kt_proj))*kind(corr_kt_proj)
         deallocate(corr_kt_proj,stat=i_stat)
         call memocc(i_stat,i_all,'corr_kt_proj','calc_gkt_projch')
         !
         i_all=-product(shape(corr_kw_proj))*kind(corr_kw_proj)
         deallocate(corr_kw_proj,stat=i_stat)
         call memocc(i_stat,i_all,'corr_kw_proj','calc_gkt_projch')
         !
         !
         if (do_sc_proj_axis=='Y') then
            i_all=-product(shape(mavg_proj_axis))*kind(mavg_proj_axis)
            deallocate(mavg_proj_axis,stat=i_stat)
            call memocc(i_stat,i_all,'mavg_proj_axis','calc_gkt_projch')
            i_all=-product(shape(mavg_proj_rotmat)*kind(mavg_proj_rotmat))
            deallocate(mavg_proj_rotmat,stat=i_stat)
            call memocc(i_stat,i_all,'mavg_proj_rotmat','calc_gkt_projch')
         end if
         !
         i_all=-product(shape(mort_axis_p))*kind(mort_axis_p)
         deallocate(mort_axis_p,stat=i_stat)
         call memocc(i_stat,i_all,'mort_axis_p','calc_gkt_projch')
         !
         i_all=-product(shape(m_proj_p))*kind(m_proj_p)
         deallocate(m_proj_p,stat=i_stat)
         call memocc(i_stat,i_all,'m_proj_p','calc_gkt_projch')
         !
      end if

      return
      !
      10005 format (i7,3f10.6,2x,i7,2x,9es16.8)
      !
   end subroutine calc_gkt_projch

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_gkt_proj
   !> @brief Calculate \f$\mathbf{S}\left(\mathbf{q},t\right)\f$ for obtaining \f$\mathbf{S}\left(\mathbf{q},\omega\right)\f$ after FT
   !-----------------------------------------------------------------------------
   subroutine calc_gkt_proj(Natom, Mensemble, coord, simid, emomM, flag, atype, nt)
      !
      use Constants
      use BLS, only : gramms!, gramms2
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      integer, intent(inout) :: flag  !< Setup, sample, or print
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      integer, intent(in) :: NT  !< Number of types of atoms
      !
      integer :: iq,iw,step,r,j,l,i_stat,i_all,k
      integer,dimension(nt) :: ncounter
      character(len=30) :: filn
      real(dblprec), dimension(sc_nstep) :: dt
      real(dblprec) :: qdr
      real(dblprec),dimension(nt) :: nainv
      real(dblprec),dimension(3,3) :: local_rotmat
      complex(dblprec) :: i, iqfac, iwfac, tt, epowqr,epowwt
      !
      i=(0.0d0,1.0d0)
      !
      if(flag==0) then
         ! First call, allocate and clear arrays
         allocate(corr_st_proj(3,nq,nt),stat=i_stat)
         call memocc(i_stat,product(shape(corr_st_proj))*kind(corr_st_proj),'corr_st_proj','calc_gkt_proj')
         !
         allocate(corr_kt_proj(3,nq,sc_nstep+1,nt),stat=i_stat)
         call memocc(i_stat,product(shape(corr_kt_proj))*kind(corr_kt_proj),'corr_kt_proj','calc_gkt_proj')
         corr_kt_proj=0.0d0
         !
         if (do_sc_proj_axis=='Y') then
            allocate(mavg_proj_axis(3,Nt,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(mavg_proj_axis))*kind(mavg_proj_axis),'mavg_proj_axis','calc_gkt_proj')
            allocate(mavg_proj_rotmat(3,3,Nt,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(mavg_proj_rotmat))*kind(mavg_proj_rotmat),'mavg_proj_rotmat','calc_gkt_proj')
            ! Sample current moment directions and weight in for average
            mavg_proj_axis=0.0d0
            do l=1,Mensemble
               do r=1,Natom
                  k=atype(r)
                  mavg_proj_axis(:,k,l)=mavg_proj_axis(:,k,l)+emomM(:,r,l)
               end do
            end do
            call find_local_rotmat(nt*Mensemble,mavg_proj_axis,mavg_proj_rotmat)
         end if
         !
         allocate(mort_axis_p(3,3,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(mort_axis_p))*kind(mort_axis_p),'mort_axis_p','calc_gkt_proj')
         mort_axis_p=0.0d0
         !
         allocate(m_proj_p(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(m_proj_p))*kind(m_proj_p),'m_proj_p','calc_gkt_proj')
         !
         !
         flag=1
         sc_samp_done=0
         qfac=2.0d0*pi
      end if

      ncounter=0
      do k=1,nt
         do r=1,Natom
            if (k==atype(r)) then
               ncounter(k)=ncounter(k)+1
            end if
         enddo
      enddo
      nainv(:)=1.0d0/ncounter(:)

      ! Calculate g(k) for the current iteration and add to G(k,t)
      if (flag==1) then
         iqfac=i*qfac
         corr_st_proj=0.0d0
         if(do_sc_proj_axis/='Y') then
            do l=1,Mensemble
               do r=1,Natom
                  m_proj_p(1,r,l)=emomM(1,r,l)
                  m_proj_p(2,r,l)=emomM(2,r,l)
                  m_proj_p(3,r,l)=emomM(3,r,l)
               end do
            end do
         else
            !$omp parallel do default(shared) private(l,r,k,local_rotmat) schedule(static) collapse(2)
            do l=1,Mensemble
               do r=1,Natom
                  k=atype(r)
                  local_rotmat=mavg_proj_rotmat(:,:,1,l)
                  m_proj_p(1,r,l)=local_rotmat(1,1)*emomM(1,r,l)+local_rotmat(2,1)*emomM(2,r,l)+local_rotmat(3,1)*emomM(3,r,l)
                  m_proj_p(2,r,l)=local_rotmat(1,2)*emomM(1,r,l)+local_rotmat(2,2)*emomM(2,r,l)+local_rotmat(3,2)*emomM(3,r,l)
                  m_proj_p(3,r,l)=local_rotmat(1,3)*emomM(1,r,l)+local_rotmat(2,3)*emomM(2,r,l)+local_rotmat(3,3)*emomM(3,r,l)
               end do
            end do
            !$omp end parallel do
         end if

         !$omp parallel do default(shared) private(r,k,iq,l,qdr,epowqr) schedule(static)
         do iq=1,nq
            do l=1,Mensemble
               corr_st_proj(:,iq,:)=(0.0d0,0.0d0)
               do r=1,Natom
                  k=atype(r)
                  qdr=q(1,iq)*coord(1,r)+q(2,iq)*coord(2,r)+q(3,iq)*coord(3,r)
                  epowqr=exp(iqfac*qdr)*nainv(k)
                  corr_st_proj(1,iq,k)=corr_st_proj(1,iq,k)+epowqr*m_proj_p(1,r,l)
                  corr_st_proj(2,iq,k)=corr_st_proj(2,iq,k)+epowqr*m_proj_p(2,r,l)
                  corr_st_proj(3,iq,k)=corr_st_proj(3,iq,k)+epowqr*m_proj_p(3,r,l)
               end do
               do k=1,nt
                  corr_kt_proj(1,iq,sc_tidx,k)=corr_kt_proj(1,iq,sc_tidx,k) + corr_st_proj(1,iq,k)
                  corr_kt_proj(2,iq,sc_tidx,k)=corr_kt_proj(2,iq,sc_tidx,k) + corr_st_proj(2,iq,k)
                  corr_kt_proj(3,iq,sc_tidx,k)=corr_kt_proj(3,iq,sc_tidx,k) + corr_st_proj(3,iq,k)
               end do
            end do
         end do
         !$omp end parallel do
      end if

      ! Final operations, transform and print
      if (flag==2) then

         ! Allocate arrays
         i_all=-product(shape(corr_st_proj))*kind(corr_st_proj)
         deallocate(corr_st_proj,stat=i_stat)
         call memocc(i_stat,i_all,'corr_st_proj','calc_gkt_proj')
         !
         allocate(corr_kw_proj(3,nq,nw,nt),stat=i_stat)
         call memocc(i_stat,product(shape(corr_kw_proj))*kind(corr_kw_proj),'corr_kw_proj','calc_gkt_proj')
         corr_kw_proj=0.0d0

         ! Finish sampling and transform (k,t)->(k,w)
         if(sc_tidx.GT.sc_nstep) then
            sc_tidx = sc_nstep
         end if

         j=1
         do while(j.LE.sc_tidx)
            dt(j) = scstep_arr(j)*deltat_corr(j)
            j = j+1
         end do

         wfac=1.0d0
         iwfac=i*wfac
         !$omp parallel do default(shared) private(iw,iq,k,step,tt,epowwt) schedule(static)
         do iw=1,nw
            do k=1,nt
               do step=1,sc_tidx
                  tt=iwfac*(step-1)*dt(step)
                  epowwt=exp(w(iw)*tt)
                  do iq=1,nq
                     corr_kw_proj(1,iq,iw,k)=corr_kw_proj(1,iq,iw,k)+epowwt*corr_kt_proj(1,iq,step,k)
                     corr_kw_proj(2,iq,iw,k)=corr_kw_proj(2,iq,iw,k)+epowwt*corr_kt_proj(2,iq,step,k)
                     corr_kw_proj(3,iq,iw,k)=corr_kw_proj(3,iq,iw,k)+epowwt*corr_kt_proj(3,iq,step,k)
                  enddo
               enddo
            end do
         enddo
         !$omp end parallel do

         ! Write S(q,w)
         do k=1,nt
            if(k>10) then
               write (filn,'(''projsqw.'',i2,''.'',a8,''.out'')') k,simid
            else
               write (filn,'(''projsqw.'',i1,''.'',a8,''.out'')') k,simid
            end if
            open(ofileno, file=filn)
            do iq=1,Nq
               do iw=1,Nw
                  write (ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, &
                     abs(corr_kw_proj(1,iq,iw,k)),abs(corr_kw_proj(2,iq,iw,k)),abs(corr_kw_proj(3,iq,iw,k)), &
                     abs((corr_kw_proj(1,iq,iw,k)**2+corr_kw_proj(2,iq,iw,k)**2+corr_kw_proj(3,iq,iw,k)**2)**0.5d0)
               end do
            end do
            close(ofileno)
         end do

         ! Write S(q,w) decomposed in real and imag
         do k=1,nt
            if(k>10) then
               write (filn,'(''cprojsqw.'',i2,''.'',a8,''.out'')') k,simid
            else
               write (filn,'(''cprojsqw.'',i1,''.'',a8,''.out'')') k,simid
            end if
            open(ofileno, file=filn)
            do iq=1,Nq
               do iw=1,Nw
                  write (ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, &
                     abs(corr_kw_proj(1,iq,iw,k)),atan2(aimag(corr_kw_proj(1,iq,iw,k)),real(corr_kw_proj(1,iq,iw,k))), &
                     abs(corr_kw_proj(2,iq,iw,k)),atan2(aimag(corr_kw_proj(2,iq,iw,k)),real(corr_kw_proj(2,iq,iw,k))), &
                     abs(corr_kw_proj(3,iq,iw,k)),atan2(aimag(corr_kw_proj(3,iq,iw,k)),real(corr_kw_proj(3,iq,iw,k)))

               end do
            end do
            close(ofileno)
         end do

         ! Deallocate arrays
         i_all=-product(shape(corr_kt_proj))*kind(corr_kt_proj)
         deallocate(corr_kt_proj,stat=i_stat)
         call memocc(i_stat,i_all,'corr_kt_proj','calc_gkt_proj')
         !
         i_all=-product(shape(corr_kw_proj))*kind(corr_kw_proj)
         deallocate(corr_kw_proj,stat=i_stat)
         call memocc(i_stat,i_all,'corr_kw_proj','calc_gkt_proj')
         !
         !
         if (do_sc_proj_axis=='Y') then
            i_all=-product(shape(mavg_proj_axis))*kind(mavg_proj_axis)
            deallocate(mavg_proj_axis,stat=i_stat)
            call memocc(i_stat,i_all,'mavg_proj_axis','calc_gkt_proj')
            i_all=-product(shape(mavg_proj_rotmat)*kind(mavg_proj_rotmat))
            deallocate(mavg_proj_rotmat,stat=i_stat)
            call memocc(i_stat,i_all,'mavg_proj_rotmat','calc_gkt_proj')
         end if
         !
         i_all=-product(shape(mort_axis_p))*kind(mort_axis_p)
         deallocate(mort_axis_p,stat=i_stat)
         call memocc(i_stat,i_all,'mort_axis_p','calc_gkt_proj')
         !
         i_all=-product(shape(m_proj_p))*kind(m_proj_p)
         deallocate(m_proj_p,stat=i_stat)
         call memocc(i_stat,i_all,'m_proj_p','calc_gkt_proj')
         !
      end if

      return
      !
      10005 format (i7,3f10.6,2x,i7,2x,9es16.8)
      !
   end subroutine calc_gkt_proj

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_sr
   !> @brief Perform only spatial correlation in real space to be able to deal with non periodic systems
   !-----------------------------------------------------------------------------
   subroutine calc_sr(Natom, Mensemble, coord, simid, emomM, cr_flag)
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
         corr_sr=0.0d0
         cr_flag=1
         sc_samp_done_sr=0
      end if

      nainv=1.0d0/Natom

      if (cr_flag==1) then
         ! Calculate s(k) for the current iteration and add to average of G(k)

         if (do_connected=='Y') then
            !$omp parallel do default(shared) private(iatom,l) schedule(static)
            do l=1,Mensemble
               do iatom=1,Natom
                  connected(1,iatom)=connected(1,iatom)+emomM(1,iatom,l)
                  connected(2,iatom)=connected(2,iatom)+emomM(2,iatom,l)
                  connected(3,iatom)=connected(3,iatom)+emomM(3,iatom,l)
               enddo
            enddo
            !$omp end parallel do
         else
            connected(:,:)=0.d0
         endif

         !$omp parallel do default(shared) private(r,iatom,l) schedule(static)
         do l=1,Mensemble
            do iatom=1,Natom
               do r=1,Natom
                  corr_sr(1,iatom)=corr_sr(1,iatom)+emomM(1,iatom,l)*emomM(1,r,l)*nainv-connected(1,iatom)*connected(1,r)/Mensemble
                  corr_sr(2,iatom)=corr_sr(2,iatom)+emomM(2,iatom,l)*emomM(2,r,l)*nainv-connected(2,iatom)*connected(2,r)/Mensemble
                  corr_sr(3,iatom)=corr_sr(3,iatom)+emomM(3,iatom,l)*emomM(3,r,l)*nainv-connected(3,iatom)*connected(3,r)/Mensemble
               end do
            enddo
         end do
         !$omp end parallel do

         sc_samp_done_sr=sc_samp_done_sr+1
      end if

      if (cr_flag==2) then
         ! Finish sampling and write S(q)
         corr_sr=corr_sr/(sc_samp_done_sr*Mensemble)

         ! Write G(r)
         write (filn,'(''dir_sr.'',a8,''.out'')') simid
         open(ofileno,file=filn,status='replace')
         do r=1,Natom
            write(ofileno,'(i10,3f10.4,5f18.8)') r,(coord(l,r),l=1,3),(((corr_sr(l,r))),l=1,3),&
               sqrt(corr_sr(1,r)**2+corr_sr(2,r)**2+corr_sr(3,r)**2),corr_sr(1,r)+corr_sr(2,r)+corr_sr(3,r)
         end do
         close(ofileno)

         ! Write G(|r|)
         write (filn,'(''dir_sra.'',a8,''.out'')') simid
         open(ofileno,file=filn,status='replace')
         do r=1,Natom
            write(ofileno,'(7f18.8)') sqrt((coord(1,r)-coord(1,1))**2+(coord(2,r)-coord(2,1))**2+(coord(3,r)-coord(3,1))**2),&
               (((corr_sr(l,r))),l=1,3),&
               sqrt(corr_sr(1,r)**2+corr_sr(2,r)**2+corr_sr(3,r)**2),corr_sr(1,r)+corr_sr(2,r)+corr_sr(3,r)
         end do
         close(ofileno)

         ! Deallocate arrays
         i_all=-product(shape(corr_sr))*kind(corr_s)
         deallocate(corr_sr,stat=i_stat)
         call memocc(i_stat,i_all,'corr_sr','calc_sr')
      end if
      return

   end subroutine calc_sr

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_dos_conv
   !> @brief Perform convolutions of the magnon DOS
   !-----------------------------------------------------------------------------
   subroutine calc_dos_conv(corr_cv)

      use Constants

      implicit none

      complex(dblprec), dimension(3,nq,nw), intent(inout) :: corr_cv

      integer :: iq,iw,j,t, conv_range_w, conv_range_q
      integer :: u, i_stat, i_all
      complex(dblprec) :: conv_cutoff_w, facw1, facw2, sfacww ! Variables for the convolution in frequency
      complex(dblprec) :: conv_cutoff_q, facq1, facq2, sfacqq ! Variables for the convolution in qpoints
      complex(dblprec) :: qq,ww
      complex(dblprec), dimension(:,:,:), allocatable  :: corr_tmp

      allocate(corr_tmp(3,nq,nw),stat=i_stat)
      call memocc(i_stat,product(shape(corr_tmp))*kind(corr_tmp),'calc_dos_conv','calc_dos_conv')

      ! Variables for the convolutions in frequency and qpoints domain, this are optional for the magnon DOS calculations
      ! This are the variables for qpoints convolutions
      if (do_conv=='GQ') then
         ! Variables for the reciprocal space for Gaussian function
         facq1 = -1.0d0/(2.0d0*sigma_q**2)
         facq2 = 1.0d0/(sqrt(2.0d0*pi)*sigma_q)
         conv_cutoff_q=0.001d0
         conv_range_q=int(sqrt(log(conv_cutoff_q)/facq1)+0.5d0)

      else if (do_conv=='LQ') then
         ! Variables for the reciprocal space for Lorentz distribution
         facq1=-1.0d0/(LQfactor**2)
         facq2=LQfactor/(Nq*pi)
         conv_cutoff_q=0.01d0
         conv_range_q=int(sqrt(log(conv_cutoff_q)/facq1)+0.5d0)

      else if (do_conv=='GW') then
         ! Variables for the frequency convolutions for Gaussian function
         facw1 = -1.0d0/(2*sigma_w**2)
         facw2 = 1.0d0/(sqrt(2*pi)*sigma_w)
         conv_cutoff_w=0.01d0
         conv_range_w=int(sqrt(log(conv_cutoff_w)/facw1)+0.5d0)

      else if (do_conv=='LW') then
         ! Variables for the frequency convolution with Lorentz distribution
         facw1=-1.0d0/(LWfactor**2)
         facw2=LWfactor/(Nw*pi)
         conv_cutoff_w=0.01d0
         conv_range_w=int(sqrt(log(conv_cutoff_w)/facw1)+0.5d0)

      else if (do_conv=='GY') then
         ! Variables for both qpoints and frequencies convolutions for Gaussian function
         facq1 = -1.0d0/(2*sigma_q**2)
         facw1 = -1.0d0/(2*sigma_w**2)
         facq2 = 1.0d0/(sqrt(2*pi)*sigma_q)
         facw2 = 1.0d0/(sqrt(2*pi)*sigma_w)
         conv_cutoff_q=0.001d0
         conv_cutoff_w=0.001d0
         conv_range_q=int(sqrt(log(conv_cutoff_q)/facq1)+0.5d0)
         conv_range_w=int(sqrt(log(conv_cutoff_w)/facw1)+0.5d0)
         print *,'range',conv_range_q,conv_range_w

      else if (do_conv=='LY') then
         ! Variables for both qpoints and frequencies convolutions for Lorentz distribution
         facq1 = -1.0d0/(LQfactor**2)
         facw1 = -1.0d0/(LWfactor**2)
         facq2 = LQfactor/(Nq*pi)
         facw2 = LWfactor/(Nw*pi)
         conv_cutoff_q=0.01d0
         conv_cutoff_w=0.01d0
         conv_range_q=int(sqrt(log(conv_cutoff_q)/facq1)+0.5d0)
         conv_range_w=int(sqrt(log(conv_cutoff_q)/facw1)+0.5d0)

      endif

      corr_tmp=corr_cv
      corr_cv=0.0d0

      ! Do the convolution of the corr_cv to smooth out the magnon DOS
      if (do_conv=='GQ') then
         write(*,*) '-------------------------> CONV_GQ',facq1,facq2
         !$omp parallel do default(shared) private(iq,iw,u,qq,sfacqq)
         do iq=1,nq
            do iw=1,nw !/2
               ! Convolution with Gaussian resolution function
               do u = max(1,iq-conv_range_q), min(nq-1,iq+conv_range_q)
                  qq=(iq-u)
                  sfacqq=exp(facq1*qq**2)*facq2 ! Factor which controls the convolution
                  if (iw==1.and.iq==nq/2) write(100,*) u,abs(sfacqq)
                  ! Convolution over the reciprocal space
                  corr_cv(1,iq,iw)=corr_cv(1,iq,iw)+corr_tmp(1,u,iw)*sfacqq
                  corr_cv(2,iq,iw)=corr_cv(2,iq,iw)+corr_tmp(2,u,iw)*sfacqq
                  corr_cv(3,iq,iw)=corr_cv(3,iq,iw)+corr_tmp(3,u,iw)*sfacqq
               enddo
            enddo
         enddo
         !$omp end parallel do
      else if (do_conv=='GW') then
         do iq=1,Nq
            !$omp parallel do default(shared) private(iw,j,t,sfacww)
            do iw=1, nw  !/2
               if (iw==1.and.iq==1) write(*,*) '-------------------------> CONV_GW'
               ! Convolution with Gaussian resolution function
               do j = max(1,iw-conv_range_w), min(nw,iw+conv_range_w)
                  ww=(iw-j)
                  sfacww=exp(facw1*ww**2)*facw2 ! This is the parameter that controls the convolution (Should the sum be over the absolute values?)
                  if (iq==10.and.iw==200) write(100,*) j,abs(sfacww)
                  ! Convolution over frequency space with a gaussian resolution function
                  corr_cv(1,iq,iw)=corr_cv(1,iq,iw)+corr_tmp(1,iq,j)*sfacww
                  corr_cv(2,iq,iw)=corr_cv(2,iq,iw)+corr_tmp(2,iq,j)*sfacww
                  corr_cv(3,iq,iw)=corr_cv(3,iq,iw)+corr_tmp(3,iq,j)*sfacww
               enddo
            enddo
            !$omp end parallel do
         enddo
      else if (do_conv=='GY') then
         do iq=1,Nq
            !$omp parallel do default(shared) private(iw,u,qq,sfacqq)
            do iw=1, nw! /2
               if (iw==1.and.iq==1) write(*,*) '-------------------------> CONV_GY',facq1,facq2
               ! Convolution with Gaussian resolution function
               !do u = 1,Nq
               do u = max(1,iq-conv_range_q), min(nq,iq+conv_range_q)
                  qq=(iq-u)
                  sfacqq=exp(facq1*qq**2)*facq2 ! This parameter controls the convolution
                  ! Convolution over the qpoints
                  corr_cv(1,iq,iw)=corr_cv(1,iq,iw)+corr_tmp(1,u,iw)*sfacqq
                  corr_cv(2,iq,iw)=corr_cv(2,iq,iw)+corr_tmp(2,u,iw)*sfacqq
                  corr_cv(3,iq,iw)=corr_cv(3,iq,iw)+corr_tmp(3,u,iw)*sfacqq
               enddo
            enddo
            !$omp end parallel do
         enddo
         corr_tmp=corr_cv
         corr_cv=0.0d0
         do iw=1, nw
            !$omp parallel do default(shared) private(iq,j,t,sfacww)
            do iq=1,Nq
               if (iw==1.and.iq==1) write(*,*) '-------------------------> CONV_GY',facw1,facw2
               ! Convolution with Gaussian resolution function
               do j = max(1,iw-conv_range_w), min(nw-1,iw+conv_range_w)
                  ww=(j-iw)
                  sfacww=exp(facw1*ww**2)*facw2 ! This is the parameter that controls the convolution (Should the sum be over the absolute values?)
                  ! Convolution over frequency space with a gaussian resolution function
                  corr_cv(1,iq,iw)=corr_cv(1,iq,iw)+corr_tmp(1,iq,j)*sfacww
                  corr_cv(2,iq,iw)=corr_cv(2,iq,iw)+corr_tmp(2,iq,j)*sfacww
                  corr_cv(3,iq,iw)=corr_cv(3,iq,iw)+corr_tmp(3,iq,j)*sfacww
               enddo
            enddo
            !$omp end parallel do
         enddo
      else if (do_conv=='LY') then
         do iq=1,Nq
            !$omp parallel do default(shared) private(iw,u,qq,sfacqq)
            do iw=1, nw/2
               if (iw==1.and.iq==1) write(*,*) '-------------------------> CONV_LY'
               ! Convolution with Gaussian resolution function
               do u = max(1,iq-conv_range_q), min(nq-1,iq+conv_range_q)
                  qq=(iq-u)
                  sfacqq=facq1/(LQfactor**2 +qq**2) ! This parameter controls the convolution
                  ! Convolution over the qpoints
                  corr_cv(1,iq,iw)=corr_cv(1,iq,iw)+corr_tmp(1,iq,iw)*sfacqq
                  corr_cv(2,iq,iw)=corr_cv(2,iq,iw)+corr_tmp(2,iq,iw)*sfacqq
                  corr_cv(3,iq,iw)=corr_cv(3,iq,iw)+corr_tmp(3,iq,iw)*sfacqq
               enddo
            enddo
            !$omp end parallel do
         enddo
         do iq=1,Nq
            !$omp parallel do default(shared) private(iw,j,t,sfacww)
            do iw=1, Nw/2
               if (iw==1.and.iq==1) write(*,*) '--------------------------> CONV_LY'
               do j = max(1,iw-conv_range_w), min(nw-1,iw+conv_range_w)
                  t=(iw-j)
                  sfacww=facw1/(LWfactor**2+t**2)
                  corr_cv(1,iq,iw)=corr_cv(1,iq,iw)+corr_tmp(1,iq,iw)*sfacww
                  corr_cv(2,iq,iw)=corr_cv(2,iq,iw)+corr_tmp(2,iq,iw)*sfacww
                  corr_cv(3,iq,iw)=corr_cv(3,iq,iw)+corr_tmp(3,iq,iw)*sfacww
               enddo
            enddo
            !$omp end parallel do
         enddo
      else if (do_conv=='LW') then
         do iq=1,Nq
            !$omp parallel do default(shared) private(iw,j,t,sfacww)
            do iw=1, Nw/2
               if (iw==1.and.iq==1) write(*,*) '--------------------------> CONV_LW'
               do j = max(1,iw-conv_range_w), min(nw-1,iw+conv_range_w)
                  t=(iw-j)
                  sfacww=facw1/(LWfactor**2+t**2)
                  corr_cv(1,iq,iw)=corr_cv(1,iq,iw)+corr_tmp(1,iq,iw)*sfacww
                  corr_cv(2,iq,iw)=corr_cv(2,iq,iw)+corr_tmp(2,iq,iw)*sfacww
                  corr_cv(3,iq,iw)=corr_cv(3,iq,iw)+corr_tmp(3,iq,iw)*sfacww
               enddo
            enddo
            !$omp end parallel do
         enddo
      else if (do_conv=='LQ') then
         do iq=1,Nq
            !$omp parallel do default(shared) private(iw,u,qq,sfacqq)
            do iw=1, nw/2
               if (iw==1.and.iq==1) write(*,*) '-------------------------> CONV_LQ'
               ! Convolution with Gaussian resolution function
               do u = max(1,iq-conv_range_q), min(nq-1,iq+conv_range_q)
                  qq=(iq-u)
                  sfacqq=facq1/(LQfactor**2 +qq**2) ! This parameter controls the convolution
                  ! Convolution over the qpoints
                  corr_cv(1,iq,iw)=corr_cv(1,iq,iw)+corr_tmp(1,iq,iw)*sfacqq
                  corr_cv(2,iq,iw)=corr_cv(2,iq,iw)+corr_tmp(2,iq,iw)*sfacqq
                  corr_cv(3,iq,iw)=corr_cv(3,iq,iw)+corr_tmp(3,iq,iw)*sfacqq
               enddo
            enddo
            !$omp end parallel do
         enddo
      endif

      i_all=-product(shape(corr_tmp))*kind(corr_tmp)
      deallocate(corr_tmp,stat=i_stat)
      call memocc(i_stat,i_all,'corr_tmp','calc_tmp')

      return

   end subroutine calc_dos_conv

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_mavrg_vec
   !> @brief Calculate current average magnetization by components for the connected \f$\mathbf{S}\left(\mathbf{q},\omega\right)\f$
   !-----------------------------------------------------------------------------
   subroutine calc_mavrg_vec(Natom, Mensemble, emomM, mavrg_vec, mavg_axis)
      !

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM !< Current unit moment vector
      real(dblprec), dimension(3), intent(out) :: mavrg_vec !< Current average magnetization vector
      real(dblprec), dimension(3,Mensemble), intent(out) :: mavg_axis !< Current average, ensemble specific, magnetization vector

      !.. Scalar variables
      integer :: i,k

      !.. Local arrays
      real(dblprec), dimension(3,Mensemble) ::  m

      !.. Executable statements
      m=0.0d0

      do k=1,Mensemble
#if _OPENMP >= 201307 && __INTEL_COMPILER < 1800
         !$omp parallel do private(i) default(shared) schedule(static) reduction(+:m)
#endif
         do i=1, Natom
            m(:,k) = m(:,k) + emomM(:,i,k)
         end do
#if _OPENMP >= 201307 && __INTEL_COMPILER < 1800
         !$omp end parallel do
#endif
         mavg_axis(:,k)=m(:,k)/Natom
      end do
      mavrg_vec(1)=sum(m(1,:))/Mensemble/Natom
      mavrg_vec(2)=sum(m(2,:))/Mensemble/Natom
      mavrg_vec(3)=sum(m(3,:))/Mensemble/Natom

   end subroutine calc_mavrg_vec

   !-----------------------------------------------------------------------------
   ! SUBROUTINE allocate_deltatcorr
   !> @brief Allocate adaptivetime step for correlation
   !-----------------------------------------------------------------------------
   subroutine allocate_deltatcorr(allocate_flag)

      implicit none

      logical, intent(in) :: allocate_flag                            !< Allocate/deallocate

      integer :: i_stat, i_all

      if(allocate_flag) then
         if(.not.allocated(deltat_corr)) then
            allocate(deltat_corr(sc_nstep+1),stat=i_stat)
            call memocc(i_stat,product(shape(deltat_corr))*kind(deltat_corr),'deltat_corr','allocate_deltatcorr')
         end if
         if(.not.allocated(scstep_arr)) then
            allocate(scstep_arr(sc_nstep+1),stat=i_stat)
            call memocc(i_stat,product(shape(scstep_arr))*kind(scstep_arr),'scstep_arr','allocate_deltatcorr')
         end if
      else
         if(allocated(deltat_corr)) then
            i_all=-product(shape(deltat_corr))*kind(deltat_corr)
            deallocate(deltat_corr,stat=i_stat)
            call memocc(i_stat,i_all,'deltat_corr','allocate_deltatcorr')
         end if
         if(allocated(scstep_arr)) then
            i_all=-product(shape(scstep_arr))*kind(scstep_arr)
            deallocate(scstep_arr,stat=i_stat)
            call memocc(i_stat,i_all,'scstep_arr','allocate_deltatcorr')
         end if
      end if

   end subroutine allocate_deltatcorr

   !-----------------------------------------------------------------------------
   ! SUBROUTINE calc_bgkt
   !> @brief Calculate bimagnon correlation function
   !-----------------------------------------------------------------------------
   subroutine calc_bgkt(Natom, Mensemble, coord, simid, emomM, sc_tidx, sc_nstep, &
         flag,q_flag,deltat_corr,do_connected,do_conv,sigma_q,sigma_w,LQfactor,LWfactor,max_no_neigh,conf_num, &
         nlist,nlistsize,ncoup)
      !
      use Constants
      use BLS, only : gramms
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      integer, intent(inout) :: sc_tidx ! < Sampling time index
      integer, intent(in) :: sc_nstep !< Number of steps to sample
      integer, intent(inout) :: flag  !< Setup, sample, or print
      character, intent(in) :: q_flag
      character(len=2), intent(in) :: do_conv
      real(dblprec), intent(in) :: sigma_q, sigma_w
      real(dblprec), intent(in) :: LQfactor,LWfactor
      real(dblprec), dimension(sc_nstep+1), intent(in) :: deltat_corr
      character(len=1), intent(in) :: do_connected
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, intent(in) :: conf_num   !< Number of configurations for LSF
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
      integer :: iq,iw,step,r,l,i_stat,i_all,j
      character(len=30) :: filn
      complex(dblprec) :: epowqr, i, iqfac,tt, epowwt
      real(dblprec) :: qdr,nainv,mavg_norm
      real(dblprec), dimension(3) :: mavrg_vec
      real(dblprec), dimension(sc_nstep+1) :: dt
      real(dblprec), dimension(3,3) :: local_rotmat
      !
      i=(0.0d0,1.0d0)

      if(flag==0) then
         ! First call, allocate and clear arrays
         !
         allocate(corr_bkt(nq,sc_nstep+1),stat=i_stat)
         call memocc(i_stat,product(shape(corr_bkt))*kind(corr_bkt),'corr_bkt','calc_bgkt')
         corr_bkt=0.0d0
         !
         !
         if (do_sc_local_axis=='Y') then
            mavg_local_axis=0.0d0
            mavg_local_axis(1,:,:)=emomM(1,1,1)
            mavg_local_axis(2,:,:)=emomM(2,1,1)
            mavg_local_axis(3,:,:)=emomM(3,1,1)
         end if
         !
         mavg_axis=0.0d0
         mort_axis=0.0d0
         !
         allocate(bm_proj(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(bm_proj))*kind(bm_proj),'bm_proj','calc_bgkt')
         !
         flag=1
         sc_samp_done=0
         qfac=2.0d0*pi
      end if

      nainv=1.0d0/Natom

      ! Calculate b(k) for the current iteration and add to b(k,t)
      if (flag==1) then
         iqfac=i*qfac
         call calc_mavrg_vec(Natom,Mensemble,emomM,mavrg_vec,mavg_axis)
         do l=1,Mensemble
            mavg_axis(:,l)=mavg_axis(:,l)/Natom
            mavg_norm=sum(mavg_axis(:,l)*mavg_axis(:,l))**0.5d0
            if(mavg_norm>1.0d-2) then
               mavg_axis(:,l)=mavg_axis(:,l)/mavg_norm
            else
               mavg_axis(:,l)=(/0.0d0,0.0d0,1.0d0/)
            end if
         end do

         if(do_connected/='Y') then
            ! Keep the average magnetization for subtraction below if do_connected
            ! otherwise put it to zero (in order to remove a number of if-statements in-loop)
            mavrg_vec=0.0d0
         end if

         if(do_sc_local_axis/='Y') then
            !$omp parallel do default(shared) private(l,r) schedule(static) collapse(2)
            do l=1,Mensemble
               do r=1,Natom
                  bm_proj(1,r,l)=emomM(1,r,l)-mavrg_vec(1)
                  bm_proj(2,r,l)=emomM(2,r,l)-mavrg_vec(2)
                  bm_proj(3,r,l)=emomM(3,r,l)-mavrg_vec(3)
               end do
            end do
            !$omp end parallel do
         else
            !$omp parallel do default(shared) private(l,r,local_rotmat) schedule(static) collapse(2)
            do l=1,Mensemble
               do r=1,Natom
                  local_rotmat=mavg_local_rotmat(:,:,r,l)
                  bm_proj(1,r,l)=local_rotmat(1,1)*emomM(1,r,l)+local_rotmat(2,1)*emomM(2,r,l)+local_rotmat(3,1)*emomM(3,r,l)
                  bm_proj(2,r,l)=local_rotmat(1,2)*emomM(1,r,l)+local_rotmat(2,2)*emomM(2,r,l)+local_rotmat(3,2)*emomM(3,r,l)
                  bm_proj(3,r,l)=local_rotmat(1,3)*emomM(1,r,l)+local_rotmat(2,3)*emomM(2,r,l)+local_rotmat(3,3)*emomM(3,r,l)
               end do
            end do
            !$omp end parallel do
         end if

         !$omp parallel do default(shared) private(r,iq,l,j,qdr,epowqr) schedule(static)
         do iq=1,nq
            do l=1,Mensemble
               do r=1,Natom
                  qdr=q(1,iq)*coord(1,r)+q(2,iq)*coord(2,r)+q(3,iq)*coord(3,r)
                  epowqr=exp(iqfac*qdr)*nainv
                  do j=1,nlistsize(r)
                     corr_bkt(iq,sc_tidx)=corr_bkt(iq,sc_tidx)-epowqr*ncoup(j,r,1)*(bm_proj(1,r,l)*bm_proj(1,nlist(j,r),l)+ &
                        bm_proj(2,r,l)*bm_proj(2,nlist(j,r),l)+bm_proj(3,r,l)*bm_proj(3,nlist(j,r),l))
                  enddo
               end do
            end do
         end do
         !$omp end parallel do
      end if

      ! Final operations, transform and print
      if (flag==2) then

         ! Allocate arrays
         !
         allocate(corr_bkw(nq,nw),stat=i_stat)
         call memocc(i_stat,product(shape(corr_bkw))*kind(corr_bkw),'corr_bkw','calc_bgkt')
         corr_bkw=0.0d0

         ! Finish sampling and transform (k,t)->(k,w)
         if(sc_tidx.GT.sc_nstep) then
            sc_tidx = sc_nstep
         end if

         j=1
         do while(j.LE.sc_tidx)
            dt(j) = scstep_arr(j)*deltat_corr(j)
            j = j+1
         end do

         wfac=1.0d0
         corr_bkw=0.0d0

         !$omp parallel do default(shared) private(iw,iq,step,tt,epowwt) schedule(static)
         do iw=1,nw
            do step=1,sc_tidx
               tt=i*wfac*(step-1)*dt(step)
               ! Transform with windowing
               epowwt=exp(w(iw)*tt)*sc_window_fac(sc_window_fun,step,sc_tidx)
               do iq=1,nq
                  corr_bkw(iq,iw)=corr_bkw(iq,iw)+epowwt*corr_bkt(iq,step)
               enddo
            enddo
         enddo
         !$omp end parallel do


         !	! Write B(q,t)
         write (filn,'(''bsqt.'',a8,''.out'')') simid
         open(ofileno, file=filn)
         do iq=1,Nq
            do step=1,sc_tidx
               write (ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),step, &
                  real(corr_bkt(iq,step)),aimag(corr_bkt(iq,step)),abs(corr_bkt(iq,step))
            end do
         end do
         close(ofileno)

         ! Write B(q,w)
         write (filn,'(''bsqw.'',a8,''.out'')') simid
         open(ofileno, file=filn)
         write (filn,'(''bsqwx.'',a8,''.out'')') simid
         open(ofileno2, file=filn)
         do iq=1,Nq
            do iw=1,Nw
               write (ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, &
                  real(corr_bkw(iq,iw)),aimag(corr_bkw(iq,iw)),abs(corr_bkw(iq,iw))
               write (ofileno2,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, &
                  abs(real(corr_bkw(iq,iw))),abs(aimag(corr_bkw(iq,iw))),abs(corr_bkw(iq,iw))
            end do
         end do
         close(ofileno)
         close(ofileno2)

         if (do_conv.eq.'LQ' .or. do_conv.eq.'LW' .or. do_conv.eq.'LW' .or. &
            do_conv.eq.'GQ' .or. do_conv.eq.'GW' .or. do_conv.eq.'GY') then
            ! Calculate the convolution for the magnon DOS
            call calc_dos_conv(corr_bkw)
         endif


         ! Calculate and write the magnon DOS
         allocate(bimagnon_dos(nw/2),stat=i_stat)
         call memocc(i_stat,product(shape(bimagnon_dos))*kind(bimagnon_dos),'bimagnon_dos','calc_bgkt')
         !
         write (filn,'(''bswdos.'',a8,''.out'')') simid
         open(ofileno, file=filn)
         bimagnon_dos=0.0d0
         if (q_flag=='I'.or.q_flag=='B') then
            do iw=1,Nw/2
               do iq=1,Nq
                  bimagnon_dos(iw)=bimagnon_dos(iw)+abs(corr_bkw(iq,iw))*q_weight(iq) ! Put a weight for the IR BZ points
               end do
               write (ofileno,10006) hbar_mev*w(iw),bimagnon_dos(iw)
            end do
         else
            do iw=1,Nw/2
               do iq=1,Nq
                  bimagnon_dos(iw)=bimagnon_dos(iw)+abs(corr_bkw(iq,iw)) ! Put a weight for the IR BZ points
               end do
               write (ofileno,10006) hbar_mev*w(iw),bimagnon_dos(iw)
            end do
         endif
         close(ofileno)

         ! Deallocate arrays
         i_all=-product(shape(corr_bkt))*kind(corr_bkt)
         deallocate(corr_bkt,stat=i_stat)
         call memocc(i_stat,i_all,'corr_bkt','calc_bgkt')
         !
         i_all=-product(shape(corr_bkw))*kind(corr_bkw)
         deallocate(corr_bkw,stat=i_stat)
         call memocc(i_stat,i_all,'corr_bkw','calc_bgkt')
         !
         i_all=-product(shape(bimagnon_dos))*kind(bimagnon_dos)
         deallocate(bimagnon_dos,stat=i_stat)
         call memocc(i_stat,i_all,'bimagnon_dos','calc_bgkt')
         !
         i_all=-product(shape(bm_proj))*kind(bm_proj)
         deallocate(bm_proj,stat=i_stat)
         call memocc(i_stat,i_all,'bm_proj','calc_bgkt')
         !
      end if
      return
      !
      10005 format (i7,3f10.6,2x,i7,2x,9es16.8)
      10006 format (5G16.8)
      !
   end subroutine calc_bgkt

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: find_local_rotmat
   !> @brief Finding the local rotational matrix for the local quantization axis
   !-----------------------------------------------------------------------------
   subroutine find_local_rotmat(ldim,mom_in,mat_out)
      !
      !
      implicit none
      !
      integer, intent(in) :: ldim
      real(dblprec),dimension(3,ldim), intent(in) :: mom_in
      real(dblprec),dimension(3,3,ldim), intent(out) :: mat_out
      !
      integer :: i
      real(dblprec) :: theta, m_norm, k_norm
      real(dblprec),dimension(3) :: k_vec, z_vec, m_vec
      real(dblprec),dimension(3,3) :: K_mat, KK_mat, eye
      !
      z_vec(1)=0.0d0
      z_vec(2)=0.0d0
      z_vec(3)=1.0d0
      !
      eye=0.0d0;eye(1,1)=1.0d0;eye(2,2)=1.0d0;eye(3,3)=1.0d0
      !
      !$omp parallel do default(shared) private(i,m_norm,m_vec,k_vec,k_norm,theta,K_mat,KK_mat)
      do i=1,ldim
         m_norm=sum(mom_in(:,i)*mom_in(:,i))**0.5d0
         m_vec=mom_in(:,i)/m_norm
         k_vec(1)=m_vec(2)*z_vec(3)-m_vec(3)*z_vec(2)
         k_vec(2)=m_vec(3)*z_vec(1)-m_vec(1)*z_vec(3)
         k_vec(3)=m_vec(1)*z_vec(2)-m_vec(2)*z_vec(1)
         k_norm=sum(k_vec*k_vec)**0.5d0+1.0d-14
         k_vec=k_vec/k_norm
         !
         theta=-acos(m_vec(1)*z_vec(1)+m_vec(2)*z_vec(2)+m_vec(3)*z_vec(3))
         !
         K_mat(1,1)=0.0d0;K_mat(2,2)=0.0d0;K_mat(3,3)=0.0d0
         K_mat(1,2)=-k_vec(3);K_mat(2,1)= k_vec(3)
         K_mat(1,3)= k_vec(2);K_mat(3,1)=-k_vec(2)
         K_mat(2,3)=-k_vec(1);K_mat(3,2)= k_vec(1)
         !
         KK_mat=matmul(K_mat,K_mat)
         !
         mat_out(:,:,i)=eye+K_mat*sin(theta)+KK_mat*(1.0d0-cos(theta))
         !
      end do
      !$omp end parallel do
      !
      !
      return
      !
   end subroutine find_local_rotmat

   !---------------------------------------------------------------------------
   !> @brief
   !> Read input parameters.
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_parameters_correlation(ifile)
      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword, cache
      integer :: rd_len, i_err, i_errb
      logical :: comment

      do
         10     continue
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
            !> - simid
            ! This is the flags for the S(q,w)
            !> - do_sc
            !! Calculate correlation (N=no/C=static/Q=dynamic)
            !! Dynamical structure factor measures the time and space correlation function
            !! \f$C^k (\mathbf{r}-\mathbf{r'},t) = \langle m_{\mathbf{r}^k}(t) m_{\mathbf{r'}^k}(0) \rangle -
            !! \langle m_{\mathbf{r}^k}(t) \rangle \langle m_{\mathbf{r'}^k}(0) \rangle\f$
            !! ,where the angular brackets signify an ensemble average and k the Cartesian component, and its Fourier Transform, the dynamical structure factor
            !! \f$S^k(\mathbf{q},\omega) = \frac{1}{\sqrt{2\pi}N} \sum_{\mathbf{r},\mathbf{r'}} e^{i\mathbf{q}\cdot(\mathbf{r}-\mathbf{r'})}
            !! \int_{-\infty}^{\infty} e^{i \omega t} C^k (\mathbf{r}-\mathbf{r'},t) dt, \f$
            case('do_sc')
               read(ifile,*,iostat=i_err) do_sc
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               !> - do_sc_bimag
               !! Calculate bimagnon form factor (Y/N)
            case('do_sc_bimag')
               read(ifile,*,iostat=i_err) do_sc_bimag
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               !> - do_sc_local_axis
               !! Resolve SQW in transversal and longitudinal components (Y/N)
            case('do_sc_local_axis')
               read(ifile,*,iostat=i_err) do_sc_local_axis
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               !> - sc_local_axis_mix
               !! How much should the local axis be updated dynamically
            case('sc_local_axis_mix')
               read(ifile,*,iostat=i_err) sc_local_axis_mix
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               !> - do_sc_proj_axis
            case('do_sc_proj_axis')
               read(ifile,*,iostat=i_err) do_sc_proj_axis
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               !> - do_sc_dosonly
               !! Magnon density of states without full SQW (Y/N)
            case('do_sc_dosonly')
               read(ifile,*,iostat=i_err) do_sc_dosonly
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               !> - do_sc_complex
               !! Print the complex valued S(q,w)
            case('do_sc_complex')
               read(ifile,*,iostat=i_err) do_sc_complex
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_sr')
               read(ifile,*,iostat=i_err) do_sr
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('sc_window_fun')
               read(ifile,*,iostat=i_err) sc_window_fun
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

               ! Flag for doing the convolutions in frequencies and qpoints can be in Gaussian (G)
               ! or Lorentzian (L) broadening
               ! ((G/L)Q==qpoints, (G/L)W==frequencies, (G/L)Y==Both)
            case('do_conv')
               read(ifile,*,iostat=i_err) do_conv
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

               ! If the Convolution is turned on one must read the sigma parameters
            case('sigma_q')
               read(ifile,*,iostat=i_err) sigma_q ! This parameter is for the gaussian in reciprocal space
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('sigma_w')
               read(ifile,*,iostat=i_err) sigma_w ! This parameter is for the gaussian in frequency
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('lorentz_q')
               read(ifile,*,iostat=i_err) LQfactor ! This parameter is for the Lorentzian in reciprocal space
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('lorentz_w')
               read(ifile,*,iostat=i_err) LWfactor ! This parameter is for the Lorentzian in frequency
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_sc_proj')
               read(ifile,*,iostat=i_err) do_sc_proj
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_sc_projch')
               read(ifile,*,iostat=i_err) do_sc_projch
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_qt_traj')
               read(ifile,*,iostat=i_err) do_qt_traj
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('sc_mode')
               read(ifile,*,iostat=i_err) sc_mode
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               !> - sc_step
               !! Sampling frequency for SQW
            case('sc_step')
               read(ifile,*,iostat=i_err) sc_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               !> - sc_nstep
               !! Number of frequencies in SQW
            case('sc_nstep')
               read(ifile,*,iostat=i_err) sc_nstep
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               !> - sc_sep
               !! Sampling period of static G(r)
            case('sc_sep')
               read(ifile,*,iostat=i_err) sc_sep
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               !> - qpoints
               !! Specify format of qfile for correlation (F=file in cart. coord,
               !! D=file in direct coord, C=full BZ + .... )
            case('qpoints')
               read(ifile,*,iostat=i_err) qpoints
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               !> - qfile
               !! Name of qfile for correlation
            case('qfile')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               qfile=trim(adjustl(cache))

            case('do_connected')
               read(ifile,*,iostat=i_err) do_connected
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            end select
         end if

         ! End of file
         if (i_errb==20) goto 20
         ! End of row
         if (i_errb==10) goto 10
      end do

      20  continue

      rewind(ifile)
      return
   end subroutine read_parameters_correlation

   !-----------------------------------------------------------------------------
   ! FUNCTION: sc_window_fac
   !> @brief Window function factor for the different types of windowing
   !> @author Anders Bergman
   !-----------------------------------------------------------------------------
   real(dblprec) function sc_window_fac(sc_window_fun,step,nstep)
      use Constants, only : pi
      !
      implicit none
      !
      integer, intent(in)  :: sc_window_fun
      integer, intent(in)  :: step
      integer, intent(in)  :: nstep
      !
      real(dblprec) :: dum
      !
      dum=1.0d0
      select case(sc_window_fun)
         ! Hann
      case(2)
         dum= (0.50d0-0.50d0*cos(2.0d0*pi*(step-1.d0)/(nstep-1.d0)))
         ! Hamming
      case(3)
         dum= (0.54d0-0.46d0*cos(2.0d0*pi*(step-1.d0)/(nstep-1.d0)))
         ! Hamming v2
      case(32)
         dum= (0.53836d0- 0.46164d0*cos(2.0d0*pi*(step-1.d0)/(nstep-1.d0)))
         ! Blackman-Harris
      case(4)
         dum=  &
         (0.35785d0-0.48829d0*cos(2.0d0*pi*(step-1.d0)/(nstep-1.0d0))+ &
         0.14128d0*cos(4.0d0*pi*(step-1.d0)/(nstep-1.0d0))   &
         -0.01168d0*cos(6.0d0*pi*(step-1.d0)/(nstep-1.0d0)))
         ! Nuttal
      case(5)
         dum=  &
         (0.355768d0-0.478396d0*cos(2.0d0*pi*(step-1.d0)/(nstep-1.0d0))+ &
         0.144232d0*cos(4.0d0*pi*(step-1.d0)/(nstep-1.0d0))   &
         -0.012604d0*cos(6.0d0*pi*(step-1.d0)/(nstep-1.0d0)))
         ! Square windows
      case default
         dum=1.0d0
      end select
      !
      sc_window_fac=dum
      return
      !
   end function sc_window_fac

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: find_rmid
   !> @brief Finds the center of the sample
   !-----------------------------------------------------------------------------
   subroutine find_rmid(rmid,coord,Natom)
      !
      implicit none
      !
      real(dblprec), dimension(3), intent(out) :: rmid
      integer, intent(in) :: Natom
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      !
      integer i,ridx
      real(dblprec) :: rnorm2_i, rnorm2_j
      real(dblprec), dimension(3) :: rcenter
      !
      !
      !
      rcenter(1)=(maxval(coord(1,:))+minval(coord(1,:)))*0.5d0
      rcenter(2)=(maxval(coord(2,:))+minval(coord(2,:)))*0.5d0
      rcenter(3)=(maxval(coord(3,:))+minval(coord(3,:)))*0.5d0
      !
      ridx=1
      rnorm2_j=sum((coord(:,ridx)-rcenter)**2)
      do i=2, Natom
         rnorm2_i=sum((coord(:,i)-rcenter)**2)
         if(rnorm2_i<rnorm2_j) then
            ridx=i
            rnorm2_j=sum((coord(:,ridx)-rcenter)**2)
         end if
      end do
      !
      rmid=coord(:,ridx)
      !
   end subroutine find_rmid

end module Correlation
