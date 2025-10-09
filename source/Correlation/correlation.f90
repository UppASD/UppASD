!------------------------------------------------------------------------------------
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
!> @todo Automatic generation of q-points for do_sc="Y"
!------------------------------------------------------------------------------------
module Correlation
   use Parameters
   use Profiling
   use Correlation_print
   use Correlation_core
   use Qvectors
   use Omegas
   !use BiMagnonCorr
   !
   implicit none
   !

   type(corr_t)  :: sc   !< Spin-correlation 
   type(corr_t)  :: uc   !< Lattice displacement-correlation 
   type(corr_t)  :: vc   !< Velocity displacement-correlation 
   type(corr_t)  :: lc   !< Lattice angular momentum  displacement-correlation 
   !type(corr_t)  :: suc   !< Velocity displacement-correlation 
   !type(corr_t)  :: uvc   !< Velocity displacement-correlation 

   !Spin-spin correlation parameters from input file
   character(len=1) :: do_sc        !< Perform spin-correlation sampling (Y/N/C)
   character(len=1) :: do_sr        !< Perform spin-correlation sampling in real space directly (Y/N)

   character(len=1) :: do_uc        !< Perform displacement-correlation sampling (Y/N/C)
   character(len=1) :: do_ur        !< Perform displacement-correlation sampling in real space directly (Y/N)

   character(len=1) :: do_vc        !< Perform velocity-correlation sampling (Y/N/C)
   character(len=1) :: do_vr        !< Perform velocity-correlation sampling in real space directly (Y/N)

   character(len=1) :: do_lc        !< Perform L-correlation sampling (Y/N/C)
   character(len=1) :: do_lr        !< Perform L-correlation sampling in real space directly (Y/N)

   character(len=1) :: do_suc        !< Perform velocity-correlation sampling (Y/N/C)
   !character(len=1) :: do_sur        !< Perform velocity-correlation sampling in real space directly (Y/N)

   character(len=1) :: do_uvc        !< Perform velocity-correlation sampling (Y/N/C)
   !character(len=1) :: do_uvr        !< Perform velocity-correlation sampling in real space directly (Y/N)

   character(len=1) :: do_slc        !< Perform velocity-correlation sampling (Y/N/C)

   !-----------------------
   character(len=1) :: do_bimag  !< Perform spin-correlation sampling of bi-magnons (Y/N) (not implemented currently)

   ! CPU/GPU control flag
   character(len=1) :: do_gpu_correlations !< Perform correlations on GPU (Y/N)
   ! Working variables to perform the printing of the correlation

   private

   public :: correlation_init, read_parameters_correlation
   public :: read_parameters_correlation_sc, correlation_init_loc
   public :: read_parameters_correlation_uc
   public :: read_parameters_correlation_vc
   public :: read_parameters_correlation_lc
   public :: correlation_wrapper,  set_correlation_average
   public :: do_sc, do_uc, do_sr, do_ur, do_vc, do_vr, do_uvc, do_suc, do_slc, do_lc, do_lr
   public :: sc, uc, vc, lc!, uvc, suc

contains

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: correlation_wrapper
   !> @brief
   !> Driver for correlation function calculations
   !---------------------------------------------------------------------------------
   subroutine correlation_wrapper(Natom,Mensemble,coord,simid,idata,mstep,delta_t,  &
         NT,atype,Nchmax,achtype,cc, do_cc, do_cr, flag)
      !
      implicit none
      !
      integer, intent(in) :: NT           !< Number of types of atoms
      integer, intent(in) :: mstep        !< Current simulation step
      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Nchmax       !< Number of chemical types
      integer, intent(in) :: Mensemble    !< Number of ensembles
      real(dblprec), intent(in) :: delta_t      !< Time step
      character(len=8), intent(in) :: simid     !< Name of simulation
      integer, dimension(Natom), intent(in) :: atype        !< Type of atom
      integer, dimension(Natom), intent(in) :: achtype      !< Type of atom
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: idata  !< Measurable to be sampled
      type(corr_t), intent(inout) :: cc
      character(len=1) :: do_cc        !< Calculate reciprocal correlations
      character(len=1) :: do_cr        !< Calculate real-space correlations
      integer, intent(inout) :: flag  !< Setup, sample, or print

      !
      integer :: cr_flag, gk_flag, gkt_flag
      !
      !if(cc%gkt_flag==0) then
      !   cc%sc_tidx=0
      !   cc%prefix='s'
      !end if

      !if(flag==2) print *,'------------------Correlation-analysis-----------------------', cc%label
      !print *,'Correlation', cc%label, flag,do_cc, do_cr
      cr_flag = flag
      gk_flag = flag
      gkt_flag = flag

      ! Sample S(r) through S(q)
      if(do_cc=='C'.or.do_cc=='Y') then

         ! AB TODO: Reinstate projected functionality
         !if(mod(mstep-1,cc%sc_sep)==0.or.flag==2) call calc_gk2(Natom, Mensemble, cc, coord, simid, idata,  mstep, flag)
         if(mod(mstep-1,cc%sc_sep)==0.or.flag==2)   call calc_gk2(Natom, Mensemble,NT,atype,Nchmax,achtype, cc, coord, simid, idata, gk_flag)

      end if
      !else if(do_cc=='Q') then
      if(do_cc=='Q'.or.do_cc=='T'.or.do_cc=='Y') then

         ! S(k,w)
         if(mod(mstep-1,cc%sc_step)==1.and.cc%sc_tidx<cc%sc_max_nstep) then

            if(cc%gkt_flag==0) call allocate_deltatcorr(.true.,cc)

            cc%sc_tidx=cc%sc_tidx+1
            cc%deltat_corr(cc%sc_tidx) = delta_t ! Save current time step
            cc%scstep_arr(cc%sc_tidx) = cc%sc_step  ! Save current sampling period

            call calc_gkt(Natom, Mensemble,NT,atype,Nchmax,achtype, cc, coord, idata, gkt_flag)

         else if (flag==2) then

            call calc_gkt(Natom, Mensemble,NT,atype,Nchmax,achtype, cc, coord, idata, gkt_flag)

            if(do_cc=='Q'.or.do_cc=='Y')  then
               call print_gkw(NT, Nchmax, cc, cc, simid, cc%label)
            end if

            if(do_cc=='T'.or.do_cc=='Y')  then
               call print_gkt(NT, Nchmax, cc, cc, simid, cc%label)
            end if

         else if (flag==3) then

            call deallocate_gkw(cc)

         end if
      end if

      ! Sample S(r) directly for non periodic systems
      if (do_cr=='Y') then
         !print *,'call calc_sr',cr_flag
         if (mod(mstep-1,cc%sc_sep)==0.or.cr_flag==2) call calc_sr(Natom, Mensemble, cc, simid, idata, cr_flag)
      endif
      !
      if(flag==2) then
         call allocate_deltatcorr(.false.,cc)
      end if

      flag = max(cr_flag,max(gk_flag,gkt_flag))

   end subroutine correlation_wrapper

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: correlation_init
   !> @brief
   !> Driver for correlation function calculations
   !---------------------------------------------------------------------------------
   subroutine correlation_init()

      implicit none


      
      ! General control flags
      do_sc        = 'N'
      do_sr        = 'N'
      do_uc        = 'N'
      do_ur        = 'N'
      do_vc        = 'N'
      do_vr        = 'N'
      do_lc        = 'N'
      do_lr        = 'N'
      do_suc       = 'N'
      !do_sur       = 'N'
      do_uvc       = 'N'
      !do_uvr       = 'N'
      do_bimag     = 'N'
      do_conv      = 'N'
      sc_window_fun = 1

      ! General convolution variables
      sigma_q      = 0.0_dblprec
      sigma_w      = 0.0_dblprec
      LWfactor     = 0.0_dblprec
      LQfactor     = 0.0_dblprec

      ! Specific flags 
      sc%label ='s'
      sc%sc_tidx=0
      uc%label ='u'
      uc%sc_tidx=0
      vc%label ='v'
      vc%sc_tidx=0
      lc%label ='l'
      lc%sc_tidx=0

      !suc%label ='su'
      !suc%sc_tidx=0
      !uvc%label ='uv'
      !uvc%sc_tidx=0



   end subroutine correlation_init

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: correlation_init
   !> @brief
   !> Driver for correlation function calculations
   !---------------------------------------------------------------------------------
   subroutine correlation_init_loc(cc)

      implicit none

      type(corr_t), intent(inout) :: cc

      !Spin-spin correlation
      cc%sc_sep       = 100
      cc%sc_step      = 1
      cc%sc_nstep     = 1
      cc%do_sc_local_axis   = 'N'
      cc%do_sc_dosonly   = 'N'
      cc%do_sc_complex   = 'N'
      cc%do_sc_tens   = 'N'
      cc%do_sc_list   = 'N'
      cc%do_qt_traj   = 'N'
      cc%sc_average = 'N'
      cc%sc_local_axis_mix = 0.0_dblprec
      cc%do_connected = 'N'
      cc%sc_naverages = 1
      cc%gkt_flag = 0
      cc%gk_flag = 0

      ! Real-space correlation cut-off
      cc%cr_cut = 1.0

      ! Energy settings
      cc%sc_emax = -1.0_dblprec
      cc%sc_eres = -1.0_dblprec

      ! Set mode according to lenght of label
      ! (len=1 : single correlation i.e. <m*m>)
      ! (len=2 : double correlation i.e. <m*u>)
      cc%cmode = len(trim(cc%label))

   end subroutine correlation_init_loc



   !---------------------------------------------------------------------------------
   !> @brief
   !> Read input parameters.
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
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

         case('do_uc')
            read(ifile,*,iostat=i_err) do_uc
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_ur')
            read(ifile,*,iostat=i_err) do_ur
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_vc')
            read(ifile,*,iostat=i_err) do_vc
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_vr')
            read(ifile,*,iostat=i_err) do_vr
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_lc')
            read(ifile,*,iostat=i_err) do_lc
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_lr')
            read(ifile,*,iostat=i_err) do_lr
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_suc')
            read(ifile,*,iostat=i_err) do_suc
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         !case('do_sur')
         !   read(ifile,*,iostat=i_err) do_sur
         !   if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_uvc')
            read(ifile,*,iostat=i_err) do_uvc
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_slc')
            read(ifile,*,iostat=i_err) do_slc
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         !case('do_uvr')
         !   read(ifile,*,iostat=i_err) do_uvr
         !   if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !> - do_bimag
            !! Calculate bimagnon form factor (Y/N)
         case('do_bimag')
            read(ifile,*,iostat=i_err) do_bimag
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


   !---------------------------------------------------------------------------------
   !> @brief
   !> Read input parameters for spin-correlations
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   subroutine read_parameters_correlation_sc(ifile)

      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword
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
            !> - do_sc_local_axis
            !! Resolve SQW in transversal and longitudinal components (Y/N)
         case('do_sc_local_axis')
            read(ifile,*,iostat=i_err) sc%do_sc_local_axis
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - sc_local_axis_mix
            !! How much should the local axis be updated dynamically
         case('sc_local_axis_mix')
            read(ifile,*,iostat=i_err) sc%sc_local_axis_mix
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - do_sc_dosonly
            !! Magnon density of states without full SQW (Y/N)
         case('do_sc_dosonly')
            read(ifile,*,iostat=i_err) sc%do_sc_dosonly
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - do_sc_complex
            !! Print the complex valued S(q,w)
         case('do_sc_complex')
            read(ifile,*,iostat=i_err) sc%do_sc_complex
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - do_sc_tens
            !! Print the tensorial elements of S(q,w)
         case('do_sc_tens')
            read(ifile,*,iostat=i_err) sc%do_sc_tens
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !! Print the max locations of s(q,w)
         case('do_sc_list')
            read(ifile,*,iostat=i_err) sc%do_sc_list
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_sc_proj')
            read(ifile,*,iostat=i_err) sc%do_proj
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_sc_projch')
            read(ifile,*,iostat=i_err) sc%do_projch
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_qt_traj')
            read(ifile,*,iostat=i_err) sc%do_qt_traj
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !> - sc_step
            !! Sampling frequency for SQW
         case('sc_step')
            read(ifile,*,iostat=i_err) sc%sc_step
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - sc_nstep
            !! Number of frequencies in SQW
         case('sc_nstep')
            read(ifile,*,iostat=i_err) sc%sc_nstep
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !! Sampling period of static G(r)
         case('sc_sep')
            read(ifile,*,iostat=i_err) sc%sc_sep
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !! How to perform averaging of S(q,w) (F)ull,(E)ven, or (N)one.
         case('sc_average')
            read(ifile,*,iostat=i_err) sc%sc_average
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('do_sc_connected')
            read(ifile,*,iostat=i_err) sc%do_connected
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            ! Real space correlation cut-off
         case('sc_cr_cut')
            read(ifile,*,iostat=i_err) sc%cr_cut
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            ! Energy cutoff for S(q,w) (in mRy)
         case('sc_emax')
            read(ifile,*,iostat=i_err) sc%sc_emax
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            ! Energy resolution for S(q,w) (in mRy)
         case('sc_eres')
            read(ifile,*,iostat=i_err) sc%sc_eres
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
   end subroutine read_parameters_correlation_sc


   !---------------------------------------------------------------------------------
   !> @brief
   !> Read input parameters for displacement correlations
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   subroutine read_parameters_correlation_uc(ifile)

      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword
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
            !> - do_sc_local_axis
            !! Resolve SQW in transversal and longitudinal components (Y/N)
         case('do_uc_local_axis')
            read(ifile,*,iostat=i_err) uc%do_sc_local_axis
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - sc_local_axis_mix
            !! How much should the local axis be updated dynamically
         case('uc_local_axis_mix')
            read(ifile,*,iostat=i_err) uc%sc_local_axis_mix
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - do_sc_dosonly
            !! Magnon density of states without full SQW (Y/N)
         case('do_uc_dosonly')
            read(ifile,*,iostat=i_err) uc%do_sc_dosonly
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - do_sc_complex
            !! Print the complex valued S(q,w)
         case('do_uc_complex')
            read(ifile,*,iostat=i_err) uc%do_sc_complex
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - do_sc_tens
            !! Print the tensorial elements of S(q,w)
         case('do_uc_tens')
            read(ifile,*,iostat=i_err) uc%do_sc_tens
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !! Print the max locations of s(q,w)
         case('do_uc_list')
            read(ifile,*,iostat=i_err) uc%do_sc_list
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_uc_proj')
            read(ifile,*,iostat=i_err) uc%do_proj
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_uc_projch')
            read(ifile,*,iostat=i_err) uc%do_projch
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_uc_qt_traj')
            read(ifile,*,iostat=i_err) uc%do_qt_traj
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !> - sc_step
            !! Sampling frequency for SQW
         case('uc_step')
            read(ifile,*,iostat=i_err) uc%sc_step
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - sc_nstep
            !! Number of frequencies in SQW
         case('uc_nstep')
            read(ifile,*,iostat=i_err) uc%sc_nstep
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !! Sampling period of static G(r)
         case('uc_sep')
            read(ifile,*,iostat=i_err) uc%sc_sep
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !! How to perform averaging of S(q,w) (F)ull,(E)ven, or (N)one.
         case('uc_average')
            read(ifile,*,iostat=i_err) uc%sc_average
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('do_uc_connected')
            read(ifile,*,iostat=i_err) uc%do_connected
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            ! Real space correlation cut-off
         case('uc_cr_cut')
            read(ifile,*,iostat=i_err) uc%cr_cut
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
   end subroutine read_parameters_correlation_uc

   !---------------------------------------------------------------------------------
   !> @brief
   !> Read input parameters for displacement correlations
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   subroutine read_parameters_correlation_vc(ifile)

      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword
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
            !> - do_sc_local_axis
            !! Resolve SQW in transversal and longitudinal components (Y/N)
         case('do_vc_local_axis')
            read(ifile,*,iostat=i_err) vc%do_sc_local_axis
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - sc_local_axis_mix
            !! How much should the local axis be updated dynamically
         case('vc_local_axis_mix')
            read(ifile,*,iostat=i_err) vc%sc_local_axis_mix
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - do_sc_dosonly
            !! Magnon density of states without full SQW (Y/N)
         case('do_vc_dosonly')
            read(ifile,*,iostat=i_err) vc%do_sc_dosonly
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - do_sc_complex
            !! Print the complex valued S(q,w)
         case('do_vc_complex')
            read(ifile,*,iostat=i_err) vc%do_sc_complex
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - do_sc_tens
            !! Print the tensorial elements of S(q,w)
         case('do_vc_tens')
            read(ifile,*,iostat=i_err) vc%do_sc_tens
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !! Print the max locations of s(q,w)
         case('do_vc_list')
            read(ifile,*,iostat=i_err) vc%do_sc_list
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_vc_proj')
            read(ifile,*,iostat=i_err) vc%do_proj
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_vc_projch')
            read(ifile,*,iostat=i_err) vc%do_projch
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_vc_qt_traj')
            read(ifile,*,iostat=i_err) vc%do_qt_traj
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !> - sc_step
            !! Sampling frequency for SQW
         case('vc_step')
            read(ifile,*,iostat=i_err) vc%sc_step
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - sc_nstep
            !! Number of frequencies in SQW
         case('vc_nstep')
            read(ifile,*,iostat=i_err) vc%sc_nstep
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !! Sampling period of static G(r)
         case('vc_sep')
            read(ifile,*,iostat=i_err) vc%sc_sep
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !! How to perform averaging of S(q,w) (F)ull,(E)ven, or (N)one.
         case('vc_average')
            read(ifile,*,iostat=i_err) vc%sc_average
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('do_vc_connected')
            read(ifile,*,iostat=i_err) vc%do_connected
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            ! Real space correlation cut-off
         case('vc_cr_cut')
            read(ifile,*,iostat=i_err) vc%cr_cut
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
   end subroutine read_parameters_correlation_vc

   !---------------------------------------------------------------------------------
   !> @brief
   !> Read input parameters for displacement correlations
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   subroutine read_parameters_correlation_lc(ifile)

      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword
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
            !> - do_sc_local_axis
            !! Resolve SQW in transversal and longitudinal components (Y/N)
         case('do_lc_local_axis')
            read(ifile,*,iostat=i_err) lc%do_sc_local_axis
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - sc_local_axis_mix
            !! How much should the local axis be updated dynamically
         case('lc_local_axis_mix')
            read(ifile,*,iostat=i_err) lc%sc_local_axis_mix
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - do_sc_dosonly
            !! Magnon density of states without full SQW (Y/N)
         case('do_lc_dosonly')
            read(ifile,*,iostat=i_err) lc%do_sc_dosonly
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - do_sc_complex
            !! Print the complex valued S(q,w)
         case('do_lc_complex')
            read(ifile,*,iostat=i_err) lc%do_sc_complex
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - do_sc_tens
            !! Print the tensorial elements of S(q,w)
         case('do_lc_tens')
            read(ifile,*,iostat=i_err) lc%do_sc_tens
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !! Print the max locations of s(q,w)
         case('do_lc_list')
            read(ifile,*,iostat=i_err) lc%do_sc_list
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_lc_proj')
            read(ifile,*,iostat=i_err) lc%do_proj
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_lc_projch')
            read(ifile,*,iostat=i_err) lc%do_projch
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_lc_qt_traj')
            read(ifile,*,iostat=i_err) lc%do_qt_traj
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !> - sc_step
            !! Sampling frequency for SQW
         case('lc_step')
            read(ifile,*,iostat=i_err) lc%sc_step
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - sc_nstep
            !! Number of frequencies in SQW
         case('lc_nstep')
            read(ifile,*,iostat=i_err) lc%sc_nstep
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !! Sampling period of static G(r)
         case('lc_sep')
            read(ifile,*,iostat=i_err) lc%sc_sep
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !! How to perform averaging of S(q,w) (F)ull,(E)ven, or (N)one.
         case('lc_average')
            read(ifile,*,iostat=i_err) lc%sc_average
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('do_lc_connected')
            read(ifile,*,iostat=i_err) lc%do_connected
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            ! Real space correlation cut-off
         case('lc_cr_cut')
            read(ifile,*,iostat=i_err) lc%cr_cut
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_gpu_correlations')
            read(ifile,*,iostat=i_err) do_gpu_correlations
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
   end subroutine read_parameters_correlation_lc


   subroutine set_correlation_average(cc,nstep)
      !
      implicit none
      !
      type(corr_t), intent(inout) :: cc
      integer, intent(in) :: nstep
      !
      if (cc%sc_average=='E'.or.cc%sc_average=='F') then
         !sc_naverages=nstep-sc_nstep*sc_step
         cc%sc_naverages=max(1,(nstep-cc%sc_nstep*cc%sc_step)/cc%sc_step+1)
      else
         cc%sc_naverages=1
      end if
      cc%sc_max_nstep=cc%sc_nstep+cc%sc_naverages
      !
      return
   end subroutine set_correlation_average


end module Correlation