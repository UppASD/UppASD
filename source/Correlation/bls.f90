!------------------------------------------------------------------------------------
! MODULE: BLS
!> Data and routines for calculate frequency based spin correlation function \f$\mathbf{S}\left(\omega\right)\f$
!> @author
!> A. Bergman
!> @copyright
!> GNU Public License.
!------------------------------------------------------------------------------------
module BLS
   use Parameters
   use Profiling
   !
   implicit none
   !
   integer :: numq   !< Number of q points to be sampled
   integer :: nw_bls !< Number of frequencies to be studied
   integer, dimension(:), allocatable :: bls_tidx  !< Storage for measurement steps
   !
   real(dblprec) :: dt_bls
   real(dblprec), dimension(:), allocatable :: w_bls            !< Frequencies w
   real(dblprec), dimension(:,:), allocatable :: qcoord         !< q-points
   real(dblprec), dimension(:,:), allocatable :: bls_dos        !< DOS
   real(dblprec), dimension(:,:), allocatable :: mavg_axis      !< Internal array for finding the global direction of magnetization
   real(dblprec), dimension(:,:,:), allocatable :: mavg_local_axis!< Internal array for finding the local average direction of magnetization
   real(dblprec), dimension(:,:,:), allocatable :: mort_axis    !< Orthogonal (and parallel ) components of magnetization
   real(dblprec), dimension(:,:,:), allocatable :: magmom_start !< M(w)
   !
   complex(dblprec), dimension(:), allocatable :: expw_bls        !< Array for storing all exponentials needed for transform
   complex(dblprec), dimension(:,:,:,:), allocatable :: magmom_w  !< M(w) in real space
   complex(dblprec), dimension(:,:,:,:), allocatable :: magmom_qw !< M(w) in reciprocal space

   ! BLS simulation
   character(len=1) :: do_bls !< Perform frequency based spin-correlation sampling (Y/N/C)
   character(len=1) :: do_bls_local_axis !< Perform BLS along local quantization axis (Y/N)
   integer :: bls_step        !< Separation between sampling steps
   integer :: bls_nstep       !< Number of steps to sample
   integer :: bls_window      !< Number of measurement windows for the power spectra
   integer :: bls_nsamples    !< Number of samples to be studied in a window
   integer :: bls_sample_step !< Number of steps between measurements in a window
   integer :: bls_window_fun  !< Choice of FFT window function (1=box, 2=Hann, 3=Hamming, 4=Blackman-Harris)
   real(dblprec) :: bls_local_axis_mix !< Mixing for local-axis sampling

   private

   public :: gramms, setup_bls, calc_bls, deallocate_bls_data, read_parameters_bls
   ! public subroutines

contains

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: setup_bls
   !> @brief
   !> Main driver for frequency based correlation function \f$\mathbf{S}\left(\omega\right)\f$ similar to
   !> Brillouin light scattering (BLS) in experiments.
   !---------------------------------------------------------------------------------
   subroutine setup_bls()

      implicit none

      do_bls='N'
      do_bls_local_axis='N'
      bls_step=1
      bls_nstep=1
      bls_window=1
      bls_nsamples=1
      bls_sample_step=1
      bls_window_fun=1
      bls_local_axis_mix=0.0_dblprec

   end subroutine setup_bls

   !---------------------------------------------------------------------------------
   ! SUBROTUINE: calc_bls
   !> Calculate S(q,t) for obtaining S(q,w) after FT
   !---------------------------------------------------------------------------------
   subroutine calc_bls(N1,N2,N3,C1,C2,C3,Natom, Mensemble, simid, coord,emomM, mstep, delta_t, flag)

      use Constants
      use FieldData, only : beff
      !
      implicit none
      !
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      integer, intent(in) :: mstep     !< Current simulation step
      real(dblprec), intent(in) :: delta_t   !< Time step
      integer, intent(in) :: flag  !< Setup, sample, or print
      !
      !
      real(dblprec), dimension(3,Mensemble) :: tot_moment
      integer :: iatom,iens,iw,i_stat, component, isamp, shifted_mstep,iq,r
      character(len=30) :: filn
      complex(dblprec) :: i
      real(dblprec) :: mavg_norm,qfac,nainv,qdr, bls_avg_norm, bls_phase_avg_avg, bls_norm_avg_avg, mm
      real(dblprec), dimension(3) :: m_proj, bls_avg, bls_norm_avg, bls_phase_avg, bls_phase_zero
      real(dblprec), dimension(:,:,:), allocatable :: bls_norm,bls_phase
      complex(dblprec) :: epowqr, eexp
      real(dblprec) :: beff_norm
      !

      if (.not.(do_bls=='Y'.or.do_bls=='D'.or.do_bls=='Q')) return

      i=(0.0_dblprec,1.0_dblprec)

      tot_moment=0.0_dblprec

      qfac=2._dblprec*pi
      nainv=1.0_dblprec/Natom

      if(flag==0) then

         ! First call, allocate and clear arrays
         if (do_bls=='Y'.or.do_bls=='Q') then
            allocate(magmom_w(bls_nstep,3,Natom,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(magmom_w))*kind(magmom_w),'magmom_w','calc_bls')
            magmom_w=0.0_dblprec
            !
            allocate(magmom_start(3,Natom,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(magmom_start))*kind(magmom_start),'magmom_start','calc_bls')
            magmom_start=0.0_dblprec
         end if

         if (do_bls=='Q') then
            call setup_bls_qcoord(N1,N2,N3,C1,C2,C3)
            allocate(magmom_qw(bls_nstep,3,numq,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(magmom_qw))*kind(magmom_qw),'magmom_qw','calc_bls')
            magmom_qw=0.0_dblprec
         endif

         if ((do_bls=='Y'.or.do_bls=='Q').and.do_bls_local_axis=='Y') then
            allocate(mavg_local_axis(3,Natom,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(mavg_local_axis))*kind(mavg_local_axis),'mavg_local_axis','calc_bls')
            mavg_local_axis=emomM
         end if

         !
         ! Set up frequencies
         nw_bls = bls_nstep
         allocate(w_bls(nw_bls),stat=i_stat)
         call memocc(i_stat,product(shape(w_bls))*kind(w_bls),'w_bls','calc_bls')
         call set_w_bls(delta_t, bls_step, bls_nstep)
         !
         allocate(bls_tidx(bls_nsamples),stat=i_stat)
         call memocc(i_stat,product(shape(bls_tidx))*kind(bls_tidx),'bls_tidx','calc_bls')
         bls_tidx=0
         !
         allocate(expw_bls(nw_bls),stat=i_stat)
         call memocc(i_stat,product(shape(expw_bls))*kind(expw_bls),'expw_bls','calc_bls')
         expw_bls=0.0_dblprec
         !
         allocate(bls_dos(3,nw_bls),stat=i_stat)
         call memocc(i_stat,product(shape(bls_dos))*kind(bls_dos),'bls_dos','calc_bls')
         bls_dos=0.0_dblprec
         !
         allocate(mavg_axis(3,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(mavg_axis))*kind(mavg_axis),'mavg_axis','calc_bls')
         mavg_axis=0.0_dblprec
         !
         allocate(mort_axis(3,3,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(mort_axis))*kind(mort_axis),'mort_axis','calc_bls')
         mort_axis=0.0_dblprec
         !
         dt_bls=delta_t*bls_step
         !

         ! Calculate g(k) for the current iteration and add to G(k,t)
      else if (flag==1) then

         do isamp=1,bls_nsamples
            ! Exit early if no need for sampling
            shifted_mstep=mstep-(isamp-1)*bls_sample_step
            if (.not.(mod(shifted_mstep,bls_step)/=1.or.bls_tidx(isamp)>=bls_nstep)) then
               ! Increase sample counter
               bls_tidx(isamp)=bls_tidx(isamp)+1

               if(bls_tidx(isamp)==1) then
#if _OPENMP >= 201307
                  !$omp parallel do default(shared) private(iens,iatom) collapse(2) reduction(+:mavg_axis)
#endif
                  do iens=1,Mensemble
                     do iatom=1,Natom
                        magmom_start(1,iatom,iens)=emomM(1,iatom,iens)
                        magmom_start(2,iatom,iens)=emomM(2,iatom,iens)
                        magmom_start(3,iatom,iens)=emomM(3,iatom,iens)
                     end do
                  end do
#if _OPENMP >= 201307
                  !$end parallel do
#endif
               end if

               ! Find magnetization axis to obtain orthogonal components to M
               mavg_axis=0.0_dblprec
#if _OPENMP >= 201307
               !$omp parallel do default(shared) private(iens,iatom) collapse(2) reduction(+:mavg_axis)
#endif
               do iens=1,Mensemble
                  do iatom=1,Natom
                     mavg_axis(1,iens)=mavg_axis(1,iens)+emomM(1,iatom,iens)
                     mavg_axis(2,iens)=mavg_axis(2,iens)+emomM(2,iatom,iens)
                     mavg_axis(3,iens)=mavg_axis(3,iens)+emomM(3,iatom,iens)
                  end do
               end do
#if _OPENMP >= 201307
               !$end parallel do
#endif
               do iens=1,Mensemble
                  mavg_axis(:,iens)=mavg_axis(:,iens)/Natom
                  mavg_norm=sum(mavg_axis(:,iens)*mavg_axis(:,iens))**0.5_dblprec
                  if(mavg_norm>1.0d-2) then
                     mavg_axis(:,iens)=mavg_axis(:,iens)/mavg_norm
                  else
                     mavg_axis(:,iens)=(/0.0_dblprec,0.0_dblprec,1.0_dblprec/)
                  end if
               end do
               ! Without orthogonalization
               mort_axis=0.0_dblprec
               mort_axis(1,1,:)=1.0_dblprec
               mort_axis(2,2,:)=1.0_dblprec
               mort_axis(3,3,:)=1.0_dblprec
               !  Measure magnetic moments and transform on the fly to M(w)
               !  Start with precalculation of exponential/cosine lookup-table
               eexp=(bls_tidx(isamp))*dt_bls*i

               do iw=1,nw_bls
                  select case(bls_window_fun)
                  ! Hann windowing
               case(2)
                  expw_bls(iw)=exp(eexp*w_bls(iw))*(0.50_dblprec-0.50*cos(2.0_dblprec*pi*(bls_tidx(isamp)-1)/(bls_nstep-1.0)))
                  ! Hamming windowing
               case(3)
                  expw_bls(iw)=exp(eexp*w_bls(iw))*(0.54_dblprec-0.46*cos(2.0_dblprec*pi*(bls_tidx(isamp)-1.0_dblprec)/(bls_nstep-1.0_dblprec)))
                  ! Blackman-Harris
               case(4)
                  expw_bls(iw)=exp(eexp*w_bls(iw))* &
                     (0.35785_dblprec-0.48829_dblprec*cos(2.0_dblprec*pi*(bls_tidx(isamp)-1)/(1.0_dblprec*bls_nstep))+ &
                     0.14128_dblprec*cos(4.0_dblprec*pi*(bls_tidx(isamp)-1)/(1.0_dblprec*bls_nstep))   &
                     -0.01168_dblprec*cos(6.0_dblprec*pi*(bls_tidx(isamp)-1)/(1.0_dblprec*bls_nstep)))
                  ! Square windows
               case default
                  expw_bls(iw)=exp(eexp*w_bls(iw))
               end select
            end do

            if((do_bls=='Y'.or.do_bls=='Q').and.do_bls_local_axis=='N') then
               !$omp parallel do default(shared) private(iens,iatom,component,iw,m_proj) collapse(2)
               do iens=1,Mensemble
                  do iatom=1,Natom
                     !m_proj(1)=mort_axis(1,1,iens)*emomM(1,iatom,iens)+mort_axis(2,1,iens)*emomM(2,iatom,iens)+mort_axis(3,1,iens)*emomM(3,iatom,iens)
                     !m_proj(2)=mort_axis(1,2,iens)*emomM(1,iatom,iens)+mort_axis(2,2,iens)*emomM(2,iatom,iens)+mort_axis(3,2,iens)*emomM(3,iatom,iens)
                     !m_proj(3)=mort_axis(1,3,iens)*emomM(1,iatom,iens)+mort_axis(2,3,iens)*emomM(2,iatom,iens)+mort_axis(3,3,iens)*emomM(3,iatom,iens)
                     m_proj(1)=emomM(1,iatom,iens)
                     m_proj(2)=emomM(2,iatom,iens)
                     m_proj(3)=emomM(3,iatom,iens)
                     do component=1,3
                        do iw=1,nw_bls
                           magmom_w(iw,component,iatom,iens)=magmom_w(iw,component,iatom,iens)+m_proj(component)*expw_bls(iw)
                        end do
                     end do
                  end do
               end do
               !$omp end parallel do
            else if ((do_bls=='Y'.or.do_bls=='Q').and.do_bls_local_axis=='B') then
               !$omp parallel do default(shared) private(iens,iatom,component,iw,m_proj) collapse(2)
               do iens=1,Mensemble
                  do iatom=1,Natom
                     beff_norm=sqrt(sum(beff(:,iatom,iens)**2))+1.0d-12
                     m_proj(1)=beff(2,iatom,iens)*emomM(3,iatom,iens)-beff(3,iatom,iens)*emomM(2,iatom,iens)
                     m_proj(2)=beff(3,iatom,iens)*emomM(1,iatom,iens)-beff(1,iatom,iens)*emomM(3,iatom,iens)
                     m_proj(3)=beff(1,iatom,iens)*emomM(2,iatom,iens)-beff(2,iatom,iens)*emomM(1,iatom,iens)
                     m_proj=m_proj/beff_norm
                     do component=1,3
                        do iw=1,nw_bls
                           magmom_w(iw,component,iatom,iens)=magmom_w(iw,component,iatom,iens)+m_proj(component)*expw_bls(iw)
                        end do
                     end do
                  end do
               end do
               !$omp end parallel do
            else if ((do_bls=='Y'.or.do_bls=='Q').and.do_bls_local_axis=='Y') then
               ! Sample current moment directions and weight in for average
               do iens=1,Mensemble
                  do iatom=1,Natom
                     mm=sum(emomM(:,iatom,iens)**2)**0.5_dblprec
                     mavg_local_axis(:,iatom,iens)=mavg_local_axis(:,iatom,iens)*(1.0_dblprec-bls_local_axis_mix)+emomM(:,iatom,iens)/mm*bls_local_axis_mix
                     !mavg_local_axis(:,iatom,iens)=mavg_local_axis(:,iatom,iens)*0.95_dblprec+emomM(:,iatom,iens)/mm*0.05_dblprec
                     mavg_local_axis(:,iatom,iens)=mavg_local_axis(:,iatom,iens)/(sum(mavg_local_axis(:,iatom,iens)**2)**0.5+1.0e-12_dblprec)
                  end do
               end do
!              call find_local_rotmat(Natom*Mensemble,mavg_local_axis,mavg_local_rotmat)
               !$omp parallel do default(shared) private(iens,iatom,component,iw,m_proj) collapse(2)
               do iens=1,Mensemble
                  do iatom=1,Natom
                     call gramms(mavg_local_axis(1:3,iatom,iens),mort_axis(1:3,1:3,1),1)
                     m_proj(1)=mort_axis(1,1,1)*emomM(1,iatom,iens)+mort_axis(2,1,1)*emomM(2,iatom,iens)+mort_axis(3,1,1)*emomM(3,iatom,iens)
                     m_proj(2)=mort_axis(1,2,1)*emomM(1,iatom,iens)+mort_axis(2,2,1)*emomM(2,iatom,iens)+mort_axis(3,2,1)*emomM(3,iatom,iens)
                     m_proj(3)=mort_axis(1,3,1)*emomM(1,iatom,iens)+mort_axis(2,3,1)*emomM(2,iatom,iens)+mort_axis(3,3,1)*emomM(3,iatom,iens)
                     do component=1,3
                        do iw=1,nw_bls
                           magmom_w(iw,component,iatom,iens)=magmom_w(iw,component,iatom,iens)+m_proj(component)*expw_bls(iw)
                        end do
                     end do
                  end do
               end do
               !$omp end parallel do
            end if
         end if
      end do
      return
      ! Final operations, sum and print
   else if (flag==2) then
      magmom_w(1,:,:,:)=(0.0_dblprec,0.0_dblprec)
      ! Fourier transform to reciprocal space
      if (do_bls=='Q') then
         !$omp parallel do default(shared) private(iw,iq,iens,r,qdr,epowqr)
         do iw=1,(nw_bls-1)/2
            do iq=1,numq
               do iens=1, Mensemble
                  do r=1,Natom
                     qdr=qcoord(1,iq)*coord(1,r)+qcoord(2,iq)*coord(2,r)+qcoord(3,iq)*coord(3,r)
                     epowqr=exp(i*qfac*qdr)*nainv
                     magmom_qw(iw,1,iq,iens) = magmom_qw(iw,1,iq,iens) + magmom_w(iw,1,r,iens)*epowqr
                     magmom_qw(iw,2,iq,iens) = magmom_qw(iw,2,iq,iens) + magmom_w(iw,2,r,iens)*epowqr
                     magmom_qw(iw,3,iq,iens) = magmom_qw(iw,3,iq,iens) + magmom_w(iw,3,r,iens)*epowqr
                  enddo
               enddo
            enddo
         enddo
         !$omp end parallel do
      endif

      ! Write the site dependent power spectra
      if(do_bls=='Y'.or.do_bls=='Q') then
         write (filn,'(''bls.'',a8,''.out'')') simid
         open(ofileno, file=filn,position='append')
         do iw=1,(nw_bls-1)/2
            do iatom=1,Natom
               bls_avg=0.0_dblprec
               do iens=1,Mensemble
                  bls_avg=bls_avg+abs(magmom_w(iw,1:3,iatom,iens))**2/bls_nsamples
               end do
               bls_avg_norm=sqrt(sum(bls_avg))
               if (bls_window.gt.1) then
                  write (ofileno,10005) mstep,iw,iatom,w_bls(iw), bls_avg,bls_avg_norm
               else
                  write (ofileno,10004) iw,iatom, w_bls(iw), bls_avg,bls_avg_norm
               endif
            end do
         end do
         close(ofileno)
         ! also writing complex bls in order to get phases (probably does not work with ensemble averaging)
         write (filn,'(''cbls.'',a8,''.out'')') simid
         open(ofileno, file=filn,position='append')
         allocate(bls_norm(3,Natom,Mensemble))
         allocate(bls_phase(3,Natom,Mensemble))
         do iw=1,(nw_bls-1)/2
            do iens=1,Mensemble
               do iatom=1,Natom
                  do component=1,3
                     bls_norm(component,iatom,iens)=abs(magmom_w(iw,component,iatom,iens))**2
                     bls_phase(component,iatom,iens)=atan2(aimag(magmom_w(iw,component,iatom,iens)),real(magmom_w(iw,component,iatom,iens)))
                  end do
               end do
               bls_phase_zero=bls_phase(:,2,iens)
               do iatom=1,Natom
                  bls_phase(:,iatom,iens)=bls_phase(:,iatom,iens)-bls_phase_zero(:)
               end do
            end do
            do iatom=1,Natom
               bls_norm_avg=0.0_dblprec
               bls_phase_avg=0.0_dblprec
               do iens=1,Mensemble
                  bls_phase_avg=bls_phase_avg+bls_phase(:,iatom,iens)
                  bls_norm_avg=bls_norm_avg+bls_norm(:,iatom,iens)
               end do
               bls_phase_avg=bls_phase_avg/Mensemble
               bls_norm_avg=bls_norm_avg/Mensemble
               bls_norm_avg_avg=sqrt(sum(bls_norm_avg(1:2)**2))
               bls_phase_avg_avg=sum(bls_phase_avg(1:2))/2
               if (bls_window.gt.1) then
                  write (ofileno,10015) mstep,iw,iatom, w_bls(iw), bls_norm_avg,bls_phase_avg,bls_norm_avg_avg,bls_phase_avg_avg
               else
                  write (ofileno,10014) iw,iatom,  w_bls(iw), bls_norm_avg,bls_phase_avg ,bls_norm_avg_avg,bls_phase_avg_avg
               endif
            end do
         end do
         close(ofileno)
         deallocate(bls_norm)
         deallocate(bls_phase)
      end if

      ! Write the site dependent power spectra in reciprocal space
      if (do_bls=='Q') then
         write (filn,'(''blsqres.'',a8,''.out'')') simid
         open(ofileno, file=filn,position='append')
         do iw=1,(nw_bls-1)/2
            do iq=1,numq
               do iens=1, Mensemble
                  if (bls_window.gt.1) then
                     write(ofileno,10005) mstep,iw,iq,iens,w_bls(iw),               &
                        (abs(magmom_qw(iw,:,iq,iens)))**2/bls_nsamples,             &
                        abs(sqrt(sum(magmom_qw(iw,:,iatom,iens)**2)))**2/bls_nsamples
                  else
                     write(ofileno,10004) iw,iq,iens,w_bls(iw),                     &
                        (abs(magmom_qw(iw,:,iq,iens)))**2/bls_nsamples,             &
                        abs(sqrt(sum(magmom_qw(iw,:,iatom,iens)**2)))**2/bls_nsamples
                  endif
               enddo
            enddo
         enddo
         close(ofileno)
      endif

      ! Calculate and write the bls DOS
      bls_dos=0.0_dblprec
      write (filn,'(''blsdos.'',a8,''.out'')') simid
      open(ofileno, file=filn,position='append')
      if ((mstep/(bls_step*bls_nstep))==1) then
         if (bls_window.gt.1) then
            write(ofileno,'(a7,5a16)') "#time","E(meV)","DOS_x","DOS_y","DOS_z","DOS_tot"
         else
            write(ofileno,'(5a16)') "#    E(meV)","DOS_x","DOS_y","DOS_z","DOS_tot"
         endif
      endif
      !do iw=1,nw_bls
      do iw=1,(nw_bls)/2
         do iens=1, Mensemble
            do iatom=1, Natom
               bls_dos(1,iw)=bls_dos(1,iw)+abs(magmom_w(iw,1,iatom,iens))
               bls_dos(2,iw)=bls_dos(2,iw)+abs(magmom_w(iw,2,iatom,iens))
               bls_dos(3,iw)=bls_dos(3,iw)+abs(magmom_w(iw,3,iatom,iens))
            enddo
         enddo
         bls_dos(1,iw)= abs(bls_dos(1,iw))**2/bls_nstep/bls_nsamples
         bls_dos(2,iw)= abs(bls_dos(2,iw))**2/bls_nstep/bls_nsamples
         bls_dos(3,iw)= abs(bls_dos(3,iw))**2/bls_nstep/bls_nsamples
         ! Taking care of the time dependent bls
         if(do_bls_local_axis=='N') then
            if (bls_window.gt.1) then
               write (ofileno,10006) mstep,hbar_mev*w_bls(iw),bls_dos(1,iw),        &
               bls_dos(2,iw), bls_dos(3,iw), norm2(bls_dos(:,iw))
            else
               write (ofileno,10007) hbar_mev*w_bls(iw),bls_dos(1,iw),bls_dos(2,iw),&
               bls_dos(3,iw), norm2(bls_dos(:,iw))
            endif
         else
            if (bls_window.gt.1) then
               write (ofileno,10006) mstep,hbar_mev*w_bls(iw),bls_dos(1,iw),        &
               bls_dos(2,iw),bls_dos(3,iw),norm2(bls_dos(:,iw))
            else
               write (ofileno,10007) hbar_mev*w_bls(iw),bls_dos(1,iw),bls_dos(2,iw),&
               bls_dos(3,iw), norm2(bls_dos(:,iw))
            endif
         end if
      end do
      close(ofileno)

      bls_dos=0.0_dblprec
      bls_tidx=0
      magmom_w=0.0_dblprec
      expw_bls=0.0_dblprec
      mavg_axis=0.0_dblprec
      if (do_bls=='Q') magmom_qw=0.0_dblprec
      magmom_start=0.0_dblprec

   end if
   return
   !
   10004 format (2i7,g18.8,3g16.8,2x,g16.8)
   10005 format (3i7,g18.8,3g16.8,2x,g16.8)
   10006 format (i7,5G16.8)
   10007 format (5G16.8)
   10014 format (2i7,g18.8,3g16.8,3g16.8,2x,2g16.8)
   10015 format (3i7,g18.8,3g16.8,3g16.8,2x,2g16.8)
   !
end subroutine calc_bls

!------------------------------------------------------------------------------------
! SUBROUTINE: set_w_bls
!> Calculate suitable values of frequencies for S(q,t) -> S(q,w) transform
!> @todo Change setup to smarter algorithm wrt the exchange strength of the
!> system
!------------------------------------------------------------------------------------
subroutine set_w_bls(delta_t, bls_step, bls_nstep)
   !
   use Constants, only : pi, hbar_mev
   !
   implicit none
   !
   real(dblprec), intent(in) :: delta_t !< Time step
   integer, intent(in) :: bls_step      !< Separation between sampling steps
   integer, intent(in) :: bls_nstep     !< Number of steps to sample
   !
   integer :: j
   real(dblprec) :: dt !< Time step
   real(dblprec) :: dww
   real(dblprec) :: emin,emax
   !
   dt=bls_step*delta_t
   nw_bls = bls_nstep
   dww = 2*pi/(nw_bls*dt)

   do j=0,nw_bls-1
      w_bls(j+1)=j*dww
   end do

   emin=abs(hbar_mev*(w_bls(2)-w_bls(1)))
   emax=abs(hbar_mev*w_bls(nw_bls))

end subroutine set_w_bls


!------------------------------------------------------------------------------------
! SUBROUTINE: gramms
!> @brief Gramm-Schmidt orthogonalization
!------------------------------------------------------------------------------------
subroutine gramms(v_in,m_out,ldim)
   !
   integer, intent(in) :: ldim
   real(dblprec), dimension(3,ldim), intent(in) :: v_in
   real(dblprec), dimension(3,3,ldim), intent(out) :: m_out
   !
   real(dblprec), dimension(3) :: v_trial
   real(dblprec)  :: v_in_norm
   integer :: i
   !
   do i=1,ldim
      v_in_norm=sqrt(v_in(1,i)*v_in(1,i)+v_in(2,i)*v_in(2,i)+v_in(3,i)*v_in(3,i))
      ! Set trial vector to x
      v_trial=(/1.0_dblprec, 0.0_dblprec, 0.0_dblprec/)
      v_trial=v_trial/sqrt(sum(v_trial*v_trial))
      ! If input vector is close to x, set trial vector to y..
      if(abs(v_in(1,i))/v_in_norm>0.9995_dblprec) then
         v_trial=(/0.0_dblprec, 1.0_dblprec, 0.0_dblprec/)
         v_trial=v_trial/sqrt(sum(v_trial*v_trial))
      end if
      ! Keep input vector as third output vector
      m_out(1:3,3,i)=v_in(1:3,i)/v_in_norm
      ! Get first orthogonal vector through cross product
      m_out(1,1,i)=m_out(2,3,i)*v_trial(3)-m_out(3,3,i)*v_trial(2)
      m_out(2,1,i)=m_out(3,3,i)*v_trial(1)-m_out(1,3,i)*v_trial(3)
      m_out(3,1,i)=m_out(1,3,i)*v_trial(2)-m_out(2,3,i)*v_trial(1)
      ! Normalize
      m_out(:,1,i)=m_out(:,1,i) / sqrt(sum(m_out(:,1,i)*m_out(:,1,i)))
      ! Get second orthogonal vector by cross product
      m_out(1,2,i)=m_out(2,3,i)*m_out(3,1,i)-m_out(3,3,i)*m_out(2,1,i)
      m_out(2,2,i)=m_out(3,3,i)*m_out(1,1,i)-m_out(1,3,i)*m_out(3,1,i)
      m_out(3,2,i)=m_out(1,3,i)*m_out(2,1,i)-m_out(2,3,i)*m_out(1,1,i)
      ! Normalize for safety..
      m_out(:,2,i)=m_out(:,2,i) / sqrt(sum(m_out(:,2,i)*m_out(:,2,i)))
      ! Get second orthogonal vector by cross product
      m_out(1,1,i)=m_out(2,2,i)*m_out(3,3,i)-m_out(3,2,i)*m_out(2,3,i)
      m_out(2,1,i)=m_out(3,2,i)*m_out(1,3,i)-m_out(1,2,i)*m_out(3,3,i)
      m_out(3,1,i)=m_out(1,2,i)*m_out(2,3,i)-m_out(2,2,i)*m_out(1,3,i)
      ! Normalize for safety..
      m_out(:,1,i)=m_out(:,1,i) / sqrt(sum(m_out(:,1,i)*m_out(:,1,i)))
   end do
   !
end subroutine gramms


!------------------------------------------------------------------------------------
! SUBROUTINE: gramms2
!> @brief Gramm-Schmidt orthogonalization alternative
!------------------------------------------------------------------------------------
subroutine gramms2(v_in,v_qt,m_out,ldim)
   !
   integer, intent(in) :: ldim
   real(dblprec), dimension(3,ldim), intent(in) :: v_in
   real(dblprec), dimension(3,ldim), intent(in) :: v_qt
   real(dblprec), dimension(3,ldim), intent(out) :: m_out
   !
   real(dblprec), dimension(3) :: v_perp1, v_perp2
   integer :: i
   !
   do i=1,ldim
      ! Get first orthogonal vector by cross product
      v_perp1(1)=v_qt(2,i)*v_in(3,i)-v_qt(3,i)*v_in(2,i)
      v_perp1(2)=v_qt(3,i)*v_in(1,i)-v_qt(2,i)*v_in(3,i)
      v_perp1(3)=v_qt(1,i)*v_in(2,i)-v_qt(1,i)*v_in(1,i)
      ! Get second orthogonal vector by cross product
      v_perp2(1)=v_perp1(2)*v_in(3,i)-v_perp1(3)*v_in(2,i)
      v_perp2(2)=v_perp1(3)*v_in(1,i)-v_perp1(2)*v_in(3,i)
      v_perp2(3)=v_perp1(1)*v_in(2,i)-v_perp1(1)*v_in(1,i)
      ! Calculate projections
      m_out(1,i)=sum(v_perp1*v_in(:,i))
      m_out(2,i)=sum(v_perp2*v_in(:,i))
      m_out(3,i)=sum(v_qt(:,i)*v_in(:,i))
   end do
   !
end subroutine gramms2

!------------------------------------------------------------------------------------
! SUBROUTINE: setup_bls_qcoord
!> @brief Setup the q-coordinates for the BLS analysis
!------------------------------------------------------------------------------------
subroutine setup_bls_qcoord(N1,N2,N3,C1,C2,C3)

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
   real(dblprec), dimension(:), allocatable :: dqcoord
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
   if(allocated(qcoord)) then
      i_all=-product(shape(qcoord))*kind(qcoord)
      deallocate(qcoord,stat=i_stat)
      call memocc(i_stat,i_all,'qcoord','setup_bls_qcoord')
   end if

   numq=0
   numq=N1*N2*N3
   allocate(ia(numq),stat=i_stat)
   call memocc(i_stat,product(shape(ia))*kind(ia),'ia','setup_bls_qcoord')
   allocate(qcoord(3,numq),stat=i_stat)
   call memocc(i_stat,product(shape(qcoord))*kind(qcoord),'qcoord','setup_bls_qcoord')
   allocate(dqcoord(numq),stat=i_stat)
   call memocc(i_stat,product(shape(dqcoord))*kind(dqcoord),'dqcoord','setup_bls_qcoord')
   iq=0
   do zq=-(N3-1)/2,N3/2
      do yq=-(N2-1)/2,N2/2
         do xq=-(N1-1)/2,N1/2
            iq=iq+1
            qcoord(:,iq)=xq/(1.0_dblprec*N1)*b1+yq/(1.0_dblprec*N2)*b2+zq/(1.0_dblprec*N3)*b3
            dqcoord(iq)=qcoord(1,iq)**2+qcoord(2,iq)**2+qcoord(3,iq)**2
         end do
      end do
   end do

   if(allocated(dqcoord)) then
      i_all=-product(shape(dqcoord))*kind(dqcoord)
      deallocate(dqcoord,stat=i_stat)
      call memocc(i_stat,i_all,'dqcoord','setup_bls_qcoord')
   end if

   open(ofileno,file='blsqpoints.out',status='replace')
   do iq=1,numq
      write(ofileno,'(i6,3f14.6)') iq,qcoord(1,iq),qcoord(2,iq),qcoord(3,iq)
   end do
   close(ofileno)

end subroutine setup_bls_qcoord

!------------------------------------------------------------------------------------
! SUBROUTINE: deallocate_bls_data
!> @brief Deallocation of the arrays for the bls analysis
!------------------------------------------------------------------------------------
subroutine deallocate_bls_data()

   implicit none

   integer :: i_all, i_stat
   if (do_bls=='Y'.or.do_bls=='D'.or.do_bls=='Q') then
      i_all=-product(shape(magmom_start))*kind(magmom_start)
      deallocate(magmom_start,stat=i_stat)
      call memocc(i_stat,i_all,'magmom_start','calc_bls')
      !
      i_all=-product(shape(magmom_w))*kind(magmom_w)
      deallocate(magmom_w,stat=i_stat)
      call memocc(i_stat,i_all,'magmom_w','calc_bls')
      !
      i_all=-product(shape(w_bls))*kind(w_bls)
      deallocate(w_bls,stat=i_stat)
      call memocc(i_stat,i_all,'w_bls','calc_bls')
      !
      i_all=-product(shape(expw_bls))*kind(expw_bls)
      deallocate(expw_bls,stat=i_stat)
      call memocc(i_stat,i_all,'expw_bls','calc_bls')
      !
      i_all=-product(shape(bls_tidx))*kind(bls_tidx)
      deallocate(bls_tidx,stat=i_stat)
      call memocc(i_stat,i_all,'bls_tidx','calc_bls')
      !
      i_all=-product(shape(bls_dos))*kind(bls_dos)
      deallocate(bls_dos,stat=i_stat)
      call memocc(i_stat,i_all,'bls_dos','calc_bls')
      !

      if (do_bls=='Q') then
         i_all=-product(shape(magmom_qw))*kind(magmom_qw)
         deallocate(magmom_qw,stat=i_stat)
         call memocc(i_stat,i_all,'magmom_qw','calc_bls')

         i_all=-product(shape(qcoord))*kind(qcoord)
         deallocate(qcoord,stat=i_stat)
         call memocc(i_stat,i_all,'qcoord','calc_bls')

      endif

      if (do_bls=='Y'.and.do_bls_local_axis=='Y') then
         i_all=-product(shape(mavg_local_axis))*kind(mavg_local_axis)
         deallocate(mavg_local_axis,stat=i_stat)
         call memocc(i_stat,i_all,'mavg_local_axis','calc_bls')

      end if

   end if
end subroutine deallocate_bls_data


!------------------------------------------------------------------------------------
!> @brief
!> Read input parameters.
!
!> @author
!> Anders Bergman
!------------------------------------------------------------------------------------
subroutine read_parameters_bls(ifile)
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
         case('do_bls')
            read(ifile,*,iostat=i_err) do_bls
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_bls_local_axis')
            read(ifile,*,iostat=i_err) do_bls_local_axis
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('bls_step')
            read(ifile,*,iostat=i_err) bls_step
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('bls_nstep')
            read(ifile,*,iostat=i_err) bls_nstep
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('bls_nsamples')
            read(ifile,*,iostat=i_err) bls_nsamples
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('bls_sample_step')
            read(ifile,*,iostat=i_err) bls_sample_step
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('bls_window')
            read(ifile,*,iostat=i_err) bls_window
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('bls_local_axis_mix')
            read(ifile,*,iostat=i_err) bls_local_axis_mix
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('bls_window_fun')
            read(ifile,*,iostat=i_err) bls_window_fun
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
end subroutine read_parameters_bls

end module BLS
