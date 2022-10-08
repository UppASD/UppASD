module BiMagnonCorr
   use Parameters
   use Profiling
   use Correlation_utils
   use Qvectors
   use Omegas
   !
   implicit none
   !

   !Spin-spin correlation parameters from input file
   !character(len=1) :: do_bimag  !< Perform spin-correlation sampling of bi-magnons (Y/N)
   real(dblprec), dimension(:), allocatable :: bimagnon_dos !< Bi-magnon density of states

   complex(dblprec), dimension(:,:), allocatable :: corr_bkt        ! Correlation for B(k,t)
   complex(dblprec), dimension(:,:), allocatable :: corr_bkw        ! Correlation for B(k,w)

   real(dblprec), dimension(:,:,:), allocatable ::  bm_proj      !< Array containing bi-moments projected to parallel/perpendicular directions
   real(dblprec) :: bm_qfac
   integer :: sc_samp_done_bm               !< Flag to keep track of if S(q) sampling is done (for do_sc='B')
   real(dblprec), dimension(:), allocatable :: bmstep_arr  ! Array storing sc_step for each sample

   character(len=1) :: do_bm_local_axis   !< Perform BQW along local quantization axis (Y/N)
   real(dblprec), dimension(:,:), allocatable :: bm_mavg_axis      !< Internal array for finding the global direction of magnetization

   real(dblprec), dimension(:,:,:), allocatable :: bm_mavg_local_axis!< Internal array for finding the local average direction of magnetization
   real(dblprec), dimension(:,:,:,:), allocatable :: bm_mavg_local_rotmat !< Internal array for finding rotatons to project the local average direction of magnetization
   real(dblprec), dimension(:,:,:), allocatable :: bm_mort_axis    !< Orthogonal (and parallel ) components of magnetization
   real(dblprec), dimension(:,:,:), allocatable ::  bm_m_proj      !< Array containing moments projected to parallel/perpendicular directions

   public


contains

   !---------------------------------------------------------------------------------
   ! SUBROUTINE calc_bgkt
   !> @brief Calculate bimagnon correlation function
   !---------------------------------------------------------------------------------
   subroutine calc_bgkt(Natom,Mensemble,nHam,coord,simid, emomM, sc_tidx, sc_nstep, &
         flag,q_flag,deltat_corr,do_connected,do_conv,sigma_q,sigma_w,LQfactor,        &
         LWfactor,max_no_neigh,conf_num,nlist,nlistsize,ncoup,aham)
      !
      use Constants
      use Math_functions, only : gramms
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: nHam !< Number of atoms in Hamiltonian
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
      integer, dimension(nHam),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(max_no_neigh,nHam,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
      integer, dimension(Natom),intent(in) :: aham !< Hamiltonian look-up table

      integer :: iq,iw,step,r,l,i_stat,i_all,j, ih
      character(len=30) :: filn
      complex(dblprec) :: epowqr, i, iqfac,tt, epowwt
      real(dblprec) :: qdr,nainv,bm_mavg_norm, wfac
      real(dblprec), dimension(3) :: mavrg_vec
      real(dblprec), dimension(sc_nstep+1) :: dt
      real(dblprec), dimension(3,3) :: local_rotmat
      !
      i=(0.0_dblprec,1.0_dblprec)

      if(flag==0) then
         ! First call, allocate and clear arrays
         !
         allocate(corr_bkt(nq,sc_nstep+1),stat=i_stat)
         call memocc(i_stat,product(shape(corr_bkt))*kind(corr_bkt),'corr_bkt','calc_bgkt')
         corr_bkt=0.0_dblprec
         !
         if (do_bm_local_axis=='Y') then
            bm_mavg_local_axis=0.0_dblprec
            bm_mavg_local_axis(1,:,:)=emomM(1,1,1)
            bm_mavg_local_axis(2,:,:)=emomM(2,1,1)
            bm_mavg_local_axis(3,:,:)=emomM(3,1,1)
         end if
         bm_mavg_axis=0.0_dblprec
         bm_mort_axis=0.0_dblprec
         !
         allocate(bm_proj(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(bm_proj))*kind(bm_proj),'bm_proj','calc_bgkt')
         !
         flag=1
         sc_samp_done_bm=0
         bm_qfac=2.0_dblprec*pi
      end if

      nainv=1.0_dblprec/Natom

      ! Calculate b(k) for the current iteration and add to b(k,t)
      if (flag==1) then
         iqfac=i*bm_qfac
         call calc_mavrg_vec(Natom,Mensemble,emomM,mavrg_vec,bm_mavg_axis)
         do l=1,Mensemble
            bm_mavg_axis(:,l)=bm_mavg_axis(:,l)/Natom
            bm_mavg_norm=sum(bm_mavg_axis(:,l)*bm_mavg_axis(:,l))**0.5_dblprec
            if(bm_mavg_norm>1.0d-2) then
               bm_mavg_axis(:,l)=bm_mavg_axis(:,l)/bm_mavg_norm
            else
               bm_mavg_axis(:,l)=(/0.0_dblprec,0.0_dblprec,1.0_dblprec/)
            end if
         end do

         if(do_connected/='Y') then
            ! Keep the average magnetization for subtraction below if do_connected
            ! otherwise put it to zero (in order to remove a number of if-statements in-loop)
            mavrg_vec=0.0_dblprec
         end if

         if(do_bm_local_axis/='Y') then
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
                  local_rotmat=bm_mavg_local_rotmat(:,:,r,l)
                  bm_proj(1,r,l)=local_rotmat(1,1)*emomM(1,r,l)+local_rotmat(2,1)*emomM(2,r,l)+local_rotmat(3,1)*emomM(3,r,l)
                  bm_proj(2,r,l)=local_rotmat(1,2)*emomM(1,r,l)+local_rotmat(2,2)*emomM(2,r,l)+local_rotmat(3,2)*emomM(3,r,l)
                  bm_proj(3,r,l)=local_rotmat(1,3)*emomM(1,r,l)+local_rotmat(2,3)*emomM(2,r,l)+local_rotmat(3,3)*emomM(3,r,l)
               end do
            end do
            !$omp end parallel do
         end if

         !$omp parallel do default(shared) private(r,iq,ih,l,j,qdr,epowqr) schedule(static)
         do iq=1,nq
            do l=1,Mensemble
               do r=1,Natom
                  ih=aham(r)
                  qdr=q(1,iq)*coord(1,r)+q(2,iq)*coord(2,r)+q(3,iq)*coord(3,r)
                  epowqr=exp(iqfac*qdr)*nainv
                  do j=1,nlistsize(ih)
                     corr_bkt(iq,sc_tidx)=corr_bkt(iq,sc_tidx)-epowqr*ncoup(j,ih,1)*(bm_proj(1,r,l)*bm_proj(1,nlist(j,r),l)+ &
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
         corr_bkw=0.0_dblprec

         ! Finish sampling and transform (k,t)->(k,w)
         if(sc_tidx.GT.sc_nstep) then
            sc_tidx = sc_nstep
         end if

         j=1
         do while(j.LE.sc_tidx)
            dt(j) = bmstep_arr(j)*deltat_corr(j)
            j = j+1
         end do

         wfac=1.0_dblprec
         corr_bkw=0.0_dblprec

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
         write (filn,'(''bsqt.'',a,''.out'')') trim(simid)
         open(ofileno, file=filn)
         do iq=1,Nq
            do step=1,sc_tidx
               write (ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),step, &
                  real(corr_bkt(iq,step)),aimag(corr_bkt(iq,step)),abs(corr_bkt(iq,step))
            end do
         end do
         close(ofileno)

         ! Write B(q,w)
         write (filn,'(''bsqw.'',a,''.out'')') trim(simid)
         open(ofileno, file=filn)
         write (filn,'(''bsqwx.'',a,''.out'')') trim(simid)
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

            ! Calculate the convolution for the magnon DOS
            call calc_sqw_conv(3,corr_bkw)

         ! Calculate and write the magnon DOS
         allocate(bimagnon_dos(nw/2),stat=i_stat)
         call memocc(i_stat,product(shape(bimagnon_dos))*kind(bimagnon_dos),'bimagnon_dos','calc_bgkt')
         !
         write (filn,'(''bswdos.'',a,''.out'')') trim(simid)
         open(ofileno, file=filn)
         bimagnon_dos=0.0_dblprec
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

end module BiMagnonCorr
