module Correlation_kernels
   use Parameters
   use Profiling
   use Correlation_type
   use Correlation_utils
   use Qvectors
   use Omegas
   !
   implicit none
   !

   public

contains

   subroutine corr_kernel_dual(Natom,Mensemble,nq,coord,r_mid,iqfac,nainv,SA,SB,SC)

      implicit none

      integer, intent(in) :: Natom
      integer, intent(in) :: Mensemble
      integer, intent(in) :: nq
      !
      real(dblprec), dimension(3,Natom), intent(in)  :: coord
      real(dblprec), dimension(3), intent(in)        :: r_mid
      !
      complex(dblprec), intent(in)    ::  iqfac
      real(dblprec), intent(in)       ::  nainv
      !
      real(dblprec), dimension(3,Natom,Mensemble), intent(in)  :: SA
      real(dblprec), dimension(3,Natom,Mensemble), intent(in)  :: SB
      real(dblprec), dimension(3,nq), intent(out)              :: SC

      integer ::  iq,l,r
      complex(dblprec) :: epowqr
      real(dblprec) :: qdr
      complex(dblprec), dimension(3)                        :: wA
      complex(dblprec), dimension(3)                        :: wB


      !$omp parallel do default(shared) private(wA,wB,r,iq,l,qdr,epowqr) schedule(static)
      do iq=1,nq
         do l=1,Mensemble
            wA=0.0_dblprec
            wB=0.0_dblprec
            do r=1,Natom
               qdr=q(1,iq)*(coord(1,r)-r_mid(1))+q(2,iq)*(coord(2,r)-r_mid(2))+q(3,iq)*(coord(3,r)-r_mid(3))
               epowqr=exp(iqfac*qdr)*nainv
               wA=wA+epowqr*SA(:,r,l)
               wB=wB+epowqr*SB(:,r,l)
            end do
            SC(:,iq)=SC(:,iq) + real(conjg(wA(:))*wB(:))
            !! wait with kt0
            !cc%corr_kt0(:,iq)=cc%corr_kt0(:,iq) + conjg(cc%corr_sA(:,iq))*cc%corr_sB(:,iq)
         end do
      end do
      !$omp end parallel do

   end subroutine corr_kernel_dual


   subroutine corr_kernel_equal(Natom,Mensemble,nq,coord,r_mid,iqfac,nainv,SA,SC)

      implicit none

      integer, intent(in) :: Natom
      integer, intent(in) :: Mensemble
      integer, intent(in) :: nq
      !
      real(dblprec), dimension(3,Natom), intent(in)  :: coord
      real(dblprec), dimension(3), intent(in)        :: r_mid
      !
      complex(dblprec), intent(in)    ::  iqfac
      real(dblprec), intent(in)       ::  nainv
      !
      real(dblprec), dimension(3,Natom,Mensemble), intent(in)  :: SA
      real(dblprec), dimension(3,nq), intent(out)              :: SC

      integer  :: iq,l,r
      complex(dblprec) :: epowqr
      real(dblprec) :: qdr
      complex(dblprec), dimension(3)                        :: wA


      !$omp parallel do default(shared) private(wA,r,iq,l,qdr,epowqr) schedule(static)
      do iq=1,nq
         do l=1,Mensemble
            wA=0.0_dblprec
            do r=1,Natom
               qdr=q(1,iq)*(coord(1,r)-r_mid(1))+q(2,iq)*(coord(2,r)-r_mid(2))+q(3,iq)*(coord(3,r)-r_mid(3))
               epowqr=exp(iqfac*qdr)*nainv
               wA=wA+epowqr*SA(:,r,l)
            end do
            SC(:,iq)=SC(:,iq) + real(conjg(wA(:))*wA(:))
         end do
      end do
      !$omp end parallel do

   end subroutine corr_kernel_equal

   subroutine corr_kernel_single(Natom,Mensemble,nq,coord,r_mid,iqfac,nainv,SA,SC)

      implicit none

      integer, intent(in) :: Natom
      integer, intent(in) :: Mensemble
      integer, intent(in) :: nq
      !
      real(dblprec), dimension(3,Natom), intent(in)  :: coord
      real(dblprec), dimension(3), intent(in)        :: r_mid
      !
      complex(dblprec), intent(in)    ::  iqfac
      real(dblprec), intent(in)       ::  nainv
      !
      real(dblprec), dimension(3,Natom,Mensemble), intent(in)  :: SA
      complex(dblprec), dimension(3,nq), intent(inout)              :: SC

      integer  :: iq,l,r
      complex(dblprec) :: epowqr
      real(dblprec) :: qdr
      complex(dblprec), dimension(3)                        :: wA


      !$omp parallel do default(shared) private(wA,r,iq,l,qdr,epowqr) schedule(static)
      do iq=1,nq
         do l=1,Mensemble
            wA=0.0_dblprec
            do r=1,Natom
               qdr=q(1,iq)*(coord(1,r)-r_mid(1))+q(2,iq)*(coord(2,r)-r_mid(2))+q(3,iq)*(coord(3,r)-r_mid(3))
               epowqr=exp(iqfac*qdr)*nainv!*sc_window_fac(sc_window_fun,iq,nq)
               wA=wA+epowqr*SA(:,r,l)
            end do
            SC(:,iq)=SC(:,iq) + wA(:)
         end do
      end do
      !$omp end parallel do

   end subroutine corr_kernel_single

   subroutine corr_kernel_proj(Natom,nproj,Mensemble,nq,coord,r_mid,aproj,iqfac,nainv,SA,SC)
      !
      implicit none

      integer, intent(in) :: Natom
      integer, intent(in) :: nproj
      integer, intent(in) :: Mensemble
      integer, intent(in) :: nq
      !
      real(dblprec), dimension(3,Natom), intent(in)  :: coord
      real(dblprec), dimension(3), intent(in)        :: r_mid
      integer, dimension(Natom), intent(in)              :: aproj
      !
      complex(dblprec), intent(in)    ::  iqfac
      real(dblprec), intent(in)       ::  nainv
      
      real(dblprec), dimension(3,Natom,Mensemble), intent(in)  :: SA
      complex(dblprec), dimension(3,nproj,nq), intent(out)              :: SC

      integer  :: iq,l,r,it
      complex(dblprec) :: epowqr
      real(dblprec) :: qdr
      complex(dblprec), dimension(3,nproj)                        :: wA


      !$omp parallel do default(shared) private(wA,r,iq,l,qdr,epowqr,it) schedule(static)
      do iq=1,nq
         do l=1,Mensemble
            wA=0.0_dblprec
            do r=1,Natom
               qdr=q(1,iq)*(coord(1,r)-r_mid(1))+q(2,iq)*(coord(2,r)-r_mid(2))+q(3,iq)*(coord(3,r)-r_mid(3))
               epowqr=exp(iqfac*qdr)*nainv
               it=aproj(r)
               wA(:,it)=wA(:,it)+epowqr*SA(:,r,l)
            end do
            SC(:,:,iq)=SC(:,:,iq) + wA(:,:)
         end do
      end do
      !$omp end parallel do

   end subroutine corr_kernel_proj



   subroutine corr_kernel_time(nq,nw,nelem,cc,dt,sc_window_fun,corr_kt,corr_kw)
      !
      implicit none
      !
      integer, intent(in) :: nq
      integer, intent(in) :: nw
      integer, intent(in) :: nelem
      type(corr_t), intent(in) :: cc
      real(dblprec), dimension(cc%sc_max_nstep), intent(in) :: dt
      integer, intent(in)             ::  sc_window_fun

      complex(dblprec), dimension(nelem*nq,cc%sc_max_nstep), intent(in) :: corr_kt
      complex(dblprec), dimension(nelem*nq,nw), intent(out) :: corr_kw

      integer :: iw, step, tidx
      complex(dblprec) :: i, tt, epowwt, wfac


      wfac=1.0_dblprec
      corr_kw=0.0_dblprec
      i=(0.0_dblprec,1.0_dblprec)
      tidx=cc%sc_max_nstep


      !$omp parallel do default(shared) private(iw,step,tt,epowwt) schedule(static)
      do iw=1,nw
         !DIR$ vector aligned
         !do step=1,min(cc%sc_tidx,cc%sc_nstep)
         do step=1,tidx
            tt=i*wfac*(step-1)*dt(step)
            ! Apply window function as determined by 'sc_window_fun'
            epowwt=exp(cc%w(iw)*tt)*sc_window_fac(sc_window_fun,step,tidx)
            !epowwt=exp(cc%w(iw)*tt)*sc_window_fac(sc_window_fun,step,min(cc%sc_tidx,cc%sc_nstep))
            !epowwt=exp(w(iw)*tt)*sc_window_fac(sc_window_fun,step,sc_tidx)
            !
            corr_kw(:,iw)=corr_kw(:,iw)+epowwt*corr_kt(:,step)
         enddo
         !corr_kw(:,iw)=corr_kw(:,iw)*conjg(corr_kw(:,iw))
      enddo
      !$omp end parallel do

   end subroutine corr_kernel_time

   subroutine combine_corr_scalar(nq, nelem, nw, corr_A, corr_B, corr_C)
      !
      implicit none
      !
      integer, intent(in) :: nq
      integer, intent(in) :: nelem
      integer, intent(in) :: nw

      complex(dblprec), dimension(nelem,nq,nw), intent(in) :: corr_A
      complex(dblprec), dimension(nelem,nq,nw), intent(in) :: corr_B
      complex(dblprec), dimension(nelem,nq,nw), intent(inout) :: corr_C

      integer :: iq, iw

      !$omp parallel do default(shared) private(iq,iw) collapse(2)
      do iw=1,nw
         do iq=1,nq
            corr_C(:,iq,iw) = corr_C(:,iq,iw) + corr_A(:,iq,iw)*conjg(corr_B(:,iq,iw))
         enddo
      enddo
      !$omp end parallel do

   end subroutine combine_corr_scalar

   subroutine combine_corr_tensor(nq, nelem, nw, corr_A, corr_B, corr_C)
      !
      implicit none
      !
      integer, intent(in) :: nq
      integer, intent(in) :: nelem

      integer, intent(in) :: nw

      complex(dblprec), dimension(nelem,nq,nw), intent(in) :: corr_A
      complex(dblprec), dimension(nelem,nq,nw), intent(in) :: corr_B
      complex(dblprec), dimension(nelem,nelem,nq,nw), intent(inout) :: corr_C

      integer :: iq, iw, ielem

      !$omp parallel do default(shared) private(iq,iw,ielem)
      do iw=1,nw
         do iq=1,nq
            do ielem=1,nelem
               corr_C(:,ielem,iq,iw) = corr_C(:,ielem,iq,iw) + corr_A(ielem,iq,iw)*conjg(corr_B(:,iq,iw))
            end do
         enddo
      enddo
      !$omp end parallel do


   end subroutine combine_corr_tensor


   subroutine combine_corr_proj_scalar(nt, nq, nelem, nw, corr_A, corr_B, corr_C)
      !     combine_corr_proj_scalar(nt, nq, 3, cc, cc%m_kt, cc%m_kt_proj, corr_kt_proj)
      !
      implicit none
      !
      integer, intent(in) :: nt
      integer, intent(in) :: nq
      integer, intent(in) :: nelem

      integer, intent(in) :: nw

      complex(dblprec), dimension(nelem,nt,nq,nw), intent(in) :: corr_A
      complex(dblprec), dimension(nelem,nt,nq,nw), intent(in) :: corr_B
      complex(dblprec), dimension(nelem,nt,nt,nq,nw), intent(out) :: corr_C

      integer :: iw,iq
      integer :: it, jt

      !$omp parallel do default(shared) private(iw,iq,it,jt)
      do iw=1,nw
         ! Change to corr_SA and corr_SB for two-component correlation
         do iq=1,nq
            do it=1,nt
               do jt=1,nt
                  corr_C(:,it,jt,iq,iw) = corr_C(:,it,jt,iq,iw) + corr_A(:,it,iq,iw)*conjg(corr_B(:,jt,iq,iw))
               end do
            end do
         end do
      enddo
      !$omp end parallel do

   end subroutine combine_corr_proj_scalar

   subroutine combine_corr_proj_tensor(nt, nq, nelem, nw, corr_A, corr_B, corr_C)
      !
      implicit none
      !
      integer, intent(in) :: nt
      integer, intent(in) :: nq
      integer, intent(in) :: nelem

      integer, intent(in) :: nw

      complex(dblprec), dimension(nelem,nt,nq,nw), intent(in) :: corr_A
      complex(dblprec), dimension(nelem,nt,nq,nw), intent(in) :: corr_B
      complex(dblprec), dimension(nelem,nelem,nt,nt,nq,nw), intent(out) :: corr_C

      integer :: iw,iq
      integer :: it, jt, ielem

      !$omp parallel do default(shared) private(iw,iq,it,jt)
      do iw=1,nw
         ! Change to corr_SA and corr_SB for two-component correlation
         do iq=1,nq
            do it=1,nt
               do jt=1,nt
                  do ielem=1,nelem
                     corr_C(:,ielem,it,jt,iq,iw) = corr_C(:,ielem,it,jt,iq,iw) + corr_A(ielem,it,iq,iw)*conjg(corr_B(:,jt,iq,iw))
                  end do
               end do
            end do
         end do
      enddo
      !$omp end parallel do

   end subroutine combine_corr_proj_tensor

end module Correlation_kernels
