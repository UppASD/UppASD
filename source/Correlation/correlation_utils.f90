module Correlation_utils
   use Parameters
   use Profiling
   use Correlation_type
   use Qvectors
   use Omegas
   !
   implicit none
   !
   character(len=2) :: do_conv      !< Flag for convolution of Magnon DOS Lorentzian(LY/LG/LW/N) Gaussian(GY/GW/GQ/N)
   real(dblprec) :: sigma_q         !< Sigma parameter in Q for the Gaussian convolution
   real(dblprec) :: sigma_w         !< Sigma parameter in W for the Gaussian convolution
   real(dblprec) :: sigma_t         !< Sigma parameter in t for the Gaussian convolution
   real(dblprec) :: LWfactor        !< Gamma parameter in W for Lorentzian convolution
   real(dblprec) :: LQfactor        !< Gamma parameter in Q for Lorentzian convolution

   public

contains

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_sqw_conv
   !> @brief Perform convolutions of the reciprocal correlation function sqw
   !---------------------------------------------------------------------------------
   subroutine calc_sqw_conv(nelem,nw,corr_cv)

      use Constants, only : pi

      implicit none

      integer, intent(in) :: nelem
      integer, intent(in) :: nw
      complex(dblprec), dimension(nelem,nq,nw), intent(inout) :: corr_cv

      integer :: iq,iw,j,t, conv_range_w, conv_range_q
      integer :: u, i_stat, i_all
      complex(dblprec) :: conv_cutoff_w, facw1, facw2, sfacww ! Variables for the convolution in frequency
      complex(dblprec) :: conv_cutoff_q, facq1, facq2, sfacqq ! Variables for the convolution in qpoints
      complex(dblprec) :: qq,ww
      complex(dblprec), dimension(:,:,:), allocatable  :: corr_tmp


      facw1 = 1.0_dblprec; facw2 = 1.0_dblprec
      facq1 = 1.0_dblprec; facq2 = 1.0_dblprec
      conv_range_w=0;conv_range_q=0;
      ! Variables for the convolutions in frequency and qpoints domain, this are optional for the magnon DOS calculations
      ! This are the variables for qpoints convolutions
      if (do_conv=='GQ') then
         ! Variables for the reciprocal space for Gaussian function
         facq1 = -1.0_dblprec/(2.0_dblprec*sigma_q**2)
         facq2 = 1.0_dblprec/(sqrt(2.0_dblprec*pi)*sigma_q)
         conv_cutoff_q=0.001_dblprec
         conv_range_q=int(sqrt(log(conv_cutoff_q)/facq1)+0.5_dblprec)

      else if (do_conv=='LQ') then
         ! Variables for the reciprocal space for Lorentz distribution
         facq1=-1.0_dblprec/(LQfactor**2)
         facq2=LQfactor/(Nq*pi)
         conv_cutoff_q=0.01_dblprec
         conv_range_q=int(sqrt(log(conv_cutoff_q)/facq1)+0.5_dblprec)

      else if (do_conv=='GW') then
         ! Variables for the frequency convolutions for Gaussian function
         facw1 = -1.0_dblprec/(2*sigma_w**2)
         facw2 = 1.0_dblprec/(sqrt(2*pi)*sigma_w)
         conv_cutoff_w=0.01_dblprec
         conv_range_w=int(sqrt(log(conv_cutoff_w)/facw1)+0.5_dblprec)

      else if (do_conv=='LW') then
         ! Variables for the frequency convolution with Lorentz distribution
         facw1=-1.0_dblprec/(LWfactor**2)
         facw2=LWfactor/(Nw*pi)
         conv_cutoff_w=0.01_dblprec
         conv_range_w=int(sqrt(log(conv_cutoff_w)/facw1)+0.5_dblprec)

      else if (do_conv=='GY') then
         ! Variables for both qpoints and frequencies convolutions for Gaussian function
         facq1 = -1.0_dblprec/(2*sigma_q**2)
         facw1 = -1.0_dblprec/(2*sigma_w**2)
         facq2 = 1.0_dblprec/(sqrt(2*pi)*sigma_q)
         facw2 = 1.0_dblprec/(sqrt(2*pi)*sigma_w)
         conv_cutoff_q=0.001_dblprec
         conv_cutoff_w=0.001_dblprec
         conv_range_q=int(sqrt(log(conv_cutoff_q)/facq1)+0.5_dblprec)
         conv_range_w=int(sqrt(log(conv_cutoff_w)/facw1)+0.5_dblprec)

      else if (do_conv=='LY') then
         ! Variables for both qpoints and frequencies convolutions for Lorentz distribution
         facq1 = -1.0_dblprec/(LQfactor**2)
         facw1 = -1.0_dblprec/(LWfactor**2)
         facq2 = LQfactor/(Nq*pi)
         facw2 = LWfactor/(Nw*pi)
         conv_cutoff_q=0.01_dblprec
         conv_cutoff_w=0.01_dblprec
         conv_range_q=int(sqrt(log(conv_cutoff_q)/facq1)+0.5_dblprec)
         conv_range_w=int(sqrt(log(conv_cutoff_q)/facw1)+0.5_dblprec)
      else

         return

      endif

      allocate(corr_tmp(nelem,nq,nw),stat=i_stat)
      call memocc(i_stat,product(shape(corr_tmp))*kind(corr_tmp),'calc_sqw_conv','calc_sqw_conv')

      corr_tmp=corr_cv
      corr_cv=0.0_dblprec

      ! Do the convolution of the corr_cv to smooth out the magnon DOS
      if (do_conv=='GQ') then
         !write(*,*) '-------------------------> CONV_GQ',facq1,facq2
         !$omp parallel do default(shared) private(iq,iw,u,qq,sfacqq)
         do iq=1,nq
            do iw=1,nw
               ! Convolution with Gaussian resolution function
               do u = max(1,iq-conv_range_q), min(nq-1,iq+conv_range_q)
                  qq=(iq-u)
                  sfacqq=exp(facq1*qq**2)*facq2 ! Factor which controls the convolution
                  !if (iw==1.and.iq==nq/2) write(100,*) u,abs(sfacqq)
                  ! Convolution over the reciprocal space
                  corr_cv(:,iq,iw)=corr_cv(:,iq,iw)+corr_tmp(:,u,iw)*sfacqq
               enddo
            enddo
         enddo
         !$omp end parallel do
      else if (do_conv=='GW') then
         do iq=1,Nq
            !$omp parallel do default(shared) private(iw,j,t,sfacww)
            do iw=1, nw
               !if (iw==1.and.iq==1) write(*,*) '-------------------------> CONV_GW'
               ! Convolution with Gaussian resolution function
               do j = max(1,iw-conv_range_w), min(nw,iw+conv_range_w)
                  ww=(iw-j)
                  sfacww=exp(facw1*ww**2)*facw2 ! This is the parameter that controls the convolution (Should the sum be over the absolute values?)
                  !if (iq==10.and.iw==200) write(100,*) j,abs(sfacww)
                  ! Convolution over frequency space with a gaussian resolution function
                  corr_cv(:,iq,iw)=corr_cv(:,iq,iw)+corr_tmp(:,iq,j)*sfacww
               enddo
            enddo
            !$omp end parallel do
         enddo
      else if (do_conv=='GY') then
         do iq=1,Nq
            !$omp parallel do default(shared) private(iw,u,qq,sfacqq)
            do iw=1, nw
               !if (iw==1.and.iq==1) write(*,*) '-------------------------> CONV_GY',facq1,facq2
               ! Convolution with Gaussian resolution function
               do u = max(1,iq-conv_range_q), min(nq,iq+conv_range_q)
                  qq=(iq-u)
                  sfacqq=exp(facq1*qq**2)*facq2 ! This parameter controls the convolution
                  ! Convolution over the qpoints
                  corr_cv(:,iq,iw)=corr_cv(:,iq,iw)+corr_tmp(:,u,iw)*sfacqq
               enddo
            enddo
            !$omp end parallel do
         enddo
         corr_tmp=corr_cv
         corr_cv=0.0_dblprec
         do iw=1, nw
            !$omp parallel do default(shared) private(iq,j,t,sfacww)
            do iq=1,Nq
               !if (iw==1.and.iq==1) write(*,*) '-------------------------> CONV_GY',facw1,facw2
               ! Convolution with Gaussian resolution function
               do j = max(1,iw-conv_range_w), min(nw-1,iw+conv_range_w)
                  ww=(j-iw)
                  sfacww=exp(facw1*ww**2)*facw2 ! This is the parameter that controls the convolution (Should the sum be over the absolute values?)
                  ! Convolution over frequency space with a gaussian resolution function
                  corr_cv(:,iq,iw)=corr_cv(:,iq,iw)+corr_tmp(:,iq,j)*sfacww
               enddo
            enddo
            !$omp end parallel do
         enddo
      else if (do_conv=='LY') then
         do iq=1,Nq
            !$omp parallel do default(shared) private(iw,u,qq,sfacqq)
            do iw=1, nw/2
               !if (iw==1.and.iq==1) write(*,*) '-------------------------> CONV_LY'
               ! Convolution with Gaussian resolution function
               do u = max(1,iq-conv_range_q), min(nq-1,iq+conv_range_q)
                  qq=(iq-u)
                  sfacqq=facq1/(LQfactor**2 +qq**2) ! This parameter controls the convolution
                  ! Convolution over the qpoints
                  corr_cv(:,iq,iw)=corr_cv(:,iq,iw)+corr_tmp(:,iq,iw)*sfacqq
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
                  corr_cv(:,iq,iw)=corr_cv(:,iq,iw)+corr_tmp(:,iq,iw)*sfacww
               enddo
            enddo
            !$omp end parallel do
         enddo
      else if (do_conv=='LW') then
         do iq=1,Nq
            !$omp parallel do default(shared) private(iw,j,t,sfacww)
            do iw=1, Nw/2
               !if (iw==1.and.iq==1) write(*,*) '--------------------------> CONV_LW'
               do j = max(1,iw-conv_range_w), min(nw-1,iw+conv_range_w)
                  t=(iw-j)
                  sfacww=facw1/(LWfactor**2+t**2)
                  corr_cv(:,iq,iw)=corr_cv(:,iq,iw)+corr_tmp(:,iq,iw)*sfacww
               enddo
            enddo
            !$omp end parallel do
         enddo
      else if (do_conv=='LQ') then
         do iq=1,Nq
            !$omp parallel do default(shared) private(iw,u,qq,sfacqq)
            do iw=1, nw/2
               !if (iw==1.and.iq==1) write(*,*) '-------------------------> CONV_LQ'
               ! Convolution with Gaussian resolution function
               do u = max(1,iq-conv_range_q), min(nq-1,iq+conv_range_q)
                  qq=(iq-u)
                  sfacqq=facq1/(LQfactor**2 +qq**2) ! This parameter controls the convolution
                  ! Convolution over the qpoints
                  corr_cv(:,iq,iw)=corr_cv(:,iq,iw)+corr_tmp(:,iq,iw)*sfacqq
               enddo
            enddo
            !$omp end parallel do
         enddo
      endif

      i_all=-product(shape(corr_tmp))*kind(corr_tmp)
      deallocate(corr_tmp,stat=i_stat)
      call memocc(i_stat,i_all,'corr_tmp','calc_tmp')

      return

   end subroutine calc_sqw_conv

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_sqt_conv
   !> @brief Perform convolutions of the reciprocal correlation function sqt
   !---------------------------------------------------------------------------------
   subroutine calc_sqt_conv(nelem,cc,dt,corr_cv)

      use Constants, only : pi

      implicit none

      integer, intent(in) :: nelem
      type(corr_t), intent(in) :: cc
      real(dblprec), dimension(cc%sc_max_nstep), intent(in) :: dt
      complex(dblprec), dimension(nelem*nq,cc%sc_max_nstep), intent(inout) :: corr_cv

      complex(dblprec), dimension(:,:,:),     allocatable :: a_kt             ! Correlation for G(k,t)
      real(dblprec) :: sigma_t, epowwt, tt
      integer :: step, i_stat, i_all

      ! Convolution over frequencies. Done in real-time space by means of the convolution theorem
      if (do_conv .eq. 'G') then

         ! Allocate arrays
         allocate(a_kt(3,nq,cc%sc_max_nstep),stat=i_stat)
         call memocc(i_stat,product(shape(a_kt))*kind(a_kt),'a_kt','calc_sqt_conv')

         sigma_t=sqrt(1.0_dblprec/(2.0_dblprec*pi*pi*sigma_w*sigma_w))

         !$omp parallel do default(shared) private(step,tt,epowwt) schedule(static)
         do step=1,min(cc%sc_tidx,cc%sc_nstep)
            tt=dt(step)*(step-min(cc%sc_tidx,cc%sc_nstep)*0.5_dblprec)*sigma_w
            epowwt=exp(-tt*tt/2.0_dblprec)*sqrt(2.0_dblprec*pi)
            !!tt=dt(step)*(step-min(sc_tidx,sc_nstep)*0.5_dblprec)/sigma_t
            !!epowwt=exp(-tt*tt/2.0_dblprec)*2.0_dblprec
            corr_cv(:,step)=epowwt*corr_cv(:,step)
         enddo
         !$omp end parallel do

         ! Deallocate arrays
         i_all=-product(shape(a_kt))*kind(a_kt)
         deallocate(a_kt,stat=i_stat)
         call memocc(i_stat,i_all,'a_kt','calc_sqt_conv')

      end if

      return

   end subroutine calc_sqt_conv

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_mavrg_vec
   !> @brief Calculate current average magnetization by components for the connected \f$\mathbf{S}\left(\mathbf{q},\omega\right)\f$
   !---------------------------------------------------------------------------------
   subroutine calc_mavrg_vec(Natom, Mensemble, emomM, mavrg_vec, mavg_axis)

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
      m=0.0_dblprec

      do k=1,Mensemble
         !$omp parallel do private(i) default(shared) schedule(static) reduction(+:m)
         do i=1, Natom
            m(:,k) = m(:,k) + emomM(:,i,k)
         end do
         !$omp end parallel do
         mavg_axis(:,k)=m(:,k)/Natom
      end do
      mavrg_vec(1)=sum(m(1,:))/Mensemble/Natom
      mavrg_vec(2)=sum(m(2,:))/Mensemble/Natom
      mavrg_vec(3)=sum(m(3,:))/Mensemble/Natom

   end subroutine calc_mavrg_vec


   !---------------------------------------------------------------------------------
   ! SUBROUTINE: find_local_rotmat
   !> @brief Finding the local rotational matrix for the local quantization axis
   !---------------------------------------------------------------------------------
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
      z_vec(1)=0.0_dblprec
      z_vec(2)=0.0_dblprec
      z_vec(3)=1.0_dblprec
      !
      eye=0.0_dblprec;eye(1,1)=1.0_dblprec;eye(2,2)=1.0_dblprec;eye(3,3)=1.0_dblprec
      !
      !$omp parallel do default(shared) private(i,m_norm,m_vec,k_vec,k_norm,theta,K_mat,KK_mat)
      do i=1,ldim
         m_norm=sum(mom_in(:,i)*mom_in(:,i))**0.5_dblprec
         m_vec=mom_in(:,i)/m_norm
         k_vec(1)=m_vec(2)*z_vec(3)-m_vec(3)*z_vec(2)
         k_vec(2)=m_vec(3)*z_vec(1)-m_vec(1)*z_vec(3)
         k_vec(3)=m_vec(1)*z_vec(2)-m_vec(2)*z_vec(1)
         k_norm=sum(k_vec*k_vec)**0.5_dblprec+1.0d-14
         k_vec=k_vec/k_norm
         !
         theta=-acos(m_vec(1)*z_vec(1)+m_vec(2)*z_vec(2)+m_vec(3)*z_vec(3))
         !
         K_mat(1,1)=0.0_dblprec;K_mat(2,2)=0.0_dblprec;K_mat(3,3)=0.0_dblprec
         K_mat(1,2)=-k_vec(3);K_mat(2,1)= k_vec(3)
         K_mat(1,3)= k_vec(2);K_mat(3,1)=-k_vec(2)
         K_mat(2,3)=-k_vec(1);K_mat(3,2)= k_vec(1)
         !
         KK_mat=matmul(K_mat,K_mat)
         !
         mat_out(:,:,i)=eye+K_mat*sin(theta)+KK_mat*(1.0_dblprec-cos(theta))
         !
      end do
      !$omp end parallel do
      !
      return
      !
   end subroutine find_local_rotmat

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
      rcenter(1)=(maxval(coord(1,:))+minval(coord(1,:)))*0.5_dblprec
      rcenter(2)=(maxval(coord(2,:))+minval(coord(2,:)))*0.5_dblprec
      rcenter(3)=(maxval(coord(3,:))+minval(coord(3,:)))*0.5_dblprec
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

   subroutine project_local_axis(Natom,Mensemble,cc,SA)
      !
      use FieldData, only : beff
      !
      implicit none
      !
      integer, intent(in) :: Natom
      integer, intent(in) :: Mensemble
      !
      type(corr_t),intent(inout) :: cc
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: SA
      !
      integer :: l,r
      real(dblprec) :: beff_norm
      real(dblprec), dimension(3) :: beff_hat
      real(dblprec), dimension(3,3) :: local_rotmat

      if(cc%do_sc_local_axis/='Y'.and.cc%do_sc_local_axis/='B') then
         !$omp parallel do default(shared) private(l,r) schedule(static) collapse(2)
         do l=1,Mensemble
            do r=1,Natom
               cc%m_loc(:,r,l)=SA(:,r,l)
            end do
         end do
         !$omp end parallel do
      else if (cc%do_sc_local_axis=='B') then
         !$omp parallel do default(shared) private(l,r,beff_hat) schedule(static) collapse(2)
         do l=1,Mensemble
            do r=1,Natom
               beff_norm=sqrt(sum(beff(:,r,l)**2))+1.0d-12
               beff_hat=beff(:,r,l)/beff_norm
               cc%m_loc(1,r,l)=beff_hat(2)*SA(3,r,l)-beff_hat(3)*SA(2,r,l)
               cc%m_loc(2,r,l)=beff_hat(3)*SA(1,r,l)-beff_hat(1)*SA(3,r,l)
               cc%m_loc(3,r,l)=beff_hat(1)*SA(2,r,l)-beff_hat(2)*SA(1,r,l)
            end do
         end do
         !$omp end parallel do
      else ! do_sc_local_axis==Y
         !! Sample current moment directions and weight in for average 
         if(cc%sc_local_axis_mix>0.0_dblprec) then
            !$omp parallel do schedule(static) default(shared) private(l,r,local_rotmat) collapse(2)
            do l=1,Mensemble
               do r=1,Natom
                  cc%mavg_local_axis(:,r,l)=cc%mavg_local_axis(:,r,l)*(1.0_dblprec-cc%sc_local_axis_mix)+ &
                     SA(:,r,l)*cc%sc_local_axis_mix
               end do
            end do
            !$omp end parallel do

            ! Find the proper rotation axes (local and individual) 
            call find_local_rotmat(Natom*Mensemble,cc%mavg_local_axis,cc%mavg_local_rotmat)
         end if

         !$omp parallel do default(shared) private(l,r,local_rotmat) schedule(static) collapse(2)
         do l=1,Mensemble
            do r=1,Natom
               local_rotmat=cc%mavg_local_rotmat(:,:,r,l)
               !call gramms(mavg_local_axis(1:3,r,l),local_rotmat,1)
               cc%m_loc(1,r,l)=local_rotmat(1,1)*SA(1,r,l)+local_rotmat(2,1)*SA(2,r,l)+local_rotmat(3,1)*SA(3,r,l)
               cc%m_loc(2,r,l)=local_rotmat(1,2)*SA(1,r,l)+local_rotmat(2,2)*SA(2,r,l)+local_rotmat(3,2)*SA(3,r,l)
               cc%m_loc(3,r,l)=local_rotmat(1,3)*SA(1,r,l)+local_rotmat(2,3)*SA(2,r,l)+local_rotmat(3,3)*SA(3,r,l)
            end do
         end do
         !$omp end parallel do
      end if

      if (cc%do_connected=='Y') then
         call calc_mavrg_vec(Natom,Mensemble,cc%m_loc,cc%SA_avrg,cc%SA_axis)
         !$omp parallel do default(shared) private(l,r,beff_hat) schedule(static) collapse(2)
         do l=1,Mensemble
            do r=1,Natom
               cc%m_loc(:,r,l)=cc%m_loc(:,r,l)-cc%SA_axis(:,l)
            end do
         end do
         !$omp end parallel do
      else
         cc%SA_avrg=0.0_dblprec
      end if

   end subroutine project_local_axis


   subroutine extract_dos(nq,nw,nelem,w,corr_out,filn)
      !
      use Constants
      !
      implicit none
      !
      integer, intent(in) :: nq
      integer, intent(in) :: nw
      integer, intent(in) :: nelem
      real(dblprec), dimension(nw) :: w
      complex(dblprec), dimension(nelem,nq,nw), intent(in) :: corr_out
      character(len=30), intent(in) :: filn !< Name of simulation
      !
      real(dblprec), dimension(:,:), allocatable :: magnon_dos !< Magnon density of states
      integer :: iq,iw, ielem, i_stat, i_all
      real(dblprec) :: deltae
      real(dblprec), dimension(nelem) :: snorm, sint
      real(dblprec) :: sintt
      !
      !
      ! Calculate and write the magnon DOS
      allocate(magnon_dos(nelem+1,max(1,nw/2)),stat=i_stat)
      call memocc(i_stat,product(shape(magnon_dos))*kind(magnon_dos),'magnon_dos','calc_gkt')
      magnon_dos=0.0_dblprec
      !
      if (qpoints=='I'.or.qpoints=='B') then
         do iq=1,Nq
            do ielem=1,nelem
               snorm(ielem)=1.0_dblprec/maxval(abs(corr_out(ielem,iq,:)))
            end do
            !              snorm(2)=1.0_dblprec/maxval(abs(corr_out(2,iq,:)))
            !              snorm(3)=1.0_dblprec/maxval(abs(corr_out(3,iq,:)))
            do iw=1,max(1,Nw/2)
               magnon_dos(1:nelem,iw)=magnon_dos(1:nelem,iw) &
                  +abs(corr_out(1:nelem,iq,iw))*q_weight(iq)*snorm(:) ! Put a weight for the IR BZ points
            end do
         end do
      else
         do iq=1,Nq
            do ielem=1,nelem
               snorm(ielem)=1.0_dblprec/maxval(abs(corr_out(ielem,iq,:)))
            end do
            !              snorm(2)=1.0_dblprec/maxval(abs(corr_out(2,iq,:)))
            !              snorm(3)=1.0_dblprec/maxval(abs(corr_out(3,iq,:)))
            do iw=1,max(1,Nw/2)
               magnon_dos(1:nelem,iw)=magnon_dos(1:nelem,iw) &
                  +abs(corr_out(1:nelem,iq,iw))*snorm(:) ! Put a weight for the IR BZ points
            end do
         end do
      endif

      ! Resetting the baseline to reduce effect of noise
      do ielem=1,nelem
         magnon_dos(ielem,:)=magnon_dos(ielem,:)-minval(magnon_dos(ielem,:))
      end do



      !if (cc%do_sc_tens=='Y') then
      !   do iw=1,Nw/2
      !         !magnon_dos(nelem+1,iw)=magnon_dos(nelem+1,iw)+abs(sum(sqwintensity(:,iw)))
      !   end do
      !else
      do iw=1,max(1,Nw/2)
         magnon_dos(nelem+1,iw)=sqrt(sum(magnon_dos(1:nelem,iw)**2))
      end do
      !end if

      deltae=hbar_mev*(w(2)-w(1))
      ! Renormalize DOS for energy in meV 
      sint=0.0_dblprec
      do iw=1,max(1,Nw/2)
         sint=sint+deltae*magnon_dos(1:nelem,iw)
      end do
      sintt=0.0_dblprec
      do iw=1,max(1,Nw/2)
         sintt=sintt+deltae*magnon_dos(nelem+1,iw)
      end do

      open(ofileno, file=filn)
      if (nelem==3) then
         write(ofileno,'(a)') &
               "   E(mev)          D(S(q,E)_x)     D(S(q,E)_y)     D(S(q,E)_z)     D(|S(q,E)|)"
      else
         write(ofileno,'(a)') &
               "   E(mev)          D(S(q,E))       D(|S(q,E)|)"
      end if
      do iw=1,max(1,Nw/2)
         write (ofileno,10006) hbar_mev*w(iw),magnon_dos(1:nelem,iw)/sint,magnon_dos(nelem+1,iw)/sintt
      end do
      close(ofileno)
      !
      i_all=-product(shape(magnon_dos))*kind(magnon_dos)
      deallocate(magnon_dos,stat=i_stat)
      call memocc(i_stat,i_all,'magnon_dos','calc_gkt')

      10006 format (14G16.8)

   end subroutine extract_dos

   !---------------------------------------------------------------------------------
   ! SUBROUTINE allocate_deltatcorr
   !> @brief Allocate adaptivetime step for correlation
   !---------------------------------------------------------------------------------
   subroutine allocate_deltatcorr(allocate_flag,cc)

      implicit none

      logical, intent(in) :: allocate_flag !< Allocate/deallocate
      type(corr_t) :: cc 
      !integer, intent(in) :: Nstep !< Total number of simulation steps

      integer :: i_stat, i_all

      if(allocate_flag) then
         if(.not.allocated(cc%deltat_corr)) then
            allocate(cc%deltat_corr(cc%sc_nstep+cc%sc_naverages+1),stat=i_stat)
            call memocc(i_stat,product(shape(cc%deltat_corr))*kind(cc%deltat_corr),'deltat_corr','allocate_deltatcorr')
         end if
         if(.not.allocated(cc%scstep_arr)) then
            allocate(cc%scstep_arr(cc%sc_nstep+cc%sc_naverages+1),stat=i_stat)
            call memocc(i_stat,product(shape(cc%scstep_arr))*kind(cc%scstep_arr),'scstep_arr','allocate_deltatcorr')
         end if
      else
         if(allocated(cc%deltat_corr)) then
            i_all=-product(shape(cc%deltat_corr))*kind(cc%deltat_corr)
            deallocate(cc%deltat_corr,stat=i_stat)
            call memocc(i_stat,i_all,'deltat_corr','allocate_deltatcorr')
         end if
         if(allocated(cc%scstep_arr)) then
            i_all=-product(shape(cc%scstep_arr))*kind(cc%scstep_arr)
            deallocate(cc%scstep_arr,stat=i_stat)
            call memocc(i_stat,i_all,'scstep_arr','allocate_deltatcorr')
         end if
      end if

   end subroutine allocate_deltatcorr

   subroutine deallocate_gkw(cc)
      !
      implicit none
      !
      type(corr_t), intent(inout) :: cc !< Derived type for correlation data
      !
      integer :: i_all, i_stat
      !


      i_all=-product(shape(cc%w))*kind(cc%w)
      deallocate(cc%w,stat=i_stat)
      call memocc(i_stat,i_all,'cc%w','deallocate_gkw')

      i_all=-product(shape(cc%m_kw))*kind(cc%m_kw)
      deallocate(cc%m_kw,stat=i_stat)
      call memocc(i_stat,i_all,'cc%m_kw','deallocate_gkw')

      i_all=-product(shape(cc%dt))*kind(cc%dt)
      deallocate(cc%dt,stat=i_stat)
      call memocc(i_stat,i_all,'cc%dt','deallocate_gkw')

      i_all=-product(shape(cc%time))*kind(cc%time)
      deallocate(cc%time,stat=i_stat)
      call memocc(i_stat,i_all,'cc%time','deallocate_gkw')

      i_all=-product(shape(cc%m_kt))*kind(cc%m_kt)
      deallocate(cc%m_kt,stat=i_stat)
      call memocc(i_stat,i_all,'cc%m_kt','calc_gkt')

      if(cc%do_proj=='Q'.or.cc%do_proj=='T'.or.cc%do_proj=='Y') then
         i_all=-product(shape(cc%m_kt_proj))*kind(cc%m_kt_proj)
         deallocate(cc%m_kt_proj,stat=i_stat)
         call memocc(i_stat,i_all,'cc%m_kt_proj','calc_gkt')
      end if

      if(cc%do_projch=='Q'.or.cc%do_projch=='T'.or.cc%do_projch=='Y') then
         i_all=-product(shape(cc%m_kt_projch))*kind(cc%m_kt_projch)
         deallocate(cc%m_kt_projch,stat=i_stat)
         call memocc(i_stat,i_all,'cc%m_kt_projch','calc_gkt')
      end if

      if(cc%do_proj=='Q'.or.cc%do_proj=='Y') then

         i_all=-product(shape(cc%m_kw_proj))*kind(cc%m_kw_proj)
         deallocate(cc%m_kw_proj,stat=i_stat)
         call memocc(i_stat,i_all,'cc%m_kw_proj','deallocate_gkw')

      end if

      if(cc%do_projch=='Q'.or.cc%do_projch=='Y') then

         i_all=-product(shape(cc%m_kw_projch))*kind(cc%m_kw_projch)
         deallocate(cc%m_kw_projch,stat=i_stat)
         call memocc(i_stat,i_all,'cc%m_kw_projch','deallocate_gkw')

      end if


      !
   end subroutine deallocate_gkw

   subroutine allocate_corr(Natom,Mensemble,cc,flag)
      !
      implicit none
      !
      integer, intent(in) :: Natom
      integer, intent(in) :: Mensemble
      type(corr_t), intent(inout) :: cc
      integer, intent(in) :: flag

      !
      integer :: i_stat, i_all
      !
      if (flag>0) then
         allocate(cc%SA_axis(3,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(cc%SA_axis))*kind(cc%SA_axis),'cc%SA_axis','allocate_corr')
         cc%SA_axis=0.0_dblprec
         !
         !allocate(cc%mort_axis(3,3,Mensemble),stat=i_stat)
         !call memocc(i_stat,product(shape(cc%mort_axis))*kind(cc%mort_axis),'mort_axis','allocate_corr')
         !cc%mort_axis=0.0_dblprec
         !
         if (cc%do_sc_local_axis=='Y'.or.cc%do_sc_local_axis=='B') then

            allocate(cc%mavg_local_axis(3,Natom,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(cc%mavg_local_axis))*kind(cc%mavg_local_axis),'cc%mavg_local_axis','allocate_corr')

            allocate(cc%mavg_local_rotmat(3,3,Natom,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(cc%mavg_local_rotmat))*kind(cc%mavg_local_rotmat), &
               'cc%mavg_local_rotmat','allocate_corr')

         end if
         !
         allocate(cc%m_loc(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(cc%m_loc))*kind(cc%m_loc),'cc%m_loc','allocate_corr')

      else
         i_all=-product(shape(cc%SA_axis))*kind(cc%SA_axis)
         deallocate(cc%SA_axis,stat=i_stat)
         call memocc(i_stat,i_all,'cc%SA_axis','allocate_corr')

         !
         if (cc%do_sc_local_axis=='Y'.or.cc%do_sc_local_axis=='B') then

            i_all=-product(shape(cc%mavg_local_axis))*kind(cc%mavg_local_axis)
            deallocate(cc%mavg_local_axis,stat=i_stat)
            call memocc(i_stat,i_all,'cc%mavg_local_axis','allocate_corr')

            i_all=-product(shape(cc%mavg_local_rotmat))*kind(cc%mavg_local_rotmat)
            deallocate(cc%mavg_local_rotmat,stat=i_stat)
            call memocc(i_stat,i_all,'cc%mavg_local_rotmat','allocate_corr')

         end if
         !
         !i_all=-product(shape(cc%mort_axis))*kind(cc%mort_axis)
         !deallocate(cc%mort_axis,stat=i_stat)
         !call memocc(i_stat,i_all,'mort_axis','allocate_corr')
         !
         i_all=-product(shape(cc%m_loc))*kind(cc%m_loc)
         deallocate(cc%m_loc,stat=i_stat)
         call memocc(i_stat,i_all,'cc%m_loc','allocate_corr')
      end if

   end subroutine allocate_corr

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_sr_init
   !> @brief Initialize real-space sampled correlation functions 
   !---------------------------------------------------------------------------------
   subroutine corr_sr_init(Natom, cc, coord, cr_flag)
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      type(corr_t) :: cc !< Correlation struct
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      integer, intent(in) :: cr_flag  !< Allocate or deallocate (1/-1)
      !
      integer :: i_stat,i_all

      if(cr_flag>=0) then
         ! First call, allocate and clear arrays
         
         ! Create neighbour list for truncated S(r) correlations
         call corr_get_nnlist(Natom,coord,cc%cr_list,cc%cr_vect,cc%cr_nnmax,1,cc%cr_cut)
         ! Identify identical r-vectors and create look-up table
         call corr_get_nnhist(Natom,cc%cr_list,cc%cr_vect,cc%cr_nnmax,cc%cr_nhist,cc%cr_uniq_vect,cc%cr_lut)
        
      else 

         i_all=-product(shape(cc%cr_vect))*kind(cc%cr_vect)
         deallocate(cc%cr_vect,stat=i_stat)
         call memocc(i_stat,i_all,'cr_vect','corr_sr_init')

         i_all=-product(shape(cc%cr_uniq_vect))*kind(cc%cr_uniq_vect)
         deallocate(cc%cr_uniq_vect,stat=i_stat)
         call memocc(i_stat,i_all,'cr_uniq_vect','corr_sr_init')

         i_all=-product(shape(cc%cr_list))*kind(cc%cr_list)
         deallocate(cc%cr_list,stat=i_stat)
         call memocc(i_stat,i_all,'cr_list','corr_sr_init')

         i_all=-product(shape(cc%cr_lut))*kind(cc%cr_lut)
         deallocate(cc%cr_lut,stat=i_stat)
         call memocc(i_stat,i_all,'cr_lut','corr_sr_init')

!         ! Real-space correlation
!         integer, dimension(:,:), allocatable :: cr_list
!         real(dblprec), dimension(:,:,:), allocatable :: cr_vect
!         real(dblprec) :: cr_cut
!         integer :: cr_nnmax
!         integer :: cr_nhist
!         integer :: cr_flag
!         real(dblprec), dimension(:), allocatable :: cr_hist
!         real(dblprec), dimension(:,:), allocatable :: cr_uniq_vect
!         integer, dimension(:,:), allocatable :: cr_lut

      end if
      return

   end subroutine corr_sr_init

   subroutine corr_get_nnlist(Natom,coord,cr_list,cr_vect,nnmax,flag,rcut)
      !
      use Parameters
      use Profiling
      use Math_functions, only : f_get_periodic_shifts

      implicit none
      !
      integer, intent(in) :: Natom
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      integer, dimension(:,:), allocatable :: cr_list
      real(dblprec), dimension(:,:,:), allocatable :: cr_vect
      integer, intent(out) :: nnmax
      real(dblprec), intent(in) :: rcut
      integer, intent(in) :: flag

      !
      integer :: nmax, amax, hmax
      integer :: i_atom,j_atom
      real(dblprec), dimension(3) :: rij
      real(dblprec), dimension(3,27) :: cell_shifts
      integer :: nshifts, ishift

      integer :: i_stat, i_all


      if (flag>0) then

         call f_get_periodic_shifts(nshifts,cell_shifts)
         nmax=0
         do i_atom=1,Natom
            amax=0
            do j_atom=1,Natom
               do ishift=1,nshifts
                  rij=coord(:,i_atom)-coord(:,j_atom)-cell_shifts(:,ishift)
                  if(norm2(rij)<=rcut) amax=amax+1
               end do
            end do
            nmax=max(nmax,amax)
         end do

         allocate(cr_list(nmax+1,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(cr_list))*kind(cr_list),'cr_list','corr_get_nnlist')
         allocate(cr_vect(3,nmax,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(cr_vect))*kind(cr_vect),'cr_vect','corr_get_nnlist')

         nnmax=nmax
         cr_list=0


         !$omp parallel do default(shared) private(i_atom,j_atom,amax,ishift,rij)
         do i_atom=1,Natom
            amax=0
            do j_atom=1,Natom
               do ishift=1,nshifts
                  rij=coord(:,i_atom)-coord(:,j_atom)-cell_shifts(:,ishift)
                  if(norm2(rij)<=rcut) then
                     amax=amax+1
                     cr_list(amax+1,i_atom)=j_atom
                     cr_vect(:,amax,i_atom)=rij
                  end if
               end do
            end do
            cr_list(1,i_atom)=amax
         end do
         !$omp end parallel do
         hmax=sum(cr_list(1,:))

      else if (flag<0) then
         i_all=-product(shape(cr_list))*kind(cr_list)
         deallocate(cr_list,stat=i_stat)
         call memocc(i_stat,i_all,'cr_list','corr_get_nnlist')

         i_all=-product(shape(cr_vect))*kind(cr_vect)
         deallocate(cr_vect,stat=i_stat)
         call memocc(i_stat,i_all,'cr_vect','corr_get_nnlist')
      end if

   end subroutine corr_get_nnlist

   subroutine corr_get_nnhist(Natom,cr_list,cr_vect,nnmax,cr_nhist,cr_uniq_vect,cr_lut)
      !  call corr_get_nnhist(Natom,cc%cr_list,cc%cr_vect,cc%cr_nnmax,cc%cr_nhist,cc%cr_uniq_vect,cc%cr_lut,1)
      ! sed -i.bak "s/cr_lut/cr_lut/g" Correlation/correlation_utils.f90
      ! sed -i.bak "s/cr_uniq_vect/cr_uniq_vect/g" Correlation/correlation_utils.f90
      ! sed -i.bak "s/cr_nhist/cr_nhist/g"  Correlation/correlation_utils.f90
      !
      use Parameters
      use Profiling

      implicit none
      !
      integer, intent(in) :: Natom
      integer, dimension(nnmax+1,Natom) :: cr_list
      real(dblprec), dimension(3,nnmax,Natom) :: cr_vect
      integer, intent(in) :: nnmax
      integer, intent(out) :: cr_nhist
      real(dblprec), dimension(:,:), allocatable :: cr_uniq_vect
      integer,  dimension(:,:), allocatable :: cr_lut

      !
      integer :: i_atom,n_idx, h_idx

      integer :: i_stat
      real(dblprec) :: eps = 1.0e-12_dblprec


      cr_nhist=nnmax*Natom
      call corr_unique_vec(cr_vect,3,cr_nhist,cr_uniq_vect)

      allocate(cr_lut(nnmax,Natom),stat=i_stat)
      call memocc(i_stat,product(shape(cr_lut))*kind(cr_lut),'nnlist_lut','corr_get_nnhist')
      cr_lut=0

      do i_atom=1,Natom
         n_idx=0
         do n_idx=1,cr_list(1,i_atom)
         
            h_idx=1
            do while(norm2(cr_vect(:,n_idx,i_atom)-cr_uniq_vect(:,h_idx))>eps)
               h_idx=h_idx+1
            end do
            cr_lut(n_idx,i_atom)=h_idx

         end do
      
      end do

   end subroutine corr_get_nnhist

   subroutine corr_unique_vec(cr_vect,x_dim,y_dim,cr_uniq_vect)
        !call corr_unique_vec(cr_vect,3,cr_nhist,cr_uniq_vect)
      !
      use Parameters
      use Profiling
      !
      implicit none
      !
      integer, intent(in) :: x_dim
      integer, intent(inout) :: y_dim
      real(dblprec), dimension(x_dim,y_dim), intent(in) :: cr_vect
      real(dblprec), dimension(:,:),allocatable,  intent(out) :: cr_uniq_vect
      !
      integer :: y_idx, u_max, u_idx
      real(dblprec), dimension(:), allocatable :: x_trial
      real(dblprec) :: eps = 1.0e-12_dblprec
      real(dblprec), dimension(:,:),allocatable :: tmp_arr
      logical :: unique

      integer :: i_stat
      !
      allocate(x_trial(x_dim))

      allocate(tmp_arr(x_dim,y_dim))

      tmp_arr=0.0_dblprec

      u_max=0

      do y_idx=1,y_dim
         
         x_trial=cr_vect(:,y_idx)
         unique=.true.
         do u_idx=1,u_max
            if ( norm2(x_trial-tmp_arr(:,u_idx))<eps)  unique=.false.
         end do

         if(unique) then
            u_max=u_max+1
            tmp_arr(:,u_max)=x_trial
         end if
      end do
      deallocate(x_trial)

      allocate(cr_uniq_vect(x_dim,u_max),stat=i_stat)
      call memocc(i_stat,product(shape(cr_uniq_vect))*kind(cr_uniq_vect),'cr_uniq_vect','corr_unique_vec')


      do u_idx=1,u_max
         cr_uniq_vect(:,u_idx)=tmp_arr(:,u_idx)
      end do
      deallocate(tmp_arr)

      y_dim=u_max

      return

   end subroutine corr_unique_vec

end module Correlation_utils
