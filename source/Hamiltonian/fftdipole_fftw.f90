!-------------------------------------------------------------------------------
! MODULE: fftdipole_fftw
!> @brief Routine for the calculation of the dipole-dipole interaction via the
!> FFT convolution method.
!> @details This routine calculates the dipole-dipole interacction making use of
!> the convolution theorem and the discrete FFT. The implementation is based in the
!> work by Moritz Sallerman https://github.com/MSallermann/fft-dipoles-proof-of-concept.
!> This module makes use of the FFTW library.
!> Notice that the FFTW makes use of 3D arrays which properly describe the
!> topology of the underlaying 3D lattice. Hence one must trasform between the
!> 1D array which contains all the data, and a temporal 3D work array.
!> @todo See the posibility of adding strides, this would eliminate the need of
!> passing between arrays of different shapes.
!> @todo See how to add the bindings and linking for the FFTW MKL interface
!> @author Jonathan Chico
!-------------------------------------------------------------------------------
module fftdipole_fftw

   use, intrinsic :: iso_c_binding
   use Parameters
   use Profiling
   use DipoleCommon, only : paddingShift, dipoleMatrix

   implicit none

#ifdef USE_FFTW
   include 'fftw3.f03'
#else
   type(C_PTR),external :: fftw_plan_dft_3d
   integer :: FFTW_FORWARD,FFTW_MEASURE, FFTW_BACKWARD
#endif

   integer :: Npadding
   integer(C_INT) :: N1_pad
   integer(C_INT) :: N2_pad
   integer(C_INT) :: N3_pad
   complex(C_DOUBLE_COMPLEX), dimension(:,:,:,:), allocatable :: bfield_ifft
   complex(C_DOUBLE_COMPLEX), dimension(:,:,:,:), allocatable :: emomM_trans
   complex(C_DOUBLE_COMPLEX), dimension(:,:,:,:), allocatable :: emomM_padded
   complex(C_DOUBLE_COMPLEX), dimension(:,:,:,:,:), allocatable :: DipMat_trans
   complex(C_DOUBLE_COMPLEX), dimension(:,:,:,:,:), allocatable :: DipMat_padded

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: geometrical_DFFT_wrapper
   !> @brief Wrapper for the calculation and initialization of the geometrical part
   !> of the FFT
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine geometrical_DFFT_wrapper(NA,N1,N2,N3,Natom,Mensemble,alat,C1,C2,C3,Bas)

      implicit none

      integer, intent(in) :: NA
      integer, intent(in) :: N1
      integer, intent(in) :: N2
      integer, intent(in) :: N3
      integer, intent(in) :: Natom
      integer, intent(in) :: Mensemble
      real(dblprec), intent(in) :: alat
      real(dblprec), dimension(3), intent(in) :: C1
      real(dblprec), dimension(3), intent(in) :: C2
      real(dblprec), dimension(3), intent(in) :: C3
      real(dblprec), dimension(3,NA), intent(in) :: Bas

      N1_pad=int(2*N1-1,C_INT)
      N2_pad=int(2*N2-1,C_INT)
      N3_pad=int(2*N3-1,C_INT)

      Npadding=N1_pad*N2_pad*N3_pad
      ! Allocate the arrays needed for the calculation of the dipole energy and
      ! field from the FFT convolution method
      call allocate_fft_dipole_structure(1,NA,Npadding,Mensemble)
      ! Calculation of the FFT of the geometrical dipole matrix
      call dipole_matrix_transform(N1,N2,N3,NA,alat,C1,C2,C3,Bas)

   end subroutine geometrical_DFFT_wrapper

   !----------------------------------------------------------------------------
   ! SUBROUTINE: calculate_fft_dipolar_field
   !> @brief Calculation of the dipolar field using the FTT convolution approach
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine calculate_fft_dipolar_field(NA,N1,N2,N3,Natom,Mensemble,emomM,bfield)

      implicit none

      integer, intent(in) :: NA
      integer, intent(in) :: N1
      integer, intent(in) :: N2
      integer, intent(in) :: N3
      integer, intent(in) :: Natom
      integer, intent(in) :: Mensemble
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: bfield

      integer :: kk,ii,jj,I3,I2,I1,I0,IA,nu
      real(dblprec) :: scale_fact,ene
      type(C_PTR) :: plan
      complex(C_DOUBLE_COMPLEX), dimension(3,Npadding) :: tmp_vec
      complex(C_DOUBLE_COMPLEX), dimension(N1_pad,N2_pad,N3_pad) :: tmp

      ! Do the FFTW of the magnetic moments
      call moment_transform(NA,N1,N2,N3,Natom,Npadding,Mensemble,emomM,emomM_padded,&
         emomM_trans)
      bfield_ifft=0.0_dblprec
      do kk=1, Mensemble
         do I0=1,NA
            do IA=1,NA
               ! Calculate the field in reciprocal space thanks to the convolution
               ! theorem
               tmp_vec=(0.0_dblprec,0.0_dblprec)
               !$omp parallel do default(shared), private(I3,I2,I1,ii)
               do I3=0,N3_pad-1
                  do I2=0,N2_pad-1
                     do I1=0,N1_pad-1
                        ii=1+I1+I2*N1_pad+I3*N2_pad*N3_pad
                        tmp_vec(:,ii)=matmul(DipMat_trans(:,:,I0,IA,ii),emomM_trans(:,IA,ii,kk))
                     enddo
                  enddo
               enddo
               !$omp end parallel do
               ! Loop over the components of the vector
               do nu=1,3
                  tmp  = (0.0_dblprec,0.0_dblprec)
                  ! Create the plan to perform an inverse FFT
                  !$OMP CRITICAL
                  plan = fftw_plan_dft_3d(N3_pad,N2_pad,N1_pad,tmp,tmp,FFTW_BACKWARD,FFTW_MEASURE)
                  !$OMP END CRITICAL
                  tmp  = (0.0_dblprec,0.0_dblprec)
                  ! Write the array to a 3D format that has the appropriate topology
                  !$omp parallel do default(shared), private(I3,I2,I1,ii)
                  do I3=0,N3_pad-1
                     do I2=0,N2_pad-1
                        do I1=0,N1_pad-1
                           ii=1+I1+I2*N1_pad+I3*N1_pad*N2_pad
                           tmp(I1+1,I2+1,I3+1)=tmp_vec(nu,ii)
                        enddo
                     enddo
                  enddo
                  !$omp end parallel do
                  ! Calcualte the inverse FFT
                  !$OMP CRITICAL
                  call fftw_execute_dft(plan,tmp,tmp)
                  !$OMP END CRITICAL
                  ! The output must be normalized by the size of the array
                  scale_fact=1.0_dblprec/dble(Npadding)
                  tmp=tmp*scale_fact
                  ! Store the field in a simple 1D array
                  !$omp parallel do default(shared), private(I3,I2,I1,ii)
                  do I3=0,N3_pad-1
                     do I2=0,N2_pad-1
                        do I1=0,N1_pad-1
                           ii=1+I1+I2*N1_pad+I3*N1_pad*N2_pad
                           bfield_ifft(nu,I0,ii,kk)=bfield_ifft(nu,I0,ii,kk)+tmp(I1+1,I2+1,I3+1)
                        enddo
                     enddo
                  enddo
                  !$omp end parallel do
                  ! Free the plan for the FFT
                  !$OMP CRITICAL
                  call fftw_destroy_plan(plan)
                  !$OMP END CRITICAL
               enddo
            enddo
         enddo
      enddo
      ! Store the data in a 1D unpadded real array
      !$omp parallel do default(shared), private(I0,I1,I2,I3,ii,jj)
      do I1=0,N1-1
         do I2=0,N2-1
            do I3=0,N3-1
               do I0=1,NA
                  ii=I0+I3*NA+I2*N3*NA+I1*N2*N3*NA
                  jj=1+I3+I2*N3_pad+I1*N2_pad*N3_pad
                  bfield(1:3,ii,1:Mensemble)=real(bfield_ifft(1:3,I0,jj,1:Mensemble))
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do

   end subroutine calculate_fft_dipolar_field

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_fft_dipole_structure
   !> @brief Allocation of the arrays for the dipole FFT convolution method
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine allocate_fft_dipole_structure(flag,NA,Npadding,Mensemble)

      implicit none

      integer, intent(in) :: flag
      integer, intent(in), optional :: NA
      integer, intent(in), optional :: Npadding
      integer, intent(in), optional :: Mensemble

      integer :: i_stat, i_all

      if (flag>0) then
         allocate(DipMat_padded(3,3,NA,NA,Npadding),stat=i_stat)
         call memocc(i_stat,product(shape(DipMat_padded))*kind(DipMat_padded),'DipMat_padded','allocate_fft_dipole_structure')
         DipMat_padded=(0.0_dblprec,0.0_dblprec)
         allocate(DipMat_trans(3,3,NA,NA,Npadding),stat=i_stat)
         call memocc(i_stat,product(shape(DipMat_trans))*kind(DipMat_trans),'DipMat_trans','allocate_fft_dipole_structure')
         DipMat_trans=(0.0_dblprec,0.0_dblprec)
         allocate(emomM_trans(3,NA,Npadding,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(emomM_trans))*kind(emomM_trans),'emomM_trans','allocate_fft_dipole_structure')
         emomM_trans=(0.0_dblprec,0.0_dblprec)
         allocate(emomM_padded(3,NA,Npadding,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(emomM_padded))*kind(emomM_padded),'emomM_padded','allocate_fft_dipole_structure')
         emomM_padded=(0.0_dblprec,0.0_dblprec)
         allocate(bfield_ifft(3,NA,Npadding,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(bfield_ifft))*kind(bfield_ifft),'bfield_ifft','allocate_fft_dipole_structure')
         bfield_ifft=(0.0_dblprec,0.0_dblprec)
      else
         i_all=-product(shape(emomM_trans))*kind(emomM_trans)
         deallocate(emomM_trans,stat=i_stat)
         call memocc(i_stat,i_all,'emomM_trans','allocate_fft_dipole_structure')
         i_all=-product(shape(emomM_padded))*kind(emomM_padded)
         deallocate(emomM_padded,stat=i_stat)
         call memocc(i_stat,i_all,'emomM_padded','allocate_fft_dipole_structure')
         i_all=-product(shape(DipMat_trans))*kind(DipMat_trans)
         deallocate(DipMat_trans,stat=i_stat)
         call memocc(i_stat,i_all,'DipMat_trans','allocate_fft_dipole_structure')
         i_all=-product(shape(bfield_ifft))*kind(bfield_ifft)
         deallocate(bfield_ifft,stat=i_stat)
         call memocc(i_stat,i_all,'bfield_ifft','allocate_fft_dipole_structure')
      endif

   end subroutine allocate_fft_dipole_structure

   !----------------------------------------------------------------------------
   ! SUBROUTINE: moment_padding
   !> @brief Creation of the padded moments array for the convolution
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine moment_padding(NA,N1,N2,N3,Natom,Npadding,Mensemble,emomM,emomM_padded)

      implicit none

      integer, intent(in) :: NA
      integer, intent(in) :: N1
      integer, intent(in) :: N2
      integer, intent(in) :: N3
      integer, intent(in) :: Natom
      integer, intent(in) :: Npadding
      integer, intent(in) :: Mensemble
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM
      complex(C_DOUBLE_COMPLEX), dimension(3,NA,Npadding,Mensemble), intent(out) :: emomM_padded

      integer :: I1,I2,I3,I0,ii,jj

      emomM_padded=(0.0_dblprec,0.0_dblprec)
      ! Store the magnetization in a padded array
      !$omp parallel do default(shared), private(I0,I1,I2,I3,ii,jj)
      do I0=1,NA
         do I3=0, N3-1
            do I2=0, N2-1
               do I1=0, N1-1
                  ii=I1*NA+N1*NA*I2+N1*N2*NA*I3+I0
                  jj=1+I1+I2*N1_pad+I3*N2_pad*N1_pad
                  emomM_padded(1:3,I0,jj,1:Mensemble)=emomM(1:3,ii,1:Mensemble)
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do

   end subroutine moment_padding

   !----------------------------------------------------------------------------
   ! SUBROUTINE: moment_transform
   !> @brief Fast Fourier Transform of the magnetic moments
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine moment_transform(NA,N1,N2,N3,Natom,Npadding,Mensemble,emomM,          &
      emomM_padded,emomM_trans)

      implicit none

      integer, intent(in) :: NA
      integer, intent(in) :: N1
      integer, intent(in) :: N2
      integer, intent(in) :: N3
      integer, intent(in) :: Natom
      integer, intent(in) :: Npadding
      integer, intent(in) :: Mensemble
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM
      complex(C_DOUBLE_COMPLEX), dimension(3,NA,Npadding,Mensemble), intent(inout) :: emomM_padded
      complex(C_DOUBLE_COMPLEX), dimension(3,NA,Npadding,Mensemble), intent(out) :: emomM_trans

      integer :: I0,kk,nu,ii,I1,I2,I3
      type(C_PTR) :: plan
      complex(C_DOUBLE_COMPLEX), dimension(N1_pad,N2_pad,N3_pad) :: tmp

      ! Padded the magnetic moment vector
      call moment_padding(NA,N1,N2,N3,Natom,Npadding,Mensemble,emomM,emomM_padded)

      ! Perform the FFT of the padded magnetic moment via the FFTW
      do kk=1,Mensemble
         do I0=1,NA
            do nu=1,3
               tmp  = (0.0_dblprec,0.0_dblprec)
               ! Create a plan for a 3D array that has the topology of the padded sample
               !$OMP CRITICAL
               plan = fftw_plan_dft_3d(N3_pad,N2_pad,N1_pad,tmp,tmp,FFTW_FORWARD,FFTW_MEASURE)
               !$OMP END CRITICAL
               tmp  = (0.0_dblprec,0.0_dblprec)
               ! Store the data in a 3D array that conserves the topology of the
               ! underlaying sample
               !$omp parallel do default(shared), private(I3,I2,I1,ii)
               do I3=0,N3_pad-1
                  do I2=0,N2_pad-1
                     do I1=0,N1_pad-1
                        ii=1+I1+I2*N1_pad+I3*N1_pad*N2_pad
                        tmp(I1+1,I2+1,I3+1)=emomM_padded(nu,I0,ii,kk)
                     enddo
                  enddo
               enddo
               !$omp end parallel do
               ! Calculate the FFT
               !$OMP CRITICAL
               call fftw_execute_dft(plan,tmp,tmp)
               !$OMP END CRITICAL
               ! Store the transformed data in a 1D array
               !$omp parallel do default(shared), private(I3,I2,I1,ii)
               do I3=0,N3_pad-1
                  do I2=0,N2_pad-1
                     do I1=0,N1_pad-1
                        ii=1+I1+I2*N1_pad+I3*N1_pad*N2_pad
                        emomM_trans(nu,I0,ii,kk)=tmp(I1+1,I2+1,I3+1)
                     enddo
                  enddo
               enddo
               !$omp end parallel do
               ! Free the information for the plan
               !$OMP CRITICAL
               call fftw_destroy_plan(plan)
               !$OMP END CRITICAL
            enddo
         enddo
      enddo
   end subroutine moment_transform

   !----------------------------------------------------------------------------
   ! SUBROUTINE: dipole_matrix_transform
   !> @brief Fast Fourier Transform of the Geometrical part of the dipole matrix
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine dipole_matrix_transform(N1,N2,N3,NA,alat,C1,C2,C3,Bas)

      implicit none

      integer, intent(in) :: N1
      integer, intent(in) :: N2
      integer, intent(in) :: N3
      integer, intent(in) :: NA
      real(dblprec), intent(in) :: alat
      real(dblprec), dimension(3), intent(in) :: C1
      real(dblprec), dimension(3), intent(in) :: C2
      real(dblprec), dimension(3), intent(in) :: C3
      real(dblprec), dimension(3,NA), intent(in) :: Bas

      type(C_PTR) :: plan
      integer :: i_stat,i_all
      integer :: I1,I2,I3,I0,IA,ii,mu,nu,jj
      integer :: I1_indx,I2_indx,I3_indx
      real(dblprec), dimension(3) :: Rij
      complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable :: tmp

      ! Allocate temporal arrays for the calculation of the FFT of the dipole matrix
      allocate(tmp(N1_pad,N2_pad,N2_pad),stat=i_stat)
      call memocc(i_stat,product(shape(tmp))*kind(tmp),'tmp','dipole_matrix_transform')
      tmp=(0.0_dblprec,0.0_dblprec)

      ! Calculate the padded dipole-dipole matrix
      Rij=0.0_dblprec
      do I0=1,NA
         do IA=1,NA
            ! Calculate the spatial part of the dipole-dipole matrix
            !$omp parallel do default(shared), private(I3,I2,I1,ii,I3_indx,I2_indx,I1_indx,mu,Rij)
            do I3=0,N3_pad-1
               do I2=0,N2_pad-1
                  do I1=0,N1_pad-1
                     ii=1+I1+I2*N1_pad+I3*N1_pad*N2_pad
                     I3_indx=paddingShift(N3,I3)
                     I2_indx=paddingShift(N2,I2)
                     I1_indx=paddingShift(N1,I1)
                     do mu=1,3
                        Rij(mu)=I1_indx*C1(mu)+I2_indx*C2(mu)+I3_indx*C3(mu)-Bas(mu,I0)+Bas(mu,IA)
                     enddo
                     DipMat_padded(:,:,I0,IA,ii)=dipoleMatrix(alat,Rij)
                  enddo
               enddo
            enddo
            !$omp end parallel do
            ! Calculate the FFT of the dipole matrix using the MKL FFT library
            do mu=1,3
               do nu=1,3
                  tmp  = (0.0_dblprec,0.0_dblprec)
                  !$OMP CRITICAL
                  ! Create a plan for a 3D array that has the topology of the padded sample
                  plan = fftw_plan_dft_3d(N3_pad,N2_pad,N1_pad,tmp,tmp,FFTW_FORWARD,FFTW_MEASURE)
                  !$OMP END CRITICAL
                  tmp  = (0.0_dblprec,0.0_dblprec)
                  ! Store the data in a 3D array that conserves the topology of the
                  ! underlaying sample
                  !$omp parallel do default(shared), private(I3,I2,I1,ii)
                  do I3=0,N3_pad-1
                     do I2=0,N2_pad-1
                        do I1=0,N1_pad-1
                           ii=1+I1+I2*N1_pad+I3*N1_pad*N2_pad
                           tmp(I1+1,I2+1,I3+1)=DipMat_padded(mu,nu,I0,IA,ii)
                        enddo
                     enddo
                  enddo
                  !$omp end parallel do
                  ! Calculate the FFT of the dipole-dipole matrix
                  !$OMP CRITICAL
                  call fftw_execute_dft(plan,tmp,tmp)
                  !$OMP END CRITICAL
                  ! Store the information in a simple 1D array
                  !$omp parallel do default(shared), private(I3,I2,I1,ii)
                  do I3=0,N3_pad-1
                     do I2=0,N2_pad-1
                        do I1=0,N1_pad-1
                           ii=1+I1+I2*N1_pad+I3*N1_pad*N2_pad
                           DipMat_trans(mu,nu,I0,IA,ii)=tmp(I1+1,I2+1,I3+1)
                        enddo
                     enddo
                  enddo
                  !$omp end parallel do
                  ! Free the information of the plan
                  !$OMP CRITICAL
                  call fftw_destroy_plan(plan)
                  !$OMP END CRITICAL
               enddo
            enddo
         enddo
      enddo

      ! Deallocate the temporal arrays for the calculation of the FFT of the dipole matrix
      i_all=-product(shape(DipMat_padded))*kind(DipMat_padded)
      deallocate(DipMat_padded,stat=i_stat)
      call memocc(i_stat,i_all,'DipMat_padded','dipole_matrix_transform')
      i_all=-product(shape(tmp))*kind(tmp)
      deallocate(tmp,stat=i_stat)
      call memocc(i_stat,i_all,'tmp','dipole_matrix_transform')

   end subroutine dipole_matrix_transform

end module fftdipole_fftw
