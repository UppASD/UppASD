!-------------------------------------------------------------------------------
! MODULE: fftdipole_mkl
!> @brief Routine for the calculation of the dipole-dipole interaction via the
!> FFT convolution method.
!> @details This routine calculates the dipole-dipole interacction making use of
!> the convolution theorem and the discrete FFT. The implementation is based in the
!> work by Moritz Sallerman https://github.com/MSallermann/fft-dipoles-proof-of-concept.
!> This module makes use of the intel MKL library.
!> Notice that the implementation makes use of the *strides* capabilities, in which
!> one can pass the topology of a 3D array to a 1D array which is our data container
!> @author Jonathan Chico
!-------------------------------------------------------------------------------
#ifdef USE_MKL_FFT
   include "mkl_dfti.f90"
#endif

module fftdipole_mkl

   use Parameters
   use Profiling
   use DipoleCommon, only : paddingShift, dipoleMatrix

#ifdef USE_MKL_FFT
   use mkl_dfti
#else
   integer, external :: DftiCreateDescriptor,DftiSetValue,DftiCommitDescriptor
   integer, external :: DftiComputeForward,DftiFreeDescriptor,DftiComputeBackward
   integer :: DFTI_COMPLEX,DFTI_DOUBLE,DFTI_FORWARD_SCALE,DFTI_INPLACE
   integer :: DFTI_PLACEMENT,DFTI_BACKWARD_SCALE,DFTI_INPUT_STRIDES,DFTI_OUTPUT_STRIDES
   type DFTI_DESCRIPTOR
   end type DFTI_DESCRIPTOR
#endif

   integer :: Npadding
   integer :: N1_pad
   integer :: N2_pad
   integer :: N3_pad
   complex(dblprec), dimension(:,:,:,:), allocatable :: bfield_ifft
   complex(dblprec), dimension(:,:,:,:), allocatable :: emomM_trans
   complex(dblprec), dimension(:,:,:,:), allocatable :: emomM_padded
   complex(dblprec), dimension(:,:,:,:,:), allocatable :: DipMat_trans
   complex(dblprec), dimension(:,:,:,:,:), allocatable :: DipMat_padded

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: geometrical_DFFT_wrapper
   !> @brief Wrapper for the calculation and initialization of the geometrical part
   !> of the FFT
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine geometrical_DFFT_wrapper_MKL(NA,N1,N2,N3,Natom,Mensemble,alat,C1,C2,C3,Bas)

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

      N1_pad=int(2*N1-1)
      N2_pad=int(2*N2-1)
      N3_pad=int(2*N3-1)
      Npadding=N1_pad*N2_pad*N3_pad
      ! Allocate the arrays needed for the calculation of the dipole energy and
      ! field from the FFT convolution method
      call allocate_fft_dipole_structure_MKL(1,NA,Npadding,Mensemble)
      ! Calculation of the FFT of the geometrical dipole matrix
      call dipole_matrix_transform_MKL(N1,N2,N3,NA,alat,C1,C2,C3,Bas)

   end subroutine geometrical_DFFT_wrapper_MKL

   !----------------------------------------------------------------------------
   ! SUBROUTINE: calculate_fft_dipolar_field
   !> @brief Calculation of the dipolar field using the FTT convolution approach
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine calculate_fft_dipolar_field_MKL(NA,N1,N2,N3,Natom,Mensemble,emomM,bfield)

      implicit none

      integer, intent(in) :: NA
      integer, intent(in) :: N1
      integer, intent(in) :: N2
      integer, intent(in) :: N3
      integer, intent(in) :: Natom
      integer, intent(in) :: Mensemble
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: bfield

      integer :: status
      integer :: kk,ii,jj,I3,I2,I1,I0,IA,nu
      type(DFTI_DESCRIPTOR), pointer :: handle
      integer, dimension(3) :: Ndim
      integer, dimension(4) :: strides
      complex(dblprec), dimension(Npadding) :: tmp
      complex(dblprec), dimension(3,Npadding) :: tmp_vec

      Ndim(1)=N1_pad
      Ndim(2)=N2_pad
      Ndim(3)=N3_pad
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! The strides contain the information of the topology of the 3D array and
      ! passes it to the 1D array. This makes sure that the FFT considers the
      ! periodicity/padding properly
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      strides(1)=0
      strides(2)=1
      strides(3)=Ndim(1)
      strides(4)=Ndim(1)*Ndim(2)

      ! Do the FFTW of the magnetic moments
      call moment_transform_MKL(NA,N1,N2,N3,Natom,Npadding,Mensemble,emomM,emomM_padded,&
         emomM_trans)

      bfield_ifft=0.0_dblprec
      do kk=1, Mensemble
         do I0=1,NA
            do IA=1,NA
               ! Calculate the field in reciprocal space thanks to the convolution
               ! theorem
               tmp_vec=(0.0_dblprec,0.0_dblprec)
               !$omp parallel do default(shared), private(I1,I2,I3,ii)
               do I3=0,Ndim(3)-1
                  do I2=0,Ndim(2)-1
                     do I1=0,Ndim(1)-1
                        ii=1+I1+I2*Ndim(1)+I3*Ndim(1)*Ndim(2)
                        tmp_vec(:,ii)=matmul(DipMat_trans(:,:,I0,IA,ii),emomM_trans(:,IA,ii,kk))
                     enddo
                  enddo
               enddo
               !$omp end parallel do
               do nu=1,3
                  ! Initialize the arrays
                  tmp = (0.0_dblprec,0.0_dblprec)
                  ! Perform the calculation for each component of the magnetization
                  tmp = tmp_vec(nu,:)
                  ! Create the pointer for the MKL calls
                  status  = DftiCreateDescriptor(handle,DFTI_DOUBLE,DFTI_COMPLEX,3,Ndim)
                  ! Set parameters to ensure that the in/out arrays are the same
                  status  = DftiSetValue(handle,DFTI_PLACEMENT,DFTI_INPLACE)
                  ! Make sure that the normalization is properly considered
                  status  = DftiSetValue(handle,DFTI_BACKWARD_SCALE,1.0_dblprec/dble(Npadding))
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ! Set the parameters for the strides, i.e. give the information
                  ! of the underlaying 3D topology which corresponds to the 1D array
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  status  = DftiSetValue(handle,DFTI_INPUT_STRIDES,strides)
                  ! Commit the information of the pointer for the FFT
                  status  = DftiCommitDescriptor(handle)
                  ! Actually perform the inverse Fourier transform
                  status  = DftiComputeBackward(handle,tmp)
                  ! Free the memory used
                  status  = DftiFreeDescriptor(handle)
                  ! Store the information in the output array
                  bfield_ifft(nu,I0,1:Npadding,kk)=bfield_ifft(nu,I0,1:Npadding,kk)+tmp
               enddo
            enddo
         enddo
      enddo
      !$omp parallel do default(shared), private(I0,I1,I2,I3,ii,jj)
      do I3=0,N3-1
         do I2=0,N2-1
            do I1=0,N1-1
               do I0=1,NA
                  ii=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA
                  jj=1+I1+I2*Ndim(1)+I3*Ndim(1)*Ndim(2)
                  bfield(1:3,ii,1:Mensemble)=real(bfield_ifft(1:3,I0,jj,1:Mensemble))
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do

   end subroutine calculate_fft_dipolar_field_MKL

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_fft_dipole_structure_MKL
   !> @brief Allocation of the arrays for the dipole FFT convolution method
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine allocate_fft_dipole_structure_MKL(flag,NA,Npadding,Mensemble)

      implicit none

      integer, intent(in) :: flag
      integer, intent(in), optional :: NA
      integer, intent(in), optional :: Npadding
      integer, intent(in), optional :: Mensemble

      integer :: i_stat, i_all

      if (flag>0) then
         allocate(DipMat_padded(3,3,NA,NA,Npadding),stat=i_stat)
         call memocc(i_stat,product(shape(DipMat_padded))*kind(DipMat_padded),'DipMat_padded','allocate_fft_dipole_structure_MKL')
         DipMat_padded=(0.0_dblprec,0.0_dblprec)
         allocate(DipMat_trans(3,3,NA,NA,Npadding),stat=i_stat)
         call memocc(i_stat,product(shape(DipMat_trans))*kind(DipMat_trans),'DipMat_trans','allocate_fft_dipole_structure_MKL')
         DipMat_trans=(0.0_dblprec,0.0_dblprec)
         allocate(emomM_trans(3,NA,Npadding,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(emomM_trans))*kind(emomM_trans),'emomM_trans','allocate_fft_dipole_structure_MKL')
         emomM_trans=(0.0_dblprec,0.0_dblprec)
         allocate(emomM_padded(3,NA,Npadding,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(emomM_padded))*kind(emomM_padded),'emomM_padded','allocate_fft_dipole_structure_MKL')
         emomM_padded=(0.0_dblprec,0.0_dblprec)
         allocate(bfield_ifft(3,NA,Npadding,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(bfield_ifft))*kind(bfield_ifft),'bfield_ifft','allocate_fft_dipole_structure_MKL')
         bfield_ifft=(0.0_dblprec,0.0_dblprec)
      else
         i_all=-product(shape(emomM_trans))*kind(emomM_trans)
         deallocate(emomM_trans,stat=i_stat)
         call memocc(i_stat,i_all,'emomM_trans','allocate_fft_dipole_structure_MKL')
         i_all=-product(shape(emomM_padded))*kind(emomM_padded)
         deallocate(emomM_padded,stat=i_stat)
         call memocc(i_stat,i_all,'emomM_padded','allocate_fft_dipole_structure_MKL')
         i_all=-product(shape(DipMat_trans))*kind(DipMat_trans)
         deallocate(DipMat_trans,stat=i_stat)
         call memocc(i_stat,i_all,'DipMat_trans','allocate_fft_dipole_structure_MKL')
         i_all=-product(shape(bfield_ifft))*kind(bfield_ifft)
         deallocate(bfield_ifft,stat=i_stat)
         call memocc(i_stat,i_all,'bfield_ifft','allocate_fft_dipole_structure_MKL')
      endif

   end subroutine allocate_fft_dipole_structure_MKL

   !----------------------------------------------------------------------------
   ! SUBROUTINE: moment_padding
   !> @brief Creation of the padded moments array for the convolution
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine moment_padding_MKL(NA,N1,N2,N3,Natom,Npadding,Mensemble,emomM,emomM_padded)

      implicit none

      integer, intent(in) :: NA
      integer, intent(in) :: N1
      integer, intent(in) :: N2
      integer, intent(in) :: N3
      integer, intent(in) :: Natom
      integer, intent(in) :: Npadding
      integer, intent(in) :: Mensemble
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM
      complex(dblprec), dimension(3,NA,Npadding,Mensemble), intent(out) :: emomM_padded

      integer :: I1,I2,I3,I0,ii,jj

      emomM_padded=(0.0_dblprec,0.0_dblprec)
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

   end subroutine moment_padding_MKL

   !----------------------------------------------------------------------------
   ! SUBROUTINE: moment_transform
   !> @brief Fast Fourier Transform of the magnetic moments
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine moment_transform_MKL(NA,N1,N2,N3,Natom,Npadding,Mensemble,emomM,          &
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
      complex(dblprec), dimension(3,NA,Npadding,Mensemble), intent(inout) :: emomM_padded
      complex(dblprec), dimension(3,NA,Npadding,Mensemble), intent(out) :: emomM_trans

      integer :: status
      integer :: I0,kk,nu,I1,I2,I3,ii
      integer, dimension(3) :: Ndim
      integer, dimension(4) :: strides
      complex(dblprec), dimension(Npadding) :: tmp
      type(DFTI_DESCRIPTOR), pointer :: handle

      Ndim(1)=N1_pad
      Ndim(2)=N2_pad
      Ndim(3)=N3_pad
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! The strides contain the information of the topology of the 3D array and
      ! passes it to the 1D array. This makes sure that the FFT considers the
      ! periodicity/padding properly
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      strides(1)=0
      strides(2)=1
      strides(3)=Ndim(1)
      strides(4)=Ndim(1)*Ndim(2)

      ! Padded the magnetic moment vector
      call moment_padding_MKL(NA,N1,N2,N3,Natom,Npadding,Mensemble,emomM,emomM_padded)

      ! Perform the FFT of the padded magnetic moment via the MKL FFT
      do kk=1,Mensemble
         do I0=1,NA
            do nu=1,3
               tmp = (0.0_dblprec,0.0_dblprec)
               ! Perform the calculation for each component of the magnetization
               tmp = emomM_padded(nu,I0,1:Npadding,kk)
               ! Create the pointer for the MKL calls
               status  = DftiCreateDescriptor(handle,DFTI_DOUBLE,DFTI_COMPLEX,3,Ndim)
               ! Set parameters to ensure that the in/out arrays are the same
               status  = DftiSetValue(handle,DFTI_PLACEMENT,DFTI_INPLACE)
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Set the parameters for the strides, i.e. give the information
               ! of the underlaying 3D topology which corresponds to the 1D array
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               status  = DftiSetValue(handle,DFTI_INPUT_STRIDES,strides)
               ! Commit the information of the pointer for the FFT
               status  = DftiCommitDescriptor(handle)
               ! Actually perform the Fourier transform
               status  = DftiComputeForward(handle,tmp)
               ! Free the memory used
               status  = DftiFreeDescriptor(handle)
               ! Store the information in the output array
               emomM_trans(nu,I0,1:Npadding,kk)=tmp
            enddo
         enddo
      enddo
   end subroutine moment_transform_MKL

   !----------------------------------------------------------------------------
   ! SUBROUTINE: dipole_matrix_transform
   !> @brief Fast Fourier Transform of the Geometrical part of the dipole matrix
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine dipole_matrix_transform_MKL(N1,N2,N3,NA,alat,C1,C2,C3,Bas)

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

      type(DFTI_DESCRIPTOR), pointer :: handle
      integer :: status
      integer :: i_stat,i_all
      integer :: I1,I2,I3,I0,IA,ii,nu,mu
      integer :: I1_indx,I2_indx,I3_indx
      integer, dimension(3) :: Ndim
      integer, dimension(4):: strides
      real(dblprec), dimension(3) :: Rij
      complex(dblprec), dimension(:), allocatable :: tmp

      Ndim(1)=N1_pad
      Ndim(2)=N2_pad
      Ndim(3)=N3_pad
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! The strides contain the information of the topology of the 3D array and
      ! passes it to the 1D array. This makes sure that the FFT considers the
      ! periodicity/padding properly
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      strides(1)=0
      strides(2)=1
      strides(3)=Ndim(1)
      strides(4)=Ndim(1)*Ndim(2)

      allocate(tmp(Npadding),stat=i_stat)
      call memocc(i_stat,product(shape(tmp))*kind(tmp),'tmp','dipole_matrix_transform_MKL')
      tmp=(0.0_dblprec,0.0_dblprec)
      ! Calculate the padded dipole-dipole matrix
      Rij=0.0_dblprec
      do I0=1,NA
         do IA=1,NA
            ! Calculate the spatial part of the dipole-dipole matrix
            do I3=0,Ndim(3)-1
               do I2=0,Ndim(2)-1
                  do I1=0,Ndim(1)-1
                     ii=1+I1+I2*Ndim(1)+I3*Ndim(1)*Ndim(2)
                     I3_indx=paddingShift(N3,I3)
                     I2_indx=paddingShift(N2,I2)
                     I1_indx=paddingShift(N1,I1)
                     do nu=1,3
                        Rij(nu)=I1_indx*C1(nu)+I2_indx*C2(nu)+I3_indx*C3(nu)-Bas(nu,I0)+Bas(nu,IA)
                     enddo
                     DipMat_padded(:,:,I0,IA,ii)=dipoleMatrix(alat,Rij)
                  enddo
               enddo
            enddo
            ! Calculate the FFT of the dipole matrix using the MKL FFT library
            do mu=1,3
               do nu=1,3
                  tmp = (0.0_dblprec,0.0_dblprec)
                  ! Perform the calculation for each component of the magnetization
                  tmp = DipMat_padded(mu,nu,I0,IA,1:Npadding)
                  ! Create the pointer for the MKL calls
                  status  = DftiCreateDescriptor(handle,DFTI_DOUBLE,DFTI_COMPLEX,3,Ndim)
                  ! Set parameters to ensure that the in/out arrays are the same
                  status  = DftiSetValue(handle,DFTI_PLACEMENT,DFTI_INPLACE)
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ! Set the parameters for the strides, i.e. give the information
                  ! of the underlaying 3D topology which corresponds to the 1D array
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  status  = DftiSetValue(handle,DFTI_INPUT_STRIDES,strides)
                  ! Commit the information of the pointer for the FFT
                  status  = DftiCommitDescriptor(handle)
                  ! Actually perform the Fourier transform
                  status  = DftiComputeForward(handle,tmp)
                  ! Free the memory used
                  status  = DftiFreeDescriptor(handle)
                  ! Store the information in the output array
                  DipMat_trans(mu,nu,I0,IA,1:Npadding)=tmp
               enddo
            enddo
         enddo
      enddo

      ! Deallocate the temporal arrays for the calculation of the FFT of the dipole matrix
      i_all=-product(shape(DipMat_padded))*kind(DipMat_padded)
      deallocate(DipMat_padded,stat=i_stat)
      call memocc(i_stat,i_all,'DipMat_padded','dipole_matrix_transform_MKL')
      i_all=-product(shape(tmp))*kind(tmp)
      deallocate(tmp,stat=i_stat)
      call memocc(i_stat,i_all,'tmp','dipole_matrix_transform_MKL')

   end subroutine dipole_matrix_transform_MKL

end module fftdipole_mkl
