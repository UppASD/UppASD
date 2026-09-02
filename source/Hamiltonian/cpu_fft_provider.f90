!-------------------------------------------------------------------------------
! MODULE: CPUFFTProvider
!> @brief Minimal provider boundary for CPU real-to-complex FFTs.
!>
!> The convolution reference owns the data layout and spectral algebra.  This
!> module owns the provider-specific plan handles and calls, keeping FFTW (and
!> a future MKL provider) out of the convolution and Hamiltonian modules.
!-------------------------------------------------------------------------------
module CPUFFTProvider

   ! fftw3.f03 refers to the complete ISO_C_BINDING namespace, not only the
   ! symbols used by the small provider entry points below.
   use, intrinsic :: iso_c_binding

   implicit none
   private

   include 'fftw3.f03'

   type, public :: cpu_fft_plan_t
      type(C_PTR) :: handle=C_NULL_PTR
   end type cpu_fft_plan_t

   public :: cpu_fft_plan_many_r2c
   public :: cpu_fft_plan_many_c2r
   public :: cpu_fft_execute_r2c
   public :: cpu_fft_execute_c2r
   public :: cpu_fft_plan_destroy

contains

   logical function cpu_fft_plan_many_r2c(plan,n1,n2,n3,howmany,real_data,complex_data)
      type(cpu_fft_plan_t), intent(inout) :: plan
      integer, intent(in) :: n1,n2,n3,howmany
      real(C_DOUBLE), intent(inout), target :: real_data(:)
      complex(C_DOUBLE_COMPLEX), intent(inout), target :: complex_data(:)
      integer(C_INT) :: n(3), inembed(3), onembed(3)

      call cpu_fft_plan_destroy(plan)
      cpu_fft_plan_many_r2c=.false.
      if (n1 < 1 .or. n2 < 1 .or. n3 < 1 .or. howmany < 1) return

      ! FFTW's C dimensions are reversed so the Fortran first dimension is
      ! the contiguous x/cell dimension.  This matches the UppASD cell map:
      ! cell = x + n1*(y + n2*z).
      n=(/int(n3,C_INT),int(n2,C_INT),int(n1,C_INT)/)
      inembed=n
      onembed=(/int(n3,C_INT),int(n2,C_INT),int(n1/2+1,C_INT)/)
      plan%handle=fftw_plan_many_dft_r2c(3,n,int(howmany,C_INT),real_data, &
         inembed,1,int(n1*n2*n3,C_INT),complex_data,onembed,1, &
         int((n1/2+1)*n2*n3,C_INT),FFTW_ESTIMATE)
      cpu_fft_plan_many_r2c=c_associated(plan%handle)
   end function cpu_fft_plan_many_r2c


   logical function cpu_fft_plan_many_c2r(plan,n1,n2,n3,howmany,complex_data,real_data)
      type(cpu_fft_plan_t), intent(inout) :: plan
      integer, intent(in) :: n1,n2,n3,howmany
      complex(C_DOUBLE_COMPLEX), intent(inout), target :: complex_data(:)
      real(C_DOUBLE), intent(inout), target :: real_data(:)
      integer(C_INT) :: n(3), inembed(3), onembed(3)

      call cpu_fft_plan_destroy(plan)
      cpu_fft_plan_many_c2r=.false.
      if (n1 < 1 .or. n2 < 1 .or. n3 < 1 .or. howmany < 1) return

      n=(/int(n3,C_INT),int(n2,C_INT),int(n1,C_INT)/)
      inembed=(/int(n3,C_INT),int(n2,C_INT),int(n1/2+1,C_INT)/)
      onembed=n
      plan%handle=fftw_plan_many_dft_c2r(3,n,int(howmany,C_INT),complex_data, &
         inembed,1,int((n1/2+1)*n2*n3,C_INT),real_data,onembed,1, &
         int(n1*n2*n3,C_INT),FFTW_ESTIMATE)
      cpu_fft_plan_many_c2r=c_associated(plan%handle)
   end function cpu_fft_plan_many_c2r


   subroutine cpu_fft_execute_r2c(plan,real_data,complex_data)
      type(cpu_fft_plan_t), intent(in) :: plan
      real(C_DOUBLE), intent(inout), target :: real_data(:)
      complex(C_DOUBLE_COMPLEX), intent(inout), target :: complex_data(:)

      call fftw_execute_dft_r2c(plan%handle,real_data,complex_data)
   end subroutine cpu_fft_execute_r2c


   subroutine cpu_fft_execute_c2r(plan,complex_data,real_data)
      type(cpu_fft_plan_t), intent(in) :: plan
      complex(C_DOUBLE_COMPLEX), intent(inout), target :: complex_data(:)
      real(C_DOUBLE), intent(inout), target :: real_data(:)

      call fftw_execute_dft_c2r(plan%handle,complex_data,real_data)
   end subroutine cpu_fft_execute_c2r


   subroutine cpu_fft_plan_destroy(plan)
      type(cpu_fft_plan_t), intent(inout) :: plan

      if (c_associated(plan%handle)) call fftw_destroy_plan(plan%handle)
      plan%handle=C_NULL_PTR
   end subroutine cpu_fft_plan_destroy

end module CPUFFTProvider
