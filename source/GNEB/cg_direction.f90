module cg_direction
  !!----------------------------------------------------------------------------
  !! PURPOSE:
  !!   Provide a robust Conjugate-Gradient search direction update suitable for
  !!   OSO (orthogonal spin optimization). Implements Fletcher–Reeves (FR) and
  !!   Polak–Ribiére+ (PRP+) with numerical safeguards and automatic restarts.
  !!
  !! USAGE:
  !!   call get_cg_direction(g_curr, g_prev, p_prev, p_curr, method, iter, restart_period)
  !!
  !! INPUTS:
  !!   g_curr(:)         - current gradient vector at x_k  (length = ndof)
  !!   g_prev(:)         - previous gradient vector at x_{k-1} (same length)
  !!   p_prev(:)         - previous search direction
  !!   method            - character(*)  'FR' or 'PRP'  (PRP means PRP+ truncation)
  !!   iter              - integer iteration index (k). Use iter=0 for first step.
  !!   restart_period    - positive integer; e.g.  (max(10, ndof)) or simply 20–50
  !!
  !! OUTPUT:
  !!   p_curr(:)         - new search direction
  !!
  !! NOTES:
  !!   - We apply standard rules:
  !!       p_k = -g_k + beta_k * p_{k-1}
  !!       beta_FR  = (g_k·g_k) / (g_{k-1}·g_{k-1})
  !!       beta_PRP = max( (g_k·(g_k - g_{k-1})) / (g_{k-1}·g_{k-1}), 0 )
  !!   - Safeguards:
  !!       * restart if p_k is not a descent direction (p_k·g_k >= 0)
  !!       * restart every 'restart_period' iterations (numerical hygiene)
  !!       * restart if ||g_{k-1}|| is tiny (avoid division blow-ups)
  !!
  !!----------------------------------------------------------------------------

  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: get_cg_direction

contains

  subroutine get_cg_direction(g_curr, g_prev, p_prev, p_curr, method_in, iter, restart_period)
    real(dp), intent(in)  :: g_curr(:)
    real(dp), intent(in)  :: g_prev(:)
    real(dp), intent(in)  :: p_prev(:)
    real(dp), intent(out) :: p_curr(:)
  character(*), intent(in) :: method_in   ! 'FR' or 'PRP'
    integer, intent(in) :: iter, restart_period

    real(dp) :: gg_curr, gg_prev, beta, pg
    real(dp), parameter :: tiny = 1.0d-300
    logical :: do_restart

    if (size(g_curr) /= size(g_prev) .or. size(g_curr) /= size(p_prev) .or. &
        size(g_curr) /= size(p_curr)) then
       error stop 'get_cg_direction: size mismatch.'
    end if

    ! First iteration or forced periodic restart → steepest descent.
    do_restart = (iter == 0)
    if (restart_period > 0) do_restart = do_restart .or. (mod(iter, restart_period) == 0)

    if (do_restart) then
       p_curr = -g_curr
       return
    end if

    ! Compute norms and beta.
    gg_curr = dot_product(g_curr, g_curr)
    gg_prev = max(dot_product(g_prev, g_prev), tiny)

  select case (adjustl(upper_trim(method_in)))
    case ('FR')
       beta = gg_curr / gg_prev
    case ('PRP')   ! Polak-Ribiere with truncation (a.k.a. PRP+)
       beta = dot_product(g_curr, g_curr - g_prev) / gg_prev
       if (beta < 0.0_dp) beta = 0.0_dp
    case default
       error stop 'get_cg_direction: method must be ''FR'' or ''PRP''.'
    end select

    p_curr = -g_curr + beta * p_prev

    ! Ensure descent direction; if not, restart.
    pg = dot_product(p_curr, g_curr)
    if (pg >= 0.0_dp .or. .not.(is_finite_vec(p_curr))) then
       p_curr = -g_curr
    end if
  end subroutine get_cg_direction


  pure logical function is_finite_vec(x) result(ok)
    real(dp), intent(in) :: x(:)
    integer :: i
    real(dp) :: xi
    ok = .true.
    do i = 1, size(x)
      xi = x(i)
      ! NaN check: NaN is not equal to itself. Infinity check: magnitude greater than huge()
      if (xi /= xi .or. abs(xi) > huge(xi)) then
        ok = .false.; return
      end if
    end do
  end function is_finite_vec


  pure function upper_trim(s) result(t)
    character(*), intent(in) :: s
    character(len=:), allocatable :: t
    integer :: n,i
    character(len=:), allocatable :: ss
    character(len=1) :: ch
    integer, parameter :: diff = iachar('a') - iachar('A')

    ! create a left-adjusted temporary and determine the trimmed length
    n = len_trim(adjustl(s))
    if (n <= 0) then
      allocate(character(len=0) :: t)
      return
    end if

    allocate(character(len=len(s)) :: ss)
    ss = adjustl(s)
    allocate(character(len=n) :: t)
    do i = 1, n
      ch = ss(i:i)
      if (iachar(ch) >= iachar('a') .and. iachar(ch) <= iachar('z')) then
        t(i:i) = achar(iachar(ch) - diff)
      else
        t(i:i) = ch
      end if
    end do
  end function upper_trim

end module cg_direction

! !------------------------------------------------------------------------------
! ! Conjugate-gradient direction (FR / PRP+ with safe restarts)
! !------------------------------------------------------------------------------
! module cg_direction_mod
!   use oso_kinds_mod, only: dp
!   use, intrinsic :: ieee_arithmetic
!   implicit none
!   private
!   public :: get_cg_direction
! 
! contains
! 
!   subroutine get_cg_direction(g_curr, g_prev, p_prev, p_curr, method, iter, restart_period)
!     real(dp), intent(in)  :: g_curr(:), g_prev(:), p_prev(:)
!     real(dp), intent(out) :: p_curr(:)
!     character(*), intent(in) :: method   ! 'FR' or 'PRP' (case-insensitive)
!     integer, intent(in) :: iter, restart_period
! 
!     real(dp) :: gg_curr, gg_prev, beta, pg
!     real(dp), parameter :: tiny = 1.0d-300
!     logical :: do_restart
!     character(len=:), allocatable :: m
! 
!     if (any([size(g_curr),size(g_prev),size(p_prev),size(p_curr)] /= size(g_curr))) &
!       error stop 'get_cg_direction: size mismatch.'
! 
!     m = to_upper(trim(method))
!     do_restart = (iter == 0)
!     if (restart_period > 0) do_restart = do_restart .or. (mod(iter,restart_period) == 0)
! 
!     if (do_restart) then
!       p_curr = -g_curr
!       return
!     end if
! 
!     gg_curr = dot_product(g_curr,g_curr)
!     gg_prev = max(dot_product(g_prev,g_prev), tiny)
! 
!     select case (m)
!     case ('FR')
!       beta = gg_curr / gg_prev
!     case ('PRP')
!       beta = dot_product(g_curr, g_curr - g_prev) / gg_prev
!       if (beta < 0.0_dp) beta = 0.0_dp
!     case default
!       error stop 'get_cg_direction: method must be FR or PRP.'
!     end select
! 
!     p_curr = -g_curr + beta * p_prev
!     pg = dot_product(p_curr, g_curr)
! 
!     if (pg >= 0.0_dp .or. .not.all_finite(p_curr)) p_curr = -g_curr
!   end subroutine get_cg_direction
! 
! end module cg_direction_mod