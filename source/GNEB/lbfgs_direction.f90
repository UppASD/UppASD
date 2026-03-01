module lbfgs_direction_mod
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: lbfgs_state, lbfgs_init, lbfgs_reset, lbfgs_push_pair, lbfgs_direction

  type :: lbfgs_state
    integer :: ndof = 0
    integer :: m    = 10        ! memory (you can change at init)
    integer :: count = 0        ! how many pairs currently stored (<= m)
    integer :: head  = 0        ! index of the most recent pair (1..m), 0 means empty
    real(dp) :: H0 = 1.0_dp     ! initial Hessian scale (γ)
    real(dp), allocatable :: S(:,:)  ! (ndof, m) columns are s_k
    real(dp), allocatable :: Y(:,:)  ! (ndof, m) columns are y_k
    real(dp), allocatable :: rho(:)  ! (m)
    real(dp), allocatable :: alpha(:)! (m) scratch for two-loop
  end type lbfgs_state

contains

  subroutine lbfgs_init(st, ndof, m_mem)
    type(lbfgs_state), intent(inout) :: st
    integer, intent(in) :: ndof, m_mem
    st%ndof = ndof
    st%m    = max(1, m_mem)
    st%count = 0
    st%head  = 0
    st%H0    = 1.0_dp
    allocate(st%S(ndof, st%m), st%Y(ndof, st%m), st%rho(st%m), st%alpha(st%m))
    st%S = 0.0_dp; st%Y = 0.0_dp; st%rho = 0.0_dp; st%alpha = 0.0_dp
  end subroutine lbfgs_init

  subroutine lbfgs_reset(st)
    type(lbfgs_state), intent(inout) :: st
    st%count = 0
    st%head  = 0
    st%H0    = 1.0_dp
  end subroutine lbfgs_reset

  pure integer function prev_idx(st, k) result(i)
    type(lbfgs_state), intent(in) :: st
    integer, intent(in) :: k
    ! k = 0 => most recent, k = st%count-1 => oldest
    integer :: t
    if (st%count == 0) then
      i = 0; return
    end if
    t = st%head - k
    do while (t <= 0) ; t = t + st%m ; end do
    i = t
  end function prev_idx

  subroutine lbfgs_push_pair(st, s, y)
    type(lbfgs_state), intent(inout) :: st
    real(dp), intent(in) :: s(:), y(:)
    real(dp) :: ys, yy, ns, ny
    integer :: loc

    ! Curvature check: y^T s must be sufficiently positive
    ys = dot_product(y, s)
    ns = sqrt(dot_product(s, s)); ny = sqrt(dot_product(y, y))
    if (ys <= 1.0d-12 * ns * ny) then
      return   ! skip this pair (fallback handled by history already present)
    end if

    ! Move head forward in a circular buffer
    if (st%head == 0) then
      st%head = 1
    else
      st%head = st%head + 1
      if (st%head > st%m) st%head = 1
    end if
    loc = st%head

    st%S(:,loc) = s
    st%Y(:,loc) = y
    st%rho(loc) = 1.0_dp / ys

    if (st%count < st%m) st%count = st%count + 1

    ! Scale H0 with γ = (s^T y)/(y^T y)
    yy = dot_product(y, y)
    if (yy > 0.0_dp) st%H0 = ys / yy
  end subroutine lbfgs_push_pair

  subroutine lbfgs_direction(st, g, p)
    type(lbfgs_state), intent(inout) :: st
    real(dp), intent(in)  :: g(:)
    real(dp), intent(out) :: p(:)
    real(dp), allocatable :: q(:), r(:)
    real(dp) :: beta
    integer :: k, idx

    if (size(g) /= st%ndof) error stop 'lbfgs_direction: size mismatch'

    if (st%count == 0) then
      p = -g
      return
    end if

    allocate(q(st%ndof), r(st%ndof))
    q = g

    ! First loop: newest -> oldest
    do k = 0, st%count-1
      idx = prev_idx(st, k)
      st%alpha(idx) = st%rho(idx) * dot_product(st%S(:,idx), q)
      q = q - st%alpha(idx) * st%Y(:,idx)
    end do

    ! Apply initial Hessian H0 ≈ γ I
    r = st%H0 * q

    ! Second loop: oldest -> newest
    do k = st%count-1, 0, -1
      idx = prev_idx(st, k)
      beta = st%rho(idx) * dot_product(st%Y(:,idx), r)
      r = r + st%S(:,idx) * (st%alpha(idx) - beta)
    end do

    p = -r   ! descent direction candidate
    if (dot_product(p, g) >= 0.0_dp) then
      ! Numerical safeguard: fall back to steepest descent if not a descent dir
      p = -g
      call lbfgs_reset(st)
    end if

    deallocate(q, r)
  end subroutine lbfgs_direction

end module lbfgs_direction_mod
