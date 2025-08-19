!==============================================================================
! OSO Conjugate-Gradient Minimization with Strong-Wolfe Line Search
!==============================================================================

module oso_kinds_mod
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: dp
end module oso_kinds_mod

!------------------------------------------------------------------------------
! Axis/angle extraction and Rodrigues rotation (batch over all sites/ensembles)
!------------------------------------------------------------------------------
module oso_rotation_mod
  use oso_kinds_mod, only: dp
  implicit none
  private
  public :: angles_axes_from_a, rodrigues_apply_batch

contains

  subroutine angles_axes_from_a(Natom, Mensemble, avec, theta, axis)
    !! Map OSO "a" parameters -> rotation angle θ and unit axis r for each site.
    !! avec = (3,Natom,Mensemble) packs (a12,a13,a23) per site.
    integer, intent(in)  :: Natom, Mensemble
    real(dp), intent(in) :: avec(3, Natom, Mensemble)   ! (a12,a13,a23)
    real(dp), intent(out):: theta(Natom, Mensemble)     ! rotation angles
    real(dp), intent(out):: axis(3, Natom, Mensemble)   ! unit axes r
    integer :: i, j
    real(dp) :: a12, a13, a23, th, invth
    do j = 1, Mensemble
      do i = 1, Natom
        a12 = avec(1,i,j);  a13 = avec(2,i,j);  a23 = avec(3,i,j)
        th  = sqrt(a12*a12 + a13*a13 + a23*a23)
        theta(i,j) = th
        if (th > 1.0d-14) then
          invth = 1.0d0 / th
          axis(1,i,j) = -a23 * invth
          axis(2,i,j) =  a13 * invth
          axis(3,i,j) = -a12 * invth
        else
          axis(:,i,j) = 0.0d0
          axis(3,i,j) = 1.0d0  ! arbitrary axis when θ≈0
        end if
      end do
    end do
  end subroutine angles_axes_from_a


  subroutine rodrigues_apply_batch(Natom, Mensemble, theta, axis, emom_in, emom_out)
    !! Apply rotation defined by (theta, axis) to all spins in emom_in.
    integer, intent(in)  :: Natom, Mensemble
    real(dp), intent(in) :: theta(Natom, Mensemble), axis(3, Natom, Mensemble)
    real(dp), intent(in) :: emom_in(3, Natom, Mensemble)
    real(dp), intent(out):: emom_out(3, Natom, Mensemble)
    integer :: i, j
    real(dp) :: th, ct, st, omt, rdote
    real(dp) :: r1, r2, r3, e1, e2, e3, cx1, cx2, cx3, normv
    do j = 1, Mensemble
      do i = 1, Natom
        th = theta(i,j)
        e1 = emom_in(1,i,j); e2 = emom_in(2,i,j); e3 = emom_in(3,i,j)
        if (th < 1.0d-14) then
          emom_out(1,i,j) = e1
          emom_out(2,i,j) = e2
          emom_out(3,i,j) = e3
        else
          r1 = axis(1,i,j);  r2 = axis(2,i,j);  r3 = axis(3,i,j)
          ct = cos(th);  st = sin(th);  omt = 1.0d0 - ct
          ! r × e
          cx1 = r2*e3 - r3*e2
          cx2 = r3*e1 - r1*e3
          cx3 = r1*e2 - r2*e1
          rdote = r1*e1 + r2*e2 + r3*e3
          emom_out(1,i,j) = ct*e1 + st*cx1 + omt*rdote*r1
          emom_out(2,i,j) = ct*e2 + st*cx2 + omt*rdote*r2
          emom_out(3,i,j) = ct*e3 + st*cx3 + omt*rdote*r3
          ! normalize to unit length (cheap safeguard)
          normv = sqrt(emom_out(1,i,j)**2 + emom_out(2,i,j)**2 + emom_out(3,i,j)**2)
          if (normv > 0.0d0) then
            emom_out(:,i,j) = emom_out(:,i,j) / normv
          end if
        end if
      end do
    end do
  end subroutine rodrigues_apply_batch

end module oso_rotation_mod

!------------------------------------------------------------------------------
! Pack/unpack helpers: OSO gradient ↔ torque, and vector <-> (3,N,M)
!------------------------------------------------------------------------------
module oso_pack_mod
  use oso_kinds_mod, only: dp
  implicit none
  private
  public :: torque_to_grad, pack3_to_vec, vec_to_pack3, step_rms_angle

contains

  subroutine torque_to_grad(Natom, Mensemble, torque, grad_vec)
    !! Map torques to OSO gradient components g=(tz,-ty,tx) per site.
    integer, intent(in)  :: Natom, Mensemble
    real(dp), intent(in) :: torque(3, Natom, Mensemble)
    real(dp), intent(out):: grad_vec(3*Natom*Mensemble)
    integer :: i, j, k, idx
    idx = 0
    do j = 1, Mensemble
      do i = 1, Natom
        k = idx*0  ! dummy to silence some compilers about unused vars
        grad_vec(idx+1) = torque(3,i,j)   ! tz  -> a12
        grad_vec(idx+2) = -torque(2,i,j)  ! -ty -> a13
        grad_vec(idx+3) = torque(1,i,j)   ! tx  -> a23
        idx = idx + 3
      end do
    end do
  end subroutine torque_to_grad


  subroutine pack3_to_vec(Natom, Mensemble, pack3, vec)
    !! Flatten (3,N,M) to vec length 3*N*M in order (a12,a13,a23) site-major.
    integer, intent(in)  :: Natom, Mensemble
    real(dp), intent(in) :: pack3(3, Natom, Mensemble)
    real(dp), intent(out):: vec(3*Natom*Mensemble)
    integer :: i, j, idx
    idx = 0
    do j = 1, Mensemble
      do i = 1, Natom
        vec(idx+1:idx+3) = pack3(:,i,j)
        idx = idx + 3
      end do
    end do
  end subroutine pack3_to_vec


  subroutine vec_to_pack3(Natom, Mensemble, vec, pack3)
    !! Unflatten vec -> (3,N,M).
    integer, intent(in)  :: Natom, Mensemble
    real(dp), intent(in) :: vec(3*Natom*Mensemble)
    real(dp), intent(out):: pack3(3, Natom, Mensemble)
    integer :: i, j, idx
    idx = 0
    do j = 1, Mensemble
      do i = 1, Natom
        pack3(:,i,j) = vec(idx+1:idx+3)
        idx = idx + 3
      end do
    end do
  end subroutine vec_to_pack3


  pure function step_rms_angle(Natom, Mensemble, step_pack3) result(theta_rms)
    !! RMS rotation angle for a proposed OSO step "a" (per-site norm).
    integer, intent(in)  :: Natom, Mensemble
    real(dp), intent(in) :: step_pack3(3, Natom, Mensemble)
    real(dp) :: theta_rms
    integer :: i, j
    real(dp) :: s
    s = 0.0d0
    do j = 1, Mensemble
      do i = 1, Natom
        s = s + step_pack3(1,i,j)**2 + step_pack3(2,i,j)**2 + step_pack3(3,i,j)**2
      end do
    end do
    theta_rms = sqrt( s / real(Natom*Mensemble,dp) )
  end function step_rms_angle

end module oso_pack_mod

!------------------------------------------------------------------------------
! Conjugate-gradient direction (FR / PRP+ with safe restarts)
!------------------------------------------------------------------------------
module cg_direction_mod
  use oso_kinds_mod, only: dp
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: get_cg_direction

contains

  subroutine get_cg_direction(g_curr, g_prev, p_prev, p_curr, method, iter, restart_period)
    real(dp), intent(in)  :: g_curr(:), g_prev(:), p_prev(:)
    real(dp), intent(out) :: p_curr(:)
    character(*), intent(in) :: method   ! 'FR' or 'PRP' (case-insensitive)
    integer, intent(in) :: iter, restart_period

    real(dp) :: gg_curr, gg_prev, beta, pg
    real(dp), parameter :: tiny = 1.0d-300
    logical :: do_restart
    character(len=:), allocatable :: m

    if (any([size(g_curr),size(g_prev),size(p_prev),size(p_curr)] /= size(g_curr))) &
      error stop 'get_cg_direction: size mismatch.'

    m = to_upper(trim(method))
    do_restart = (iter == 0)
    if (restart_period > 0) do_restart = do_restart .or. (mod(iter,restart_period) == 0)

    if (do_restart) then
      p_curr = -g_curr
      return
    end if

    gg_curr = dot_product(g_curr,g_curr)
    gg_prev = max(dot_product(g_prev,g_prev), tiny)

    select case (m)
    case ('FR')
      beta = gg_curr / gg_prev
    case ('PRP')
      beta = dot_product(g_curr, g_curr - g_prev) / gg_prev
      if (beta < 0.0_dp) beta = 0.0_dp
    case default
      error stop 'get_cg_direction: method must be FR or PRP.'
    end select

    p_curr = -g_curr + beta * p_prev
    pg = dot_product(p_curr, g_curr)

    if (pg >= 0.0_dp .or. .not.all_finite(p_curr)) p_curr = -g_curr
  end subroutine get_cg_direction


  pure logical function all_finite(x)
    real(dp), intent(in) :: x(:)
    integer :: i
    all_finite = .true.
    do i=1,size(x)
      if (.not. ieee_is_finite(x(i))) then
        all_finite = .false.; return
      end if
    end do
  end function all_finite

  pure function to_upper(s) result(t)
    character(*), intent(in) :: s
    character(len=:), allocatable :: t
    integer :: i, n, c
    integer, parameter :: diff = iachar('a') - iachar('A')

    n = len_trim(s)
    if (n <= 0) then
      allocate(character(len=0) :: t)
      return
    end if

    allocate(character(len=n) :: t)
    do i=1,n
      c = iachar(s(i:i))
      if (c >= iachar('a') .and. c <= iachar('z')) then
        t(i:i) = achar(c - diff)
      else
        t(i:i) = s(i:i)
      end if
    end do
  end function to_upper

end module cg_direction_mod

!------------------------------------------------------------------------------
! Strong-Wolfe line search (auto-bracketing + zoom; cubic/quad + bisection)
!------------------------------------------------------------------------------
module line_search_wolfe_mod
  use oso_kinds_mod, only: dp
  implicit none
  private
  public :: line_search_strong_wolfe

  abstract interface
    subroutine eval_phi_iface(alpha, phi, dphi)
      import dp
      real(dp), intent(in)  :: alpha
      real(dp), intent(out) :: phi
      real(dp), intent(out) :: dphi
    end subroutine eval_phi_iface
  end interface

contains

  subroutine line_search_strong_wolfe(eval_phi, c1, c2, alpha_init, alpha_max, &
                                      max_evals, phi0, dphi0, alpha_star, status)
    procedure(eval_phi_iface) :: eval_phi
    real(dp), intent(in)  :: c1, c2, alpha_init, alpha_max
    integer, intent(in)   :: max_evals
    real(dp), intent(in)  :: phi0, dphi0             ! at alpha=0
    real(dp), intent(out) :: alpha_star
    integer, intent(out)  :: status                  ! 0 ok; >0 suboptimal

    real(dp) :: a0, a1, phi1, dphi1
    real(dp) :: ap, phip, dphip
    integer  :: ne

    if (dphi0 >= 0.0_dp) then
      alpha_star = 0.0_dp; status = 3; return   ! not a descent direction
    end if

    a0 = 0.0_dp
    ap = a0; phip = phi0; dphip = dphi0
    a1 = min(alpha_init, alpha_max)
    ne = 0

    do
      call eval_phi(a1, phi1, dphi1); ne = ne + 1

      if ((phi1 > phi0 + c1*a1*dphi0) .or. ((ne > 1) .and. (phi1 >= phip))) then
        call zoom(eval_phi, c1, c2, ap, a1, phip, dphip, phi0, dphi0,          &
                  alpha_star, status, max_evals-ne)
        return
      end if

      if (abs(dphi1) <= c2*abs(dphi0)) then
        alpha_star = a1; status = 0; return
      end if

      if (dphi1 >= 0.0_dp) then
        call zoom(eval_phi, c1, c2, a1, ap, phi1, dphi1, phi0, dphi0,          &
                  alpha_star, status, max_evals-ne)
        return
      end if

      ap   = a1; phip = phi1; dphip = dphi1

      if (a1 >= alpha_max*0.999_dp .or. ne >= max_evals) then
        alpha_star = a1; status = 1; return
      end if

      a1 = min(2.0_dp*a1, alpha_max)
    end do
  end subroutine line_search_strong_wolfe


  subroutine zoom(eval_phi, c1, c2, alo, ahi, philo, dphilo, phi0, dphi0,     &
                  alpha_star, status, budget)
    procedure(eval_phi_iface) :: eval_phi
    real(dp), intent(in)  :: c1, c2, alo, ahi, philo, dphilo, phi0, dphi0
    real(dp), intent(out) :: alpha_star
    integer, intent(out)  :: status
    integer, intent(in)   :: budget

    real(dp) :: a_lo, a_hi, phi_lo, dphi_lo, phi_hi, dphi_hi
    real(dp) :: aj, phi_j, dphi_j, width, aj_min, aj_max
    integer :: ne
    real(dp), parameter :: tau = 0.1_dp

    a_lo = min(alo, ahi); a_hi = max(alo, ahi)
    phi_lo = philo; dphi_lo = dphilo
    call eval_phi(a_hi, phi_hi, dphi_hi); ne = 1

    do
      width = a_hi - a_lo
      if (width <= 1.0d-12*a_hi) then
        alpha_star = 0.5_dp*(a_lo+a_hi); status = 5; return
      end if

      ! Safeguarded interpolation (cubic/quad) -> [a_lo+τw, a_hi-τw]
      aj = cubic_or_quadratic(a_lo, a_hi, phi_lo, dphi_lo, phi_hi, dphi_hi)
      aj_min = a_lo + tau*width; aj_max = a_hi - tau*width
      if ( (aj <= aj_min) .or. (aj >= aj_max) ) aj = 0.5_dp*(a_lo+a_hi)

      call eval_phi(aj, phi_j, dphi_j); ne = ne + 1

      if ((phi_j > phi0 + c1*aj*dphi0) .or. (phi_j >= phi_lo)) then
        a_hi = aj; phi_hi = phi_j; dphi_hi = dphi_j
      else
        if (abs(dphi_j) <= c2*abs(dphi0)) then
          alpha_star = aj; status = 0; return
        end if
        if (dphi_j*(a_hi - a_lo) >= 0.0_dp) then
          a_hi = a_lo; phi_hi = phi_lo; dphi_hi = dphi_lo
        end if
        a_lo = aj; phi_lo = phi_j; dphi_lo = dphi_j
      end if

      if (ne >= max(10,budget)) then
        alpha_star = aj; status = 4; return
      end if
    end do
  end subroutine zoom


  pure real(dp) function cubic_or_quadratic(a, b, fa, dfa, fb, dfb) result(x)
    !! Return an interpolated point inside (a,b) using cubic if possible,
    !! otherwise quadratic based on (a,fa,dfa) and (b,fb,dfb).
    real(dp), intent(in) :: a, b, fa, dfa, fb, dfb
    real(dp) :: d1, d2, Acoef, Bcoef, Ccoef, Dcoef, denom, rad, xm
    real(dp), parameter :: eps = 1.0d-18
    xm = 0.5_dp*(a+b)

    ! Try cubic: p(t)=A t^3 + B t^2 + C t + D with t in {0,1} mapped to {a,b}
    ! Solve constraints: p(0)=fa, p'(0)=dfa*(b-a), p(1)=fb, p'(1)=dfb*(b-a).
    d1 = dfa*(b-a); d2 = dfb*(b-a)
    Acoef = 2.0d0*(fa - fb) + d1 + d2
    Bcoef = -3.0d0*(fa - fb) - 2.0d0*d1 - d2
    Ccoef = d1
    Dcoef = fa
    ! derivative p'(t)=3At^2 + 2Bt + C ; set to 0 for stationary point(s).
    rad = Bcoef*Bcoef - 3.0d0*Acoef*Ccoef
    if (abs(Acoef) > eps .and. rad >= 0.0_dp) then
      rad = sqrt(max(0.0_dp, rad))
      ! two candidates in t-space
      x = a + (b-a) * max( min( (-2.0d0*Bcoef + 2.0d0*rad) / (6.0d0*Acoef) , 0.9d0 ), 0.1d0 )
      return
    end if

    ! Quadratic fallback using (a,fa,dfa)
    denom = 2.0d0*(fa - fb) + (dfa + dfb)*(b-a)
    if (abs(denom) > eps) then
      x = a + (b-a) * (dfa*(b-a)) / denom
      x = min( max(x, a+0.1d0*(b-a)), b-0.1d0*(b-a) )
    else
      x = xm
    end if
  end function cubic_or_quadratic

end module line_search_wolfe_mod

!------------------------------------------------------------------------------
! OSO problem definition and the CG driver with Wolfe line search
!------------------------------------------------------------------------------
module oso_optimize_mod
  use oso_kinds_mod,        only: dp
  use oso_rotation_mod,     only: angles_axes_from_a, rodrigues_apply_batch
  use oso_pack_mod,         only: torque_to_grad, pack3_to_vec, vec_to_pack3, step_rms_angle
  use cg_direction_mod,     only: get_cg_direction
  use line_search_wolfe_mod,only: line_search_strong_wolfe
  implicit none
  private
  public :: calculate_torques_and_energy, oso_problem, oso_cg_minimize
  public :: oso_run

! In module oso_optimize_mod, replace the current abstract interface block:
abstract interface
  subroutine energy_torque_iface(Natom, Mensemble, emom, energy, torque)
    import dp
    integer, intent(in) :: Natom, Mensemble
    real(dp), intent(in)  :: emom(3, Natom, Mensemble)       ! unit directions from optimizer
    real(dp), intent(out) :: energy
    real(dp), intent(out) :: torque(3, Natom, Mensemble)     ! τ_i = - m_i × Beff_i  (see note)
  end subroutine energy_torque_iface
end interface

type :: oso_problem
  integer :: Natom = 0, Mensemble = 1
  procedure(energy_torque_iface), pointer, nopass :: eval => null()
end type oso_problem

! Module local array for magnetic fields
real(dp), dimension(:,:,:), allocatable :: bfield

contains

    subroutine calculate_torques_and_energy(Natom, Mensemble, moments, energy, torques)
        use math_functions, only : f_cross_product
        use HamiltonianActions, only : effective_field_standalone
          implicit none

          integer, intent(in) :: Natom        !< Number of atoms in system
          integer, intent(in) :: Mensemble    !< Number of ensembles
          real(dp), dimension(3,Natom,Mensemble), intent(in) :: moments !! Magnetic torque
          real(dp), dimension(3,Natom,Mensemble), intent(out) :: torques !! Magnetic torque
          real(dp), intent(out) :: energy !! Total energy

          integer :: iatom, iens

          ! Get effective field as bfield (module-local array, needs to be allocated)
          if (.not. allocated(bfield)) allocate(bfield(3,Natom,Mensemble))

          call effective_field_standalone(moments, bfield, energy)
          !print *,'Calculating effective field and energy', energy

          ! Loop to get torques
          !$omp parallel do default(shared) private(iens, iatom) collapse(2)
          do iens = 1, Mensemble
                do iatom = 1, Natom
                      torques(:, iatom, iens) = f_cross_product(moments(:, iatom, iens), bfield(:, iatom, iens))
                end do
          end do
          !$omp end parallel do

          return

    end subroutine calculate_torques_and_energy

subroutine oso_cg_minimize(prob, emom, maxiter, tol_grad,                    &
                           c1, c2, theta_rms_cap, reanchor_period,           &
                           method, energy_out, iters_done)
  use oso_kinds_mod,         only : dp
  use oso_rotation_mod,      only : angles_axes_from_a, rodrigues_apply_batch
  use oso_pack_mod,          only : torque_to_grad, vec_to_pack3, step_rms_angle
  use cg_direction_mod,      only : get_cg_direction
  use line_search_wolfe_mod, only : line_search_strong_wolfe
  use lbfgs_direction_mod,   only : lbfgs_state, lbfgs_init, lbfgs_reset,     &
                                    lbfgs_push_pair, lbfgs_direction
  implicit none
  ! --- inputs/outputs ---
  type(oso_problem), intent(in)        :: prob
  real(dp), intent(inout)              :: emom(3, prob%Natom, prob%Mensemble) ! UPDATED in place
  integer,  intent(in)                 :: maxiter
  real(dp), intent(in)                 :: tol_grad
  real(dp), intent(in)                 :: c1, c2           ! Wolfe; e.g. 1e-4, 0.9
  real(dp), intent(in)                 :: theta_rms_cap    ! e.g. 0.35–0.40 (radians)
  integer,  intent(in)                 :: reanchor_period  ! e.g. 50; 0 means never
  character(*), intent(in)             :: method           ! 'PRP' | 'FR' | 'LBFGS'
  real(dp), intent(out)                :: energy_out
  integer,  intent(out)                :: iters_done

  ! --- locals ---
  integer :: N, M, ndof, iter
  real(dp), allocatable :: a_total(:), g(:), g_old(:), p(:), p_old(:)
  real(dp), allocatable :: a_step_pack(:,:,:), a_trial_pack(:,:,:)
  real(dp), allocatable :: theta(:,:), axis(:,:,:)
  real(dp), allocatable :: emom_ref(:,:,:), emom_trial(:,:,:), torque(:,:,:)
  real(dp), allocatable :: g_tmp(:), s_vec(:), y_vec(:)
  real(dp) :: energy, phi0, dphi0, alpha_init, alpha_max, alpha_star
  real(dp) :: gnorm_inf, alpha_prev
  integer  :: status_ls
  type(lbfgs_state) :: lb
  integer :: lbfgs_m
  real(dp) :: theta_rms_cap_0

  N = prob%Natom;  M = prob%Mensemble;  ndof = 3*N*M ; theta_rms_cap_0 = theta_rms_cap
  if (.not. associated(prob%eval)) error stop 'oso_cg_minimize: prob%eval not set.'

  allocate(a_total(ndof), g(ndof), g_old(ndof), p(ndof), p_old(ndof))
  allocate(a_step_pack(3,N,M), a_trial_pack(3,N,M))
  allocate(theta(N,M), axis(3,N,M))
  allocate(emom_ref(3,N,M), emom_trial(3,N,M), torque(3,N,M))
  allocate(g_tmp(ndof), s_vec(ndof), y_vec(ndof))

  ! L-BFGS history (active only if method='LBFGS')
  lbfgs_m = 10
  call lbfgs_init(lb, ndof, lbfgs_m)

  ! Reference configuration and zero initial OSO parameters
  emom_ref = emom
  a_total  = 0.0_dp

  ! Initial energy/gradient
  call prob%eval(N, M, emom, energy, torque)
  call torque_to_grad(N, M, torque, g)
  energy_out = energy

  gnorm_inf = maxval(abs(g))
  if (gnorm_inf <= tol_grad) then
    iters_done = 0
    return
  end if

  ! Initial search direction (steepest) and histories
  p      = -g
  p_old  = p
  g_old  = g
  alpha_prev = 1.0_dp

  do iter = 1, maxiter
    ! --- Line-search setup on current (g, p) ---
    phi0  = energy
    dphi0 = dot_product(g, p)

    ! Ensure descent direction; if not, restart (and reset LBFGS history)
    if (dphi0 >= 0.0_dp) then
      p = -g
      dphi0 = -dot_product(g, g)
      if (is_lbfgs(method)) call lbfgs_reset(lb)
    end if

    ! α_max from θ_rms cap for α=1 (scale if needed); α_init from previous step
    call vec_to_pack3(N, M, p, a_step_pack)
    alpha_max = 1.0d6
    if (theta_rms_cap > 0.0_dp) then
      alpha_max = min(alpha_max, theta_rms_cap / max(step_rms_angle(N,M,a_step_pack), 1.0e-16_dp))
    end if
    alpha_init = min( max(0.25_dp, 1.5_dp*alpha_prev), alpha_max )

    call line_search_strong_wolfe( eval_phi, c1, c2, alpha_init, alpha_max,   &
                                   40, phi0, dphi0, alpha_star, status_ls )

    ! --- Accept step (update a_total and emom) ---
    a_total = a_total + alpha_star * p
    call vec_to_pack3(N, M, a_total, a_trial_pack)
    call angles_axes_from_a(N,M, a_trial_pack, theta, axis)
    call rodrigues_apply_batch(N,M, theta, axis, emom_ref, emom)

    ! New energy/gradient at accepted point
    call prob%eval(N, M, emom, energy, torque)
    call torque_to_grad(N,M, torque, g)
    energy_out = energy
    alpha_prev = alpha_star

    ! Convergence test
    gnorm_inf = maxval(abs(g))
    if (gnorm_inf <= tol_grad) then
      iters_done = iter
      exit
    end if

    print *, 'Iteration', iter, ': energy =', energy/N, ', |g|_inf =', gnorm_inf
    write(200,'(A,i7, A,F12.8,A,F10.4)') 'Iteration', iter, ' : energy =', energy/N, ' , |g|_inf =', gnorm_inf

    ! --- Build step/curvature pair for LBFGS (s = α*p_old, y = g - g_old) ---
    s_vec = alpha_star * p
    y_vec = g - g_old
    if (is_lbfgs(method)) call lbfgs_push_pair(lb, s_vec, y_vec)

    ! --- Choose next direction ---
    if (is_lbfgs(method)) then
      call lbfgs_direction(lb, g, p)
    else
      call get_cg_direction(g, g_old, p_old, p, method, iter, 50)
    end if

    ! --- Shift histories for next iteration ---
    g_old = g
    p_old = p

    ! --- Optional re-anchoring of the reference frame ---
    if (reanchor_period > 0 .and. mod(iter, reanchor_period) == 0) then
      emom_ref = emom
      a_total  = 0.0_dp
      if (is_lbfgs(method)) call lbfgs_reset(lb)
    end if
  end do

contains
  logical function is_lbfgs(meth) result(ans)
    character(*), intent(in) :: meth
    ans = (meth == 'LBFGS' .or. meth == 'lbfgs' .or. meth == 'Lbfgs')
  end function is_lbfgs

  subroutine eval_phi(alpha, phi, dphi)
    !! φ(α) = E(e_ref rotated by a_total + α p);  φ'(α) = g(α)·p
    real(dp), intent(in)  :: alpha
    real(dp), intent(out) :: phi, dphi

    ! Trial a and trial spins
    call vec_to_pack3(N, M, a_total + alpha*p, a_trial_pack)
    call angles_axes_from_a(N,M, a_trial_pack, theta, axis)
    call rodrigues_apply_batch(N,M, theta, axis, emom_ref, emom_trial)

    ! Energy and gradient at trial
    call prob%eval(N, M, emom_trial, phi, torque)
    call torque_to_grad(N, M, torque, g_tmp)
    dphi = dot_product(g_tmp, p)
  end subroutine eval_phi
end subroutine oso_cg_minimize

!!!   subroutine oso_cg_minimize(prob, moments, maxiter, tol_grad,                    &
!!!                              c1, c2, theta_rms_cap, reanchor_period,           &
!!!                              method, energy_out, iters_done)
!!!     !! Main optimizer.
!!!     type(oso_problem), intent(in)        :: prob
!!!     real(dp), intent(inout)              :: moments(3, prob%Natom, prob%Mensemble) ! spins updated in place
!!!     integer,  intent(in)                 :: maxiter
!!!     real(dp), intent(in)                 :: tol_grad
!!!     real(dp), intent(in)                 :: c1, c2           ! Wolfe; typical 1e-4, 0.9
!!!     real(dp), intent(in)                 :: theta_rms_cap    ! e.g. 0.35 (radians)
!!!     integer,  intent(in)                 :: reanchor_period  ! e.g. 50; 0 means never
!!!     character(*), intent(in)             :: method           ! 'PRP' or 'FR'
!!!     real(dp), intent(out)                :: energy_out
!!!     integer,  intent(out)                :: iters_done
!!! 
!!!     integer :: N, M, ndof, iter
!!!     real(dp), allocatable :: a_total(:), g(:), g_prev(:), p(:), p_prev(:)
!!!     real(dp), allocatable :: a_step_pack(:,:,:), a_trial_pack(:,:,:)
!!!     real(dp), allocatable :: theta(:,:), axis(:,:,:)
!!!     real(dp) :: energy, phi0, dphi0, alpha_init, alpha_max, alpha_star, alpha_prev
!!!     real(dp) :: gnorm_inf, status_dummy
!!!     integer  :: status_ls
!!!     real(dp), allocatable :: emom_ref(:,:,:), emom_trial(:,:,:), torque(:,:,:)
!!! 
!!!     N = prob%Natom
!!!     M = prob%Mensemble
!!!     ndof = 3*N*M
!!! 
!!!     alpha_prev = 1.0_dp
!!! 
!!!     if (.not. associated(prob%eval)) error stop 'oso_cg_minimize: prob%eval not set.'
!!! 
!!!     allocate(a_total(ndof), g(ndof), g_prev(ndof), p(ndof), p_prev(ndof))
!!!     allocate(a_step_pack(3,N,M), a_trial_pack(3,N,M))
!!!     allocate(theta(N,M), axis(3,N,M))
!!!     allocate(emom_ref(3,N,M), emom_trial(3,N,M), torque(3,N,M))
!!! 
!!!     ! Reference spins and zero initial "a".
!!!     emom_ref = moments
!!!     a_total  = 0.0d0
!!! 
!!!     call prob%eval(N, M, moments, energy, torque)
!!!     call torque_to_grad(N, M, torque, g)
!!!     energy_out = energy
!!! 
!!!     gnorm_inf = maxval(abs(g))
!!!     if (gnorm_inf <= tol_grad) then
!!!       iters_done = 0; return
!!!     end if
!!! 
!!!     ! Initial direction: steepest descent
!!!     p = -g
!!!     g_prev = g
!!!     p_prev = p
!!! 
!!! 
!!!     do iter = 1, maxiter
!!! 
!!!       ! ---- Line-search callback φ(α) and φ'(α) = g(α)·p ----
!!!       phi0  = energy
!!!       dphi0 = dot_product(g, p)
!!! 
!!!       ! Cap alpha via RMS rotation for α=1:  α_max = min(user cap, θ_cap / θ_rms(α=1))
!!!       call vec_to_pack3(N, M, p, a_step_pack)
!!!       alpha_max = 1.0d6
!!!       if (theta_rms_cap > 0.0d0) then
!!!   alpha_max = min(alpha_max, theta_rms_cap / max(step_rms_angle(N,M,a_step_pack), 1.0e-16_dp))
!!!       end if
!!!       alpha_init = min(1.0d0, alpha_max)
!!!       ! alpha_init = min( max(0.25d0, 1.5d0*alpha_prev), alpha_max )
!!! 
!!! 
!!!       call line_search_strong_wolfe( eval_phi, c1, c2, alpha_init, alpha_max, &
!!!                                      40, phi0, dphi0, alpha_star, status_ls )
!!! 
!!!       ! ---- Accept step ----
!!!       a_total = a_total + alpha_star * p
!!!       alpha_prev = alpha_star
!!! 
!!!       call vec_to_pack3(N, M, a_total, a_trial_pack)
!!!       call angles_axes_from_a(N,M, a_trial_pack, theta, axis)
!!!       call rodrigues_apply_batch(N,M, theta, axis, emom_ref, moments)
!!! 
!!!       call prob%eval(N, M, moments, energy, torque)
!!!       call torque_to_grad(N,M, torque, g)
!!!       energy_out = energy
!!! 
!!!     gnorm_inf = maxval(abs(g))
!!!     print *, 'Iteration', iter, ': energy =', energy/N, ', |g|_inf =', gnorm_inf
!!!     write(200,'(A,i7, A,F12.8,A,F10.4)') 'Iteration', iter, ' : energy =', energy/N, ' , |g|_inf =', gnorm_inf
!!!       if (gnorm_inf <= tol_grad) then
!!!         iters_done = iter; exit
!!!       end if
!!! 
!!!       ! ---- CG update ----
!!!         call get_cg_direction(g, g_prev, p_prev, p, method, iter, 50)
!!! 
!!!         ! Shift state for next iteration
!!!         g_prev = g
!!!         p_prev = p
!!! 
!!!       ! ---- CG update ---- (wrong order)
!!!       ! g_prev = g
!!!       ! p_prev = p
!!!       ! call get_cg_direction(g, g_prev, p_prev, p, method, iter, 50)
!!! 
!!!       ! ---- Optional re-anchoring ----
!!!       if (reanchor_period > 0 .and. mod(iter, reanchor_period) == 0) then
!!!         emom_ref = moments
!!!         a_total  = 0.0d0
!!!       end if
!!! 
!!!     end do
!!! 
!!!   contains
!!! 
!!!     subroutine eval_phi(alpha, phi, dphi)
!!!       !! Evaluate φ(α) = E(e_ref rotated by a_total + α p), and φ'(α)=g·p.
!!!       real(dp), intent(in)  :: alpha
!!!       real(dp), intent(out) :: phi, dphi
!!! 
!!!       ! Make trial a = a_total + alpha * p
!!!       a_trial_pack = 0.0d0
!!!       call vec_to_pack3(N, M, a_total + alpha*p, a_trial_pack)
!!! 
!!!       call angles_axes_from_a(N,M, a_trial_pack, theta, axis)
!!!       call rodrigues_apply_batch(N,M, theta, axis, emom_ref, emom_trial)
!!! 
!!!       call prob%eval(N, M, emom_trial, phi, torque)
!!!       call torque_to_grad(N, M, torque, g_prev)  ! reuse g_prev as temp
!!!       dphi = dot_product(g_prev, p)
!!!     end subroutine eval_phi
!!! 
!!!   end subroutine oso_cg_minimize


subroutine oso_run(N, M, emom, eval_proc, maxiter, tol_grad, c1, c2,          &
                   theta_rms_cap, reanchor_period, method, energy_out, iters_done)
  use oso_kinds_mod,     only : dp
!  use oso_optimize_mod,  only : oso_problem, oso_cg_minimize, energy_torque_iface
  implicit none
  integer,  intent(in)    :: N, M
  real(dp), intent(inout) :: emom(3,N,M)
  procedure(energy_torque_iface) :: eval_proc
  integer,  intent(in)    :: maxiter, reanchor_period
  real(dp), intent(in)    :: tol_grad, c1, c2, theta_rms_cap
  character(*), intent(in):: method
  real(dp), intent(out)   :: energy_out
  integer,  intent(out)   :: iters_done

  type(oso_problem) :: prob
  integer :: i, j
  real(dp) :: nrm
  character(len=:), allocatable :: meth

  ! 0) light sanity on method (oso_cg_minimize is case-insensitive anyway)
  meth = adjustl(method)

  ! 1) ensure unit directions (cheap safeguard)
  do j=1,M
    do i=1,N
      nrm = sqrt(emom(1,i,j)**2 + emom(2,i,j)**2 + emom(3,i,j)**2)
      if (nrm > 0.0_dp) emom(:,i,j) = emom(:,i,j) / nrm
    end do
  end do

  ! 2) bind problem + callback and dispatch to the driver
  prob%Natom      = N
  prob%Mensemble  = M
  prob%eval       => eval_proc

  call oso_cg_minimize(prob, emom, maxiter, tol_grad, c1, c2,                  &
                       theta_rms_cap, reanchor_period, meth, energy_out, iters_done)
end subroutine oso_run


end module oso_optimize_mod

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
