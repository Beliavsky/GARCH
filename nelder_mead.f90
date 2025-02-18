module nelder_mead_mod
use kind_mod, only: dp
use obj_fun_mod, only: obj_fun
implicit none
private
public :: nelder_mead
contains
  !------------------------------------------------------------
  ! Subroutine implementing the Nelder–Mead (downhill simplex)
  ! algorithm for unconstrained minimization.
  !
  ! Input:
  !   n        - Dimension of parameter space.
  !   x0       - Initial guess (vector of length n).
  !   max_iter - Maximum number of iterations.
  !   tol      - Tolerance (stopping criterion based on function values).
  !
  ! Output:
  !   xopt     - Estimated minimizer (vector of length n).
  !   fopt     - Function value at xopt.
  !   info     - Convergence flag (0 if converged, 1 if max_iter reached).
  !------------------------------------------------------------
  subroutine nelder_mead(x0, max_iter, tol, xopt, fopt, info, &
                         niter)
    integer, intent(in) :: max_iter
    real(kind=dp), intent(in) :: x0(:), tol
    real(kind=dp), intent(out) :: xopt(:), fopt
    integer, intent(out) :: info
    integer, optional, intent(out) :: niter
    integer :: i, j, k, iter, n, nvertices
    real(kind=dp), allocatable :: simplex(:,:)
    real(kind=dp), allocatable :: fvals(:)
    real(kind=dp), allocatable :: x_r(:), x_e(:), x_c(:), x_centroid(:)
    real(kind=dp) :: alpha, gamma, rho, sigma
    real(kind=dp) :: f_r, f_e, f_c
    real(kind=dp) :: f_diff
    integer, allocatable :: idx(:)
    logical, parameter :: debug = .false.
    n = size(x0)
    if (size(xopt) /= n) error stop "in nelder_mead, xopt has wrong size"
    ! Coefficients.
    alpha = 1.0_dp    ! reflection
    gamma = 2.0_dp    ! expansion
    rho   = 0.5_dp    ! contraction
    sigma = 0.5_dp    ! shrink

    nvertices = n + 1
    allocate(simplex(n, nvertices))
    allocate(fvals(nvertices))
    allocate(idx(nvertices))
    allocate(x_r(n), x_e(n), x_c(n), x_centroid(n))

    !--- Initialize the simplex.
    simplex(:,1) = x0
    do i = 2, nvertices
       simplex(:,i) = x0
       ! Perturb the (i-1)th coordinate.
       simplex(i-1, i) = simplex(i-1, i) + 0.05_dp*(abs(x0(i-1)) + 0.001_dp)
    end do
    do i = 1, nvertices
       fvals(i) = obj_fun(simplex(:,i))
    end do

    iter = 0
    do while (iter < max_iter)
       iter = iter + 1
       ! Sort the simplex vertices according to fvals.
       do i = 1, nvertices
          idx(i) = i
       end do
       ! Simple (bubble) sort.
       do i = 1, nvertices-1
          do j = i+1, nvertices
             if (fvals(idx(i)) > fvals(idx(j))) then
                k = idx(i)
                idx(i) = idx(j)
                idx(j) = k
             end if
          end do
       end do

       ! Best vertex is at idx(1), worst at idx(nvertices).
       ! Compute centroid of all points except the worst.
       x_centroid = 0.0_dp
       do i = 1, nvertices-1
          x_centroid = x_centroid + simplex(:, idx(i))
       end do
       x_centroid = x_centroid / real(n, dp)

       ! Reflection.
       x_r = x_centroid + alpha*(x_centroid - simplex(:, idx(nvertices)))
       f_r = obj_fun(x_r)

       if (f_r < fvals(idx(1))) then
          ! Expansion.
          x_e = x_centroid + gamma*(x_r - x_centroid)
          f_e = obj_fun(x_e)
          if (f_e < f_r) then
             simplex(:, idx(nvertices)) = x_e
             fvals(idx(nvertices)) = f_e
          else
             simplex(:, idx(nvertices)) = x_r
             fvals(idx(nvertices)) = f_r
          end if
       else if (f_r < fvals(idx(nvertices-1))) then
          simplex(:, idx(nvertices)) = x_r
          fvals(idx(nvertices)) = f_r
       else
          ! Contraction.
          if (f_r < fvals(idx(nvertices))) then
             ! Outside contraction.
             x_c = x_centroid + rho*(x_r - x_centroid)
             f_c = obj_fun(x_c)
          else
             ! Inside contraction.
             x_c = x_centroid - rho*(x_centroid - simplex(:, idx(nvertices)))
             f_c = obj_fun(x_c)
          end if
          if (f_c < fvals(idx(nvertices))) then
             simplex(:, idx(nvertices)) = x_c
             fvals(idx(nvertices)) = f_c
          else
             ! Shrink the simplex.
             do i = 2, nvertices
                simplex(:, idx(i)) = simplex(:, idx(1)) + sigma*(simplex(:, idx(i)) - simplex(:, idx(1)))
                fvals(idx(i)) = obj_fun(simplex(:, idx(i)))
             end do
          end if
       end if

       ! Convergence test.
       f_diff = abs(fvals(idx(nvertices)) - fvals(idx(1)))
       if (debug) then
          print "(/,a,i0,*(f12.4))", "iter, f_diff, fvals = ", iter, f_diff, fvals
          print "(a,*(f12.4))", "fopt, simplex(:, idx(1)) =", fvals(idx(1)), simplex(:, idx(1))
       end if
       if (f_diff < tol*(1.0_dp + abs(fvals(idx(1))))) exit
    end do

    ! Return the best vertex.
    xopt = simplex(:, idx(1))
    fopt = fvals(idx(1))
    info = merge(1, 0, iter >= max_iter)
    if (present(niter)) niter = iter
    deallocate(simplex, fvals, idx, x_r, x_e, x_c, x_centroid)
  end subroutine nelder_mead

end module nelder_mead_mod