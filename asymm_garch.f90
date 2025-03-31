module asymm_garch_mod
use        kind_mod, only: dp
use optim_methods_mod, only: nelder_mead_str, uobyqa_str
use nelder_mead_mod, only: nelder_mead
use     obj_fun_mod, only: xret, xdist, gjr_garch_sigma, obj_sub
use basic_stats_mod, only: variance
use      random_mod, only: random_normal
use   constants_mod, only: pi
use        util_mod, only: assert_equal, default
use      uobyqa_mod, only: uobyqa
implicit none
private
public :: fit_gjr_garch, simulate_gjr_garch, nparam_gjr_garch
integer, parameter :: nparam_gjr_garch = 5
contains
!------------------------------------------------------------
! Subroutine to fit the GJR-GARCH(1,1) model.
!
! Input:
! ret      -- Array of returns.
! dist     -- Distribution flag (such as "normal" or "laplace").
! max_iter -- (optional) max # of Nelder-Mead iterations
! tol      -- (optional) convergence criterion for value of likelihood
! alpha0 -- (optional) guess for weight on past squared return
! gamma0 -- (optional) guess for weight on past squared return when return is negative
! beta0  -- (optional) guess for weight on previous variance
! opt_method -- (optional) optimization method
!
! Output:
! par_out -- Estimated parameters: [mu, ω, α, γ, β].
! logL    -- Final (maximized) log likelihood.
! info    -- Convergence flag (0 = converged, 1 = max_iter reached).
! niter   -- # of iterations in Nelder-Mead
!------------------------------------------------------------
subroutine fit_gjr_garch(ret, dist, par_out, logL, info, niter, max_iter, &
   tol, sigma, alpha0, gamma0, beta0, opt_method, rhobeg, rhoend)
character(len=*), intent(in) :: dist ! name of conditional probability distribution
real(kind=dp)   , intent(in) :: ret(:) ! returns
real(kind=dp)   , intent(out) :: par_out(:), logL
integer         , intent(out) :: info ! convergence flag (0 if converged)
integer, intent(out), optional :: niter ! # of iterations used
integer, intent(in) , optional :: max_iter
real(kind=dp), intent(in), optional :: tol
real(kind=dp), intent(out), optional :: sigma(:)
real(kind=dp), intent(in), optional :: alpha0, gamma0, beta0
character (len=*), intent(in), optional :: opt_method ! optimization method
real(kind=dp), intent(in), optional :: rhobeg, rhoend ! initial and final values of trust region radius used in uobyqa
real(kind=dp) :: x0(nparam_gjr_garch)
real(kind=dp) :: tol_
integer :: n, max_iter_
logical, parameter :: print_guess_ = .false.
character (len=*), parameter :: msg = "in fit_gjr_garch, "
character (len=20) :: opt_method_
integer :: iprint_
iprint_ = 0
opt_method_ = default(uobyqa_str, opt_method)
call assert_equal(size(par_out), nparam_gjr_garch, &
   msg // "size(par_out)")
max_iter_ = default(1000, max_iter)
n = size(ret)
if (present(sigma)) call assert_equal(size(sigma), n, msg // "size(sigma)")
!--- Set the module variables for the objective function.
xret = ret
xdist = dist
tol_ = default(1.0e-10_dp, tol)

!--- Set up an initial guess.
x0(1) = sum(ret) / n ! mu
x0(2) = 0.1_dp * variance(ret) ! ω
x0(3) = default(0.05_dp, alpha0) ! α
x0(4) = default(0.05_dp, gamma0) ! γ
x0(5) = default(0.9_dp, beta0)  ! β
if (print_guess_) print "(/,'in fit_gjr_garch, x0:', *(1x,f12.6))", x0
!--- Run the Nelder–Mead optimizer.
if (opt_method_ == nelder_mead_str) then
   call nelder_mead(x0, max_iter_, tol_, par_out, logL, info, &
         niter) ! out: par_out, logL
else if (opt_method_ == uobyqa_str) then
   call uobyqa(size(x0), x0, default(0.01_dp, rhobeg), default(1.0d-8, rhoend), iprint_, max_iter_, &
      obj_sub, fval=logL, nfun=niter)
   info = 0
   par_out = x0
else
   error stop "opt_method " // trim(opt_method_) // " not recognized, must be one of " &
              // nelder_mead_str // " " // uobyqa_str
end if
! Since we minimized the negative log likelihood, the maximized
! log–likelihood is its negative.
logL = -logL
if (present(sigma)) sigma = gjr_garch_sigma(ret, par_out)
end subroutine fit_gjr_garch

!------------------------------------------------------------
! Subroutine to simulate a GJR-GARCH(1,1) process.
!
! Input:
! mu    -- mean return
! omega -- variance constant
! alpha -- symmetric weight on past squared return
! gamma -- weight on past squared return when return is negative
! beta  -- weight on previous variance
!
! Output:
! ret -- Simulated return series (allocated inside).
!------------------------------------------------------------
subroutine simulate_gjr_garch(mu, omega, alpha, gamma, beta, ret, sigma)
real(kind=dp), intent(in) :: mu, omega, alpha, gamma, beta
real(kind=dp), intent(out) :: ret(:)
real(kind=dp), intent(out), optional :: sigma(:)
real(kind=dp) :: uc_var, sigma2, noise
integer :: i, n
real(kind=dp) :: z, sigma_i
n = size(ret)
if (present(sigma)) call assert_equal(size(sigma), n, &
   "in simulate_gjr_garch, size(sigma)")
uc_var = omega / (1.0_dp - alpha - 0.5_dp*gamma - beta)
sigma2 = uc_var
do i = 1, n
   z = random_normal()
   sigma_i = sqrt(sigma2)
   noise = sigma_i * z
   if (present(sigma)) sigma(i) = sigma_i
   ret(i) = mu + noise
   sigma2 = (alpha + merge(gamma,0.0_dp,noise < 0.0_dp)) * noise**2 &
            + beta*sigma2 + omega
end do
end subroutine simulate_gjr_garch
end module asymm_garch_mod