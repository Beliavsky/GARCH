module garch_mod
use        kind_mod, only: dp
use nelder_mead_mod, only: nelder_mead
use     obj_fun_mod, only: xret, xdist, gjr_garch_sigma, garch_sigma, &
                           garch_model, nparam_gjr_garch, nparam_garch, &
                           garch_str, gjr_garch_str, nparam_ew_arch, &
                           ew_arch_str, ew_arch_sigma, arch_str, arch_sigma
use basic_stats_mod, only: variance, basic_stats_names, basic_stats, kurtosis, &
                           acf
use      random_mod, only: random_normal
use   constants_mod, only: pi
use        util_mod, only: assert_equal, default, write_merge
use iso_fortran_env, only: output_unit
implicit none
private
public :: fit_gjr_garch, simulate_gjr_garch, nparam_gjr_garch, &
   fit_garch, simulate_garch, nparam_garch, print_garch_results, &
   garch_str, gjr_garch_str, print_garch_lik, fit_ew_arch, &
   nparam_ew_arch, ew_arch_str, simulate_ew_arch, simulate_arch, &
   fit_arch
contains
!------------------------------------------------------------
! Subroutine to fit the GJR-GARCH(1,1) model.
!
! Input:
! ret      -- Array of returns.
! dist     -- Distribution flag ("normal" or "laplace").
! max_iter -- (optional) max # of Nelder-Mead iterations
! tol      -- (optional) convergence criterion for value of likelihood
! alpha0 -- (optional) guess for weight on past squared return
! gamma0 -- (optional) guess for weight on past squared return when return is negative
! beta0  -- (optional) guess for weight on previous variance
!
! Output:
! par_out -- Estimated parameters: [mu, ω, α, γ, β].
! logL    -- Final (maximized) log likelihood.
! info    -- Convergence flag (0 = converged, 1 = max_iter reached).
! niter   -- # of iterations in Nelder-Mead
!------------------------------------------------------------
subroutine fit_gjr_garch(ret, dist, par_out, logL, info, niter, max_iter, &
   tol, sigma, alpha0, gamma0, beta0)
character(len=*), intent(in) :: dist
real(kind=dp)   , intent(in) :: ret(:)
real(kind=dp)   , intent(out) :: par_out(:), logL
integer         , intent(out) :: info
integer, intent(out), optional :: niter
integer, intent(in) , optional :: max_iter
real(kind=dp), intent(in), optional :: tol
real(kind=dp), intent(out), optional :: sigma(:)
real(kind=dp), intent(in), optional :: alpha0, gamma0, beta0
real(kind=dp) :: x0(nparam_gjr_garch)
real(kind=dp) :: tol_
integer :: n, max_iter_
logical, parameter :: print_guess_ = .false.
character (len=*), parameter :: msg = "in fit_gjr_garch, "
garch_model = gjr_garch_str
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
call nelder_mead(x0, max_iter_, tol_, par_out, logL, info, &
         niter) ! output: par_out, logL
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

!------------------------------------------------------------
! Subroutine to fit the GARCH(1,1) model.
!
! Input:
! ret      -- Array of returns.
! dist     -- Distribution flag ("normal" or "laplace").
! max_iter -- (optional) max # of Nelder-Mead iterations
! tol      -- (optional) convergence criterion for value of likelihood
! alpha0 -- (optional) guess for weight on past squared return
! beta0  -- (optional) guess for weight on previous variance
!
! Output:
! par_out -- Estimated parameters: [mu, ω, α, β].
! logL    -- Final (maximized) log likelihood
! info    -- Convergence flag (0 = converged, 1 = max_iter reached).
! niter   -- # of iterations in Nelder-Mead
!------------------------------------------------------------
subroutine fit_garch(ret, dist, par_out, logL, info, niter, max_iter, &
   tol, sigma, alpha0, beta0)
character(len=*), intent(in) :: dist
real(kind=dp)   , intent(in) :: ret(:)
real(kind=dp)   , intent(out) :: par_out(:), logL
integer         , intent(out) :: info
integer, intent(out), optional :: niter
integer, intent(in) , optional :: max_iter
real(kind=dp), intent(in), optional :: tol
real(kind=dp), intent(out), optional :: sigma(:)
real(kind=dp), intent(in), optional :: alpha0, beta0
real(kind=dp) :: x0(nparam_garch)
real(kind=dp) :: tol_
integer :: n, max_iter_
logical, parameter :: print_guess_ = .false.
character (len=*), parameter :: msg = "in fit_garch, "
garch_model = garch_str
call assert_equal(size(par_out), nparam_garch, &
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
x0(4) = default(0.9_dp, beta0)  ! β
if (print_guess_) print "(/,'in fit_gjr_garch, x0:', *(1x,f12.6))", x0
!--- Run the Nelder–Mead optimizer.
call nelder_mead(x0, max_iter_, tol_, par_out, logL, info, &
         niter) ! output: par_out, logL
! Since we minimized the negative log likelihood, the maximized
! log–likelihood is its negative.
logL = -logL
if (present(sigma)) sigma = garch_sigma(ret, par_out)
end subroutine fit_garch

!------------------------------------------------------------
! Subroutine to simulate a GARCH(1,1) process.
!
! Input:
! mu    -- mean return
! omega -- variance constant
! alpha -- symmetric weight on past squared return
! beta  -- weight on previous variance
!
! Output:
! ret -- Simulated return series (allocated inside).
!------------------------------------------------------------
subroutine simulate_garch(mu, omega, alpha, beta, ret, sigma)
real(kind=dp), intent(in) :: mu, omega, alpha, beta
real(kind=dp), intent(out) :: ret(:)
real(kind=dp), intent(out), optional :: sigma(:)
real(kind=dp) :: uc_var, sigma2, noise
integer :: i, n
real(kind=dp) :: z, sigma_i
n = size(ret)
if (present(sigma)) call assert_equal(size(sigma), n, &
   "in simulate_garch, size(sigma)")
uc_var = omega / (1.0_dp - alpha - beta)
sigma2 = uc_var
do i = 1, n
   z = random_normal()
   sigma_i = sqrt(sigma2)
   noise = sigma_i * z
   if (present(sigma)) sigma(i) = sigma_i
   ret(i) = mu + noise
   sigma2 = alpha*noise**2 + beta*sigma2 + omega
end do
end subroutine simulate_garch

subroutine print_garch_results(garch_model, dist, est_param, logl, &
   niter, info, sigma_est, ret, nacfsq, nacfabs, fmt_header)
! print estimated GARCH parameters and properties of conditional standard
! deviations and normalized residuals
character (len=*), intent(in) :: garch_model  ! type of GARCH model
character (len=*), intent(in) :: dist         ! conditional distribution
real(kind=dp)    , intent(in) :: est_param(:) ! estimated parameters
real(kind=dp)    , intent(in) :: logl         ! log-likelihood
integer                       :: outu_        ! output unit
integer, intent(in), optional :: niter        ! # of iterations used by optimizer
integer, intent(in), optional :: info         ! whether optimizer converted
real(kind=dp)    , intent(in), optional :: sigma_est(:) ! conditional standard deviation
real(kind=dp)    , intent(in), optional :: ret(:)       ! returns
integer          , intent(in), optional :: nacfsq       ! # of autocorrelations of squared returns to print
integer          , intent(in), optional :: nacfabs      ! # of autocorrelations of absolute returns to print
character (len=*), intent(in), optional :: fmt_header
character (len=*), parameter  :: fmt_cr = "(a25,':',*(1x,f20.6))", &
   fmt_cc = "(a25,':',*(1x,a))", fmt_par = "(a25,':',*(f10.6))", &
   fmt_ci = "(a25,':',*(1x,i20))", fmt_acf = "(a25,':',*(f7.3))", &
   fmt_labels = "(/,a25,':',*(a10))"
outu_ = output_unit
if (present(fmt_header)) write (outu_, fmt_header)
write (outu_, fmt_cc) "GARCH model", garch_model
write (outu_, fmt_cc) "Conditional distribution", trim(dist)
select case (garch_model)
   case (ew_arch_str)
      write (outu_, fmt_labels) "parameters", "mu", "omega", "alpha"
   case (garch_str)
      write (outu_, fmt_labels) "parameters", "mu", "omega", "alpha", "beta"
   case (gjr_garch_str)
      write (outu_, fmt_labels) "parameters", "mu", "omega", "alpha", "gamma", "beta"
end select
write (outu_, fmt_par) "estimated", est_param
write (outu_, fmt_cr) "Log-likelihood", logL
if (present(niter)) write (outu_, fmt_ci) "#iterations", niter
if (present(info)) call write_merge(info==0, "Convergence achieved.", &
        "Maximum iterations reached without full convergence.")
if (present(sigma_est)) then
   write (outu_, "(/,a25,':',*(a10))") "stat", basic_stats_names
   write (outu_, fmt_par) "sigma_est", basic_stats(sigma_est)
end if
if (present(sigma_est) .and. present(ret)) then
   write (outu_, "(/, a, *(1x,f10.4))") "kurtosis of ret, ret/sigma_est", &
      kurtosis(ret), kurtosis(ret/sigma_est)
   if (nacfsq > 0) then
      write (outu_,fmt_acf) "ACF(ret^2)", acf(ret**2, nacfsq)
      write (outu_,fmt_acf) "ACF((ret/sigma_est)^2)", acf((ret/sigma_est)**2, nacfsq)
   end if
   if (nacfabs > 0) then
      write (outu_,fmt_acf) "ACF(|ret|)", acf(abs(ret), nacfabs)
      write (outu_,fmt_acf) "ACF(|ret/sigma_est|)", acf(abs(ret/sigma_est), nacfabs)
   end if
end if
end subroutine print_garch_results

subroutine print_garch_lik(garch_model, sym, dist, lik)
! Print the log-likelihoods for a GARCH model
character(len=*), intent(in) :: garch_model  ! Name of the GARCH model (e.g., "garch" or "gjr_garch")
character(len=*), intent(in) :: sym(:)   ! stock symbols
character(len=*), intent(in) :: dist(:)  ! distribution names
real(kind=dp)   , intent(in) :: lik(:, :) ! Log-likelihoods (nsym x ndist)
integer :: isym, idist, idist_best
print "(/, a, ' log likelihoods', /, *(a10))", trim(garch_model), "symbol", "best_dist", &
      (trim(dist(idist)), idist=1, size(dist))
do isym = 1, size(lik, 1)
  idist_best = maxloc(lik(isym, :), dim=1)  ! Find the best distribution for this symbol
  print "(2a10, *(f10.1))", trim(sym(isym)), trim(dist(idist_best)), lik(isym, :)
end do
end subroutine print_garch_lik

!------------------------------------------------------------
! Subroutine to fit the equal-weight ARCH model.
!
! Input:
! ret      -- Array of returns.
! dist     -- Distribution flag ("normal" or "laplace").
! max_iter -- (optional) max # of Nelder-Mead iterations
! tol      -- (optional) convergence criterion for value of likelihood
! alpha0 -- (optional) guess for weight on past squared return
!
! Output:
! par_out -- Estimated parameters: [mu, ω, α].
! logL    -- Final (maximized) log likelihood
! info    -- Convergence flag (0 = converged, 1 = max_iter reached).
! niter   -- # of iterations in Nelder-Mead
!------------------------------------------------------------
subroutine fit_ew_arch(ret, dist, par_out, logL, info, niter, max_iter, &
   tol, sigma, alpha0)
character(len=*), intent(in) :: dist
real(kind=dp)   , intent(in) :: ret(:)
real(kind=dp)   , intent(out) :: par_out(:), logL
integer         , intent(out) :: info
integer, intent(out), optional :: niter
integer, intent(in) , optional :: max_iter
real(kind=dp), intent(in), optional :: tol
real(kind=dp), intent(out), optional :: sigma(:)
real(kind=dp), intent(in), optional :: alpha0
real(kind=dp) :: x0(nparam_ew_arch)
real(kind=dp) :: tol_
integer :: n, max_iter_
logical, parameter :: print_guess_ = .false.
character (len=*), parameter :: msg = "in fit_ew_arch, "
garch_model = ew_arch_str
call assert_equal(size(par_out), nparam_ew_arch, &
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
x0(3) = default(0.9_dp, alpha0) ! α
if (print_guess_) print "(/,'in fit_ew_garch, x0:', *(1x,f12.6))", x0
!--- Run the Nelder–Mead optimizer.
call nelder_mead(x0, max_iter_, tol_, par_out, logL, info, &
         niter) ! output: par_out, logL
! Since we minimized the negative log likelihood, the maximized
! log–likelihood is its negative.
logL = -logL
if (present(sigma)) sigma = ew_arch_sigma(ret, par_out)
end subroutine fit_ew_arch

!------------------------------------------------------------
! Subroutine to simulate an equal-weight ARCH process.
!
! Input:
! mu    -- mean return
! omega -- variance constant
! alpha -- symmetric weight on past squared return
!
! Output:
! ret -- Simulated return series (allocated inside).
!------------------------------------------------------------
subroutine simulate_ew_arch(mu, omega, alpha, narch, ret, sigma)
real(kind=dp), intent(in) :: mu, omega, alpha
integer, intent(in) :: narch
real(kind=dp), intent(out) :: ret(:)
real(kind=dp), intent(out), optional :: sigma(:)
real(kind=dp) :: noise, noise_sq(size(ret)), uc_var, &
   sigma2, noise_sq_sum, noise_sq_avg
integer :: i, n
real(kind=dp) :: z, sigma_i
n = size(ret)
if (present(sigma)) call assert_equal(size(sigma), n, &
   "in simulate_ew_arch, size(sigma)")
uc_var = omega / (1.0_dp - narch*alpha)
sigma2 = uc_var
do i = 1, n
   z = random_normal()
   sigma_i = sqrt(sigma2)
   noise = sigma_i * z
   if (present(sigma)) sigma(i) = sigma_i
   ret(i) = mu + noise
   noise_sq(i) = noise**2
   if (i == 1) then
      noise_sq_sum = noise_sq(1)
      noise_sq_avg = noise_sq_sum
   else if (i <= narch) then
      noise_sq_sum = noise_sq_sum + noise_sq(i)
      noise_sq_avg = noise_sq_sum/i
   else
      noise_sq_sum = noise_sq_sum + noise_sq(i) - noise_sq(i-narch)
      noise_sq_avg = noise_sq_sum/narch
   end if
   sigma2 = min(i,narch)*alpha*noise_sq_avg + omega
end do
end subroutine simulate_ew_arch

!------------------------------------------------------------
! Subroutine to fit the equal-weight ARCH model.
!
! Input:
! ret      -- Array of returns.
! dist     -- Distribution flag ("normal" or "laplace").
! max_iter -- (optional) max # of Nelder-Mead iterations
! tol      -- (optional) convergence criterion for value of likelihood
!
! Output:
! par_out -- Estimated parameters: [mu, ω, α].
! logL    -- Final (maximized) log likelihood
! info    -- Convergence flag (0 = converged, 1 = max_iter reached).
! niter   -- # of iterations in Nelder-Mead
!------------------------------------------------------------
subroutine fit_arch(ret, dist, par_out, logL, info, niter, max_iter, &
   tol, sigma, alpha0)
character(len=*), intent(in) :: dist
real(kind=dp)   , intent(in) :: ret(:)
real(kind=dp)   , intent(out) :: par_out(:), logL
integer         , intent(out) :: info
integer, intent(out), optional :: niter
integer, intent(in) , optional :: max_iter
real(kind=dp), intent(in), optional :: tol, alpha0
real(kind=dp), intent(out), optional :: sigma(:)
real(kind=dp) :: x0(size(par_out))
real(kind=dp) :: tol_
integer :: n, max_iter_
logical, parameter :: print_guess_ = .false.
character (len=*), parameter :: msg = "in fit_arch, "
garch_model = arch_str
max_iter_ = default(1000, max_iter)
n = size(ret)
if (present(sigma)) call assert_equal(size(sigma), n, msg // "size(sigma)")
!--- Set the module variables for the objective function.
xret = ret
xdist = dist
tol_ = default(1.0e-10_dp, tol)

!--- Set up an initial guess.
x0(1) = sum(ret) / n ! mu
x0(2) = variance(ret) ! ω
x0(3:) = default(0.0_dp, alpha0) ! 0.1_dp ! 0.0_dp
if (print_guess_) print "(/,'in fit_arch, x0:', *(1x,f12.6))", x0
!--- Run the Nelder–Mead optimizer.
call nelder_mead(x0, max_iter_, tol_, par_out, logL, info, &
         niter) ! output: par_out, logL
! Since we minimized the negative log likelihood, the maximized
! log–likelihood is its negative.
logL = -logL
if (present(sigma)) sigma = arch_sigma(ret, par_out)
end subroutine fit_arch

!------------------------------------------------------------
! Subroutine to simulate a general ARCH process.
!
! Input:
! mu    -- mean return
! omega -- variance constant
! alpha -- array of ARCH coefficients (alpha(1) for lag 1, etc.)
!
! Output:
! ret   -- Simulated return series (allocated inside).
! sigma -- (optional) Simulated conditional standard deviation.
!------------------------------------------------------------
subroutine simulate_arch(mu, omega, alpha, ret, sigma)
real(kind=dp), intent(in) :: mu, omega
real(kind=dp), intent(in) :: alpha(:)   ! ARCH coefficients (alpha(1) for lag 1, etc.)
real(kind=dp), intent(out) :: ret(:)
real(kind=dp), intent(out), optional :: sigma(:)
integer :: n, narch, i, j
real(kind=dp) :: sigma2, sigma_i, z, noise, uc_var
real(kind=dp), allocatable :: noise_sq(:)
n = size(ret)
narch = size(alpha)
if (present(sigma)) &
   call assert_equal(size(sigma), n, "simulate_arch: size(sigma) mismatch")
allocate(noise_sq(n))
! Use the unconditional variance as the starting value.
uc_var = omega / (1.0_dp - sum(alpha))
sigma2 = uc_var
do i = 1, n
   sigma_i = sqrt(sigma2)
   z = random_normal()
   noise = sigma_i * z
   ret(i) = mu + noise
   if (present(sigma)) sigma(i) = sigma_i
   noise_sq(i) = noise**2
   ! Update sigma2 for the next time step
   sigma2 = omega
   if (i < narch) then
      do j = 1, i
         sigma2 = sigma2 + alpha(j) * noise_sq(i + 1 - j)
      end do
   else
      do j = 1, narch
         sigma2 = sigma2 + alpha(j) * noise_sq(i + 1 - j)
      end do
   end if
end do
end subroutine simulate_arch

end module garch_mod
