program xew_arch_sim
! simulate an equal-weight ARCH process and fit a GARCH model, an
! equal-weight ARCH model, and an unconstrained ARCH model to it
use kind_mod, only: dp
use garch_mod, only: simulate_garch, nparam_garch, fit_garch, &
   simulate_ew_arch, fit_ew_arch, nparam_ew_arch, simulate_arch, &
   fit_arch
use obj_fun_mod, only: narch
use random_mod, only: random_seed_init
use util_mod, only: write_merge
use basic_stats_mod, only: mean, sd, kurtosis, basic_stats, &
   basic_stats_names, acf, print_basic_stats, print_corr_mat
implicit none
! Define format strings for integer, real and character output.
character(len=*), parameter :: fmt_ci = "(a30,':',*(1x,i20))", &
   fmt_cr = "(a30,':',*(1x,f20.6))", fmt_cc = "(a30,':',*(1x,a))",&
   fmt_par = "(a30,':',*(f10.6))"
integer, parameter :: nacf=50, iseed = 456, n = 10**4 ! # of observations
integer :: info, niter_nm_garch, niter_nm_arch
real(kind=dp) :: ret(n), sigma(n), sigma_est_garch(n), sigma_est_arch(n), &
   sigma_est_ew_arch(n)
! True parameters for simulation declared as named constants.
real(kind=dp), parameter :: true_mu = 0.0_dp, true_omega = 0.4_dp, &
   true_alpha = 0.19_dp, alpha0 = 0.0_dp
real(kind=dp) :: logl_garch, logl_arch ! log-likelihood of garch and arch models
real(kind=dp) :: par_out_garch(nparam_garch), par_out_ew_arch(nparam_ew_arch) ! estimated parameters
real(kind=dp), allocatable :: par_out_arch(:)
! Conditional return distribution
character(len=*), parameter :: dist = "normal"
narch = 5
allocate (par_out_arch(narch+2))
print fmt_ci, "iseed", iseed
print fmt_ci, "#obs", n
print fmt_ci, "#lags ARCH", narch
if (iseed /= 0) call random_seed_init(iseed, print_seed=.true.)
! garch_model = garch_str
! Simulate the equal-weight ARCH process
call simulate_ew_arch(true_mu, true_omega, true_alpha, narch, ret, sigma=sigma)
! Fit the GARCH model to the simulated returns.
call fit_garch(ret, dist, par_out_garch, logl_garch, info, niter=niter_nm_garch, &
   sigma=sigma_est_garch)
call fit_ew_arch(ret, dist, par_out_ew_arch, logl_arch, info, niter=niter_nm_arch, &
   sigma=sigma_est_ew_arch)
! write (*, fmt_cc) "GARCH model", garch_model
write (*, fmt_cc) "Conditional distribution", trim(dist)
write (*, fmt_ci) "Number of observations", n
write (*, "(/,a30,':',*(a10))") "parameters", "mu", "omega", "alpha", "beta"
write (*, fmt_par) "true", true_mu, true_omega, true_alpha
write (*, fmt_par) "estimated", par_out_garch
write (*, fmt_cr) "Log-likelihood garch", logl_garch
write (*, fmt_cr) "Log-likelihood ew_arch", logl_arch
write (*, fmt_ci) "#Nelder-Mead iterations", niter_nm_garch
call write_merge(info==0, "Convergence achieved.", &
     "Maximum iterations reached without full convergence.")
write (*,"(/,a30,':',*(a10))") "stat", basic_stats_names
write (*,fmt_par) "sigma", basic_stats(sigma)
write (*,fmt_par) "sigma_est_garch", basic_stats(sigma_est_garch)
write (*,fmt_par) "sigma_est_ew_arch", basic_stats(sigma_est_ew_arch)
write (*, "(/, a, *(1x,f10.4))") "kurtosis of ret, ret/sigma, ret/sigma_est_garch, ret/sigma_est_ew_arch", &
   kurtosis(ret), kurtosis(ret/sigma), kurtosis(ret/sigma_est_garch), kurtosis(ret/sigma_est_ew_arch)
write (*,fmt_par) "ACF(ret^2)", acf(ret**2, nacf)
write (*,fmt_par) "ACF((ret/sigma_est_garch)^2)", acf((ret/sigma_est_garch)**2, nacf)
write (*,fmt_par) "ACF((ret/sigma_est_ew_arch)^2)", acf((ret/sigma_est_ew_arch)**2, nacf)
write (*,fmt_par) "true EW ARCH", true_mu, true_omega, true_alpha
write (*,fmt_par) "estimated EW ARCH", par_out_ew_arch
print*,"calling fit_arch"
call fit_arch(ret, dist, par_out_arch, logl_arch, info, niter=niter_nm_arch, &
   sigma=sigma_est_arch, max_iter=10**5, alpha0=alpha0)
write (*,fmt_par) "estimated ARCH", par_out_arch
write (*, fmt_cr) "Log-likelihood", logl_arch
call print_corr_mat(reshape([sigma, sigma_est_garch, sigma_est_ew_arch, &
   sigma_est_arch], [n, 4]), col_names=[character(len=10) :: &
   "true", "garch", "ew_arch", "arch"], fmt_header="(/,'correlations of sigma')", &
   fmt_trailer="()")
print*,"(5) finished xew_arch_sim"
end program xew_arch_sim
