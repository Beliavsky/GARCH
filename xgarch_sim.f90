program xgarch_sim
! simulate a GARCH process and fit the GARCH model to it
use kind_mod, only: dp
use garch_mod, only: simulate_garch, nparam_garch, garch_str, fit_garch
use obj_fun_mod, only: garch_model, neg_loglik_garch
use random_mod, only: random_seed_init
use util_mod, only: write_merge
use basic_stats_mod, only: mean, sd, kurtosis, basic_stats, &
   basic_stats_names, correl, acf
implicit none
! Define format strings for integer, real and character output.
character(len=*), parameter :: fmt_ci = "(a25,':',*(1x,i20))", &
   fmt_cr = "(a25,':',*(1x,f20.6))", fmt_cc = "(a25,':',*(1x,a))",&
   fmt_par = "(a25,':',*(f10.6))"
integer, parameter :: nacf=10, iseed = 123, n = 10**4 ! # of observations
integer :: info, niter_nm
real(kind=dp) :: ret(n), sigma(n), sigma_est(n)
! True parameters for simulation declared as named constants.
real(kind=dp), parameter :: true_mu = 0.0_dp, true_omega = 0.1_dp, &
   true_alpha = 0.15_dp, true_beta  = 0.8_dp
real(kind=dp) :: logL ! log-likelihood
real(kind=dp) :: par_out(nparam_garch) ! estimated parameters
! Conditional return distribution
character(len=*), parameter :: dist = "normal" 
print fmt_ci, "iseed", iseed
print fmt_ci, "#obs", n
if (iseed /= 0) call random_seed_init(iseed, print_seed=.true.)
garch_model = garch_str
! Simulate the GARCH process.
call simulate_garch(true_mu, true_omega, true_alpha,true_beta, ret, sigma=sigma)
! Fit the GARCH model to the simulated returns.
call fit_garch(ret, dist, par_out, logL, info, niter=niter_nm, &
   sigma=sigma_est)
! Print output using the defined format strings.
write (*, fmt_cc) "GARCH model", garch_model
write (*, fmt_cc) "Conditional distribution", trim(dist)
write (*, fmt_ci) "Number of observations", n
write (*, "(/,a25,':',*(a10))") "parameters", "mu", "omega", "alpha", "beta"
write (*, fmt_par) "true", true_mu, true_omega, true_alpha, true_beta
write (*, fmt_par) "estimated", par_out
write (*, fmt_cr) "Log-likelihood", logL
write (*, fmt_cr) "loglik_garch", -neg_loglik_garch(par_out)
write (*, fmt_ci) "#Nelder-Mead iterations", niter_nm
call write_merge(info==0, "Convergence achieved.", &
     "Maximum iterations reached without full convergence.")
write (*,"(/,a25,':',*(a10))") "stat", basic_stats_names
write (*,fmt_par) "sigma", basic_stats(sigma)
write (*,fmt_par) "sigma_est", basic_stats(sigma_est)
write (*, "(/, a, *(1x,f10.4))") "kurtosis of ret, ret/sigma, ret/sigma_est", &
   kurtosis(ret), kurtosis(ret/sigma), kurtosis(ret/sigma_est)
write (*,fmt_cr) "corr(sigma, sigma_est)", correl(sigma, sigma_est)
write (*,fmt_par) "ACF(ret^2)", acf(ret**2, nacf)
write (*,fmt_par) "ACF((ret/sigma_est)^2)", acf((ret/sigma_est)**2, nacf)
print*,"(2) finished xgarch_sim"
end program xgarch_sim
