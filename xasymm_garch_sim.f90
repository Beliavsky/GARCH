program xasymm_garch_sim
! simulates from asymmetric GARCH model and prints stats on returns
! and conditional volatilities
use kind_mod           , only: dp
use asymm_garch_sim_mod, only: simulate_shift_twist_garch
use random_mod         , only: random_normal
use basic_stats_mod    , only: print_basic_stats
implicit none
integer :: iter
integer, parameter :: n = 10**6, niter = 20
real(kind=dp) :: z(n), ret(n), sigma(n), mu, omega, alpha, beta, gamma, &
   c_shift
mu = 0.0_dp
omega = 0.1_dp
beta = 0.7_dp
alpha = 0.1_dp
gamma = 0.3_dp
c_shift = 0.0_dp
print "(*(a10))",  "#obs", "mu", "omega", "beta", "alpha", "gamma", "c_shift"
print "(i10,*(f10.4))", n,  mu ,  omega ,  beta ,  alpha ,  gamma ,  c_shift 
print*
do iter=1,niter
   z = random_normal(n)
   call simulate_shift_twist_garch(z, ret, mu, omega, alpha, gamma, beta, &
      c_shift, sigma)
   call print_basic_stats(reshape([z, ret, sigma], [n, 3]), &
      labels=[character (len=5) :: "z", "ret", "sigma"], &
      fmt_cr="(a10,*(f10.4))")
   print "('first, last sigma:',2(1x,f10.4), /)", sigma([1, n])
end do
end program xasymm_garch_sim