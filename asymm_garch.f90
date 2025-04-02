module asymm_garch_mod
use        kind_mod, only: dp
! use optim_methods_mod, only: nelder_mead_str, uobyqa_str
! use nelder_mead_mod, only: nelder_mead
! use     obj_fun_mod, only: xret, xdist, gjr_garch_sigma, obj_sub
! use basic_stats_mod, only: variance
! use      random_mod, only: random_normal
! use   constants_mod, only: pi
! use        util_mod, only: assert_equal, default
! use      uobyqa_mod, only: uobyqa
implicit none
private
public :: simulate_shift_twist_garch
integer, parameter :: nparam_gjr_garch = 5
contains

subroutine simulate_shift_twist_garch(z, ret, mu, omega, alpha, &
   gamma, beta, c_shift, sigma)
! Subroutine to simulate a GJR-GARCH(1,1) process.
!
! Input:
! z     -- standardized noise 
! mu    -- mean return
! omega -- variance constant
! alpha -- symmetric weight on past squared return
! gamma -- weight on past squared return when return is negative
! beta  -- weight on previous variance
! c_shift -- shift of standardized noise in variance equation
!
! Output:
! ret -- Simulated return series
! sigma (optional) -- conditional standard deviation
real(kind=dp), intent(in)  :: mu, omega, alpha, gamma, beta, c_shift
real(kind=dp), intent(in)  :: z(:) ! noise with 0 mean and unit variance
real(kind=dp), intent(out) :: ret(:)
real(kind=dp), intent(out), optional :: sigma(:)
real(kind=dp) :: uc_var, sigma2, noise, noise_shift
integer :: i, n
real(kind=dp) :: sigma_i
n = size(ret)
if (present(sigma)) call assert_equal(size(sigma), n, &
   "in simulate_shift_twist_garch, size(sigma)")
uc_var = omega / (1.0_dp - (1.0_dp + c_shift**2) * (alpha + 0.5_dp*gamma) - &
   beta)
sigma2 = uc_var
do i = 1, n
   sigma_i = sqrt(sigma2)
   noise = sigma_i * z(i)
   noise_shift = sigma_i * (z(i) - c_shift)
   if (present(sigma)) sigma(i) = sigma_i
   ret(i) = mu + noise
   sigma2 = (alpha + merge(gamma,0.0_dp,noise_shift < 0.0_dp)) * &
              noise_shift**2 + beta*sigma2 + omega
end do
end subroutine simulate_shift_twist_garch
end module asymm_garch_mod
