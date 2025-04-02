module asymm_garch_sim_mod
use        kind_mod, only: dp
use        util_mod, only: assert_equal, default
implicit none
private
public :: simulate_shift_twist_garch
integer, parameter :: nparam_gjr_garch = 5
contains

subroutine simulate_shift_twist_garch(z, ret, mu, omega, alpha, &
   gamma, beta, c_shift, sigma, iprint)
! simulate an asymmetric GARCH(1,1) process with shift and twist parameters
!
! Input:
! z     -- standardized noise 
! mu    -- mean return
! omega -- variance constant
! alpha -- symmetric weight on past squared return
! gamma -- weight on past squared return when return is negative
! beta  -- weight on previous variance
! c_shift -- shift of standardized noise in variance equation
! iprint
!
! Output:
! ret -- Simulated return series
! sigma (optional) -- conditional standard deviation
real(kind=dp), intent(in)  :: mu, omega, alpha, gamma, beta, c_shift
real(kind=dp), intent(in)  :: z(:) ! noise with 0 mean and unit variance
real(kind=dp), intent(out) :: ret(:)
real(kind=dp), intent(out), optional :: sigma(:)
integer, intent(in), optional :: iprint
real(kind=dp) :: uc_var, sigma2, noise, noise_shift
integer :: i, iprint_, n
real(kind=dp) :: sigma_i
iprint_ = default(0, iprint)
n = size(ret)
if (present(sigma)) call assert_equal(size(sigma), n, &
   "in simulate_shift_twist_garch, size(sigma)")
uc_var = omega / (1.0_dp - (1.0_dp + c_shift**2) * (alpha + 0.5_dp*gamma) - &
   beta)
if (iprint_ > 0) print*,"uc_var:", uc_var
sigma2 = uc_var
do i = 1, n
   sigma_i = sqrt(sigma2)
   noise = sigma_i * z(i)
   noise_shift = sigma_i * (z(i) - c_shift)
   if (present(sigma)) sigma(i) = sigma_i
   ret(i) = mu + noise
   if (iprint_ > 1) print "(i8, *(f10.6))", i,sigma_i, z(i), noise, ret(i)
   sigma2 = (alpha + merge(gamma,0.0_dp,noise_shift < 0.0_dp)) * &
              noise_shift**2 + beta*sigma2 + omega
end do
end subroutine simulate_shift_twist_garch
end module asymm_garch_sim_mod
