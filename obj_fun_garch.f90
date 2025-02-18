module obj_fun_mod
use        kind_mod, only: dp
use   constants_mod, only: pi, log_two, log_two_pi, sqrt_two
use basic_stats_mod, only: variance, moving_average
use     density_mod, only: minus_log_density
implicit none
private
public :: obj_fun, xret, xdist, garch_model, neg_loglik_garch, &
   neg_loglik_gjr_garch, gjr_garch_sigma, garch_sigma, narch, &
   nparam_gjr_garch, nparam_garch, neg_loglik_ew_arch, &
   garch_str, gjr_garch_str, nparam_ew_arch, ew_arch_str, &
   ew_arch_sigma, neg_loglik_arch, arch_sigma, arch_str
integer, parameter :: nparam_gjr_garch = 5, nparam_garch=4, &
   nparam_ew_arch=3
character (len=*), parameter :: garch_str="garch", gjr_garch_str="gjr_garch", &
   ew_arch_str="ew_arch", arch_str="arch"
logical, parameter :: call_log_density = .true.
real(kind=dp), allocatable :: xret(:)  ! returns
character (len=10)         :: xdist="normal"    ! conditional distribution ("normal" or something else)
character (len=10)         :: garch_model
real(kind=dp), parameter   :: bad_nll = 1.0e10_dp
integer                    :: narch
contains

function obj_fun(x) result(y)
! Selects the appropriate negative log likelihood
! function depending on the chosen model.
real(kind=dp), intent(in) :: x(:)
real(kind=dp)             :: y

select case (garch_model)
   case (gjr_garch_str) ; y = neg_loglik_gjr_garch(x)
   case (garch_str)     ; y = neg_loglik_garch(x)
   case (ew_arch_str)   ; y = neg_loglik_ew_arch(x)
   case (arch_str)      ; y = neg_loglik_arch(x)
   case default         ; error stop "in obj_fun, invalid garch_model '" &
                              // trim(garch_model) // "'"
end select
end function obj_fun

function neg_loglik_garch(par) result(nll)
!------------------------------------------------------------
! Function to compute the negative log likelihood for the
! symmetric GARCH(1,1) model.
!
! Model:
!    r_t = mu + ε_t,
!    σ_t² = ω + α (r_{t-1} - mu)² + β σ_{t-1}².
!
! The density is chosen via the character variable XDIST ("normal" or some other distribution).
! A penalty (returning a huge value) is applied if any constraints are violated:
!    ω > 0,  α ≥ 0,  β ≥ 0, and (α + β) < 1.
!------------------------------------------------------------
real(kind=dp), intent(in) :: par(nparam_garch)
real(kind=dp)             :: nll
real(kind=dp)             :: mu, omega, alpha, beta, r
real(kind=dp), allocatable :: sig2(:)
integer                   :: n, t

n = size(xret)

! Unpack parameters.
mu    = par(1)
omega = par(2)
alpha = par(3)
beta  = par(4)

! Impose parameter constraints.
if (omega <= 0.0_dp .or. alpha < 0.0_dp .or. beta < 0.0_dp .or. (alpha + beta) >= 1.0_dp) then
   nll = bad_nll
   return
end if

allocate(sig2(n))
! Use the sample variance as the starting value for σ².
sig2(1) = variance(xret)
if (sig2(1) <= 0.0_dp) sig2(1) = 1.0e-6_dp

nll = 0.0_dp
do t = 1, n
   if (t > 1) then
      r = xret(t-1) - mu
      sig2(t) = omega + alpha*r**2 + beta*sig2(t-1)
      if (sig2(t) <= 0.0_dp) then
         nll = bad_nll
         return
      end if
   end if
   nll = nll + minus_log_density(xret(t) - mu, xdist, sqrt(sig2(t)))
end do
end function neg_loglik_garch

function neg_loglik_gjr_garch(par) result(nll)
!------------------------------------------------------------
! Function to compute the negative log likelihood for the
! GJR-GARCH(1,1) model.
!
! Model:
!    r_t = mu + ε_t,
!    σ_t² = ω + α ε_{t-1}² + γ I(ε_{t-1}<0) ε_{t-1}² + β σ_{t-1}².
!
! The density is chosen via the character variable XDIST ("normal" or some other distribution).
! A penalty (returning a huge value) is applied if any constraints are violated:
!    ω > 0,  α ≥ 0,  γ ≥ 0,  β ≥ 0, and (α + 0.5*γ + β) < 1.
!------------------------------------------------------------
real(kind=dp), intent(in) :: par(nparam_gjr_garch)
real(kind=dp)             :: nll
real(kind=dp)             :: mu, omega, alpha, gamma, beta, r
real(kind=dp), allocatable :: sig2(:)
integer                   :: n, t

n = size(xret)

! Unpack parameters.
mu    = par(1)
omega = par(2)
alpha = par(3)
gamma = par(4)
beta  = par(5)

! Impose parameter constraints.
if (omega <= 0.0_dp .or. alpha < 0.0_dp .or. gamma < 0.0_dp .or. beta < 0.0_dp .or. &
(alpha + 0.5_dp*gamma + beta) >= 1.0_dp) then
nll = bad_nll
return
end if

allocate(sig2(n))
! Use the sample variance as the starting value for σ².
sig2(1) = variance(xret)
if (sig2(1) <= 0.0_dp) sig2(1) = 1.0e-6_dp

nll = 0.0_dp
do t = 1, n
   if (t > 1) then
      r = xret(t-1) - mu
      sig2(t) = omega + alpha*r**2 + gamma*merge(r**2, 0.0_dp, (r < 0.0_dp)) + beta*sig2(t-1)
      if (sig2(t) <= 0.0_dp) then
         nll = bad_nll
         return
      end if
   end if
   nll = nll + minus_log_density(xret(t) - mu, xdist, sqrt(sig2(t)))
end do
end function neg_loglik_gjr_garch

function gjr_garch_sigma(xxret, par) result(sigma)
!> Computes the conditional standard deviation (σ_t) for the GJR-GARCH(1,1) model.
!! Given the parameter vector `par` and a series of returns `xxret`, it iteratively computes
!! the time-varying variance using:
!!    σ_t² = ω + α ε_{t-1}² + γ I(ε_{t-1} < 0) ε_{t-1}² + β σ_{t-1}²,
!! and returns the square root of the variance (σ_t) at each time step.
!! If any variance value becomes negative, the corresponding σ_t is set to -1.0.
real(kind=dp), intent(in) :: par(nparam_gjr_garch)
real(kind=dp), intent(in) :: xxret(:)
real(kind=dp)             :: sigma(size(xxret))
real(kind=dp)             :: mu, omega, alpha, gamma, beta, r
real(kind=dp), allocatable :: sig2(:)
integer                   :: n, t
n = size(xxret)
! Unpack parameters.
mu    = par(1)
omega = par(2)
alpha = par(3)
gamma = par(4)
beta  = par(5)
allocate(sig2(n))
! Use the sample variance as the starting value for σ².
sig2(1) = variance(xxret)
if (sig2(1) <= 0.0_dp) sig2(1) = 1.0e-6_dp
do t = 1, n
   if (t > 1) then
      r = xxret(t-1) - mu
      sig2(t) = omega + alpha*r**2 + gamma*merge(r**2, 0.0_dp, (r < 0.0_dp)) + beta*sig2(t-1)
      if (sig2(t) < 0.0_dp) then
         sigma(t) = -1.0_dp
      else
         sigma(t) = sqrt(sig2(t))
      end if
   end if
   r = xxret(t) - mu
end do
where (sig2 >= 0.0_dp)
   sigma = sqrt(sig2)
elsewhere
   sigma = -1.0_dp
end where
end function gjr_garch_sigma

function garch_sigma(xxret, par) result(sigma)
!> Computes the conditional standard deviation (σ_t) for the GARCH(1,1) model.
!! Given the parameter vector `par` and a series of returns `xxret`, it iteratively computes
!! the time-varying variance using:
!!    σ_t² = ω + α ε_{t-1}² +  β σ_{t-1}²,
!! and returns the square root of the variance (σ_t) at each time step.
!! If any variance value becomes negative, the corresponding σ_t is set to -1.0.
real(kind=dp), intent(in) :: par(nparam_garch)
real(kind=dp), intent(in) :: xxret(:)
real(kind=dp)             :: sigma(size(xxret))
real(kind=dp)             :: mu, omega, alpha, beta, r
real(kind=dp), allocatable :: sig2(:)
integer                   :: n, t
n = size(xxret)
! Unpack parameters.
mu    = par(1)
omega = par(2)
alpha = par(3)
beta  = par(4)
allocate(sig2(n))
! Use the sample variance as the starting value for σ².
sig2(1) = variance(xxret)
if (sig2(1) <= 0.0_dp) sig2(1) = 1.0e-6_dp
do t = 1, n
   if (t > 1) then
      r = xxret(t-1) - mu
      sig2(t) = omega + alpha*r**2 + beta*sig2(t-1)
      if (sig2(t) < 0.0_dp) then
         sigma(t) = -1.0_dp
      else
         sigma(t) = sqrt(sig2(t))
      end if
   end if
   r = xxret(t) - mu
end do
where (sig2 >= 0.0_dp)
   sigma = sqrt(sig2)
elsewhere
   sigma = -1.0_dp
end where
end function garch_sigma

function neg_loglik_ew_arch(par) result(nll)
!------------------------------------------------------------
! Function to compute the negative log likelihood for the
! symmetric equal-weight ARCH(narch) model.
!
! The density is chosen via the character variable XDIST
! A penalty (returning a huge value) is applied if any constraints are violated:
!    ω > 0,  α ≥ 0
!------------------------------------------------------------
real(kind=dp), intent(in)  :: par(nparam_ew_arch)
real(kind=dp)              :: nll
real(kind=dp)              :: mu, omega, alpha
real(kind=dp)              :: sig2
integer                    :: n, t
real(kind=dp), allocatable :: var_ret(:)
logical, parameter :: debug = .false.
n = size(xret)
! Unpack parameters.
mu    = par(1)
omega = par(2)
alpha = par(3)
var_ret = moving_average((xret-mu)**2, narch)
! Impose parameter constraints.
if (narch < 0 .or. omega <= 0.0_dp .or. alpha < 0.0_dp) then
   nll = bad_nll
   return
end if
nll = 0.0_dp
if (debug) then
   print*,"in neg_loglik_ew_arch, narch, var_ret =", narch, var_ret ! debug
   print*,"size(xret), xret =", size(xret), xret ! debug
end if
do t = max(1,narch)+1, n
   sig2 = omega + narch*alpha*var_ret(t-1)
   nll = nll + minus_log_density(xret(t) - mu, xdist, sqrt(sig2))
   if (debug) print*,"t, sig2, nll =", t, sig2, nll ! debug
   if (sig2 <= 0.0_dp) error stop "sig2 must be positive"
end do
end function neg_loglik_ew_arch

function ew_arch_sigma(xxret, par) result(sigma)
!> Computes the conditional standard deviation (σ_t) for the equal-weight ARCH model.
!! Given the parameter vector `par` and a series of returns `xxret`, it computes
!! the time-varying variance
!! and returns the square root of the variance (σ_t) at each time step.
!! If any variance value becomes negative, the corresponding σ_t is set to -1.0.
real(kind=dp), intent(in)  :: par(nparam_ew_arch)
real(kind=dp), intent(in)  :: xxret(:)
real(kind=dp)              :: sigma(size(xxret))
real(kind=dp)              :: mu, omega, alpha
real(kind=dp), allocatable :: sig2(:)
integer                    :: n, t
real(kind=dp), allocatable :: var_ret(:)
n = size(xxret)
! Unpack parameters.
mu    = par(1)
omega = par(2)
alpha = par(3)
allocate(sig2(n))
var_ret = moving_average((xret-mu)**2, narch)
! Use the sample variance as the starting value for σ².
sig2(1) = variance(xxret)
do t = 2, n
   sig2(t) = omega + narch*alpha*var_ret(t-1)
end do
where (sig2 >= 0.0_dp)
   sigma = sqrt(sig2)
elsewhere
   sigma = -1.0_dp
end where
end function ew_arch_sigma

function neg_loglik_arch(par) result(nll)
!------------------------------------------------------------
! Function to compute the negative log likelihood for the
! symmetric ARCH(narch) model.
!
! The density is chosen via the character variable XDIST
! A penalty (returning a huge value) is applied if any constraints are violated:
!    ω > 0,  α ≥ 0
!------------------------------------------------------------
real(kind=dp), intent(in)  :: par(:)
real(kind=dp)              :: nll
real(kind=dp)              :: mu, omega
real(kind=dp)              :: sig2
integer                    :: n, t, nlags
real(kind=dp), allocatable :: alpha(:), ret_sq(:)
logical, parameter :: debug = .false.
n = size(xret)
! Unpack parameters.
mu    = par(1)
omega = par(2)
alpha = par(3:)
nlags = size(alpha)
ret_sq = (xret-mu)**2
! Impose parameter constraints.
if (omega <= 0.0_dp .or. any(alpha < 0.0_dp)) then
   nll = bad_nll
   return
end if
nll = 0.0_dp
if (debug) then
   print*,"in neg_loglik_arch, par =", par
   print*,"alpha =", alpha ! debug
end if
do t = nlags+1, n
   sig2 = omega + sum(alpha(nlags:1:-1)*ret_sq(t-nlags:t-1))
   nll = nll + minus_log_density(xret(t) - mu, xdist, sqrt(sig2))
   if (debug) print*,"t, sig2, nll =", t, sig2, nll ! debug
   if (sig2 <= 0.0_dp) error stop "sig2 must be positive"
end do
end function neg_loglik_arch

function arch_sigma(xxret, par) result(sigma)
!> Computes the conditional standard deviation (σ_t) for the ARCH model.
!! Given the parameter vector `par` and a series of returns `xxret`, it computes
!! the time-varying variance
!! and returns the square root of the variance (σ_t) at each time step.
!! If any variance value becomes negative, the corresponding σ_t is set to -1.0.
! real(kind=dp), intent(in)  :: par(nparam_ew_arch)
real(kind=dp), intent(in)  :: par(:)
real(kind=dp), intent(in)  :: xxret(:)
real(kind=dp)              :: sigma(size(xxret))
real(kind=dp)              :: mu, omega
real(kind=dp), allocatable :: sig2(:), alpha(:)
integer                    :: n, t, nlags
real(kind=dp), allocatable :: ret_sq(:)
n = size(xxret)
! Unpack parameters.
mu    = par(1)
omega = par(2)
alpha = par(3:)
nlags = size(alpha)
allocate(sig2(n))
ret_sq = (xret-mu)**2
! Use the sample variance as the starting value for σ².
sig2(1:nlags) = variance(xxret)
do t = nlags+1, n
   sig2(t) = omega + sum(alpha(nlags:1:-1)*ret_sq(t-nlags:t-1))
end do
where (sig2 >= 0.0_dp)
   sigma = sqrt(sig2)
elsewhere
   sigma = -1.0_dp
end where
end function arch_sigma

end module obj_fun_mod
