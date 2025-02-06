module basic_stats_mod
use iso_fortran_env, only: output_unit
use kind_mod, only: dp
use util_mod, only: default
implicit none
private
public :: mean, variance, sd, mean_and_sd, kurtosis, basic_stats, &
   print_basic_stats, basic_stats_names, correl, acf
integer, parameter :: nbasic_stats = 6
character (len=*), parameter :: basic_stats_names(nbasic_stats) = &
   [character(len=4) :: "mean", "sd", "skew", "kurt", "min", "max"]
contains

function mean(x) result(xmean)
! return the mean of x(:)
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: xmean
xmean = sum(x)/max(1,size(x))
end function mean

function sd(x) result(xsd)
! return the standard deviation of x(:)
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: xsd
real(kind=dp) :: m, var
integer :: n
n = size(x)
m = sum(x) / n
var = sum((x - m)**2) / (n-1)
xsd = sqrt(max(0.0_dp, var))
end function sd

function mean_and_sd(x) result(res)
! return the mean and standard deviation of x(:)
real(kind=dp), intent(in) :: x(:)
real(kind=dp)             :: res(2)
real(kind=dp)             :: var
integer :: n
n = size(x)
res(1) = sum(x) / n
var = sum((x - res(1))**2) / (n-1)
res(2) = sqrt(max(0.0_dp, var))
end function mean_and_sd

function variance(x) result(var)
! return the variance of x(:)
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: var, m
integer :: n
n = size(x)
m = sum(x) / n
var = sum((x - m)**2) / (n-1)
end function variance

function skew(x) result(skew_val)
! return the skewness of x
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: skew_val
real(kind=dp) :: mean_x, sd_x
integer :: n
n = size(x)
mean_x = mean(x)
sd_x = sd(x)
skew_val = sum(((x - mean_x) / sd_x)**3) / n
end function skew

function kurtosis(x) result(kurtosis_val)
! return the kurtosis of x
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: kurtosis_val
real(kind=dp) :: mean_x, sd_x
integer :: n
n = size(x)
mean_x = mean(x)
sd_x = sd(x)
kurtosis_val = sum(((x - mean_x) / sd_x)**4) / n - 3.0_dp
end function kurtosis

function basic_stats(x) result(stats)
real(kind=dp), intent(in) :: x(:)
real(kind=dp)             :: stats(nbasic_stats)
stats = [mean(x), sd(x), skew(x), kurtosis(x), minval(x), maxval(x)]
end function basic_stats

subroutine print_basic_stats(x, outu)
real(kind=dp), intent(in) :: x(:)
integer, intent(in), optional :: outu
integer :: outu_
outu_ = default(output_unit, outu)
write (outu_, "(*(f10.6))") basic_stats(x)
end subroutine print_basic_stats

pure function correl(x, y) result(corr_xy)
! Returns the linear Pearson correlation of x(:) and y(:)
! Returns a correlation < -1.0_dp to signal an error
real(kind=dp), intent(in) :: x(:), y(:)
real(kind=dp) :: corr_xy
real(kind=dp) :: x_mean, y_mean, cov_xy, var_x, var_y
integer :: n
n = size(x)
if (n /= size(y) .or. n == 0) then
   corr_xy = -2.0_dp
   return
end if
x_mean = sum(x) / n
y_mean = sum(y) / n
cov_xy = sum((x - x_mean) * (y - y_mean))
var_x  = sum((x - x_mean)**2)
var_y  = sum((y - y_mean)**2)
if (var_x <= 0.0_dp .or. var_y <= 0.0_dp) then
   corr_xy = -3.0_dp
else
   corr_xy = cov_xy / sqrt(var_x * var_y)
end if
end function correl

function acf(x, nacf) result(xacf)
! return the autocorrelations at lags 1 through nacf
real(kind=dp), intent(in) :: x(:)         ! Input array
integer, intent(in) :: nacf               ! Number of autocorrelations to compute
real(kind=dp) :: xacf(nacf)               ! Output array for autocorrelations
real(kind=dp) :: denom
real(kind=dp), allocatable :: xdm(:)      ! Demeaned version of x
integer :: n, lag
n = size(x)
xdm = x - mean(x)                          ! Compute demeaned x
denom = sum(xdm**2)
! Compute autocorrelation for each lag from 1 to nacf
do lag = 1, nacf
   xacf(lag) = sum(xdm(1:n-lag) * xdm(lag+1:n)) / denom
end do
end function acf
end module basic_stats_mod
