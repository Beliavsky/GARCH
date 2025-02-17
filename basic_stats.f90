module basic_stats_mod
use iso_fortran_env, only: output_unit
use kind_mod, only: dp
use util_mod, only: default
implicit none
private
public :: mean, variance, sd, mean_and_sd, kurtosis, basic_stats, &
   print_basic_stats, basic_stats_names, correl, acf, nbasic_stats, &
   stat, stats, corr_mat, rms, moving_sum, moving_average
integer, parameter :: nbasic_stats = 6
character (len=*), parameter :: basic_stats_names(nbasic_stats) = &
   [character(len=4) :: "mean", "sd", "skew", "kurt", "min", "max"]
real(kind=dp), parameter :: bad_value = -huge(1.0d0)
interface stats
   module procedure stats_many_vec, stats_many_mat
end interface stats
contains

function stats_many_vec(funcs, x) result(y)
! return statistics on x(:)
character (len=*), intent(in) :: funcs(:)
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: y(size(funcs))
integer :: i
do i=1,size(funcs)
   y(i) = stat(funcs(i), x)
end do
end function stats_many_vec

function stats_many_mat(funcs, x) result(y)
! return a matrix of statistics on each column of x(:,:)
character (len=*), intent(in) :: funcs(:)
real(kind=dp), intent(in) :: x(:,:)
real(kind=dp) :: y(size(funcs), size(x,2))
integer :: i
do i=1,size(x,2)
   y(:,i) = stats_many_vec(funcs, x(:,i))
end do
end function stats_many_mat

function stat(func, x) result(y)
! return a statistic on x(:)
character (len=*), intent(in) :: func
real(kind=dp), intent(in) :: x(:)
real(kind=dp)             :: y
select case(func)
   case ("mean")    ; y = mean(x)
   case ("sd")      ; y = sd(x)
   case ("variance"); y = variance(x)
   case ("skew")    ; y = skew(x)
   case ("kurt")    ; y = kurtosis(x)
   case ("min")     ; y = minval(x)
   case ("max")     ; y = maxval(x)
   case ("first")
      if (size(x) > 0) then
         y = x(1)
      else
         y = bad_value
      end if 
   case ("last")
      if (size(x) > 0) then
         y = x(size(x))
      else
         y = bad_value
      end if 
   case default ; y = -huge(x)
end select
end function stat

pure function mean(x) result(xmean)
! return the mean of x(:)
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: xmean
xmean = sum(x)/max(1,size(x))
end function mean

pure function sd(x) result(xsd)
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

pure function rms(x) result(xrms)
! return the root-mean-square of x(:)
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: xrms
xrms = sqrt(sum(x**2)/size(x))
end function rms

pure function mean_and_sd(x) result(res)
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

pure function variance(x) result(var)
! return the variance of x(:)
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: var, m
integer :: n
n = size(x)
m = sum(x) / n
var = sum((x - m)**2) / (n-1)
end function variance

pure function skew(x) result(skew_val)
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

pure function kurtosis(x) result(kurtosis_val)
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

pure function basic_stats(x) result(stats)
real(kind=dp), intent(in) :: x(:)
real(kind=dp)             :: stats(nbasic_stats)
stats = [mean(x), sd(x), skew(x), kurtosis(x), minval(x), maxval(x)]
end function basic_stats

subroutine print_basic_stats(x, outu, fmt_header)
real(kind=dp), intent(in) :: x(:)
integer, intent(in), optional :: outu
character (len=*), intent(in), optional :: fmt_header
integer :: i, outu_
outu_ = default(output_unit, outu)
if (present(fmt_header)) write (outu_, fmt_header)
write (outu_, "(*(a10))") (trim(basic_stats_names(i)), i=1,nbasic_stats)
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

pure function acf(x, nacf) result(xacf)
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

pure function corr_mat(x) result(cor)
    ! return the correlation matrix of the columns of x(:,:)
    real(kind=dp), intent(in) :: x(:,:)
    real(kind=dp)             :: cor(size(x,2), size(x,2))
    real(kind=dp)             :: mean_vec(size(x,2)), std_vec(size(x,2))
    real(kind=dp)             :: centered_x(size(x,1), size(x,2))
    integer                   :: n, p

    n = size(x, 1)  ! Number of rows
    p = size(x, 2)  ! Number of columns

    ! Compute the mean of each column
    mean_vec = sum(x, dim=1) / n

    ! Center the matrix by subtracting the mean of each column
    centered_x = x - spread(mean_vec, dim=1, ncopies=n)

    ! Compute the standard deviation of each column
    std_vec = sqrt(sum(centered_x**2, dim=1) / (n - 1))

    cor = matmul(transpose(centered_x), centered_x) / (n - 1)
    cor = cor / spread(std_vec, dim=1, ncopies=p)
    cor = cor / spread(std_vec, dim=2, ncopies=p)
end function corr_mat

pure function moving_sum(x, k) result(xsum)
! return a moving sum of x(:) with k terms, using fewer terms for i < k
real(kind=dp), intent(in) :: x(:)
integer      , intent(in) :: k
real(kind=dp)             :: xsum(size(x))
integer                   :: i, n
n = size(x)
if (n < 1) return
if (k < 1) then
   xsum = 0.0_dp
   return
end if
xsum(1) = x(1)
do i=2,min(k, n)
   xsum(i) = xsum(i-1) + x(i)
end do
do i=k+1, n
   xsum(i) = xsum(i-1) + x(i) - x(i-k)
end do
end function moving_sum

pure function moving_average(x, k) result(xma)
! return a moving average of x(:) with k terms, using fewer terms for i < k
real(kind=dp), intent(in) :: x(:)
integer      , intent(in) :: k
real(kind=dp)             :: xma(size(x))
integer                   :: i, n
real(kind=dp)             :: xsum(size(x))
n = size(x)
if (k < 1) then
   xma = 0.0_dp
   return
end if
xsum = moving_sum(x, k)
do i=1,min(k, n)
   xma(i) = xsum(i)/i
end do
do i=k+1,n
   xma(i) = xsum(i)/k
end do
end function moving_average

end module basic_stats_mod
