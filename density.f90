module density_mod
  use kind_mod, only: dp
  use constants_mod, only: log_two, log_two_pi, sqrt_two, pi
  implicit none
  private
  public :: minus_log_density, log_density, normal_str, laplace_str, &
             t3_str, t4_str, t5_str, t6_str, t7_str, t8_str, t9_str, t10_str, &
             t20_str, t30_str, moment_stand_4

  ! Distribution string parameters
  character(len=*), parameter :: normal_str = "normal", &
                               laplace_str = "laplace", &
                               t3_str      = "t3", &
                               t4_str      = "t4", &
                               t5_str      = "t5", &
                               t6_str      = "t6", &
                               t7_str      = "t7", &
                               t8_str      = "t8", &
                               t9_str      = "t9", &
                               t10_str     = "t10", &
                               t20_str     = "t20", &
                               t30_str     = "t30"

contains

  elemental function minus_log_density(x, dist, sigma) result(y)
    ! Return -log(density(x)) for the probability distribution 'dist'
    ! with scale parameter sigma at x.
    real(kind=dp),    intent(in) :: x
    character(len=*), intent(in) :: dist
    real(kind=dp),    intent(in) :: sigma
    real(kind=dp)                :: y
    real(kind=dp)                :: nu, x_, factor

    select case(dist)
    case(normal_str)
      y = 0.5_dp * log_two_pi + log(sigma) + 0.5_dp * ((x/sigma)**2)
    case(laplace_str)
      y = 0.5_dp * log_two + log(sigma) + sqrt_two * abs(x)/sigma
    case(t3_str)
      nu = 3.0_dp
      factor = sqrt(nu/(nu-2.0_dp))
      x_ = factor * x
      y = log(sigma/factor) + 0.5_dp*log(nu*pi) &
          - log_gamma((nu+1.0_dp)/2.0_dp) + log_gamma(nu/2.0_dp) &
          + 0.5_dp*(nu+1.0_dp)*log(1.0_dp + (x_/sigma)**2/nu)
    case(t4_str)
      nu = 4.0_dp
      factor = sqrt(nu/(nu-2.0_dp))
      x_ = factor * x
      y = log(sigma/factor) + 0.5_dp*log(nu*pi) &
          - log_gamma((nu+1.0_dp)/2.0_dp) + log_gamma(nu/2.0_dp) &
          + 0.5_dp*(nu+1.0_dp)*log(1.0_dp + (x_/sigma)**2/nu)
    case(t5_str)
      nu = 5.0_dp
      factor = sqrt(nu/(nu-2.0_dp))
      x_ = factor * x
      y = log(sigma/factor) + 0.5_dp*log(nu*pi) &
          - log_gamma((nu+1.0_dp)/2.0_dp) + log_gamma(nu/2.0_dp) &
          + 0.5_dp*(nu+1.0_dp)*log(1.0_dp + (x_/sigma)**2/nu)
    case(t6_str)
      nu = 6.0_dp
      factor = sqrt(nu/(nu-2.0_dp))
      x_ = factor * x
      y = log(sigma/factor) + 0.5_dp*log(nu*pi) &
          - log_gamma((nu+1.0_dp)/2.0_dp) + log_gamma(nu/2.0_dp) &
          + 0.5_dp*(nu+1.0_dp)*log(1.0_dp + (x_/sigma)**2/nu)
    case(t7_str)
      nu = 7.0_dp
      factor = sqrt(nu/(nu-2.0_dp))
      x_ = factor * x
      y = log(sigma/factor) + 0.5_dp*log(nu*pi) &
          - log_gamma((nu+1.0_dp)/2.0_dp) + log_gamma(nu/2.0_dp) &
          + 0.5_dp*(nu+1.0_dp)*log(1.0_dp + (x_/sigma)**2/nu)
    case(t8_str)
      nu = 8.0_dp
      factor = sqrt(nu/(nu-2.0_dp))
      x_ = factor * x
      y = log(sigma/factor) + 0.5_dp*log(nu*pi) &
          - log_gamma((nu+1.0_dp)/2.0_dp) + log_gamma(nu/2.0_dp) &
          + 0.5_dp*(nu+1.0_dp)*log(1.0_dp + (x_/sigma)**2/nu)
    case(t9_str)
      nu = 9.0_dp
      factor = sqrt(nu/(nu-2.0_dp))
      x_ = factor * x
      y = log(sigma/factor) + 0.5_dp*log(nu*pi) &
          - log_gamma((nu+1.0_dp)/2.0_dp) + log_gamma(nu/2.0_dp) &
          + 0.5_dp*(nu+1.0_dp)*log(1.0_dp + (x_/sigma)**2/nu)
    case(t10_str)
      nu = 10.0_dp
      factor = sqrt(nu/(nu-2.0_dp))
      x_ = factor * x
      y = log(sigma/factor) + 0.5_dp*log(nu*pi) &
          - log_gamma((nu+1.0_dp)/2.0_dp) + log_gamma(nu/2.0_dp) &
          + 0.5_dp*(nu+1.0_dp)*log(1.0_dp + (x_/sigma)**2/nu)
    case(t20_str)
      nu = 20.0_dp
      factor = sqrt(nu/(nu-2.0_dp))
      x_ = factor * x
      y = log(sigma/factor) + 0.5_dp*log(nu*pi) &
          - log_gamma((nu+1.0_dp)/2.0_dp) + log_gamma(nu/2.0_dp) &
          + 0.5_dp*(nu+1.0_dp)*log(1.0_dp + (x_/sigma)**2/nu)
    case(t30_str)
      nu = 30.0_dp
      factor = sqrt(nu/(nu-2.0_dp))
      x_ = factor * x
      y = log(sigma/factor) + 0.5_dp*log(nu*pi) &
          - log_gamma((nu+1.0_dp)/2.0_dp) + log_gamma(nu/2.0_dp) &
          + 0.5_dp*(nu+1.0_dp)*log(1.0_dp + (x_/sigma)**2/nu)
    case default
      error stop "invalid dist '" // trim(dist) // "' in minus_log_density"
    end select
  end function minus_log_density

  elemental function log_density(x, dist, sigma) result(y)
    ! Return the density at x (i.e. exp(-minus_log_density)).
    real(kind=dp),    intent(in) :: x
    character(len=*), intent(in) :: dist
    real(kind=dp),    intent(in) :: sigma
    real(kind=dp)                 :: y
    y = exp(-minus_log_density(x, dist, sigma))
  end function log_density

  elemental function moment_stand_4(dist) result(y)
    ! Return the 4th central moment of the distribution with variance = 1.
    ! For Student t distributions the fourth moment exists only if df > 4.
    character(len=*), intent(in) :: dist
    real(kind=dp)                 :: y
    select case(dist)
    case(normal_str)
      y = 3.0_dp
    case(laplace_str)
      y = 6.0_dp
    case(t3_str, t4_str)
      y = -1.0_dp ! moment does not exist
    case(t5_str)
      ! For df=5, standardized 4th moment = 3*(5-2)/(5-4)=9.
      y = 9.0_dp
    case(t6_str)
      ! For df=6, standardized 4th moment = 3*(6-2)/(6-4)=6.
      y = 6.0_dp
    case(t7_str)
      ! For df=7, standardized 4th moment = 3*(7-2)/(7-4)=5.
      y = 5.0_dp
    case(t8_str)
      ! For df=8, standardized 4th moment = 3*(8-2)/(8-4)=4.5.
      y = 4.5_dp
    case(t9_str)
      ! For df=9, standardized 4th moment = 3*(9-2)/(9-4)=4.2.
      y = 4.2_dp
    case(t10_str)
      ! For df=10, standardized 4th moment = 3*(10-2)/(10-4)=4.
      y = 4.0_dp
    case(t20_str)
      ! For df=20, standardized 4th moment = 3*(20-2)/(20-4)=3.375.
      y = 3.375_dp
    case(t30_str)
      ! For df=30, standardized 4th moment = 3*(30-2)/(30-4)=3.176.
      y = 3.176_dp
    case default
      error stop "invalid dist '" // trim(dist) // "' in moment_stand_4"
    end select
  end function moment_stand_4

end module density_mod
