module gradient_mod
use kind_mod, only: dp
implicit none
private
public :: gradient
! Default step size (square root of machine epsilon for double precision)
real(dp), parameter :: DEFAULT_DX = sqrt(epsilon(1.0_dp))
character (len=*), parameter :: central_str = "central", &
   forward_str = "forward", backward_str = "backward"
contains

function gradient(func, x, dx, method) result(grad)
!> Computes the gradient of func at x using central finite differences
!! Arguments:
!!   func - function to differentiate (real scalar function of real vector)
!!   x    - input point where gradient is computed
!!   dx   - optional step size (default: sqrt(epsilon) scaled by |x|)
!! Returns: gradient vector (same size as x)
!! Accuracy: O(dx^2)
! Interface for the input function
interface
   function func(x) result(y)
   import dp
   real(kind=dp), intent(in) :: x(:)
   real(kind=dp)       :: y
end function func
end interface
      
   ! Arguments
real(kind=dp), intent(in) :: x(:)! Input array
real(kind=dp), intent(in), optional :: dx! Step size (optional)
character (len=*), intent(in), optional :: method
character (len=20) :: method_
! Result
real(kind=dp) :: grad(size(x))   ! Gradient array

! Local variables
real(kind=dp) :: x_plus(size(x)) ! x + dx
real(kind=dp) :: x_minus(size(x))! x - dx
real(dp) :: h ! Actual step size used
integer :: i, n
if (present(method)) then
   method_ = method
else
   method_ = central_str
end if
! Get size of input array
n = size(x)

! Set step size: use provided dx if present, otherwise default
if (present(dx)) then
   h = dx
else
! Scale default step size by magnitude of x if non-zero
   if (any(abs(x) > 0.0_dp)) then
     h = DEFAULT_DX * maxval(abs(x))
   else
     h = DEFAULT_DX
   end if
end if    
! Ensure h is positive and reasonable
h = max(abs(h), DEFAULT_DX)
! Compute gradient for each component using central difference
do i = 1, n
   if (method_ == central_str .or. method_ == forward_str) then
      x_plus = x
      x_plus(i) = x(i) + h  ! Perturb the i-th component
   end if
   if (method_ == central_str .or. method_ == backward_str) then
      x_minus = x
      x_minus(i) = x(i) - h ! Perturb the i-th component
   end if
   
   select case (method_)
      case(central_str)   ! central difference formula: (f(x+dx) - f(x-dx)) / (2*dx)
         grad(i) = (func(x_plus) - func(x_minus)) / (2.0_dp * h)
      case (forward_str)  ! forward difference formula: (f(x+dx) - f(x)) / dx
         grad(i) = (func(x_plus) - func(x)) / h
      case (backward_str) ! backward difference formula: (f(x) - f(x-dx)) / dx
         grad(i) = (func(x) - func(x_minus)) / h
      case default
         error stop "finite difference method " // trim(method_) // " not implemented"
   end select
end do
end function gradient

end module gradient_mod
