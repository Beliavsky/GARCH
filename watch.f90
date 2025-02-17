module watch_mod
use kind_mod, only: dp
implicit none
private
public              :: watch,print_elapsed_times
integer, parameter  :: max_times = 1000
real(kind=dp)       :: cpu_times(max_times)
integer             :: wall_times(max_times)
integer, save       :: ntimes = 0
character (len=100) :: time_labels(max_times)
contains

subroutine watch(label)
! store the CPU and wall times after the computation 
! described in label was done
character (len=*), intent(in), optional :: label
ntimes = ntimes + 1
if (ntimes <= max_times) then
   call cpu_time(cpu_times(ntimes))
   call system_clock(wall_times(ntimes))
   if (present(label)) then
      time_labels(ntimes) = label
   else
      time_labels(ntimes) = ""
   end if
end if
end subroutine watch
!
subroutine print_elapsed_times()
! print elapsed CPU and wall times based on the module 
! variables time_labels, cpu_times, wall_times
integer :: i, itick
call system_clock(count_rate = itick)
write (*,"(/,a40,2a12)") "task","cpu_time","wall_time"
do i=2,ntimes
   write (*,"(a40,2f12.6)") trim(time_labels(i)),cpu_times(i)-cpu_times(i-1), &
                           (wall_times(i)-wall_times(i-1))/real(itick,kind=dp)
   if (i == ntimes) write (*,"(a40,2f12.6)") "TOTAL",cpu_times(i)-cpu_times(1), &
                           (wall_times(i)-wall_times(1))/real(itick,kind=dp)
end do
end subroutine print_elapsed_times
end module watch_mod
