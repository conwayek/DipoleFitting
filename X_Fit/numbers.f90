program numbers
implicit none

integer :: i

open(1,file='numbers.txt',action='write',status='unknown')

do i=1, 665
write(1,"(i3)") i
end do


end program numbers
