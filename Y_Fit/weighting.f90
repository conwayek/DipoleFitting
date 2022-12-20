program weighting

implicit none

integer :: i, ierr, lines
double precision :: dummy, emax, etop, alpha, test
double precision,allocatable :: r1(:), r2(:), theta(:), energy(:),s(:), weight(:),dipole(:)

open(1,file='matched.6z.dipoles.txt',status='old',action='read')
open(2,file='weighted.6z.dipx.txt',status='unknown',action='write')

lines = 0
do i =1, 18000
lines = lines + 1
read(1,*,iostat=ierr) dummy
if(ierr/= 0) exit
end do

lines = lines - 1
rewind 1

allocate(r1(lines), r2(lines), theta(lines), energy(lines),s(lines),weight(lines),dipole(lines))

etop = 30000
alpha = 6*(10**(-4))

do i =1,lines
read(1,*) dummy, r1(i), r2(i), theta(i), energy(i), dipole(i)

dummy = 25000

emax = max(dummy,energy(i))

s(i) = ( TANH( -0.006*(  energy(i) - etop  )    ) + 1.002002002 )/2.002002002

weight(i) = dummy*s(i)/emax




write(2,"(f8.4,2x,f8.4,2x,f10.5,4x,f12.9,4x,f20.10,4x,f20.10)") r1(i), r2(i), theta(i), dipole(i), weight(i), energy(i)

end do







end program weighting
