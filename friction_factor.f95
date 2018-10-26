subroutine friction_factor(f, m, node, loop)
use data
implicit none

integer :: node, loop
real(8) :: m, Re, f, f1, f2, my_ax

my_ax = Ax(node)
if (loop .eq. 1 .and. (node .ge. 7 .and. node .le. 13) ) then
  my_ax = Ax(node)*tubes_pluged
endif

f = 0.0
Re = ((abs(m)/my_ax)*(De(node)/12))/mu

if (Re .eq. 0) then
  f = 0.1
else if (Re .lt. 2300) then
  f = 64/Re
else if (Re .gt. 4200 .and. Re .lt. 30000) then
  f = 0.3164*Re**(-0.25)
else if (Re .gt. 30000) then
  f = 0.184*Re**(-0.2)
else
  f1 = 64/2300
  f2 = 0.3164*4200**(-0.25)
  f  = ((Re-4200)/(2300-4200))*(f1-f2)+f2
endif

end subroutine friction_factor
