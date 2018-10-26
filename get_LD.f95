subroutine get_LD(ld, node)
use data
implicit none

integer :: node
real(8) :: ld

if (node .eq. 6) then
  ld = 20.0
else if (node .eq. 10) then
  ld = 55.0
else if (node .eq. 14) then
  ld = 18.0
else
  ld = length(node)/(De(node)/12) !length -> ft, De -> in
endif

end subroutine get_LD
