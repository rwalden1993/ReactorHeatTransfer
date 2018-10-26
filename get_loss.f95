subroutine get_loss(loss, node, m, loop)
use data
implicit none

integer :: node, loop
real(8) :: loss, k, k_plus, k_minus, m, my_ax

loss = 0.0
k    = 0.0

my_ax = Ax(node)
if (loop .eq. 1 .and. (node .ge. 7 .and. node .le. 13) ) then
  my_ax = Ax(node)*tubes_pluged
endif

if  (loop .eq. 2) then
  ! regular
  k_plus = 1.0e-12
  k_minus = 1.0e12
else
  !regular
  k_plus = 1.0e-12
  k_minus = 1.0e12

  !lock rotor
  ! k_plus = 1.0e12
  ! k_minus = 1.0e12

  !broken shaft
  ! k_plus = 0.0
  ! k_minus = 0.0
endif

if (node .eq. 1) then
  k = 5.25
else if (node .eq. 2) then
  k = 1.0
elseif ( node .eq. 3 ) then
  k = 1.0
elseif ( node .eq. 4 ) then
  k = 5.25
elseif ( node .eq. 6 ) then
  k = 1.5
elseif ( node .eq. 7 ) then
  k = 0.5
elseif ( node .eq. 13 ) then
  k = 1.0
elseif ( node .eq. 14 ) then
  k = 5.1 + 0.5*(k_plus + k_minus) + 0.5*(abs(m)/m)*(k_plus -  k_minus)
endif

loss = k/(2*rho_const*gc)*(1/my_ax)**2

end subroutine get_loss
