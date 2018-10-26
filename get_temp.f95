subroutine get_temp(u, T)
use data
implicit none

integer :: i
real(8), intent(inout) :: u(nnodes), T(nnodes)

do i = 1, nnodes
  T(i) = -1.8768e-6*u(i)**3 + 2.4635e-3*u(i)**2 - 0.20891*u(i) + 241.35
end do

end subroutine get_temp
