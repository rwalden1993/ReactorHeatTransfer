subroutine  get_rho(T, rho)
use data
implicit none

integer :: i
real(8), intent(inout) :: T(nnodes), rho(nnodes)

do i=1, nnodes
  rho(i) = -1.4985e-6*T(i)**3 + 2.3749e-3*T(i)**2 - 1.3174*T(i) + 302.33
enddo

end subroutine get_rho
