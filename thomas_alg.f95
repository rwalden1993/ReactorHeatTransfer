subroutine thomas_alg(n, A, b_cols, B)
implicit none

integer, intent(in) :: n, b_cols
real(8), intent(inout) :: A(n,n), B(n,b_cols)

integer :: i, j
real(8), allocatable :: alpha(:), g(:)

allocate(alpha(n), g(n))

do j=1, b_cols

  alpha(1) = A(1,1)
  do i=2, n
    alpha(i) = A(i,i) - A(i,i-1)*(A(i-1,i)/alpha(i-1))
  enddo

  g(1) = B(1,j)/alpha(1)
  do i=2, n
    g(i) = (B(i,j) - A(i,i-1)*g(i-1))/alpha(i)
  enddo

  B(n,j) = g(n)
  do i=n-1, 1, -1
    B(i,j) = g(i) - (A(i,i+1)/alpha(i))*B(i+1,j)
  enddo

enddo

end subroutine thomas_alg

subroutine solve_2x2(my_A, my_B)

real(8), intent(in) :: my_A(2,2)
real(8), intent(inout) :: my_B(2)
real(8) :: x, y

if(sum(abs(my_A(1,:))) .lt. 1e-4 .or. sum(abs(my_A(2,:))) .lt. 1e-4 &
  .or. sum(abs(my_A(:,1))) .lt. 1e-4 .or. sum(abs(my_A(:,2))) .lt. 1e-4) stop "Could not solve 2x2 matrix"

if(abs(my_A(2,1)) .lt. 1e-4 .and. abs(my_A(1,2)) .lt. 1e-4) then !diagonal
  !write(*,*) "diag"
  x = my_B(1)/my_A(1,1)
  y = my_B(2)/my_A(2,2)
else if(abs(my_A(1,1)) .lt. 1e-4) then  !A(1,1) = 0
  !write(*,*) "1"
  y = my_B(1)/my_A(1,2)
  x = (my_B(2)-my_A(2,2)*y)/my_A(2,1)
else if(abs(my_A(1,2)) .lt. 1e-4) then !A(1,2) = 0
  !write(*,*) "2"
  x = my_B(1)/my_A(1,1)
  y = (my_B(2)-my_A(2,1)*x)/my_A(2,2)
else if(abs(my_A(2,1)) .lt. 1e-4) then !A(2,1) = 0
  !write(*,*) "3"
  y = my_B(2)/my_A(2,2)
  x = (my_B(1)-my_A(1,2)*y)/my_A(1,1)
else if(abs(my_A(2,2)) .lt. 1e-4) then !A(2,2) = 0
  !write(*,*) "4"
  x = my_B(2)/my_A(2,1)
  y = (my_B(1) - my_A(1,1)*x)/my_A(1,2)
else ! general
  !write(*,*) "general"
  y = (my_B(2) - (my_A(1,2)*my_B(1))/my_A(1,1)) / (my_A(2,2) - (my_A(2,1)*my_A(1,2))/my_A(1,1))
  x = (my_B(1) - my_A(1,2)*y)/my_A(1,1)
endif

my_B(1) = x
my_B(2) = y

end subroutine solve_2x2
