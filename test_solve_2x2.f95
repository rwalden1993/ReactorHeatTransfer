program test_solve_2x2
implicit none

real(8) :: A(2,2), B(2)

!Solve general case
A(1,1) = 1.0d0
A(1,2) = 2.0d0
A(2,1) = 2.0d0
A(2,2) = 6.0d0
B = (/ 3.0d0, 8.0d0 /)

write(*,*) A
write(*,*) B

call solve_2x2(A, B)

write(*,*) "general case"
write(*,*) 1, B(1)
write(*,*) 1, B(2)


!solve A(1,1) = 0
A(1,1) = 0.0d0
B = (/ 3.0d0, 8.0d0 /)

call solve_2x2(A, B)

write(*,*) "A(1,1) = 0"
write(*,*) -0.5, B(1)
write(*,*) 1.5, B(2)

!solve A(1,2) = 0
A(1,1) = 1.0d0
A(1,2) = 0.0d0
B = (/ 3.0d0, 8.0d0 /)

call solve_2x2(A, B)

write(*,*) "A(1,2) = 0"
write(*,*) 3.0, B(1)
write(*,*) 0.33, B(2)

!solve A(2,1) = 0
A(1,2) = 2.0d0
A(2,1) = 0.0d0
B = (/ 3.0d0, 8.0d0 /)

call solve_2x2(A, B)

write(*,*) "A(2,1) = 0"
write(*,*) 0.33, B(1)
write(*,*) 1.33, B(2)

!solve A(2,2) = 0
A(2,1) = 2.0d0
A(2,2) = 0.0d0
B = (/ 3.0d0, 8.0d0 /)

call solve_2x2(A, B)

write(*,*) "A(2,2) = 0"
write(*,*) 4.0, B(1)
write(*,*) -0.5, B(2)

!solve diag
A = 0.0d0
A(1,1) = 2.0d0
A(2,2) = 4.0d0
B = (/ 3.0d0, 8.0d0 /)

call solve_2x2(A, B)

write(*,*) "diag"
write(*,*) 1.5, B(1)
write(*,*) 2.0, B(2)

!stop at wrong condition
A(2,1) = 0.0d0
B = (/ 3.0d0, 8.0d0 /)

call solve_2x2(A,B)

end program test_solve_2x2
