subroutine energy(m1, m2, mc, m1_prev, m2_prev, mc_prev)
use data
implicit none

real(8), intent(in) :: m1, m2, mc, m1_prev, m2_prev, mc_prev

integer :: i, n = 26
real(8),allocatable :: A(:,:), B(:), T(:,:), R(:,:), Bot(:,:), F(:,:), S1(:), S2(:)
real(8) :: my_q, vol_lp, vol_up

allocate(A(n,n), B(n), T(n-2,n-2), R(n-2,2), Bot(2, n-2), F(2,2), S1(n-2), S2(2))
A = 0.0
B = 0.0


!Core
A(1,1) = (Ax(1)*length(1)*rho2(1))/dt + abs(mc)
A(1,2) = (mc/2)*(1 - abs(mc)/mc)
A(1,n-1) = -(mc/2)*(1 + abs(mc)/mc)
call get_q(my_q, 1, mc_prev)
B(1) = (Ax(1)*length(1)*rho2(1)*u2(1))/dt + my_q
do i=2, 3
  A(i,i)   = (Ax(i)*length(i)*rho2(i))/dt + abs(mc)
  A(i,i+1) = (mc/2)*(1 - abs(mc)/mc)
  A(i,i-1) = -(mc/2)*(1 + abs(mc)/mc)
  call get_q(my_q, i, mc_prev)
  B(i) = (Ax(i)*length(i)*rho2(i)*u2(i))/dt + my_q
enddo
A(4,4)   = (Ax(4)*length(4)*rho2(4))/dt + abs(mc)
A(4,n) = (mc/2)*(1 - abs(mc)/mc)
A(4,3) = -(mc/2)*(1 + abs(mc)/mc)
call get_q(my_q, 4, mc_prev)
B(4) = (Ax(4)*length(4)*rho2(4)*u2(4))/dt + my_q

!loop 1
A(5,5) = (Ax(6)*length(6)*rho1(6))/dt + abs(m1)
A(5,6) = (m1/2)*(1 - abs(m1)/m1)
A(5,n) = -(m1/2)*(1 + abs(m1)/m1)
call get_q(my_q, 5, m1_prev)
B(5) = (Ax(6)*length(6)*rho1(6)*u1(6))/dt + my_q
do i=6, 13
  A(i,i)   = (Ax(i+1)*length(i+1)*rho1(i+1))/dt + abs(m1)
  A(i,i+1) = (m1/2)*(1 - abs(m1)/m1)
  A(i,i-1) = -(m1/2)*(1 + abs(m1)/m1)
  call get_q(my_q, i, m1_prev)
  B(i) = (Ax(i+1)*length(i+1)*rho1(i+1)*u1(i+1))/dt + my_q
enddo
A(14,14)   = (Ax(15)*length(15)*rho1(15))/dt + abs(m1)
A(14,n-1) = (m1/2)*(1 - abs(m1)/m1)
A(14,13) = -(m1/2)*(1 + abs(m1)/m1)
call get_q(my_q, 14, m1_prev)
B(14) = (Ax(15)*length(15)*rho1(15)*u1(15))/dt + my_q


!loop 2
A(15,15)   = (Ax(6)*length(6)*rho2(6))/dt + abs(m2)
A(15,16) = (m2/2)*(1 - abs(m2)/m2)
A(15,n) = -(m2/2)*(1 + abs(m2)/m2)
call get_q(my_q, 15, m2_prev)
B(15) = (Ax(6)*length(6)*rho2(6)*u2(6))/dt + my_q
do i=16, 23
  A(i,i)   = (Ax(i-9)*length(i-9)*rho2(i-9))/dt + abs(m2)
  A(i,i+1) = (m2/2)*(1 - abs(m2)/m2)
  A(i,i-1) = -(m2/2)*(1 + abs(m2)/m2)
  call get_q(my_q, i, m2_prev)
  B(i) = (Ax(i-9)*length(i-9)*rho2(i-9)*u2(i-9))/dt + my_q
enddo
A(24,24)   = (Ax(15)*length(15)*rho2(15))/dt + abs(m2)
A(24,n-1) = (m2/2)*(1 - abs(m2)/m2)
A(24,23) = -(m2/2)*(1 + abs(m2)/m2)
call get_q(my_q, 24, m2_prev)
B(24) = (Ax(15)*length(15)*rho2(15)*u2(15))/dt + my_q

!manafolds
vol_lp = 784.4 !ft^3
!vol_lp = 784.4 * (12**3) !in^3
A(25,25) = (vol_lp*rho2(16))/dt + (0.5)*(abs(m1)+abs(3*m2)+abs(mc))
A(25,14) = -(m1/2)*(1 + abs(m1)/m1)
A(25,24) = -(3*m2/2)*(1 + abs(3*m2)/(3*m2))
A(25,1)  = (mc/2)*(1 - abs(mc)/mc)
B(25)    = (vol_lp*rho2(16)*u2(16))/dt

vol_up = 1373.7 !ft^3
!vol_up = 1373.7 * (12**3) !in^3
A(26,26) = (vol_up*rho2(5))/dt + (0.5)*(abs(m1)+abs(3*m2)+abs(mc))
A(26,5)  = (m1/2)*(1 - abs(m1)/m1)
A(26,15) = (3*m2/2)*(1 - abs(3*m2)/(3*m2))
A(26,4)  = -(mc/2)*(1 + abs(mc)/mc)
B(26)    = (vol_up*rho2(5)*u2(5))/dt

! do j = 1, n
!   do i =1, n
!     if (A(j,i) .ne. 0.0) then
!       write(*,'(E11.4)',advance='no') A(j,i)
!       !write(*,'(I2)',advance='no') 1
!     else
!       write(*,'(I2)',advance='no') 0
!     endif
!   enddo
!   write(*,*)
! enddo
!
! write(*,*)
! do i =1, n
!   write(*,'(E11.3)') B(i)
! enddo
! write(*,*)
! stop

T   = A(1:n-2,1:n-2)
R   = A(1:n-2,n-1:n)
Bot = A(n-1:n,1:n-2)
F   = A(n-1:n,n-1:n)
S1  = B(1:n-2)
S2  = B(n-1:n)

call thomas_alg(n-2, T, 1, S1)
call thomas_alg(n-2, T, 2, R)

F  = F - matmul(Bot,R)
S2 = S2 - matmul(Bot,S1)

call solve_2x2(F,S2)

B(n-1:n) = S2
B(1:n-2) = S1 - matmul(R,S2)

u1(1:4) = B(1:4)
u2(1:4) = B(1:4)
u1(5)   = B(26)
u2(5)   = B(26)
u1(6:15)= B(5:14)
u2(6:15)= B(15:24)
u1(16)  = B(25)
u2(16)  = B(25)

end subroutine energy


subroutine get_q(my_q, node, m)
use data
implicit none

integer, intent(in) :: node
real(8), intent(in) :: m
real(8), intent(inout) :: my_q

real(8) :: my_UA_sg, Ax_c, De_c, D, S, hc, Re, Nu, t_co

my_q = 0.0

Ax_c = 52.74
De_c = 0.4635/12
D = 0.374/12
S = 0.496/12

!core
if (node .le. 4) then
  my_q = Q/4 * Fk(node)
  ! call clad_temp(t_co, node)
  ! Re = ((m/Ax_c) * De_c)/mu
  ! Nu = (0.042*(S/D) - 0.024)*(Re**0.8)*(Pr**(1.0/3.0))
  ! hc = (Nu * k_th)/De_c
  ! my_q = hc*(t_co - T2(node))
endif

!SG
if (node .ge. 6 .and. node .le. 12) then
  call get_UA_sg(my_UA_sg, m, node)
  my_q = my_UA_sg * (Tsat_sg - T1(node+1))
endif
if (node .ge. 16 .and. node .le. 22) then
  call get_UA_sg(my_UA_sg, m, node)
  my_q = my_UA_sg * (Tsat_sg - T2(node-9))
endif

end subroutine get_q

subroutine get_UA_sg(UA_sg, m, node)
use data
implicit none

integer :: node
real(8), intent(in) :: m
real(8), intent(inout) :: UA_sg

integer :: n_tubes1, n_tubes2
real(8) :: Re, Nu, hc, chi, Ax_sg, De_sg, SA_sg, length_sg

chi = 2.175

n_tubes1 = 6633*tubes_pluged
n_tubes2 = 6633

length_sg = 9.54
Ax_sg = 13.35
De_sg = 0.6075/12.0

if (node .le. 12) then
  SA_sg = n_tubes1*Pi*De_sg*length_sg
else
  SA_sg = n_tubes2*Pi*De_sg*length_sg
endif

Re = ((abs(m)/Ax_sg) * De_sg)/mu

Nu = 0.023*(Re**0.8)*(Pr**0.3)

hc = (Nu * k_th)/De_sg

UA_sg = SA_sg / (1/hc + chi)

end subroutine get_UA_sg
