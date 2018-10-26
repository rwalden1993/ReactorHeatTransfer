program final_project
Use data
implicit none

character(16) :: infile
integer :: i, error, ios, in_unit
real(8) :: time, max_time, dt_min, dt_max, time_o, time_s
real(8) :: m1, m1_prev, m2, m2_prev, mc, mc_prev, dP_core, dP_core_prev, t_co
real(8) :: U(3), U_prev(3), epsilon

!Intialize data
call get_data()
dt_min = 0.1
dt_max = 1.0
dt = dt_min
time = 0.0
max_time = 60
epsilon = 0.0
time_o = 1.0e9

call getarg(1, infile)
in_unit = 13
open(unit=in_unit, file=infile, status='old', action='read', iostat=ios)
if (ios.ne.0) then
  write(*,*) "An error occured when opening the input file"
endif

read(13,*) Q, mc, m1, m2
read(13,*) u1

!test mass data
call set_dP_pumps(m1, m2, time)
dP_0 = dP_pump2
!dP_pump1 = dP_pump_rated
!dP_pump2 = dP_pump1
! m1 = 36332.16/4
! m2 = m1
! mc = m1*4
dP_core = 0.0
!u1 = 531.758
u2 = u1
call get_temp(u1, T1)
!T1 = 544.65
T2 = T1
call get_rho(T1, rho1)
!rho1 = 47.1832
rho2 = rho1
!Q = 0.0
T_fuel = T1(1)

open(unit=11, file="fuel_conduction_debug.txt", iostat=ios, status="replace", action="write")
if ( ios .ne. 0 ) stop "Error opening file "

! do i=1, nnodes
!   call get_q(time, i, T1, m1)
!   write(*,*) time
! enddo
! stop

do
  m1_prev = m1
  m2_prev = m2
  mc_prev = mc
  dP_core_prev = dP_core

  ! write(11,'(F6.1, E10.3, F9.2)',advance='no') time, Q, mc
  ! !write(*,*)
  ! do i=1, fuel_nodes
  !   write(11,'(F7.1)',advance='no') T_fuel(1,i)
  ! enddo
  ! write(11,*)
  ! call fuel_heat_conduction()

  !ramp power
  ! if (Q .lt. Q_rated) then
  !   Q = Q_rated*0.05*(time/60)
  !   if (Q .gt. Q_rated) then
  !     Q = Q_rated
  !   endif
  ! endif
  if (time .ge. 2.0) then
    time_s = time - 2.0
    Q = Q_rated*0.1*((time_s + 10)**(-0.2) - (time_o + time_s + 10)**(-0.2) &
        + 0.87*(time_o + time_s + 2.0e7)**(-0.2) - 0.87*(time_s + 2.0e7)**(-0.2))
  endif

  !write out variables
  ! test 1
  ! write(*,'(F6.2, 4F9.2)') time, mc, m1, m2, dP_core

  !test 2
  write(*,'(F6.1, E10.3, 2F9.2)',advance='no') time, Q, mc, m1
  !write(*,*)
  do i=1, nnodes
    write(*,'(F7.1)',advance='no') T1(i)
  enddo
  !write(*,*)
  do i=1, 4
    call clad_temp(t_co, i, mc)
    write(*,'(F7.1)',advance='no') t_co
  enddo
  write(*,*)


  if(time .gt. max_time) exit

  !mass momentum equation
  call set_dP_pumps(m1, m2, time)
  call mass_momentum(m1, m2, mc, dP_core, error)

  !internal energy equation
  call energy(m1, m2, mc, m1_prev, m2_prev, mc_prev)

  !state equations for temp, density
  call get_temp(u1, T1)
  call get_temp(u2, T2)
  call get_rho(T1, rho1)
  call get_rho(T2, rho2)

  !Newton Time Iteration
  if(error .eq. -1) then
    if (dt .eq. dt_min) then
      write(*,*) "Did not converge at time ", time
      exit
    endif
    dt = max(dt_min, dt/2)
    m1 = m1_prev
    m2 = m2_prev
    mc = mc_prev
    dP_core = dP_core_prev
  else
    !calculate epsilon
    U_prev = (/ m1_prev, m2_prev, mc_prev /)
    U = (/ m1, m2, mc /)

    epsilon = maxval( abs(U-U_prev) / U )

    ! if (epsilon .gt. 0.01) then
    !   dt = max(dt_min, dt/2)
    !   m1 = m1_prev
    !   m2 = m2_prev
    !   mc = mc_prev
    !   dP_core = dP_core_prev
    ! else if (epsilon .lt. 0.001) then
    !   time = time + dt
    !   dt = min(dt_max, 1.2*dt)
    ! endif
    time = time  + dt
  endif
enddo

open(unit=12, file="ouput.txt", iostat=ios, status="replace", action="write")
if ( ios .ne. 0 ) stop "Error opening file "

write(12,*) Q, mc, m1, m2
write(12,*) u2

end program final_project
