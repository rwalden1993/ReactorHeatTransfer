subroutine mass_momentum(m1_in, m2_in, mc_in, dP_in, error)
use data
implicit none

real(8), intent(inout) :: m1_in, m2_in, mc_in, dP_in

integer :: i, ite, max_ite, error
real(8) :: m1, m1_prev, m2, m2_prev, mc, mc_prev, dP, dP_prev
real(8) :: a1, a2, ac, b1, b2, bc, conv
real(8) :: my_f, my_ld, my_loss, friction_sum, loss_sum, dt_sum, boun_sum, my_ax

!set intial conditions
m1_prev = m1_in
m2_prev = m2_in
mc_prev = mc_in
dP_prev = dP_in
ite = 0
max_ite = 1000
error = 0


do
  ite = ite + 1

  !a1
  dt_sum =0.0
  loss_sum = 0.0
  friction_sum = 0.0
  do i=6, 15
    my_ax = Ax(i)
    if (i .ge. 7 .and. i .le. 13) then
      my_ax = Ax(i)*tubes_pluged
    endif
    call friction_factor(my_f, m1_prev, i, 1)
    call get_LD(my_ld, i)
    call get_loss(my_loss, i, m1_prev, 1)
    dt_sum = dt_sum + (length(i)/my_ax) * (1/dt)
    friction_sum = friction_sum + (my_f * my_ld * (1/my_ax)**2)
    loss_sum = loss_sum + my_loss
  enddo

  a1 = (1/gc)*dt_sum + 2*(friction_sum*(1/(2*rho_const*gc)) + loss_sum)*abs(m1_prev)

  !b1
  dt_sum =0.0
  loss_sum = 0.0
  friction_sum = 0.0
  boun_sum = 0.0
  do i=6, 15
    my_ax = Ax(i)
    if (i .ge. 7 .and. i .le. 13) then
      my_ax = Ax(i)*tubes_pluged
    endif
    call friction_factor(my_f, m1_prev, i, 1)
    call get_LD(my_ld, i)
    call get_loss(my_loss, i, m1_prev, 1)
    dt_sum = dt_sum + (length(i)/my_ax) * (1/dt)
    friction_sum = friction_sum + my_f * my_ld * (1/my_ax)**2
    loss_sum = loss_sum + my_loss
    boun_sum = boun_sum + (rho1(i)*dH(i))
  enddo

  b1 = (1/gc)*dt_sum*m1_in + (friction_sum*(1/(2*rho_const*gc)) + loss_sum) * m1_prev*abs(m1_prev) &
        - boun_sum + dP_pump1

  !a2
  dt_sum =0.0
  loss_sum = 0.0
  friction_sum = 0.0
  do i=6, 15
    call friction_factor(my_f, m2_prev, i, 2)
    call get_LD(my_ld, i)
    call get_loss(my_loss, i, m2_prev, 2)
    dt_sum = dt_sum + (length(i)/Ax(i)) * (1/dt)
    friction_sum = friction_sum + my_f * my_ld * (1/Ax(i))**2
    loss_sum = loss_sum + my_loss
  enddo

  a2 = (1/gc)*dt_sum + 2*(friction_sum*(1/(2*rho_const*gc)) + loss_sum)*abs(m2_prev)

  !b2
  dt_sum =0.0
  loss_sum = 0.0
  friction_sum = 0.0
  boun_sum = 0.0
  do i=6, 15
    call friction_factor(my_f, m2_prev, i, 2)
    call get_LD(my_ld, i)
    call get_loss(my_loss, i, m2_prev, 2)
    dt_sum = dt_sum + (length(i)/Ax(i)) * (1/dt)
    friction_sum = friction_sum + my_f * my_ld * (1/Ax(i))**2
    loss_sum = loss_sum + my_loss
    boun_sum = boun_sum + (rho2(i)*dH(i))
  enddo

  b2 = (1/gc)*dt_sum*m2_in + (friction_sum*(1/(2*rho_const*gc)) + loss_sum) * m2_prev*abs(m2_prev) &
        - boun_sum + dP_pump2

  !ac
  dt_sum =0.0
  loss_sum = 0.0
  friction_sum = 0.0
  do i=1, 5
    call friction_factor(my_f, mc_prev, i, 2)
    call get_LD(my_ld, i)
    call get_loss(my_loss, i, mc_prev, 2)
    dt_sum = dt_sum + (length(i)/Ax(i)) * (1/dt)
    friction_sum = friction_sum + my_f * my_ld * (1/Ax(i))**2
    loss_sum = loss_sum + my_loss
  enddo
  i=16
  call friction_factor(my_f, mc_prev, i, 2)
  call get_LD(my_ld, i)
  call get_loss(my_loss, i, mc_prev, 2)
  dt_sum = dt_sum + (length(i)/Ax(i)) * (1/dt)
  friction_sum = friction_sum + my_f * my_ld * (1/Ax(i))**2
  loss_sum = loss_sum + my_loss

  ac = (1/gc)*dt_sum + 2*(friction_sum*(1/(2*rho_const*gc)) + loss_sum)*abs(mc_prev)

  !bc
  dt_sum =0.0
  loss_sum = 0.0
  friction_sum = 0.0
  boun_sum = 0.0
  do i=1, 5
    call friction_factor(my_f, mc_prev, i, 2)
    call get_LD(my_ld, i)
    call get_loss(my_loss, i, mc_prev, 2)
    dt_sum = dt_sum + (length(i)/Ax(i)) * (1/dt)
    friction_sum = friction_sum + my_f * my_ld * (1/Ax(i))**2
    loss_sum = loss_sum + my_loss
    boun_sum = boun_sum + (rho1(i)*dH(i))
  enddo
  i=16
  call friction_factor(my_f, mc_prev, i, 2)
  call get_LD(my_ld, i)
  call get_loss(my_loss, i, mc_prev, 2)
  dt_sum = dt_sum + (length(i)/Ax(i)) * (1/dt)
  friction_sum = friction_sum + my_f * my_ld * (1/Ax(i))**2
  loss_sum = loss_sum + my_loss
  boun_sum = boun_sum + (rho1(i)*dH(i))

  bc = (1/gc)*dt_sum*mc_in + (friction_sum*(1/(2*rho_const*gc)) + loss_sum) * mc_prev*abs(mc_prev) &
        - boun_sum

  !Solve for m1, m2, mc, dP_core
  m1 = (b1*a2 + 3*b1*ac - 3*ac*b2 + a2*bc) / (a1*a2 + 3*a1*ac +a2*ac)
  dP = b1 - a1*m1
  m2 = (b2 - dP) / a2
  mc = m1 + 3*m2

  conv = 0.0
  conv = max(conv, abs(m1/m1_prev - 1), abs(m2/m2_prev - 1), abs(mc/mc_prev - 1), abs(dP/dP_prev - 1))

  if (conv .lt. 0.01) then
    m1_in = m1
    m2_in = m2
    mc_in = mc
    dP_in = dP
    exit
  endif
  if (ite .gt. max_ite) then
    !throw some exit condition
    error = -1
    exit
  endif

  m1_prev = m1
  m2_prev = m2
  mc_prev = mc
  dP_prev = dP
end do

end subroutine mass_momentum

subroutine set_dP_pumps(my_m1, my_m2, time)
use data
implicit none

real(8), intent(in):: my_m1, my_m2, time

real(8) :: beta

beta = 0.425

if (time .eq. 0.0) then
  dP_pump1 = dP_pump_rated*(1.094 + 0.089*(my_m1/m_rated) - 0.183*(my_m1/m_rated)**2)
  dP_pump2 = dP_pump_rated*(1.094 + 0.089*(my_m2/m_rated) - 0.183*(my_m2/m_rated)**2)
else
  dP_pump1 = dP_pump_rated*(1.094 + 0.089*(my_m1/m_rated) - 0.183*(my_m1/m_rated)**2)
  dP_pump2 = dP_pump_rated*(1.094 + 0.089*(my_m2/m_rated) - 0.183*(my_m2/m_rated)**2)

  !pump trip
  dP_pump1 = dP_0 /(1 + time/beta)
  dP_pump2 = dP_pump1

  !locked rotor / broken shaft
  !dP_pump1 = 0.0
endif

end subroutine set_dP_pumps
