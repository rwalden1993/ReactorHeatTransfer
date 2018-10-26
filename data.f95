Module data
implicit none

integer, parameter :: nnodes = 16
real(8), parameter :: Pi = 4*atan(1.0_8)
real(8) :: length(nnodes), De(nnodes), Ax(nnodes), Ph(nnodes), dH(nnodes), Fk(4), Q
real(8) :: dt, gc, mu, rho_const, Pr, k_th, Tsat_sg, dx_fuel, rs, tubes_pluged
real(8) :: u1(nnodes), u2(nnodes), T1(nnodes), T2(nnodes), rho1(nnodes), rho2(nnodes)
real(8) :: dP_pump1, dP_pump2, dP_pump_rated, m_rated, Q_rated, dP_0
real(8), allocatable ::  T_fuel(:,:), k_fuel(:,:)
integer :: fuel_nodes

contains
  subroutine get_data()

    integer :: i

    !Constants
    gc = 32.2d0
    mu = 0.2078d0/3600
    rho_const = 44.4677d0
    Pr = 0.880968d0
    k_th = 0.317531d0
    Tsat_sg = 544.65389158152391 !544.65

    Q_rated = 11638813484.0/3600 !btu/s
    dP_pump_rated = 16517 !lbm/ft^2 ->   psi
    m_rated  = 9083 !lbm/s

    tubes_pluged = 1.0
    !tubes_pluged = 0.619177

    rs = 0.3225/24
    dx_fuel = 0.001/12
    fuel_nodes = int(rs/dx_fuel)
    dx_fuel = rs/fuel_nodes
    allocate(T_fuel(4,fuel_nodes))
    allocate(k_fuel(4,fuel_nodes))


    !Core nodes
    do i=1, 4
      length(i) = 3.0
      De(i)     = 0.4635
      Ax(i)     = 52.74
      Ph(i)     = 4988.86
      dH(i)     = 3.0
    enddo
    Fk(1) = 0.586
    Fk(2) = 1.414
    Fk(3) = Fk(2)
    Fk(4) = Fk(1)

    !Upper Plenum
    length(5) = 1.5
    De(5)     = 158
    Ax(5)     = 136.2
    Ph(5)     = 0.0
    dH(5)     = 1.5

    !Hot leg
    length(6) = 20.0
    De(6)     = 29.0
    Ax(6)     = 4.59
    Ph(6)     = 0.0
    dH(6)     = 0.0

    !Steam Gens
    do i=7, 13
      length(i) = 9.54
      De(i)     = 0.6075
      Ax(i)     = 13.35
      Ph(i)     = 1054.9
    enddo
    do i=7, 9
      dH(i)     = 9.54
    enddo
    dH(10)    = 0.0
    do i=11, 13
      dH(i)     = -9.54
    enddo

    !Cold Leg
    length(14) = 40.0
    De(14)     = 27.5
    Ax(14)     = 4.12
    Ph(14)     = 0.0
    dH(14)     = 0.0

    !Downcomer
    length(15) = 18.4
    De(15)     = 15.0
    Ax(15)     = 6.77
    Ph(15)     = 0.0
    dH(15)     = -13.5

    !Lower Plenum
    length(16) = 7.2
    De(16)     = 173
    Ax(16)     = 163.2
    Ph(16)     = 0.0
    dH(16)     = 0.0
  end subroutine get_data

end module data
