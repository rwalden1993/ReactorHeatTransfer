subroutine fuel_heat_conduction()
use data
implicit none

integer :: i, j
real(8) :: rho_fuel, cp, my_q, T_fuel_new(4,fuel_nodes)

rho_fuel = 685.46
cp = 0.1195028

do j=1, 4
  do i=1, fuel_nodes
    k_fuel(j,i) =  3978.1/(692.6 + T_fuel(j,i)) + 6.02366e-12*((T_fuel(j,i) + 460)**3)
  enddo
enddo

do i=1,4
  my_q = ((Q/4*Fk(i))/(50952*Pi*rs**2*length(i)))

  j=1
  T_fuel_new(i,j) = T_fuel(i,j) + dt/(rho_fuel*cp*dx_fuel)*((k_fuel(i,j+1)-k_fuel(i,j))/dx_fuel* &
                    (T_fuel(i,j+1)-T_fuel(i,j))/dx_fuel &
                    + my_q*dx_fuel)
  do j=2, fuel_nodes-1
    T_fuel_new(i,j) = T_fuel(i,j) + dt/(rho_fuel*cp*dx_fuel)*((k_fuel(i,j+1)-k_fuel(i,j))/dx_fuel* &
          (T_fuel(i,j+1)-T_fuel(i,j))/dx_fuel &
        - (k_fuel(i,j)-k_fuel(i,j-1))/dx_fuel*(T_fuel(i,j)-T_fuel(i,j-1))/dx_fuel + my_q*dx_fuel)
  enddo
  j=fuel_nodes
  T_fuel_new(i,j) = T_fuel(i,j) + dt/(rho_fuel*cp*dx_fuel)* &
        (-(k_fuel(i,j)-k_fuel(i,j-1))/dx_fuel*(T_fuel(i,j)-T_fuel(i,j-1))/dx_fuel + my_q*dx_fuel)
enddo

T_fuel = T_fuel_new

end subroutine fuel_heat_conduction

subroutine clad_temp(t_co, node, m)
use data
implicit none

integer, intent(in) :: node
real(8), intent(in) :: m
real(8), intent(inout) :: t_co

real(8) :: Hg, k_c, ro, ri, my_q, Ax_c, De_c, D, S, Re, Nu, hc

Hg = 1000.0
k_c = 9.6
ro = 0.374/24
ri = ro - 0.0225/12
!rs = 0.3225/24

Ax_c = 52.74
De_c = 0.4635/12
D = 0.374/12
S = 0.496/12

Re = ((abs(m)/Ax_c) * De_c)/mu
Nu = (0.042*(S/D) - 0.024)*(Re**0.8)*(Pr**(1.0/3.0))
hc = (Nu * k_th)/De_c

!my_q = (((Q/4*Fk(node))/(50952*Pi*rs**2*length(node)))*rs**2)/(2*ro)
my_q = Q/4 * Fk(node)

t_co = my_q/(hc*Ph(node)*length(node)) + T2(node)
!t_co = T_fuel(node,1) + my_q*rs*(1/(Hg*ri) + 1/k_c*log(ro/ri))

end subroutine clad_temp
