GCC = gfortran
FFLAGS = -Wall -fbacktrace -fbounds-check


all: final_project test_solve_2x2

final_project: main.o data.o mass_momentum.o friction_factor.o get_LD.o get_loss.o get_temp.o get_rho.o energy.o thomas_alg.o fuel_heat_conduction.o
	$(GCC) $(FFLAGS) main.o data.o mass_momentum.o friction_factor.o get_LD.o get_loss.o get_temp.o get_rho.o energy.o fuel_heat_conduction.o thomas_alg.o -o $@

main.o: main.f95 data.o mass_momentum.o
	$(GCC) $(FFLAGS) -c main.f95

data.o: data.f95
	$(GCC) $(FFLAGS) -c data.f95

mass_momentum.o: mass_momentum.f95 friction_factor.o get_LD.o data.o get_loss.o
	$(GCC) $(FFLAGS) -c mass_momentum.f95

friction_factor.o: friction_factor.f95 data.o
	$(GCC) $(FFLAGS) -c friction_factor.f95

get_LD.o: get_LD.f95 data.o
	$(GCC) $(FFLAGS) -c get_LD.f95

get_loss.o: get_loss.f95 data.o
	$(GCC) $(FFLAGS) -c get_loss.f95

get_temp.o: get_temp.f95 data.o
	$(GCC) $(FFLAGS) -c get_temp.f95

get_rho.o: get_rho.f95 data.o
	$(GCC) $(FFLAGS) -c get_rho.f95

energy.o: energy.f95 data.o thomas_alg.o
	$(GCC) $(FFLAGS) -c energy.f95

thomas_alg.o: thomas_alg.f95
	$(GCC) $(FFLAGS) -c thomas_alg.f95

test_solve_2x2: test_solve_2x2.f95 thomas_alg.o
	$(GCC) $(FFLAGS) test_solve_2x2.f95 thomas_alg.o -o $@

fuel_heat_conduction.o: fuel_heat_conduction.f95 data.o
	$(GCC) $(FFLAGS) -c fuel_heat_conduction.f95
