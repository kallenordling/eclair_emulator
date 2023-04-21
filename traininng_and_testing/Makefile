# should contain paths to NLOPT and GPF
MOD = ./obj   
LIB = ./lib   

F90FLAGS = -std=f2008 -Wall -Wextra -O3 $(addprefix -I,$(MOD))
F90LINKFLAGS =  $(addprefix -L,$(LIB)) -lgpf -lnlopt -llapack -lblas

gp_train_precip_day: gp_train_precip_day.f90
	$(FC) $(F90FLAGS) -o gp_train_precip_day gp_train_precip_day.f90 $(F90LINKFLAGS)

gp_train_precip_night: gp_train_precip_night.f90
	$(FC) $(F90FLAGS) -o gp_train_precip_night gp_train_precip_night.f90 $(F90LINKFLAGS)

gp_train_w_night: gp_train_w_night.f90
	$(FC) $(F90FLAGS) -o gp_train_w_night gp_train_w_night.f90 $(F90LINKFLAGS)


gp_train_w_day: gp_train_w_day.f90
	$(FC) $(F90FLAGS) -o gp_train_w_day gp_train_w_day.f90 $(F90LINKFLAGS)

all: gp_train_precip_day gp_train_precip_night gp_train_w_day gp_train_w_night test_emu

test_emu: test_emu.f90
	$(FC) $(F90FLAGS) -o test_emu test_emu.f90 $(F90LINKFLAGS)
