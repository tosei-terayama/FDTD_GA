#Makefile
OBJS = main.o fdtd_calc.o memory_allocate2d.o memory_allocate2cd.o memory_allocate3d.o \
	memory_allocate3cd.o memory_allocate4d.o memory_allocate5d.o \
	delete_2d.o delete_2cd.o delete_3d.o delete_3cd.o delete_4d.o delete_5d.o delete_PML.o\
	sigma_calc.o D_update.o D_update_pml.o E_update.o H_update.o H_update_pml.o \
	pml_class.o ny_allocate.o Ne_allocate.o geomagnetic.o surface_impe_calc.o \
	surface_H_update.o PML_field_initialize.o PML_idx_initialize.o set_matrix.o \
	set_perturbation.o geocoordinate_class.o perturbation_class.o date_class.o \
	output_profile.o output_model.o obs_initial.o \
	GA.o GA_agent.o set_parameter.o calc_score.o

HEADERS = fdtd3d.h GA_agent.h pml.h geocoordinate.h perturbation.h date.h

OPTS = -I/opt/include/eigen3 -I/opt/include -L/opt/lib -lnrlmsise -std=c++1z -O3 -Wall

main: $(OBJS)
	mpic++ -o $@ $(OBJS) $(OPTS)
%.o: %.cpp $(HEADERS)
	mpic++ -c $< $(OPTS)
