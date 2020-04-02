#Makefile
OBJS = fdtd3d.o memory_allocate2d.o memory_allocate3d.o memory_allocate3cd.o memory_allocate4d.o \
	sigma_calc.o D_update.o D_update_pml.o E_update.o H_update.o H_update_pml.o \
	pml_class.o sigma_real_calc.o make_rotate_mat.o sigma_car_calc.o E_update_iono.o coordinate_trans.o \
	Ne_allocate.o ny_allocate.o geomagnetic.o surface_impe_calc.o surface_H_update.o \
	PML_field_initialize.o PML_idx_initialize.o

HEADERS = fdtd3d.h pml.h nrlmsise-00.h
OPTS = -I/opt/include/eigen3 -std=c++1z -O3
LIBS = -L. -lnrlmsise

all: main libnrlmsise.a

.PHONY: all clean

main: $(OBJS) libnrlmsise.a
	g++ -o $@ $(OBJS) $(OPTS) $(LIBS)

%.o: %.cpp $(HEADERS)
	g++ -c $< $(OPTS)

%.o: %.c
	g++ -c $< $(OPTS)

LIBOBJS = nrlmsise-00.o nrlmsise-00_data.o
libnrlmsise.a: $(LIBOBJS)
	ar rcs $@ $(LIBOBJS)

clean:
	rm -rf *.o main *.dat
	