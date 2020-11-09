#Makefile
OBJS = D_update.o D_update_pml.o E_update.o H_update.o H_update_pml.o Ne_allocate.o Ne_generate.o \
		Ne_allocate.o PML_field_initialize.o PML_idx_initialize.o data_class.o delete_2cd.o delete_2d.o \
		delete_3cd.o delete_3d.o delete_4d.o delete_5d.o delete_PML.o fdtd_calc.o geocoordinate_class.o \
		geomagnetic.o main.o memory_allocate2cd.o memory_allocate2d.o memory_allocate3cd.o memory_allocate3d.o \
		memory_allocate4d.o memory_allocate5d.o ny_allocate.o perturbation_class.o pml_class.o set_matrix.o \
		set_perturbation.o sigma_calc.o surface_H_update.o surface_impe_calc.o

HEADERS = date.h fdtd3d.h geocoordinate.h nrlmsise-00.h perturbation.h pml.h

OPTS = -I/opt/include/eigen3 -std=c++1z -O3 -Wall
LIBS = -L. -lnrlmsise

all: main libnrlmsise.a
.PHONY: all clean

main: $(OBJS) libnrlmsise.a
	mpic++ -o $@ $(OBJS) $(OPTS) $(LIBS)

%.o: %.cpp $(HEADERS)
	mpic++ -c $< $(OPTS)

%.o: %.c
	g++ -c $< $(OPTS)

#%.o: %.for
#	gfortran -c $< -Wall -O3

LIBOBJS = nrlmsise-00.o nrlmsise-00_data.o
libnrlmsise.a: $(LIBOBJS)
	ar rcs $@ $(LIBOBJS)

clean:
	rm -rf *.o main *.dat
	
