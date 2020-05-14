FC = gfortran
FFLAGS = -O3 -Wunused


all_objects =  \
	memory_manipulation.o \
	data.o \
	graphics.o \
	lagstep_fi.o \
	write_tio.o \
	write_silo.o \
	sod.o \
	prog.o

PROD = /usr/local/prod
INCLUDES = -I$(PROD)/typhonio/include -I$(PROD)/hdf5_oldapi/include -I$(PROD)/silo-4.10.2-bsd/include
LINKS = -L$(PROD)/typhonio/lib -ltyphonio_f -ltyphonio -L$(PROD)/silo-4.10.2-bsd/lib/ -lsiloh5 -L$(PROD)/hdf5_oldapi/lib -lhdf5_fortran -lhdf5
RPATH = -Wl,-rpath $(PROD)/hdf5/lib

.SUFFIXES:
.SUFFIXES: .o .f95

.f95.o:
	$(FC) $(FFLAGS) $(INCLUDES) $< -c

chocobof: $(all_objects)
	$(FC) $(FFLAGS) $(all_objects) $(LINKS) $(RPATH) -o chocobof

clean:
	-rm chocobof
	-rm *.o
	-rm *.mod


# Dependencies
prog.o: data.o memory_manipulation.o sod.o lagstep_fi.o graphics.o

sod.o: data.o	# For input and file unit

graphics.o: write_tio.o write_silo.o

