FC = gfortran
FFLAGS = -O3


all_objects =  \
	globalconstants.o \
	data.o \
	graphics.o \
	lagstep.o \
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
	rm chocobof
	rm *.o
	rm *.mod

prog.o: data.o globalconstants.o  sod.o lagstep.o graphics.o

lagstep.o: data.o globalconstants.o  sod.o

graphics.o: data.o globalconstants.o  lagstep.o write_tio.o write_silo.o


sod.o: data.o globalconstants.o

globalconstants.o: globalconstants.f95

write_tio.o: data.o lagstep.o

write_silo.o: data.o lagstep.o
