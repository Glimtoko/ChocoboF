FC = gfortran
FFLAGS = -O3


all_objects =  \
	geom.o \
	globalconstants.o \
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

gascode: $(all_objects)
	$(FC) $(FFLAGS) $(all_objects) $(LINKS) $(RPATH) -o gascode

clean:
	rm gascode
	rm *.o
	rm *.mod

prog.o: globalconstants.o geom.o sod.o lagstep.o graphics.o

lagstep.o: globalconstants.o geom.o sod.o

graphics.o: globalconstants.o geom.o lagstep.o write_tio.o write_silo.o

geom.o: globalconstants.o

sod.o: globalconstants.o geom.o

globalconstants.o: globalconstants.f95

write_tio.o: lagstep.o geom.o

write_silo.o: lagstep.o geom.o
