# compiler
FC = gfortran

# compile flags
# Set free length and include netcdf-fortran
# CSD3:
FCFLAGS = -g3 -ffree-line-length-none -I/usr/local/software/spack/spack-views/rocky8-icelake-20220710/netcdf-fortran-4.5.4/gcc-11.3.0/intel-oneapi-mpi-2021.6.0/mscqvjsc7bwypshmhnqfc2u3zxnims3r/include
# jwa34 local:
# FCFLAGS = -g -ffree-line-length-none -I/opt/netcdf-fortran/include
# ifort options:
# -fdefault-real-8 -fbacktrace -fno-align-commons -fbounds-check 

# link flags
# Link to netcdf and netcdf-fortran
# CSD3:
LDFLAGS = -L/usr/local/software/spack/spack-views/rhel8-icelake-20211027_2/netcdf-fortran-4.5.3/gcc-11.2.0/intel-oneapi-mpi-2021.4.0/g4qjb23rucofcg5uitt4jwrkgyf7gba7/lib -lnetcdff -lnetcdf -lm
# jwa34 local:
# LDFLAGS = -L/opt/netcdf-c/lib -lnetcdf -L/opt/netcdf-fortran/lib -lnetcdff -lm


PROGRAM = test
SRC = nn_cf_net.f90 nn_convection_flux.f90 nn_interface_SAM.f90 test.f90
OBJECTS = $(SRC:.f90=.o)


all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
		$(FC) -o $@ $^ $(LDFLAGS)

%.o: %.f90
		$(FC) $(FCFLAGS) -c $<

clean:
		rm -f *.o *.mod

