# compiler
FC = gfortran

# compile flags
FCFLAGS = -g -ffree-line-length-none
# -fdefault-real-8 -fbacktrace -fno-align-commons -fbounds-check 

# link flags
LDFLAGS = -I/usr/local/software/spack/spack-views/rocky8-icelake-20220710/netcdf-fortran-4.5.4/gcc-11.3.0/intel-oneapi-mpi-2021.6.0/mscqvjsc7bwypshmhnqfc2u3zxnims3r/include

PROGRAM = test
SRC = domain.f90 grid.f90 params.f90 vars.f90 nn_convection_flux_mod.f90 test.f90
OBJECTS = $(SRC:.f90=.o)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
		$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
		$(FC) $(FCFLAGS) $(LDFLAGS) -c $<

clean:
		rm -f *.o *.mod

