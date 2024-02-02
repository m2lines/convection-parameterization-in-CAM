# Set compiler and flags - uncomment for ifort or gfortran as appropriate
# Set free length and include netcdf-fortran
# Try using nf-config tool if available, otherwise need to manually add path as per examples below

# ifort configuration
FC = ifort
FCFLAGS = $(shell nf-config --fflags)

# gfortran configuration
# FC = gfortran
# FCFLAGS = -g -ffree-line-length-none $(shell nf-config --fflags)

# link flags
# Link to netcdf and netcdf-fortran
# Try using nf-config tool if available, otherwise manually add path as per examples below
LDFLAGS = $(shell nf-config --flibs)

# If flags need adding manually this can be done here
# CSD3:
# FCFLAGS = -g -ffree-line-length-none -I/usr/local/software/spack/spack-views/rocky8-icelake-20220710/netcdf-fortran-4.5.4/gcc-11.3.0/intel-oneapi-mpi-2021.6.0/mscqvjsc7bwypshmhnqfc2u3zxnims3r/include
# LDFLAGS = -L/usr/local/software/spack/spack-views/rhel8-icelake-20211027_2/netcdf-fortran-4.5.3/gcc-11.2.0/intel-oneapi-mpi-2021.4.0/g4qjb23rucofcg5uitt4jwrkgyf7gba7/lib -lnetcdff -lnetcdf -lm
# jwa34 local:
# FCFLAGS = -g -ffree-line-length-none -I/opt/netcdf-fortran/include
# LDFLAGS = -L/opt/netcdf-c/lib -lnetcdf -L/opt/netcdf-fortran/lib -lnetcdff -lm


PROGRAM = test
SRC = precision.f90 test_utils.f90 SAM_consts.f90 nn_cf_net.f90 nn_convection_flux.f90 nn_interface_SAM.f90 test.f90
# SRC = SAM_consts.f90 nn_cf_net.f90 nn_convection_flux.f90 nn_interface_CAM.f90 test_cam_interface.f90
# PROGRAM = test_cam
# SRC = nn_cf_net.f90 nn_convection_flux.f90 nn_interface_SAM.f90 test.f90 nn_interface_CAM.f90 test_cam_interface.f90
OBJECTS = $(SRC:.f90=.o)


all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
		$(FC) -o $@ $^ $(LDFLAGS)

%.o: %.f90
		$(FC) $(FCFLAGS) -c $<

clean:
		rm -f *.o *.mod

