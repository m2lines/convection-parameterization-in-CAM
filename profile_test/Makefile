# Set compiler and flags - uncomment for ifort or gfortran as appropriate
# Set free length and include netcdf-fortran
# Try using nf-config tool if available, otherwise need to manually add path as per examples below

# ifort configuration
# FC = ifort
# FCFLAGS = $(shell nf-config --fflags)

# gfortran configuration
FC = gfortran
FCFLAGS = -g -ffree-line-length-none $(shell nf-config --fflags)

# link flags
# Link to netcdf and netcdf-fortran
# Try using nf-config tool if available, otherwise manually add path as per examples below
LDFLAGS = $(shell nf-config --flibs)

# If flags need adding manually this can be done here
# FCFLAGS = -g -I/opt/homebrew/Cellar/netcdf-fortran/4.6.1/include -I/opt/homebrew/Cellar/netcdf-fortran/4.6.1/include -I/opt/homebrew/Cellar/netcdf/4.9.2_1/include
# LDFLAGS = -L/opt/homebrew/Cellar/netcdf-fortran/4.6.1/lib -L/opt/homebrew/Cellar/netcdf/4.9.2_1/lib -lnetcdff -lnetcdf

# CAM PROFILE TESTS
PROGRAM = cam_profile
SRC = SAM_consts.f90 write_nf_data.f90 nn_cf_net.f90 nn_convection_flux.f90 nn_interface_CAM.f90 cam_profile.f90 test_cam_profile.f90

OBJECTS = $(SRC:.f90=.o)


all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
		$(FC) -o $@ $^ $(LDFLAGS)

%.o: %.f90
		$(FC) $(FCFLAGS) -c $<

clean:
		rm -f *.o *.mod

