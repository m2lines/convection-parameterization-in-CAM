# compiler
FC      := gfortran

# compile flags
FCFLAGS = -g -c -fdefault-real-8 -fbacktrace -fno-align-commons -fbounds-check -ffree-line-length-none
# link flags
FLFLAGS = -I/usr/local/software/spack/spack-views/rhel8-icelake-20211027_2/netcdf-fortran-4.5.3/gcc-11.2.0/intel-oneapi-mpi-2021.4.0/g4qjb23rucofcg5uitt4jwrkgyf7gba7/include

domain.o: domain.f90
	$(FC) $(FCFLAGS) -o $@ $^

grid.o: grid.f90
	$(FC) $(FCFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<

%.mod: %.h
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -f *.o *.mod
