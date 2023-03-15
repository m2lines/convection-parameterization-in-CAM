# Simple shell script to build a module on CSD3 icelakes

module load gcc/11.3.0/gcc/4zpip55j
module load netcdf-c/4.8.1/gcc/intel-oneapi-mpi/xhznpeda
module load netcdf-fortran/4.5.3/gcc/intel-oneapi-mpi/g4qjb23r

gfortran -ffree-line-length-none -c domain.f90
gfortran -ffree-line-length-none -c grid.f90
gfortran -ffree-line-length-none -c vars.f90
gfortran -ffree-line-length-none -c params.f90
gfortran -ffree-line-length-none -c -I/usr/local/software/spack/spack-views/rhel8-icelake-20211027_2/netcdf-fortran-4.5.3/gcc-11.2.0/intel-oneapi-mpi-2021.4.0/g4qjb23rucofcg5uitt4jwrkgyf7gba7/include nn_convection_flux.f90
