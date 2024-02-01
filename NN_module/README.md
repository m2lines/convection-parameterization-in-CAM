# NN_module

This code is a standalone fortran module containing the fortran neural net and any dependencies extracted from the [original codebase](https://github.com/yaniyuval/Neural_nework_parameterization/tree/v.1.0.3).

It has dependencies on the NetCDF and NetCDF-Fortran libraries.

It can be built using the included makefile (may need modification).

Files:

- `nn_cf_net.f90` - The core neural net routines
- `nn_convection_flux.f90` - wrapper around the neural net routines to perform physics-based operations
- `nn_interface_SAM.f90` - wrapper around the parameterisation routines to interface with the SAM model
- `test.f90` - simple smoke tests for parameterisation routines
- `NN_weights_YOG_convection.nc` - NetCDF file with weights for the neural net
- Makefile - Makefile to compile these files

## Running tests

To run the tests (requires `ifort` and `netcdf`):

    make test