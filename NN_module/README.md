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

Note that you will need to generate the SAM sounding as a NetCDF file from the data
files in `resources/` if you are using the CAM interface.\
This can be done from within the `resources/` directory as follows:
```
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python sounding_to_netcdf.py
deactivate
```

## Running tests

There are some tests for the NN_module code in the test/ subdirectory.
These require `ifort` or `gfortran` Fortran compiler, `netcdf` and
`netcdf-fortran` libraries, and CMake.

They can be built and run with CMake using the following commands:
```
cd test
mkdir build
cd build
cmake ..
cmake --build .
```
This will create an executable `yogtest` in the test subdirectory which can be run to execute the tests using:
```
./yogtest
```
with output printed to the console.
