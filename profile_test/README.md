# CAM Profile Tests

## Description

This code is a standalone test rig around the YOGconvection parameterisation routines for debugging.
It uses a series of input profiles from CAM to pass to and run the parameterisation before
returning predicted profiles.

The test-rig is contained in `test_cam_profile.f90`.
Profiles are output at relevant points in the code for analysis.

## Running

This code has dependencies on the NetCDF and NetCDF-Fortran libraries.
It requires a Fortran compiler and python.

Set up a virtual python environment and install any dependencies required for the code:
```
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Before running the code generate the SAM sounding profiles in the `YOG_convection/resources/`
directory as described in `YOG_convection/README.md`:

Then build and run the fortran code using CMake.\
```
mkdir build
cd build
cmake ..
cmake --build .

./CAM_profile
```

This will generate `NN_test_output.nc` in the `profile_test/` directory containing
column profiles for a various timesteps and intermetiate variables.

Finally, run the `plot_nn_profiles.py` python code from within the virtual environment
to produce comparison plots of the profiles at different stages in the parameterisation.
These will be output to screen using matplotlib.

```
python plot_nn_profiles.py
```
