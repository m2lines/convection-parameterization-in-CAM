# CAM Profile Tests

## Description

This code is a standalone test rig around the YOG parameterisation routines for debugging.
It uses a series of input profiles from CAM to pass to and run the parameterisation before
returning predicted profiles.

The test-rig is contained in `test_cam_profile.f90`.
Profiles are output at relevant points in the code for analysis.

## Running

This code has dependencies on the NetCDF and NetCDF-Fortran libraries.
It requires a Fortran compiler and python.

First clone this code and checkout the relevant branch before moving to this directory:
```
git clone https://github.com/m2lines/convection-parameterization-in-CAM.git
cd convection-parameterization-in-CAM
git checkout cam-profile-test
cd profile_test
```

Set up a virtual python environment and install any dependencies required for the code:
```
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Before running the code generate the SAM sounding profiles in the `resources/` directory
as follows:
```
cd resources
python sounding_to_netcdf.py
cd ../
```

Then build and run the fortran code using the supplied Makefile.\
Note, you can edit the Makefile to change the Fortran compiler.
You may also need to edit it to specify the linking flags for the NetCDF libraries on your
system if you do not have the nc-config tools installed.
```
make
./cam_profile
```

This will generate `NN_test_output.nc` containing profiles for a single timestep.
If you wish to examine a different step in the time-series edit line 49 of `test_cam_profile.f90`,
re-build, and re-run as above.

Finally, run the `plot_nn_profiles.py` python code to produce comparison plots of the
profiles at different stages in the parameterisation.
These will be output to screen using matplotlib.

```
python plot_nn_profiles.py
```
