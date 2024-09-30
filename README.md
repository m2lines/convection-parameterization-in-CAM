# Convection Parameterization in CAM

## Description
This repository contains code as part of an effort to deploy machine learning (ML) models of geophysical parameterisations into the [Community Earth System Model (CESM)](https://www.cesm.ucar.edu/).
This work is part of the [M<sup>2</sup>LInES](https://m2lines.github.io/) project aiming to improve performance of climate models using ML models for subgrid parameterizations.

A Neural Net providing a subgrid parameterization of atmospheric convection in a [single column model](https://www.arm.gov/publications/proceedings/conf04/extended_abs/randall_da.pdf) has been developed and successfully deployed as part of an atmospheric simulation.
The work is described in a [GRL paper](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020GL091363) with [accompanying code available](https://github.com/yaniyuval/Neural_nework_parameterization/tree/v.1.0.3). The repository contains the neural net and its implementation into a simple system for atmospheric modelling, [SAM](http://rossby.msrc.sunysb.edu/~marat/SAM.html).

The aims of this repository are to:
1. provide a standalone fortran module of the convection parameterisation based on this neural net that can be used elsewhere,
2. deploy the module in another atmospheric model, and
3. evaluate its performance (and provide instructions/code to do this).

We have also started to perform an investigation into interfacing the pytorch implementation of the Neural Net using the [pytorch-fortran bridging code](https://github.com/Cambridge-ICCS/fortran-pytorch-lib) developed at the [Institute of Computing for Climate Science](https://cambridge-iccs.github.io/).

The model was first deployed into the [Single Column Atmospheric Model (SCAM)](https://www.cesm.ucar.edu/models/simple/scam) - a single column version of the CESM.
We evaluate performance using Single Column CAM (SCAM) in the gateIII configuration for tropical convection in a similar manner described by the [SCAM6 pulication in JAMES](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2018MS001578), as well as the TOGA II configuration described in the same paper.
This compares model performance to data from an intense observation period (IOP) described in an [AMS publication](https://journals.ametsoc.org/view/journals/atsc/36/1/1520-0469_1979_036_0053_saposs_2_0_co_2.xml).

Long term developments of this project will seek to re-deploy more complex ML parameterizations into mode complex atmospheric models such as the [Community Atmospheric Model (CAM)](https://www.cesm.ucar.edu/models/cam) part of the CESM.

## Repository structure

```
├── YOG_convection
│   └── ...
├── torch_nets
│   └── ...
├── CAM_interface
│   └── ...
├── profile_test
│   └── ...
└── tests
    └── ...

```

### Contents

### `YOG_convection/`

This code is a standalone fortran convection parameterisation consisting of a fortran neural net and accompanying physics code and any dependencies extracted from the [original codebase](https://github.com/yaniyuval/Neural_nework_parameterization/tree/v.1.0.3).

It can be integrated as a standalone fortran convection parameterisation in a larger atmospheric model.

It has dependencies on the NetCDF and NetCDF-Fortran libraries and can be built using the included makefile (may need modification).

Files:

- `nn_cf_net.f90` - The core neural net routines
- `nn_convection_flux.f90` - wrapper around the neural net routines to perform physics-based operations
- `test.f90` - simple smoke tests for parameterisation routines
- `NN_weights_YOG_convection.nc` - NetCDF file with weights for the neural net
- Makefile - Makefile to compile these files

Note that you will need to generate the SAM sounding (an [atmospheric sounding](https://en.wikipedia.org/wiki/Atmospheric_sounding) that defines the SAM grid that the parameterisation is operating on --
we need this grid for the parameterisation, but also to interpolate any inputs from a different model onto this grid so that they can be used in the neural net) as a NetCDF file from the data files in `YOG_convection/rresources/` if you are using the CAM interface.\
This can be done from within the `YOG_convection/resources/` directory as follows:
```
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python sounding_to_netcdf.py
deactivate

```
, with a Python version >= 3.7 and < 3.11.

There are test routines associated with this code in `/tests/test_YOG_convection/`.
Guidance on running these can be found below.


### ``torch_nets/``
This directory contains the PyTorch versions of the neural networks used in the YOG convection parameterisation.

### ``CAM_interface/``
The directory contains the additional files or details to interface the YOG code with the CAM atmospheric model as part of the CESM model suite.

An implementation of this code within CAM, following the process outlined below, is available in the [m2lines/CAM-ML](https://github.com/m2lines/CAM-ML/tree/CAM-ML) repository.
This contains a full description of how to install and run the code as part of the CESM model suite.

That specific implementation is run within CESM v2.1.5 and is based of a version of CAM with the `cam_cesm2_1_rel_60` release tag.
There are notable future changes to the CAM source for CESM v2.2 and v3.0 that may
require revisions to these files in future.

In addition to the files in [/YOG_convection](/YOG_convection) containing the parameterisation,
this files here are required:

- `nn_interface_CAM.F90` - The interface for performing conversion from CAM variables and grid into the variables/grid expected by the YOG parameterisation.
- `yog_intr.F90` - The interface between the CAM model and the YOG parameterisation.

Changes also have to be made to:

- `bld/namelist_files/namelist_definition.xml` - To add the new namelist variables.
- `bld/namelist_files/namelist_defaults_cam.xml` - To set the new namelist variable defaults.
- `/src/physics/cam/physcontrol.F90` - To read in the new namelist variables and broadcast them.
- `/src/physics/cam/physpkg.F90` - To call the new parameterisation routines detailed in `yog_intr.F90`.

Full details of these can be seen in the specific implementation in CAM linked above.

There are test routines associated with this code in `/tests/test_CAM_interface/`.
Guidance on running these can be found below.

### ``profile_test/``

The CAM Profile Test code is a standalone test rig around the YOG convection parameterisation routines for debugging.

The [README in the subdirectory](https://github.com/m2lines/convection-parameterization-in-CAM/blob/main/profile_test/README.md) contains detailed information how to run these tests.


### ``tests/``

There are some tests for the NN_module and interface code in the `test/`
subdirectory.
These require a Fortran compiler, `netcdf` and
`netcdf-fortran` libraries, and CMake to build.

They can be built and run with CMake using the following commands:
```
cd tests
mkdir build
cd build
cmake ..
cmake --build .
```
This will create executables in the `build/` subdirectory which can be
run to execute the tests using:
```
./test_YOG_convection
./test_CAM_interface
```
with output printed to the console.

Note: To build only a specific subset of tests instead of all of them use:
```
cmake --build . --target <name_of_test_executable>
```


## Contributing

Contributions to the repository are welcome, particularly from anyone seeking to implement the parameterisation in another atmospheric model.
We welcome the addition of details for interfacing to other models, including code, from anyone who has used the parameterisation in another model.

Open tickets can be viewed at ['Issues'](https://github.com/m2lines/convection-parameterization-in-CAM/issues).

To contribute find a relevant issue or open a new one and assign yourself to work on it.
Then create a branch in which to add your contribution and open a pull request.
Once ready assign a reviewer and request a code review.
Merging should _only_ be performed once a code review has taken place.
