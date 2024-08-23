# CAM Interface

This directory contains/details the modifications to the CAM source code that need to be
made in order to run the NN parameterisation in CAM.

An implementation of this code within CAM, following the process outlined below, is 
available in the [m2lines/CAM-ML](https://github.com/m2lines/CAM-ML/tree/CAM-ML)
repository.
This contains a full description of how to install and run the code as part of the CESM 
model suite.

That specific implementation is run within CESM v2.1.5 and is based of a version of CAM
with the `cam_cesm2_1_rel_60` release tag.
There are notable future changes to the CAM source for CESM v2.2 and v3.0 that may
require revisions to these files in future.

In addition to the files in `[/NN_Module](/NN_Module)` containing the parameterisation,
the additional files here are required:

- `yog_intr.F90` - The interface between the CAM model and the YOG parameterisation.

Changes also have to be made to:

- `bld/namelist_files/namelist_definition.xml` - To add the new namelist variables.
- `bld/namelist_files/namelist_defaults_cam.xml` - To set the new namelist variable defaults.
- `/src/physics/cam/physcontrol.F90` - To read in the new namelist variables and broadcast them.
- `/src/physics/cam/physpkg.F90` - To call the new parameteristion routines detailed in `yog_intr.F90`.

Full details of these can be seen in the specific implementation in CAM linked above.
