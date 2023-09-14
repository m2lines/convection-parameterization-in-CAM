# CAM mods

This directory contains the modifications to the CESM/CAM source code that need to be made in order to run the NN parameterisation in CAM.

They should be placed in 
```
my_cesm_sandbox/cime/scripts/<testcase>/SourceMods/src.cam/
```
where they will be picked up during the build process for compilation and override any CAM source files with the same name.

The modifications here consist of:

- Adding a select case structure to

- `yog_intr.F90` - based on `zm_intr.F90`, interface file for parameterisation scheme.
- `yog_deep.F90` - based on `zm_deep.F90`, contains code for the parameterisation scheme.
- `physpkg.F90` - modified to provide a select case structure for `'YOG'` to camm the parameterisation scheme.


