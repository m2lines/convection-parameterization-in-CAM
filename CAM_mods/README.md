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
- `physpkg.F90` - modified to provide a select case structure for `'YOG'` to call the parameterisation scheme.


## Building and Running in CAM
After you [create a case](https://ncar.github.io/CAM/doc/build/html/CAM6.0_users_guide/building-and-running-cam.html) in CESM you need to follow the following instructions:

1. Edit `user_nl_cam`  
   This is a CAM namelist generated from `my_cesm_sandbox/components/cam/cime_config/usermods_dirs/scam_gateIII/user_nl_cam`  
   We need to add a line with `deep_scheme = 'YOG'`  
   This will be the identifier for our new convection scheme.
2. The main CESM code checks the vailidity of options for various arguments at build time, so we need to modify the CESM source file `my_cesm_sandbox/components/cam/bld/namelist_files/namelist_definition.xml` to add `YOG` as a possible option for `deep_scheme`.
3. We can now run `./case.setup` and `./case.build`.

## Details of the changes

### `physpkg.F90`

- Add `deep_scheme` as a module variable and read in during the `phys_register` subroutine through calling `phys_getopts()`
- Adapt `tphysbc` subroutine which applies physical atmospheric processes in CAM before coupling to other models.  
  Add a `select_case()` structure so that `convect_deep_tend` is now called for `'ZM'`, with new calls to our routined for `'YOG'`.


## Questions:

- Where should we place NN weights?
