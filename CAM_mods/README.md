# CAM mods

This directory contains the modifications to the CESM/CAM source code that need to be made in order to run the NN parameterisation in CAM.

They should be placed in 
```
my_cesm_sandbox/cime/scripts/<testcase>/SourceMods/src.cam/
```
where they will be picked up during the build process for compilation and override any CAM source files with the same name.

The modifications here consist of:

- Adding a select case structure to

- `physpkg.F90` - modified to provide a select case structure for `'YOG'` to call the parameterisation scheme.
- `phys_control.F90` - modified to read filenames from namelist and broadcast.
- `convect_deep.F90` - modified to account for `'YOG'` being an option for `convect_deep` but still call ZM.

## Building and Running in CAM
After you [create a case](https://ncar.github.io/CAM/doc/build/html/CAM6.0_users_guide/building-and-running-cam.html) in CESM you need to follow the following instructions:

1. Edit `user_nl_cam`  
   This is a CAM namelist generated from `my_cesm_sandbox/components/cam/cime_config/usermods_dirs/scam_gateIII/user_nl_cam`  
   We need to add the following lines:
   1. `deep_scheme = 'YOG'`  
      This will be the identifier for our new convection scheme.
   2. `nn_weights = '<PATH/TO/WEIGHTS.nc>'`  
      The path to the nn weights.
   3. `SAM_sounding = '<PATH/TO/SAM/SOUNDING.nc>'`  
      The path to the SAM sounding for the NN.  
      This file is generated using the `sounding_to_netcdf.py` script in the resources of the NN code.
2. The latter two are new variables so we need to add the following lines to `my_cesm_sandbox/components/cam/bld/namelist_files/namelist_defaults_cam.xml`:
   ```
   <nn_weights                      >NONE     </nn_weights>
   <SAM_sounding                    >NONE     </SAM_sounding>
   ```

3. The main CESM code checks the vailidity of options for various arguments at build time, so we need to modify the CESM source file `my_cesm_sandbox/components/cam/bld/namelist_files/namelist_definition.xml` to add:
   1. `YOG` as a possible option for `deep_scheme`.
   2. `nn_weights` as a `char*132`
   3. `SAM_sounding` as a `char*132`
3. We can now run `./case.setup` and `./case.build`.

An example `user_nl_cam` is provided in this repo.

**Note:**  
By default CESM will place output in `/glade/scratch/user/case/`
and logs/restart files in `/glade/scratch/user/archive/case/`.
To place all output with logs in `archive/case` switch 'short term archiving' on by
editing `env_run.xml` in the case directory to change `DOUT_S` from `FALSE` to `TRUE`.

## Details of the changes

### `physpkg.F90`

- Add `deep_scheme` as a module variable and read in during the `phys_register` subroutine through calling `phys_getopts()`
- Adapt `tphysbc` subroutine which applies physical atmospheric processes in CAM before coupling to other models.  
  Add a `select_case()` structure so that `convect_deep_tend` is now called for `'ZM'`, with new calls to our routines for `'YOG'`.

### `phys_control.F90`

- Add `nn_weights` and `SAM_sounding` as namelist variables.

### `Test case`

On cheyenne we have been running a case generated from:
```
./create_newcase --case test_scam_gateIII_clean --compset FSCAM --res T42_T42 --user-mods-dir ../../components/cam/cime_config/usermods_dirs/scam_gateIII --project NCGD0054
```


## Questions:

- Useful link: [adding variables to output](https://ncar.github.io/CESM-Tutorial/README.html)
- NB: `zm_conv_intr.F90` has details of the ZM tendencies that can be added to output.
- Useful link: [CAM phys docs](https://www2.cesm.ucar.edu/models/atm-cam/docs/phys-interface/)

- There is a call to `zm_convect_deep_tend_2` in `physpkg.F90` called when `if ( .not. deep_scheme_does_scav_trans() ) then` But the call has different arguments, so how should YOG handle this?

! jwa34 - this writes to outfield - i.e. diaganostics to file
call outfld('CAPE', cape, pcols, lchnk)        ! RBN - CAPE output
