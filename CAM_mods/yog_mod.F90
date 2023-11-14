module yog_mod
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to YOG deep convection. Currently includes:
!
!
! Author: J.W. Atkinson, Sep 2023
!
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8=>shr_kind_r8
   use ppgrid,       only: pver, pcols, pverp
   use cam_logfile,  only: iulog

   implicit none

   save
   private                         ! Make default type private to the module

! Public methods

   public ::&
      yog_register,                    &! register fields in physics buffer
      yog_init,                    &! register fields in physics buffer
   
   integer ::& ! indices for fields in the physics buffer
      yog_var_idx

!=========================================================================================
  contains 

!=========================================================================================
subroutine yog_register

!----------------------------------------
! Purpose: register fields with the physics buffer
!----------------------------------------
! Used for fields that persist across timesteps and written out to restart files
  
  use physics_buffer, only : pbuf_add_field, dtype_r8

  implicit none

  call pbuf_add_field('YOG_VAR',    'physpkg',dtype_r8,(/pcols,pver/),yog_var_idx)

end subroutine yog_register

!=========================================================================================

subroutine yog_init

!----------------------------------------
! Purpose: declare output fields, initialise variables needed for YOG
!----------------------------------------
! Used for fields that will be written out to 
! Use addfld to add to master list as an option to be written out, use
! add_default to force writing out.
  
  use physics_buffer, only : pbuf_get_index
  use cam_history,    only: addfld, horiz_only

  implicit none

  call addfld('PRECYOG',  horix_only,  'A', 'ms', 'total precipitation from YOG convection')
  call addfld('YOGDQ',    (/ 'lev' /), 'A', 'kg/kg/s', 'Q tendency - YOG convection')

end subroutine yog_init


end module yog_mod
