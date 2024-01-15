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
      yog_init, & ! register fields in physics buffer
      convect_deep_tend_yog ! 
      

! Physics buffer indices 
   integer     ::  icwmrdp_idx      = 0 
   integer     ::  cldtop_idx       = 0 
   integer     ::  cldbot_idx       = 0 

   integer     ::  pblh_idx        = 0 
   integer     ::  tpert_idx       = 0 

   integer     ::  ttend_dp_idx        = 0

   integer     ::  yog_var_idx        = 0
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

subroutine yog_init(nn_weights, SAM_sounding)

!----------------------------------------
! Purpose: declare output fields, initialise variables needed for YOG
!----------------------------------------
! Used for fields that will be written out to 
! Use addfld to add to master list as an option to be written out, use
! add_default to force writing out.
  
  use physics_buffer, only : pbuf_get_index
  use cam_history,    only: addfld, horiz_only
  use nn_interface_CAM, only: nn_convection_flux_CAM_init


  implicit none

  character(len=1024) :: nn_weights    ! location of weights for the YOG NN, set in namelist
  character(len=1024) :: SAM_sounding  ! location of SAM sounding profile for the YOG NN, set in namelist

  call addfld('PRECYOG',  horiz_only,  'A', 'ms', 'total precipitation from YOG moist convection')
  call addfld('YOGDQ',    (/ 'lev' /), 'A', 'kg/kg/s', 'Q tendency - YOG moist convection')
  call addfld('YOGDT',    (/ 'lev' /), 'A', 'K/s','T tendency - YOG moist convection')

  call nn_convection_flux_CAM_init(nn_weights, SAM_sounding)

end subroutine yog_init

!=========================================================================================

subroutine convect_deep_tend_yog( &
     mcon    ,cme     ,          &
     pflx    ,zdu      , &
     rliq    ,rice     , &
     ztodt   , &
     state   ,ptend   ,landfrac ,pbuf)


   use physics_types, only: physics_state, physics_ptend, physics_tend, physics_ptend_init
   use phys_control,  only: use_gw_convect_dp
   
   use cam_history,    only: outfld
   use constituents,   only: pcnst
   use zm_conv_intr,   only: zm_conv_tend
   use cam_history,    only: outfld
   use physconst,      only: cpair
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index

! Arguments
   type(physics_state), intent(in ) :: state   ! Physics state variables
   type(physics_ptend), intent(out) :: ptend   ! individual parameterization tendencies
   

   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8), intent(in) :: ztodt               ! 2 delta t (model time increment)
   real(r8), intent(in) :: landfrac(pcols)     ! Land fraction
      

   real(r8), intent(out) :: mcon(pcols,pverp)  ! Convective mass flux--m sub c
   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
   real(r8), intent(out) :: cme(pcols,pver)    ! cmf condensation - evaporation
   real(r8), intent(out) :: zdu(pcols,pver)    ! detraining mass flux

   real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals
   real(r8), intent(out) :: rice(pcols) ! reserved ice (not yet in cldice) for energy integrals

   real(r8), pointer :: prec(:)   ! total precipitation
   real(r8), pointer :: snow(:)   ! snow from ZM convection 

   real(r8), pointer, dimension(:) :: jctop
   real(r8), pointer, dimension(:) :: jcbot
   real(r8), pointer, dimension(:,:,:) :: cld        
   real(r8), pointer, dimension(:,:) :: ql        ! wg grid slice of cloud liquid water.
   real(r8), pointer, dimension(:,:) :: rprd      ! rain production rate
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble

   real(r8), pointer, dimension(:,:) :: evapcdp   ! Evaporation of deep convective precipitation

   real(r8), pointer :: pblh(:)                ! Planetary boundary layer height
   real(r8), pointer :: tpert(:)               ! Thermal temperature excess

   ! Temperature tendency from deep convection (pbuf pointer).
   real(r8), pointer, dimension(:,:) :: ttend_dp

   real(r8) zero(pcols, pver)

   integer i, k

   cldtop_idx = pbuf_get_index('CLDTOP')
   cldbot_idx = pbuf_get_index('CLDBOT')
   icwmrdp_idx     = pbuf_get_index('ICWMRDP')
   pblh_idx   = pbuf_get_index('pblh')
   tpert_idx  = pbuf_get_index('tpert')
   call pbuf_get_field(pbuf, cldtop_idx,  jctop )
   call pbuf_get_field(pbuf, cldbot_idx,  jcbot )
   call pbuf_get_field(pbuf, icwmrdp_idx, ql    )

   call pbuf_get_field(pbuf, pblh_idx,  pblh)
   call pbuf_get_field(pbuf, tpert_idx, tpert)

   call zm_conv_tend( pblh    ,mcon    ,cme     , &
        tpert   ,pflx    ,zdu      , &
        rliq    ,rice    , &
        ztodt   , &
        jctop, jcbot , &
        state   ,ptend   ,landfrac, pbuf)


  ! If we added temperature tendency to pbuf, set it now.

  if (use_gw_convect_dp) then
  ttend_dp_idx  = pbuf_get_index('TTEND_DP')
  if (ttend_dp_idx > 0) then
     call pbuf_get_field(pbuf, ttend_dp_idx, ttend_dp)
     ttend_dp(:state%ncol,:pver) = ptend%s(:state%ncol,:pver)/cpair
  end if
  end if

  call outfld( 'ICWMRDP ', ql  , pcols, state%lchnk )

end subroutine convect_deep_tend_yog
!=========================================================================================

end module yog_mod
