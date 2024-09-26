module yog_intr
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to the Yoval-O'Gorman deep convection scheme.
!   based off of the Zhang-McFarlane deep convection scheme
!
! Author: J.W. Atkinson
! January 2024
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8=>shr_kind_r8
   use physconst,    only: cpair                              
   use ppgrid,       only: pver, pcols, pverp, begchunk, endchunk
   use cam_abortutils,   only: endrun
   use physconst,        only: pi
   use spmd_utils,       only: masterproc
   use perf_mod
   use cam_logfile,  only: iulog
   use constituents, only: cnst_add
   
   implicit none
   private
   save

   ! Public methods

   public ::&
      yog_readnl,             &! read namelist for YOG scheme/module
      yog_init,               &! initialize YOG scheme/module
      yog_final,              &! finalize YOG scheme/module
      yog_tend                 ! return tendencies

   character(len=136) :: yog_nn_weights    ! location of weights for the YOG NN, set in namelist
   character(len=136) :: SAM_sounding  ! location of SAM sounding profile for the YOG NN, set in namelist
   
!=========================================================================================
contains
!=========================================================================================

! There is currently no need for a yog_register as nothing to add to the physics buffer
! subroutine yog_register
! 
! !----------------------------------------
! ! Purpose: register fields with the physics buffer
! !----------------------------------------
! 
!   use physics_buffer, only : pbuf_add_field, dtype_r8, dtype_i4
! 
!   implicit none
! 
! end subroutine yog_register

!=========================================================================================

subroutine yog_readnl(nlfile)

   use spmd_utils,      only: mpi_character, masterprocid, mpicom
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'yog_readnl'

   namelist /yog_params_nl/ yog_nn_weights, SAM_sounding
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'yog_params_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, yog_params_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

   end if

   ! Broadcast namelist variables
   call mpi_bcast(yog_nn_weights,   len(yog_nn_weights),   mpi_character, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("yog_readnl: FATAL: mpi_bcast: yog_nn_weights")
   call mpi_bcast(SAM_sounding, len(SAM_sounding), mpi_character, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("yog_readnl: FATAL: mpi_bcast: SAM_sounding")

end subroutine yog_readnl

!=========================================================================================

subroutine yog_init()

!----------------------------------------
! Purpose:  declare output fields, initialize variables needed by convection
!----------------------------------------

  use cam_history,    only: addfld, add_default, horiz_only
  use nn_interface_CAM,   only: nn_convection_flux_CAM_init
  
  use spmd_utils,     only: masterproc
  use phys_control,   only: phys_getopts

  implicit none

  integer istat
  logical :: history_budget ! output tendencies and state variables for
                            ! temperature, water vapor, cloud ice and cloud
                            ! liquid budgets.
  integer :: history_budget_histfile_num ! output history file number for budget fields

  ! Register fields with the output buffer
  call addfld ('YOGDT  ',   (/ 'lev' /),  'A', 'K/s','T tendency - Yuval-OGorman moist convection')
  call addfld ('YOGDQ  ',   (/ 'lev' /),  'A', 'kg/kg/s','Q tendency - Yuval-OGorman moist convection')
  call addfld ('YOGDICE',   (/ 'lev' /),  'A', 'kg/kg/s','Cloud ice tendency - Yuval-OGorman convection')
  call addfld ('YOGDLIQ',   (/ 'lev' /),  'A', 'kg/kg/s','Cloud liq tendency - Yuval-OGorman convection')
  call addfld ('YOGPREC',   horiz_only ,  'A', 'm/s','Surface preciptation - Yuval-OGorman convection')
  if (masterproc) then
     write(iulog,*)'YOG output fields added to buffer'
  end if

  call phys_getopts( history_budget_out = history_budget, &
                     history_budget_histfile_num_out = history_budget_histfile_num)

  if ( history_budget ) then
     call add_default('YOGDT  ', history_budget_histfile_num, ' ')
     call add_default('YOGDQ  ', history_budget_histfile_num, ' ')
     call add_default('YOGDICE', history_budget_histfile_num, ' ')
     call add_default('YOGDLIQ', history_budget_histfile_num, ' ')
     call add_default('YOGPREC', history_budget_histfile_num, ' ')
  end if

  call nn_convection_flux_CAM_init(yog_nn_weights, SAM_sounding)
  if (masterproc) then
     write(iulog,*)'yog_nn_weights at: ', yog_nn_weights
     write(iulog,*)'SAM_sounding at: ', SAM_sounding
     write(iulog,*)'YOG scheme initialised'
  endif

end subroutine yog_init

!=========================================================================================

subroutine yog_final()

!----------------------------------------
! Purpose:  finalization of the YOG scheme
!----------------------------------------

  use nn_interface_CAM,   only: nn_convection_flux_CAM_finalize
  
  implicit none

  call nn_convection_flux_CAM_finalize()

end subroutine yog_final

!=========================================================================================

subroutine yog_tend(ztodt, state, ptend)

!----------------------------------------
! Purpose:  tendency calculation for YOG scheme
!----------------------------------------

   use cam_history,   only: outfld
   use physics_types, only: physics_state, physics_ptend
   use physics_types, only: physics_ptend_init
   use constituents,  only: pcnst, cnst_get_ind

   use nn_interface_CAM,   only: nn_convection_flux_CAM

   ! Arguments

   type(physics_state), intent(in)          :: state          ! Physics state variables
   type(physics_ptend), intent(out)         :: ptend          ! individual parameterization tendencies

   real(r8), intent(in) :: ztodt                       ! 2 delta t (model time increment)

   ! Local variables

   integer :: i
   integer :: nstep
   integer :: ixcldice, ixcldliq      ! constituent indices for cloud liquid and ice water.
   integer :: lchnk                   ! chunk identifier
   integer :: ncol                    ! number of atmospheric columns

   real(r8) :: ftem(pcols,pver)       ! Temporary workspace for outfld variables

   logical  :: lq(pcnst)              ! Logical array corresponding to constituents present in state

   ! physics buffer fields

   real(r8), pointer, dimension(:)   :: prec         ! total precipitation
   real(r8), pointer, dimension(:)   :: snow         ! snow
   real(r8), pointer, dimension(:,:) :: cld

   real(r8) :: yog_precsfc(ncols)  ! scattered precip flux at each level

   lchnk = state%lchnk
   ncol  = state%ncol

   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   lq(:) = .true.
   call physics_ptend_init(ptend, state%psetcols, 'yogNN', ls=.true., lq=lq(:))! initialize ptend type for YOG
 
   call nn_convection_flux_CAM(state%pmid(:,pver:1:-1), state%pint(:,pverp:1:-1), state%ps, &
                               state%t(:,pver:1:-1), state%q(:,pver:1:-1,1), &
                               state%q(:,pver:1:-1,ixcldliq), state%q(:,pver:1:-1,ixcldice), &
                               cpair, &
                               ztodt, &
                               ncol, pver, &
                               yog_precsfc, &
                               ptend%q(:,pver:1:-1,ixcldice), ptend%q(:,pver:1:-1,1), &
                               ptend%q(:,pver:1:-1,ixcldliq), ptend%s(:,pver:1:-1))
 
   ftem(:ncol,:pver) = ptend%s(:ncol,:pver)/cpair
   call outfld('YOGDT   ',ftem               ,pcols   ,lchnk   )
   call outfld('YOGDQ   ',ptend%q(1,1,1) ,pcols   ,lchnk   )
   call outfld('YOGDICE ',ptend%q(1,1,ixcldice) ,pcols   ,lchnk   )
   call outfld('YOGDLIQ ',ptend%q(1,1,ixcldliq) ,pcols   ,lchnk   )
   call outfld('YOGPREC ',yog_precsfc ,pcols   ,lchnk   )
end subroutine yog_tend

!=========================================================================================

end module yog_intr
