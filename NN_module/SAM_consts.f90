module SAM_consts_mod
  !! Physical constants from the SAM Model required for variable conversions
  
  implicit none
  public

  ! From params.f90 in SAM
  ! These are constants used in equations as part of the parameterisation.
  ! Since the net is trained with them they should be left as defined here, rather than
  ! forced to match the external model

  ! Physical Constants:

  != unit (m / s^2)
  real, parameter :: ggr = 9.81
      !! Gravity acceleration, m/s2

  != unit J / kg :: lcond, lsub, lfus
  real, parameter :: lfus = 0.3336e+06
      !! Latent heat of fusion
  real, parameter :: lcond = 2.5104e+06
      !! Latent heat of condensation
  real, parameter :: lsub = 2.8440e+06
      !! Latent heat of sublimation, J/kg

  != unit (J / kg) / K :: cp
  real, parameter :: cp = 1004.
      !! Specific heat of air
  
  ! unit K :: fac_cond, fac_fus, fac_sub
  real, parameter :: fac_cond = lcond/cp
      !!
  real, parameter :: fac_fus = lfus/cp
      !!
  real, parameter :: fac_sub = lsub/cp
      !!

  ! Temperatures limits for various hydrometeors
  != unit K :: tprmin
  real, parameter :: tprmin = 268.16
      !! Minimum temperature for rain
  
  != unit K :: tprmax
  real, parameter :: tprmax = 283.16
      !! Maximum temperature for snow+graupel, K
  
  != unit K :: tbgmin
  real, parameter :: tbgmin = 253.16
      !! Minimum temperature for cloud water.
  
  != unit K :: tbgmin
  real, parameter :: tbgmax = 273.16
      !! Maximum temperature for cloud ice, K

  != unit K :: tbgmin
  real, parameter :: tgrmin = 223.16
      !! Maximum temperature for snow, K

  != unit K :: tbgmin
  real, parameter :: tgrmax = 283.16
      !! Maximum temperature for graupel, K

  ! Misc. microphysics variables
  ! != unit 1 / K :: a_pr
  real, parameter :: a_pr = 1./(tprmax-tprmin)
      !! Misc. microphysics variables
  
  ! != unit 1 / K :: a_bg
  ! real, parameter :: a_bg = 1./(tbgmax-tbgmin)
  !     !! Misc. microphysics variables

  real, parameter :: an = 1./(tbgmax-tbgmin)
  real, parameter :: bn = tbgmin * an
  real, parameter :: ap = 1./(tprmax-tprmin)
  real, parameter :: bp = tprmin * ap
  real, parameter :: fac1 = fac_cond+(1+bp)*fac_fus
  real, parameter :: fac2 = fac_fus*ap
  real, parameter :: ag = 1./(tgrmax-tgrmin)

  ! ---------------------------
  ! SAM Grid Variables
  ! ---------------------------
  integer, parameter :: input_ver_dim = 48
      !! Set to 48 in setparm.f90 of SAM. Same as nz_gl??
  
  ! Outputs from NN are supplied at lowest 30 half-model levels for sedimentation fluxes,
  ! and at 29 levels for fluxes (as flux at bottom boundary is zero).
  integer, parameter :: nrf = 30
      !! number of vertical levels the NN uses
  integer, parameter :: nrfq = nrf - 1
      !! number of vertical levels the NN uses when boundary condition is set to 0


!---------------------------------------------------------------------
! Functions and Subroutines

contains

  subroutine check(err_status)
  !! Check error status after netcdf call and print message for
  !! error codes.

  use netcdf

  integer, intent(in) :: err_status
    !! error status from nf90 function

  if(err_status /= nf90_noerr) then
     write(*, *) trim(nf90_strerror(err_status))
  end if

  end subroutine check

end module SAM_consts_mod
