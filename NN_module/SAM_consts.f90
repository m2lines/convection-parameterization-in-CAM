module SAM_consts_mod
  !! Physical constants from the SAM Model required for variable conversions
  
  implicit none
  public

  ! From params.f90 in SAM
  ! These are constants used in equations as part of the parameterisation.
  ! Since the net is trained with them they should be left as defined here, rather than
  ! forced to match the external model

  ! Physical Constants:
  
  != unit J / kg :: lfus
  real, parameter :: lfus = 0.3336e+06
      !! Latent heat of fusion
  
  != unit J / kg :: lcond
  real, parameter :: lcond = 2.5104e+06
      !! Latent heat of condensation
  
  != unit (J / kg) / K :: cp
  real, parameter :: cp = 1004.
      !! Specific heat of air
  
  ! unit K - not checkable yet
  real, parameter :: fac_cond = lcond/cp
      !!
  
  ! unit K - not checkable yet
  real, parameter :: fac_fus = lfus/cp
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



end module SAM_consts_mod
