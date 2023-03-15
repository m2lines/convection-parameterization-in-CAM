module params

use grid, only: nzm

implicit none

!   Constants:

real, parameter :: cp = 1004.             ! Specific heat of air, J/kg/K
real, parameter :: ggr = 9.81             ! Gravity acceleration, m/s2
real, parameter :: lcond = 2.5104e+06     ! Latent heat of condensation, J/kg
real, parameter :: lfus = 0.3336e+06      ! Latent heat of fusion, J/kg
real, parameter :: lsub = 2.8440e+06      ! Latent heat of sublimation, J/kg
real, parameter :: rv = 461.              ! Gas constant for water vapor, J/kg/K
real, parameter :: rgas = 287.            ! Gas constant for dry air, J/kg/K
real, parameter :: diffelq = 2.21e-05     ! Diffusivity of water vapor, m2/s
real, parameter :: therco = 2.40e-02      ! Thermal conductivity of air, J/m/s/K
real, parameter :: muelq = 1.717e-05      ! Dynamic viscosity of air

real, parameter :: fac_cond = lcond/cp 
real, parameter :: fac_fus = lfus/cp
real, parameter :: fac_sub = lsub/cp

!  Microphysics stuff:

! Densities of hydrometeors

real, parameter :: rhor = 1000. ! Density of water, kg/m3
real, parameter :: rhos = 100.  ! Density of snow, kg/m3
real, parameter :: rhog = 400.  ! Density of graupel, kg/m3
!eal, parameter :: rhog = 917.  ! hail - Lin 1983    

! Temperatures limits for various hydrometeors

real, parameter :: tbgmin = 253.16    ! Minimum temperature for cloud water., K
real, parameter :: tbgmax = 273.16    ! Maximum temperature for cloud ice, K
real, parameter :: tprmin = 268.16    ! Minimum temperature for rain, K
real, parameter :: tprmax = 283.16    ! Maximum temperature for snow+graupel, K
real, parameter :: tgrmin = 223.16    ! Minimum temperature for snow, K
real, parameter :: tgrmax = 283.16    ! Maximum temperature for graupel, K

! Terminal velocity coefficients

real, parameter :: a_rain = 842. ! Coeff.for rain term vel 
real, parameter :: b_rain = 0.8  ! Fall speed exponent for rain
real, parameter :: a_snow = 4.84 ! Coeff.for snow term vel
real, parameter :: b_snow = 0.25 ! Fall speed exponent for snow
!eal, parameter :: a_grau = 40.7! Krueger (1994) ! Coef. for graupel term vel
real, parameter :: a_grau = 94.5 ! Lin (1983) (rhog=400)
!eal, parameter :: a_grau = 127.94! Lin (1983) (rhog=917)
real, parameter :: b_grau = 0.5  ! Fall speed exponent for graupel

! Autoconversion

real, parameter :: qcw0 = 1.e-3      ! Threshold for water autoconversion, g/g  
real, parameter :: qci0 = 1.e-4     ! Threshold for ice autoconversion, g/g
real, parameter :: alphaelq = 1.e-3  ! autoconversion of cloud water rate coef
real, parameter :: betaelq = 1.e-3   ! autoconversion of cloud ice rate coef

! Accretion

real, parameter :: erccoef = 1.0   ! Rain/Cloud water collection efficiency
real, parameter :: esccoef = 1.0   ! Snow/Cloud water collection efficiency
real, parameter :: esicoef = 0.1   ! Snow/cloud ice collection efficiency
real, parameter :: egccoef = 1.0   ! Graupel/Cloud water collection efficiency
real, parameter :: egicoef = 0.1   ! Graupel/Cloud ice collection efficiency

! Interseption parameters for exponential size spectra

real, parameter :: nzeror = 8.e6   ! Intercept coeff. for rain  
real, parameter :: nzeros = 3.e6   ! Intersept coeff. for snow
real, parameter :: nzerog = 4.e6   ! Intersept coeff. for graupel
!eal, parameter :: nzerog = 4.e4   ! hail - Lin 1993 

! Cloud droplet distribution properties, used in sedimentation scheme.
real, parameter :: Nc0 = 65. ! cloud droplet concentration [cm^{-3}]

real, parameter :: qp_threshold = 1.e-8 ! minimal rain/snow water content




!  Variables:

            
real  pres0      ! Reference surface pressure, Pa
real  ug	 ! Velocity of the Domain's drift in x direction
real  vg	 ! Velocity of the Domain's drift in y direction
real  fcor       ! Coriolis parameter	
real  fcorz      ! Vertical Coriolis parameter
real  pi

real longitude0  ! latitude of the domain's center 
real latitude0   ! longitude of the domain's center 
real coszrs      ! solar zenith angle


!  Surface stuff:   	

real   tabs_s	! surface temperature,K
real   fluxt0   ! surface sensible flux, Km/s
real   fluxq0   ! surface latent flux, m/s
real   tau0     ! surface stress, m2/s2
real   z0	! roughness length
real   soil_wetness ! wetness coeff for soil (from 0 to 1.)
integer ocean_type ! type of SST forcing
 
!  Large-scale stuff

real  timelargescale ! time to start large-scale forcing
real  tauls	     ! nudging-to-large-scaler-profile time-scale


! Misc. microphysics variables

real(4) gam3       ! Gamma function of 3
real(4) gams1      ! Gamma function of (3 + b_snow)
real(4) gams2      ! Gamma function of (5 + b_snow)/2
real(4) gams3      ! Gamma function of (4 + b_snow)
real(4) gamg1      ! Gamma function of (3 + b_grau)
real(4) gamg2      ! Gamma function of (5 + b_grau)/2
real(4) gamg3      ! Gamma function of (4 + b_grau)
real(4) gamr1      ! Gamma function of (3 + b_rain)
real(4) gamr2      ! Gamma function of (5 + b_rain)/2
real(4) gamr3      ! Gamma function of (4 + b_rain)
      
real accrsc(nzm),accrsi(nzm),accrrc(nzm),coefice(nzm)
real accrgc(nzm),accrgi(nzm)
real evaps1(nzm),evaps2(nzm),evapr1(nzm),evapr2(nzm)
real evapg1(nzm),evapg2(nzm)
            
real a_bg, a_pr, a_gr 

!==========================================================================
! UW ADDITIONS

! Krueger microphysics additions
real, parameter :: ttfrz = 233.1
real, parameter :: tice = 273.16
real, parameter :: a_bgkru = 1./(tice-ttfrz)

real vconr, vcons, vcong ! fall speed constants for Kreuger microphysics
real, parameter :: vexpr = 0.2 ! exponents for Kreugers fall speed formulas
real, parameter :: vexpg = 0.125
real, parameter :: vexps = 0.0625
real, parameter :: vmaxr = 10. ! maximum precipitation velocities
real, parameter :: vmaxg = 20.

! Tracer nudging timescales.

real  itauo3	     ! inverse of nudging-to-large-scaler-O3-profile time-scale
real  itautrz	     ! inverse of the decay time of trz
real  itautrx	     ! inverse of the decay time of trx
real  itautry	     ! inverse of the decay time of try
real  itautrzz	     ! inverse of the decay time of trzz

real, parameter :: vk = 0.4 ! von Karman constant (added by cw 2/16/06)

! END UW ADDITIONS
!==========================================================================

end module params
