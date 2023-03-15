module vars

use grid

implicit none
!--------------------------------------------------------------------
! prognostic variables:

real u   (dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm) ! x-wind
real v   (dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm) ! y-wind
real w   (dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz ) ! z-wind
real t   (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! moist static energy
real q   (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! total water
real qp	 (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! rain+snow
real tke (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! SGS TKE

!--------------------------------------------------------------------
! diagnostic variables:

real p      (0:nx, (1-YES3D):ny, nzm)     ! pressure
real tabs   (nx, ny, nzm)                 ! temperature
real qn     (nx, ny, nzm)                 ! cloud water+cloud ice
! tk and tkh split into _xy and _z by cw 2/17/06
real tk_xy  (0:nxp1, (1-YES3D):nyp1, nzm) ! horizontal SGS eddyviscosity 
real tk_z   (0:nxp1, (1-YES3D):nyp1, nzm) ! vertical SGS eddyviscosity
real tkh_xy (0:nxp1, (1-YES3D):nyp1, nzm) ! horizontal SGS eddy conductivity
real tkh_z  (0:nxp1, (1-YES3D):nyp1, nzm) ! vertical SGS eddy conductivity
! following added by cw 2/23/06
real tk_z_uw(nx, ny, nzm)                 ! vertical SGS eddy viscosity
real tkh_z_uw(nx, ny, nzm)                ! vertical SGS eddy conductivity 
real cldfrc (nx, ny, nzm)                 ! Cloud fraction (from UWPBL)
real pblh   (nx, ny)                      ! PBL height (from UWPBL)
real tke_uw (nx, ny, nz)                  ! TKE (from UWPBL)
real fer    (nx, ny, nzm)                 ! lateral entrainment rate (UWCU)
real fdr    (nx, ny, nzm)                 ! lateral detrainment rate (UWCU)
real cin    (nx, ny)                      ! Convective inhibition
real cbmf   (nx, ny)                      ! cloud base mass flux
real cush   (nx, ny)                      ! convective scale height (from UWPBL)

real uwcu_pert_cldsrc(nx, ny)             ! cloud root random perturbation from UWCU scheme


real tkh_x_rf (0:nxp1, (1-YES3D):nyp1, nzm) ! JY - janniy Yani anded for RF to predict conductidvity. 
real tkh_y_rf (0:nxp1, (1-YES3D):nyp1, nzm)
real tkh_z_rf (nx,ny,nzm)

!Fields from the beginning of the time step JY
real u_i (nx,ny,nzm)
real v_i (nx,ny,nzm)
real w_i (nx,ny,nzm)
real t_i (nx,ny,nzm)
real q_i (nx,ny,nzm)
real qp_i (nx,ny,nzm)


integer uwpbl_ndiff(nx, ny)                  ! number of PBL substeps
real qtend_uwpbl(0:nxp1, 0:nyp1, nzm)
real qptend_uwpbl(0:nxp1, 0:nyp1, nzm)
real qvtend_uwpbl(0:nxp1, 0:nyp1, nzm)
real qctend_uwpbl(0:nxp1, 0:nyp1, nzm)
real qitend_uwpbl(0:nxp1, 0:nyp1, nzm)
real utend_uwpbl(0:nxp1, 0:nyp1, nzm)
real vtend_uwpbl(0:nxp1, 0:nyp1, nzm)
real wtend_uwpbl(0:nxp1, 0:nyp1, nz)
real ttend_uwpbl(0:nxp1, 0:nyp1, nzm)

real uwcu_err(nx, ny)                     ! error flag from UWCU scheme
real ttend_uwcu (nx, ny, nzm)             ! convective heating (UWCU)
real qvtend_uwcu (nx, ny, nzm)            ! convective moistening (UWCU)
real qltend_uwcu (nx, ny, nzm)            ! convective moistening (UWCU)
real qitend_uwcu (nx, ny, nzm)            ! convective moistening (UWCU)
real umf (nx, ny, nzm)                    ! convective mass flux(UWCU)
real slflx_uwcu (nx, ny, nzm)             ! flux of dry static energy(UWCU)
real qtflx_uwcu (nx, ny, nzm)             ! flux of total moisture (UWCU)
real kvh_out_flip_3d(nx, ny, nz)          ! save diffusivities from last time-
real kvm_out_flip_3d(nx, ny, nz)          ! step for UW PBL scheme

!--------------------------------------------------------------------
! time-tendencies for prognostic variables

real dudt   (nxp1, ny, nzm, 3)
real dvdt   (nx, nyp1, nzm, 3)
real dwdt   (nx, ny, nz,  3)

!----------------------------------------------------------------
! Temporary storage array:

	real misc(nx, ny, nz)
!------------------------------------------------------------------
! fluxes at the top and bottom of the domain:

real uwpbl_fluxbu (nx, ny), uwpbl_fluxbv (nx, ny), uwpbl_fluxbt (nx, ny)
real uwpbl_fluxbq (nx, ny)

real fluxbu (nx, ny), fluxbv (nx, ny), fluxbt (nx, ny)
real fluxbq (nx, ny), fluxtu (nx, ny), fluxtv (nx, ny)
real fluxtt (nx, ny), fluxtq (nx, ny), fzero  (nx, ny)
real precsfc(nx,ny) ! surface precip. rate
                
!-----------------------------------------------------------------
! profiles 

real   t0(nzm), q0(nzm), qc0(nzm), qv0(nzm),tabs0(nzm), &
       qi0(nzm),tv0(nzm), rel0(nzm), u0(nzm), v0(nzm), &
       tg0(nzm), qg0(nzm), ug0(nzm), vg0(nzm), p0(nzm), &
       tke0(nzm), t01(nzm), q01(nzm)
!----------------------------------------------------------------
! "observed" (read from snd file) surface characteristics 

real  sstobs, lhobs, shobs
!----------------------------------------------------------------
!  Domain top stuff:

real   gamt0    ! gradient of t() at the top,K/m
real   gamq0    ! gradient of q() at the top,g/g/m

!-----------------------------------------------------------------
! reference vertical profiles:
 
real   prespot(nzm)  ! (1000./pres)**R/cp
real   rho(nzm)	  ! air density at pressure levels,kg/m3 
real   rhow(nz)   ! air density at vertical velocity levels,kg/m3
real   bet(nzm)	  ! = ggr/tv0
real   betp(nzm)  ! = ggr/p0
real   gamaz(nzm) ! ggr/cp*z
real   wsub(nz)   ! Large-scale subsidence velocity,m/s
real   qtend(nzm) ! Large-scale tendency for total water
real   ttend(nzm) ! Large-scale tendency for temp.

!---------------------------------------------------------------------
! Large-scale and surface forcing:

integer nmaxlsf ! dimension of lsf forcing arrays (time)
integer nmaxrfc ! dimension of rad forcing arrays (time) 
integer nmaxsfc ! dimension of sfc forcing arrays (time) 
integer nmaxsnd ! dimension of observed sounding arrays (time) 
parameter (nmaxrfc=1000,nmaxlsf=1000,nmaxsfc=1000, nmaxsnd=1000)
integer nlsf	! number of large-scale forcing profiles
integer nrfc	! number of radiative forcing profiles
integer nsfc	! number of surface forcing profiles
integer nsnd	! number of observed soundings

real   dqls(nzm,nmaxlsf) ! Large-scale tendency for total water
real   dtls(nzm,nmaxlsf) ! Large-scale tendency for temp.
real   ugls(nzm,nmaxlsf) ! Large-scale wind in X-direction
real   vgls(nzm,nmaxlsf) ! Large-scale wind in Y-direction
real   wgls(nzm,nmaxlsf) ! Large-scale subsidence velocity,m/s
real   pres0ls(nmaxlsf) ! Surface pressure, mb
real   dayls(nmaxlsf)    ! Large-scale forcing arrays time (days) 
real   dtrfc(nzm,nmaxrfc)! Radiative tendency for pot. temp.
real   dayrfc(nmaxrfc)   ! Radiative forcing arrays time (days) 
real   sstsfc(nmaxsfc)   ! SSTs
real   hsfc(nmaxsfc) 	 ! Sensible heat flux,W/m2
real   lesfc(nmaxsfc) 	 ! Latent heat flux,W/m2
real   tausfc(nmaxsfc) 	 ! Surface drag,m2/s2
real   daysfc(nmaxsfc)   ! Surface forcing arrays time (days) 
real   usnd(nzm,nmaxsnd) ! Observed zonal wind
real   vsnd(nzm,nmaxsnd) ! Observed meriod wind
real   tsnd(nzm,nmaxsnd) ! Observed Abs. temperature
real   qsnd(nzm,nmaxsnd) ! Observed Moisture
real   daysnd(nmaxsnd)
 
!---------------------------------------------------------------------
!  Horizontally varying stuff (as a function of xy)
!
real sstxy(0:nxp1,1-YES3D:nyp1)	!  surface temperature xy-distribution
real fcory(0:ny)      !  Coriolis parameter xy-distribution !JY - Yani added 0:ny
real fcorzy(0:ny)      !  z-Coriolis parameter xy-distribution  !JY - Yani added ny+1
real latitude(nx,ny)	     ! latitude (degrees)
real longitude(nx,ny)	     ! longitude(degrees)
real prec_xy(nx,ny) ! surface precipitation rate
real shf_xy(nx,ny) ! surface sensible heat flux
real lhf_xy(nx,ny) ! surface latent heat flux
real lwns_xy(nx,ny) ! mean net lw at SFC
real swns_xy(nx,ny) ! mean net sw at SFC
real lwnsc_xy(nx,ny) ! clear-sky mean net lw at SFC
real swnsc_xy(nx,ny) ! clear-sky mean net sw at SFC
real lwnt_xy(nx,ny) ! mean net lw at TOA
real swnt_xy(nx,ny) ! mean net sw at TOA
real lwntc_xy(nx,ny) ! clear-sky mean net lw at TOA
real swntc_xy(nx,ny) ! clear-sky mean net sw at TOA
real solin_xy(nx,ny) ! solar TOA insolation
real pw_xy(nx,ny)   ! precipitable water
real cw_xy(nx,ny)   ! cloud water path
real iw_xy(nx,ny)   ! ice water path
real u200_xy(nx,ny) ! u-wind at 200 mb
real usfc_xy(nx,ny) ! u-wind at at the surface
real v200_xy(nx,ny) ! v-wind at 200 mb
real vsfc_xy(nx,ny) ! v-wind at the surface
real w500_xy(nx,ny) ! w at 500 mb
real qocean_xy(nx,ny) ! ocean cooling in W/m2
real pblh_xy(nx, ny) ! PBL height (from UW PBL scheme)
real cin_xy(nx, ny) ! CIN (from UW PBL scheme)
real cbmf_xy(nx, ny) ! cloud base mass flux (from UW PBL scheme)

!----------------------------------------------------------------------
!  Tendencies due to convective parameterization:
!

real   qtend_cup(nx,ny,nzm) ! CUP tendency for total water
real   ttend_cup(nx,ny,nzm) ! CUP tendency for temp.
real   utend_cup(nx,ny,nzm) ! CUP tendency for u
real   vtend_cup(nx,ny,nzm) ! CUP tendency for v

!----------------------------------------------------------------------
!	Vertical profiles of quantities sampled for statitistics purposes:

real &
    twle(nz), twsb(nz), qwle(nz), qwsb(nz), tkewle(nz), &
    tkewsb(nz), qpwle(nz), qpwsb(nz), precflux(nz), &
    uwle(nz), uwsb(nz), vwle(nz), vwsb(nz), &
    cloud_factor(nz), core_factor(nz), coredn_factor(nz), &
    tkeleadv(nz), tkelepress(nz), tkelediss(nz), tkelediff(nz), &
    tkesbbuoy(nz), tkesbshear(nz),tkesbdiss(nz), tkesbdiff(nz), &
    tkelebuoy(nz), radlwup(nz), radlwdn(nz), radswup(nz), radswdn(nz), &
    radqrlw(nz), radqrsw(nz), w_max, s_acld, s_acldcold, s_ar, p_conv, p_strat,&
    s_acldl, s_acldm, s_acldh, s_acldisccp, &
    s_acldlisccp, s_acldmisccp, s_acldhisccp, &
    s_flns,s_flnt,s_flnsc,s_flntc,s_flds,s_fsns, &
    s_fsnt,s_fsnsc,s_fsntc,s_fsds,s_solin, & 
    t2leadv(nz),t2legrad(nz),t2lediff(nz),t2leprec(nz),t2lediss(nz), &
    q2leadv(nz),q2legrad(nz),q2lediff(nz),q2leprec(nz),q2lediss(nz), &
    s2leadv(nz),s2legrad(nz),s2lediff(nz),s2lediss(nz), &
    twleadv(nz),twlediff(nz),twlepres(nz),twlebuoy(nz),twleprec(nz), &
    qwleadv(nz),qwlediff(nz),qwlepres(nz),qwlebuoy(nz),qwleprec(nz), &
    swleadv(nz),swlediff(nz),swlepres(nz),swlebuoy(nz), &
    momleadv(nz,3),momlepress(nz,3),momlebuoy(nz,3), &
    momlediff(nz,3),tadv(nz),tdiff(nz),tlat(nz), tlatqi(nz), qifall(nz), &
    qadv(nz),qdiff(nz),qpadv(nz),qpdiff(nz),qpsrc(nz),qpfall(nz),qpevp(nz)


! register functions:


        real esatw,esati,dtesatw,dtesati
        real qsatw,qsati,dtqsatw,dtqsati
        external esatw,esati,dtesatw,dtesati,qsatw,qsati,dtqsatw,dtqsati

	real omegan, omegap, omegag
	external omegan, omegap, omegag


        integer lenstr
        external lenstr

! energy conservation diagnostics:
 
  double precision total_water_before, total_water_after
  double precision total_water_evap, total_water_prec, total_water_ls

!===========================================================================
! UW ADDITIONS

! slab ocean stuff
real :: qoceanxy(nx,ny) ! instantaneous slab ocean heating
double precision :: sstxy_dble(nx,ny) ! double precision slab ocean temp

!kzm: o3g0 is background o3 to nudge toward, 
!     itauz is the inverse nudging timescale profile
real   o3g0(nzm) , itauz(nzm),thmode1(nzm),thmode2(nzm)

real   qlsvadv(nzm) ! Large-scale vertical advection tendency for total water
real   tlsvadv(nzm) ! Large-scale vertical advection tendency for temperature
real   qnudge(nzm) ! Nudging of horiz.-averaged total water profile
real   tnudge(nzm) ! Nudging of horiz.-averaged temperature profile

real swntm_xy(nx,ny) ! mean net incoming sw at top of model
real radqrclw(nz) ! mean clearsky longwave radiative heating
real radqrcsw(nz) ! mean clearsky shortwave radiative heating

! 200 mbar geopotential height
real h200_xy(nx,ny) ! geopotential height at 200 mb (in meters)

! 850 mbar horizontal winds and geopotential height
real u850_xy(nx,ny) ! zonal velocity at 850 mb
real v850_xy(nx,ny) ! meridional velocity at 850 mb
real h850_xy(nx,ny) ! geopotential height at 850 mb (in meters)

! 850 mbar theta_e and theta_es
real the850_xy(nx,ny)  ! equivalent potential temperature at 850 mb
real thes850_xy(nx,ny) ! saturated equivalent potential temperature at 850 mb

! 1000 mbar theta_e and theta_es
real the1000_xy(nx,ny)  ! equivalent potential temperature at 1000 mb
real thes1000_xy(nx,ny) ! saturated equivalent potential temperature at 1000 mb

! Surface pressure
real psfc_xy(nx,ny) ! pressure (in millibar) at lowest grid point

! Saturated water vapor path, useful for computing column relative humidity
real swvp_xy(nx,ny)  ! saturated water vapor path (wrt water)

! Cloud and echo top heights, and cloud top temperature (instantaneous)
real cloudtopheight(nx,ny), echotopheight(nx,ny), cloudtoptemp(nx,ny)

! w, q, theta projected onto first two gravity wave modes
real wmode1_xy(nx,ny)  ! first mode w
real wmode2_xy(nx,ny)  ! second mode w
real qmode1_xy(nx,ny)  ! first mode q
real qmode2_xy(nx,ny)  ! second mode q
real thmode1_xy(nx,ny)  ! first mode theta
real thmode2_xy(nx,ny)  ! second mode theta
real wmode1i_xy(nx,ny)  ! first mode w, instantaneous
real wmode2i_xy(nx,ny)  ! second mode w, instantaneous

! variables for nudging_hadley, kzm Aug 17,2005
real ug0x(ny_gl,nzm),vg0x(ny_gl,nzm),wg0x(ny_gl,nzm),tg0x(ny_gl,nzm),qg0x(ny_gl,nzm)

! First/second baroclinic mode energetics
double precision oldw1_xy(nx,ny), oldw2_xy(nx,ny)  ! vertical velocity
double precision wdwdt1_xy(nx,ny), wdwdt2_xy(nx,ny)
double precision olds1_xy(nx,ny), olds2_xy(nx,ny)  ! liquid-ice static energy
double precision sdsdt1_xy(nx,ny), sdsdt2_xy(nx,ny)
double precision oldh1_xy(nx,ny), oldh2_xy(nx,ny)  ! frozen moist static energy
double precision hdhdt1_xy(nx,ny), hdhdt2_xy(nx,ny)

! Storage and advection terms for column-integrated budgets.
real hstor_xy(nx,ny) ! column-avg. frozen moist static energy storage
real hadv_xy(nx,ny)  ! column-avg. frozen moist static energy advection
real sstor_xy(nx,ny) ! column-avg. liquid water ice static energy storage
real sadv_xy(nx,ny)  ! column-avg. liquid water ice static energy advection
real old_fmse(nx,ny) ! old column-avg. frozen moist static energy storage
real old_sli(nx,ny)  ! old column-avg. liquid water ice static energy storage

! Radar reflectivity slices computed in statistics.f90
real dbZe1km_xy(nx,ny) ! Equivalent Radar Reflectivity at 1km
real dbZe3km_xy(nx,ny) ! Equivalent Radar Reflectivity at 3km
real dbZe6km_xy(nx,ny) ! Equivalent Radar Reflectivity at 6km

real prec_inst(nx,ny) ! Instantaneous surface precipitation flux (mm/day)

! weighting function for the sst when driving mean sst to target
real sst_wgt(nx,ny) !  surface temperature weighting


! add variables for dokruegermicro

real, dimension(:,:,:), allocatable :: qv, qc, qi, qr, qs, qg ! [kg/kg]
real, dimension(:), allocatable :: &
     qcwle,qcwsb,qcadv,qcdiff,qcsrc,qcfall,qcevp,&
     qiwle,qiwsb,qiadv,qidiff,qisrc,qievp,&
     qrwle,qrwsb,qradv,qrdiff,qrsrc,qrfall,qrevp,&
     qgwle,qgwsb,qgadv,qgdiff,qgsrc,qgfall,qgevp,&
     qswle,qswsb,qsadv,qsdiff,qssrc,qsfall,qsevp,&
     rainfrac, snowfrac, graufrac, clifrac, clwfrac, &
     vtrz, vtsz, vtgz, vtiz, &
     rainflux, snowflux, grauflux, &
     pclwz, pvaporz, pcliz, pimltz, pihomz, &
     pidwz, prainz, prautz, pracwz, prevpz, &
     psnowz, psautz, psfwz, psfiz, praciz, &
     piacrz, psaciz, psacwz, psdepz, pssubz, &
     pracsz, psacrz, psmltz, psmltevpz, pladjz, &
     piadjz, pgraupelz, pgautz,pgfrz, pgacwz, &
     pgaciz, pgacrz, pgacsz, pgacipz, pgacrpz, &
     pgacspz, pgwetz, pdryz, pgsubz, pgdepz, &
     pgmltz, pgmltevpz, pgwordz

! add variables for hydrometeor isotope mixture fractions (future addition)

real, dimension(:,:,:), allocatable :: iso_qv, iso_qc, iso_qi, &
                                       iso_qr, iso_qs, iso_qg ! [kg/kg]
real, dimension(:), allocatable :: &
     iso_qcwle,iso_qcwsb,iso_qcadv,iso_qcdiff,iso_qcsrc,iso_qcfall,iso_qcevp,&
     iso_qiwle,iso_qiwsb,iso_qiadv,iso_qidiff,iso_qisrc,iso_qifall,iso_qievp,&
     iso_qrwle,iso_qrwsb,iso_qradv,iso_qrdiff,iso_qrsrc,iso_qrfall,iso_qrevp,&
     iso_qgwle,iso_qgwsb,iso_qgadv,iso_qgdiff,iso_qgsrc,iso_qgfall,iso_qgevp,&
     iso_qswle,iso_qswsb,iso_qsadv,iso_qsdiff,iso_qssrc,iso_qsfall,iso_qsevp

real, dimension(:,:,:), allocatable :: iso2_qv, iso2_qc, iso2_qi, &
                                       iso2_qr, iso2_qs, iso2_qg ! [kg/kg]
real, dimension(:), allocatable :: &
     iso2_qcwle,iso2_qcwsb,iso2_qcadv,iso2_qcdiff,&
       iso2_qcsrc,iso2_qcfall,iso2_qcevp,&
     iso2_qiwle,iso2_qiwsb,iso2_qiadv,iso2_qidiff,&
       iso2_qisrc,iso2_qifall,iso2_qievp,&
     iso2_qrwle,iso2_qrwsb,iso2_qradv,iso2_qrdiff,&
       iso2_qrsrc,iso2_qrfall,iso2_qrevp,&
     iso2_qgwle,iso2_qgwsb,iso2_qgadv,iso2_qgdiff,&
       iso2_qgsrc,iso2_qgfall,iso2_qgevp,&
     iso2_qswle,iso2_qswsb,iso2_qsadv,iso2_qsdiff,&
       iso2_qssrc,iso2_qsfall,iso2_qsevp

!tracers added by kzm

! variables related to tracer of x grid
real, dimension(:,:,:), allocatable :: trx
real, dimension(:,:), allocatable :: fluxbtrx
real, dimension(:), allocatable :: trxadv,trxdiff,trxwsb,trxwle, &
     trxwleadv,trxwlediff,trx2leadv,trx2lediff,trx2legrad,trx2lediss

! variables related to tracer of y grid
real, dimension(:,:,:), allocatable :: try
real, dimension(:,:), allocatable :: fluxbtry
real, dimension(:), allocatable :: tryadv,trydiff,trywsb,trywle, &
     try2leadv,trywleadv,try2legrad,try2lediff,try2lediss,trywlediff

! variables related to tracer of vertical grid
real, dimension(:,:,:), allocatable :: trz
real, dimension(:,:), allocatable :: fluxbtrz
real, dimension(:), allocatable :: trzadv,trzdiff,trzwsb,trzwle, &
     trz2leadv,trzwleadv,trz2legrad,trz2lediff,trz2lediss,trzwlediff

! variables related to tracer of vertical grid^2
real, dimension(:,:,:), allocatable :: trzz
real, dimension(:,:), allocatable :: fluxbtrzz
real, dimension(:), allocatable :: trzzadv,trzzdiff,trzzwsb,trzzwle, &
     trzz2leadv,trzzwleadv,trzz2legrad,trzz2lediff,trzz2lediss,trzzwlediff

! variables related to interactive ozone (o3)
real, dimension(:,:,:), allocatable :: tro3
real, dimension(:,:), allocatable :: fluxbtro3
real, dimension(:), allocatable :: tro3adv,tro3diff,tro3wsb,tro3wle, &
     tro32leadv,tro3wleadv,tro32legrad,tro32lediff,tro32lediss,tro3wlediff

!real trz  (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! tracer of vertical grid
!real trx  (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! tracer of x grid
!real try  (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! tracer of y grid
!real trzz  (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! tracer of vertical grid^2
!real tro3 (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! interactive o3
!
!
!real &   !kzm add tracer related variables
!    trx2leadv(nz),trxwleadv(nz),trx2legrad(nz),trx2lediff(nz),trx2lediss(nz),&
!    trxwlediff,trxwle(nz),trxadv(nz),trxdiff(nz),trxwsb(nz), &
!    try2leadv(nz),trywleadv(nz),try2legrad(nz),try2lediff(nz),try2lediss(nz),&
!    trywlediff,trywle(nz),tryadv(nz),trydiff(nz),trywsb(nz), &
!    trz2leadv(nz),trzwleadv(nz),trz2legrad(nz),trz2lediff(nz),trz2lediss(nz),&
!    trzwlediff,trzwle(nz),trzadv(nz),trzdiff(nz),trzwsb(nz), &
!    trzz2leadv(nz),trzzwleadv(nz),trzz2legrad(nz),trzz2lediff(nz),trzz2lediss(nz),&
!    trzzwlediff,trzzwle(nz),trzzadv(nz),trzzdiff(nz),trzzwsb(nz), &
!    tro32leadv(nz),tro3wleadv(nz),tro32legrad(nz),tro32lediff(nz), &
!    tro32lediss(nz),tro3wlediff,tro3wle(nz),tro3adv(nz),tro3diff(nz), &
!    tro3wsb(nz), radqrclw(nz), radqrcsw(nz)

! END UW ADDITIONS
!===========================================================================

end module vars
