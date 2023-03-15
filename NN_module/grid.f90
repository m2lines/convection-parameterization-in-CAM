module grid

use domain

implicit none
        
integer, parameter :: nx = nx_gl/nsubdomains_x
integer, parameter :: ny = ny_gl/nsubdomains_y 
integer, parameter :: nz = nz_gl+1
integer, parameter :: nzm = nz-1
        
integer, parameter :: nsubdomains = nsubdomains_x * nsubdomains_y

logical, parameter :: RUN3D = ny_gl.gt.1
logical, parameter :: RUN2D = .not.RUN3D

integer, parameter :: nxp1 = nx + 1
integer, parameter :: nyp1 = ny + 1 * YES3D
integer, parameter :: nxp2 = nx + 2
integer, parameter :: nyp2 = ny + 2 * YES3D
integer, parameter :: nxp3 = nx + 3
integer, parameter :: nyp3 = ny + 3 * YES3D
integer, parameter :: nxp4 = nx + 4
integer, parameter :: nyp4 = ny + 4 * YES3D

integer, parameter :: dimx1_u = -1
integer, parameter :: dimx2_u = nxp3
integer, parameter :: dimy1_u = 1-2*YES3D
integer, parameter :: dimy2_u = nyp2
integer, parameter :: dimx1_v = -1
integer, parameter :: dimx2_v = nxp2
integer, parameter :: dimy1_v = 1-2*YES3D
integer, parameter :: dimy2_v = nyp3
integer, parameter :: dimx1_w = -2
integer, parameter :: dimx2_w = nxp3
integer, parameter :: dimy1_w = 1-3*YES3D
integer, parameter :: dimy2_w = nyp3
!integer, parameter :: dimx1_w = -1
!integer, parameter :: dimx2_w = nxp2
!integer, parameter :: dimy1_w = 1-2*YES3D
!integer, parameter :: dimy2_w = nyp2
integer, parameter :: dimx1_s = -2
integer, parameter :: dimx2_s = nxp3
integer, parameter :: dimy1_s = 1-3*YES3D
integer, parameter :: dimy2_s = nyp3

integer, parameter :: ncols = nx*ny


real dx 	! grid spacing in x direction
real dy		! grid spacing in y direction
real dz		! grid spacing in z direction for the lowest grid layer

real z(nz)      ! height of the pressure levels above surface,m
real pres(nzm)  ! pressure,mb at scalar levels
real zi(nz)     ! height of the interface levels
real presi(nz)  ! pressure,mb at interface levels
real adz(nzm)   ! ratio of the grid spacing to dz for pressure levels
real adzw(nz)	! ratio of the grid spacing to dz for w levels
real grdf_x(nzm)! grid factor for eddy diffusion in x
real grdf_y(nzm)! grid factor for eddy diffusion in y
real grdf_z(nzm)! grid factor for eddy diffusion in z

real at, bt, ct ! coefficients for the Adams-Bashforth scheme 
real dt		! dynamical timestep
real dtn	! current dynamical timestep (can be smaller than dt)
real dt3(3) 	! dynamical timesteps for three most recent time steps
real dt_cup	! timestep for convection parameterization 
real time	! current time in sec.
real day0	! starting day (including fraction)
real day	! current day (including fraction)
real dtfactor   ! dtn/dt
        
integer nstep	! current number of performed time steps 
integer nstop   ! time step number to stop the integration
integer nelapse ! time step number to elapse before stoping
integer na, nb, nc ! indeces for swapping the rhs arrays for AB scheme
integer ncycle  ! number of subcycles over the dynamical timestep
integer icycle  ! current subcycle 
integer nadams	! the order of the AB scheme (should be kept at 3)        
integer nstat	! the interval in time steps to compute statistics
integer nstatis	! the interval between substeps to compute statistics
integer nstatfrq! frequency of computing statistics 
integer nprint 	! frequency of printing a listing (steps)
integer nrestart! switch to control starting/restarting of the model
logical restart_sep ! write separate restart files for sub-domains
integer nrad	! frequency of calling the radiation routines
logical save3Dbin ! save 3D data in binary format(no byte compression)
logical save3Dsep !write a separate file for snapshots for 2-D model  
logical save2Dsep !write a separate file for 2D horizontal fields for 3-D model  
integer nsave3D ! frequency of writting 3D fields (steps)
integer nsave3Dstart ! timestep to start writting 3D fields
integer nsave3Dend   ! timestep to end writting 3D fields
real    qnsave3D !threshold manimum cloud water(kg/kg) to save 3D fields
logical dogzip3D ! gzip compress a 3D output file   
logical dogzip2D ! gzip compress a 2D output file if save2Dsep=.true.   
logical save2Dbin !save 2D data in binary format, rather than compressed
integer nsave2D ! frequency of writting 2D fields (steps)
integer nsave2Dstart ! timestep to start writting 2D fields
integer nsave2Dend   ! timestep to end writting 2D fields
character *40 caseid! 8-symbol id-string to identify a run	
character *40 case  ! 8-symbol id-string to identify a case-name	
logical dostatis! flag to permit the gathering of statistics
logical dostatisrad ! flag to permit the gathering of radiation statistics
integer nensemble ! the number of subensemble set of perturbations
logical notopened2D ! flag to see if the 2D output datafile is opened	
logical notopened3D ! flag to see if the 3D output datafile is opened	
character *256 raddir ! path to data directory containing files for CAM radiation
character *256 rf_filename ! input file for random forest
character *256 rf_filename_tkh ! input file for random forest
character *256 rf_filename_diffusion ! input file for random forest


!   Flags:

logical CEM     ! flag for Cloud Ensemble Model
logical LES     ! flag for Large-Eddy Simulation
logical OCEAN   ! flag indicating that surface is water
logical LAND    ! flag indicating that surface is land
logical SFC_FLX_FXD  ! surface sensible flux is fixed
logical SFC_TAU_FXD! surface drag is fixed
logical SFC_TAU_MASK! surface drag is masked

!       Multitasking staff:     
          
integer rank   ! rank of the current subdomain task (default 0) 
integer ranknn ! rank of the "northern" subdomain task
integer rankss ! rank of the "southern" subdomain task
integer rankee ! rank of the "eastern"  subdomain task
integer rankww ! rank of the "western"  subdomain task
integer rankne ! rank of the "north-eastern" subdomain task
integer ranknw ! rank of the "north-western" subdomain task
integer rankse ! rank of the "south-eastern" subdomain task
integer ranksw ! rank of the "south-western" subdomain task
logical dompi  ! logical switch to do multitasking
logical masterproc ! .true. if rank.eq.0 
integer rad_lev_pred !Yani added
integer input_ver_dim
integer diffusivity_levels_rf
integer int1_dum
integer  int2_dum


!   Logical switches and flags:

logical   dodamping, doupperbound, docloud, doprecip, &
          dorandomforest, do_rf_diffusion, do_rf_q_surf_flux, rf_uses_rh, rf_uses_qp, &  !Yani
          do_wind_input_rf, do_z_diffusion_rf, do_q_T_surf_fluxes_rf, do_surf_wind_rf, do_q_surf_fluxes_out_rf, & !New parameters after matlab corrections.
          do_sedimentation_rf, do_fall_tend_rf, do_radiation_output_rf, do_flux_rf, do_hor_advection_rf, do_hor_diffusion_rf, &
          Tin_feature_rf, qin_feature_rf, Tin_z_grad_feature_rf, qin_z_grad_feature_rf, predict_tendencies_rf, do_tkh_z, do_tkh_xy, &
          do_vert_wind_input_rf, do_yin_input,  &
          dolongwave, doshortwave, dosgs, dosubsidence, &
          docoriolis, dosurface, dolargescale, doradforcing, &
          dosfcforcing, doradsimple, donudging_uv, donudging_tq, donudging_q, donudging_hadley, & 
          dosmagor, doscalar, doensemble, doxy, dowallx, dowally, docup, &
          docolumn, doperpetual, doseasons, doradhomo, dosfchomo,doradhomox, dosfchomox, &
          doisccp, dodynamicocean, dosolarconstant, donoevap, donoqpload

! For dosolarconstant simulations, allow solar constant and zenith
! angle to be specified individually
real solar_constant  ! solar constant (in W/m2)
real zenith_angle    ! zenith angle (in degrees)

!===================================================================
! UW ADDITIONS

integer nelapsemin ! number of minutes to elapse before stoping
logical dorestart_last ! only write restart file after last timestep
logical dooldrestart ! use restart file from old UW version of SAM
integer nwriterestart   ! peters write restart file every nwriterestart steps
                        ! instead of every nstat steps

logical dokruegermicro ! use Krueger microphysics with prognostic
                       ! equations for all microphysical variables
                       ! and water vapor.  

logical dokruegereffectiveradii ! use specification of effective radii from
                        ! Luo et al., Journal of the Atmospheric
                        ! Sciences, vol. 60, pp. 510--525.
                        ! cloud: 10 um, cloud ice: 25 um, snow: 75 um.

integer tracer_option ! tracer options

integer iyear	! current year (used for computing solar zenith angle)
        
! Parameters for walker circulation -- dynamic ocean
real hml        ! Ocean Mixed Layer Height in meters
real deltaS     ! Change in Ocean Mixed Layer Heating across domain (W/m^2)
real Szero      ! Baseline Ocean Mixed Layer Heating across domain (W/m^2)
double precision  sst_target ! Target horizontal-averaged SST (K)

! Parameters for walker circulation -- static ocean
real delta_sst  ! Amplitude of sinusoidal sst variation (K)

! For dynamicocean simulations, store filtered sst bias
double precision  sst_bias_filtered, sst_bias_accumulated

! Multiply CO2 by the following factor when dodoubleco2=.true.
! (Default behavior is to double CO2 levels in radiation calculation.)
real    co2factor, o3factor, nudging_q_factor

! Increase beta effect by some multiple (allows simulation of a large
! range of latitudes in a reasonably sized domain.
  real    betafactor,mpbetafactor,ravefactor,dxdzfactor
  logical dorave
! For periodic forcing
  real    day0_forcing, period_forcing !in days

! user-specified sst profile
real ssty(ny_gl)

! user-specified ocean cooling profile
real mlohflx(ny_gl)

! user-specified ocean cooling profile
real sfctaumask(ny)

! uw options for slab ocean
logical   dofixedoceancooling, resetsst, firststep

! uw options for radiation
logical   doequinox, dodoubleco2
logical   doannual           ! peters -- no diurnal, ann avg. rad
logical   doDARErad          ! peters -- use spatially varying lat/lon for
                             ! computing zenith angle.  Also, re-scale time
                             ! when calling zenith so diurnal cycle/seasons
                             ! evolve appropriately
logical   doDAREradlongitude ! peters -- make longitude periodic in x,
                             ! independent of dx, darefactor


! uw options for 2D output
logical   doclearsky, doinstantolr, doreflectivity, &      
          dosaturatedwvp, dowaveoutput, docolumnbudgets, & 
          do850mbarwinds, dooutputheights, dolowlevelthetae, &
          dorceintercomparison, doinstantoutput, dowaveenergetics, &
          docloudechoheights, docasefive

! uw options for tracers (kzm)
logical   dotrx,dotry,dotrz,dotrzz, dotro3, doactiveo3, &
	  doreset_tracer, dorestart_tracer,dotauz

! uw options for nudging, subsidence, forcing
logical   dosubsidence_tonly, domoistadiabaticnudging, doperiodic_forcing, do3Daverage

! uw options for beta plane runs (kzm)
logical   doydamping, doxdamping, dobetaplane, dowaterbubble,dobetafactor
logical daremicrophysics
logical domicroscaling,domicroscalingdiffusion,doperturb_realtime, dompinstantremoval,domicroscalingvariance,domicroscalinghydro
real  cmin,cmax
real hydrostatic_criterion, hydrostatic_width, pbl_diffusivity_enhancement, microscaling_power
integer nsmooth_nonhydro

! uw options for large-scale forcing
logical   douwforcing ! use one of the uw-developed schemes for
                      ! large-scale forcing

! uw option for specified radiative cooling in boundary layer
!              with subsidence profile prescribed to match bl depth.
logical   doblradforcing
real      blqrad ! boundary layer radiative cooling rate [K/day on input] 
real      blqmin ! q threshold [g/kg on input] below which no qrad is applied
real      blqmax ! q threshold [g/kg on input] above which blqrad is applied
                 !  - between blqmin and blqmax, qrad ramps up linearly
real      bldivg ! fixed divergence rate [1/s] used to compute wsub
                 ! in boundary layer subsidence velocity falls off to
                 ! zero as mean(q) --> blqmin 

! uw option for boundary layer simulation using interactive radiation
!              with uniform divergence in bl and uniform subsidence
!              above and nudging of moisture/temperature from 150m above bl
!              (Modeled after Wyant et al 1997) 
logical   doblforcing ! generic option for bl forcing.  
                      ! If doblforcing=.true in the prm file, implies doblsubs
logical   doblsubs ! user-specified subsidence above the inversion
real      blsubs ! fixed subsidence velocity [m/s] above boundary layer
                 ! subsidence velocity in bl increases linearly from zero
                 ! at the surface to blsubs at the inversion.

logical   dobldivg ! user-specified divergence with bl -- subsidence
                   ! above bl will vary with bl depth
real      bltauz ! relaxation timescale [s] for nudging in free troposphere.

logical   doblsubs2 ! fixed subsidence profile
real      blsubs2 ! fixed subsidence profile

! NOTE: THE FOLLOWING OPTION IS NOT CURRENTLY IMPLEMENTED.
! uw option for wtg simulations (based on Raymond & Zeng, 2005).
logical   dowtg  ! compute subsidence velocity to relax  temperature
                 !  profile back to that in snd file.
real      wtgtau ! relaxation timescale used in computation of w_wtg.

!kzm stuff having to do with linear wave in a column
 logical doevolvewavebg
 real wavenumber_factor, first_mode_speed,second_mode_speed
 integer nstartlinearwave,nsteplinearwavebg,nsteplinearwave
 real wsub_conv(nz), alfafactor,wavedampingtime,pbldampingfactor
 real tv_wavebg(nzm),t_wavebg(nzm),q_wavebg(nzm),heating_wavebg(nzm)
 real tmin_evap,evapfactor !kzm limit evaporation to temperatures above tmin_evap, and enhance evaporation by evapfactor
 real tv_wave(nzm),t_wave(nzm),q_wave(nzm),pres_wave(nzm),heating_wave(nzm),w_wave(nzm),wwave_conv(nzm)
 real ttend_wave(nzm),qtend_wave(nzm)
 real t_wave_local(nzm),q_wave_local(nzm),h_wave_local(nzm),t_target(nz)
 real t_wavebg_hist(nzm,10),q_wavebg_hist(nzm,10),heating_wavebg_hist(nzm,10)
 real heating_anom_hist(nzm,360)
 integer hanomcounter,nhwave_avg
 logical dosmoothwconv,dowavesubsidence_tonly,dosteadylinearwave,dolinearwavelid,dowavedampingt,doadvectbg,dozerowsub
 real qtendfactor
 integer wavecounter
logical doaddqvrestart,doparameterizedwave,dobulksfc
real second_mode_factor
!JAA 
logical dosubsidenceonbackground
logical dosubsidenceonobs

real    maxcfl     ! maximum cfl found in field at beginning of time step
real    maxtkz     ! maximum vertical diffusivity
real    maxndiff   ! maximum number of PBL subtimesteps
integer ncycle0    ! initial number of subcycles over the dynamical timestep
integer ncyclemax  ! maximum number of subcycles over the dynamical timestep
integer ncyclemin  ! minimum number of subcycles over the dynamical timestep
integer totalcycles ! total number of cycles

real    qvcond(nzm) ! water vapor sink due to cloud condensation [kg/kg/s]
real    tlcond(nzm) ! sensible energy source due to cloud condensation [K/s]
real    tlevp(nzm) ! sensible energy sink due to precip evaporation [K/s]

! option for cloud drop sedimentation
logical doclouddropsed ! cloud drop sedimentation in cloud() subroutine.

logical douwpbl       ! use Bretherton McCaa PBL scheme
integer uwpbl_ndiv    ! number of points to average for uwpbl scheme
logical douwcu        ! use UW shallow convection scheme
integer uwcu_ndiv     ! number of points to average for uwcu scheme
logical dowrfml       ! use WRF-like mixing lengths for Smagorinsky scheme
integer fixedndiff    ! fix NDIFF at this value (0 turns this off)

! peters
logical dosurfacefix  ! call surface three times, once at each grid?

!kzm 01/06/08 block_size: the number of points in a block used in sfchomox, block_size=1 meaning
!surface fluxes are homogeneous in x. This number needs to be smaller than the number of point in x per processor 
integer block_size

!only do homox within this range from the equator
integer hadley_lim

!large-scale mean surface wind, used in surface flux calculations
real sfcum

! END UW ADDITIONS
!===================================================================

!===================================================================

end module grid
