module nn_convection_flux_mod
    !! Code to perform the Convection Flux parameterisation
    !! and interface to the nn_cf_net Neural Net
    !! Reference: https://doi.org/10.1029/2020GL091363
    !! Also see YOG20: https://doi.org/10.1038/s41467-020-17142-3

!---------------------------------------------------------------------
! Libraries to use
use nn_cf_net_mod, only: nn_cf_net_init, net_forward, nn_cf_net_finalize
implicit none
private


!---------------------------------------------------------------------
! public interfaces
public  nn_convection_flux, nn_convection_flux_init, nn_convection_flux_finalize


!---------------------------------------------------------------------
! local/private data

! Neural Net parameters used throughout module
integer :: n_inputs, n_outputs
    !! Length of input/output vector to the NN
integer, parameter :: nrf = 30
    !! number of vertical levels the NN uses
integer, parameter :: nrfq = 29
    !! number of vertical levels the NN uses when boundary condition is set to 0
integer, parameter :: input_ver_dim = 48
    !! Set to 48 in setparm.f90 of SAM. Same as nz_gl??
logical :: do_init=.true.
    !! model initialisation is yet to be performed



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

! != unit K :: tbgmin
! real, parameter :: tbgmin = 253.16
!     !! Minimum temperature for cloud water.
!
! != unit K :: tbgmin
! real, parameter :: tbgmax = 273.16
!     !! Maximum temperature for cloud ice, K

! Misc. microphysics variables
real, parameter :: a_pr = 1./(tprmax-tprmin)
    !! Misc. microphysics variables

! real, parameter :: a_bg = 1./(tbgmax-tbgmin)
!     !! Misc. microphysics variables



!---------------------------------------------------------------------
! Functions and Subroutines

contains

    subroutine nn_convection_flux_init(nn_filename)
        !! Initialise the NN module

        character(len=1024), intent(in) :: nn_filename
            !! NetCDF filename from which to read model weights

        ! Initialise the Neural Net from file and get info
        call nn_cf_net_init(nn_filename, n_inputs, n_outputs)

        ! Set init flag as complete
        do_init = .false.

        ! These have been made module parameters
        ! Hardcoded magic numbers
        ! nrf is the fact results are supplied at lowest 30 half-model levels for sedimentation fluxes, and at 29 levels for fluxes
        ! (as flux at bottom boundary is zero).
        ! nrf = 30
        ! nrfq = 29

    end subroutine nn_convection_flux_init

    subroutine nn_convection_flux(t_i, q_i, qp_i, y_in, &
                                  rho, adz, tabs, dz, dtn, &
                                  t, q, qn, precsfc, prec_xy)
        !! Interface to the neural net that applies physical constraints and reshaping
        !! of variables.

        ! -----------------------------------
        ! Input Variables
        real, intent(in) :: y_in(:)
            !! Distance of column from equator (proxy for insolation and sfc albedo)
        ! Fields from beginning of time step used as NN inputs
        real, intent(in) :: t_i(:, :, :)
            !! Temperature
        real, intent(in) :: q_i(:, :, :)
            !! Non-precipitating water mixing ratio
        real, intent(in) :: qp_i (:, :, :)
            !! Precipitating water mixing ratio
        ! reference vertical profiles:
        != unit (kg / m**3) :: rho
        real, intent(in) :: rho(:)
            !! air density at pressure levels
        ! != unit mb :: pres
        ! real, intent(in) pres(nzm)
        !     !! pressure,mb at scalar levels
        real, intent(in) :: adz(:)
            !! ratio of the grid spacing to dz for pressure levels from grid.f90
        real, intent(in) :: tabs(:, :, :)
            !! absolute temperature
        real, intent(in) :: dz
            !! grid spacing in z direction for the lowest grid layer
        != unit s :: dtn
        real dtn
            !! current dynamical timestep (can be smaller than dt)


        ! -----------------------------------
        ! Output Variables
        ! Taken from vars.f90 in SAM
        ! Prognostic variables
        != unit J :: t
        real, intent(inout) :: t(:, :, :)
            !! moist static energy
        real, intent(inout) :: q(:, :, :)
            !! total water
        ! Diagnostic variables:
        real, intent(inout) :: qn(:, :, :)
            !! cloud water+cloud ice
        ! fluxes at the top and bottom of the domain:
        real, intent(inout) :: precsfc(:, :)
            !! surface precip. rate
        !  Horizontally varying stuff (as a function of xy)
        real, intent(inout) :: prec_xy(:, :)
            !! surface precipitation rate


        ! -----------------------------------
        ! Local Variables
        integer nx
            !! Number of x points in a subdomain
        integer ny
            !! Number of y points in a subdomain
        integer nzm
            !! Number of z points in a subdomain - 1

        real,   dimension(nrf) :: t_tendency_adv, q_tendency_adv, q_tendency_auto, &
                                  q_tendency_sed, t_tendency_auto
        real,   dimension(nrf) :: q_flux_sed, qp_flux_fall, t_tendency_sed, q_tend_tot
        real,   dimension(nrf) :: t_flux_adv, q_flux_adv, t_sed_flux, t_rad_rest_tend, &
                                  omp, fac ! Do not predict surface adv flux
        real,   dimension(size(t_i, 3)) :: qsat, irhoadz, irhoadzdz, irhoadzdz2

        ! -----------------------------------
        ! NN variables
        real(4), dimension(n_inputs) :: features
            !! Vector of input features for the NN
        real(4), dimension(n_outputs) :: outputs
            !! vector of output features from the NN

!         real :: omn
!             !! variable to store result of omegan function
        real :: rev_dz
            !! variable to store 1/dz factor
        integer  i, j, k, dim_counter, out_dim_counter

        nx = size(t_i, 1)
        ny = size(t_i, 2)
        nzm = size(t_i, 3)

        ! Check that we have initialised all of the variables.
        if (do_init) call error_mesg('NN has not yet been initialised using nn_convection_flux_init.')

        ! Define TODO Magical Mystery variables
        rev_dz = 1/dz
        do k=1,nzm
            irhoadz(k) = dtn/(rho(k)*adz(k)) ! Useful factor
            irhoadzdz(k) = irhoadz(k)/dz ! Note the time step
            irhoadzdz2(k) = irhoadz(k)/(dz*dz) ! Note the time step
        end do

        ! The NN operates on atmospheric columns, so loop over x and y coordinates in turn
        do j=1,ny
            do i=1,nx

                ! Initialize variables
                features = 0.
                outputs = 0.
                t_tendency_adv = 0.
                q_tendency_adv = 0.
                q_tendency_auto = 0.
                t_tendency_auto = 0.
                q_tendency_sed = 0.
                t_tendency_sed = 0.
                t_rad_rest_tend = 0.
                q_tend_tot = 0.
                t_flux_adv = 0.
                q_flux_adv = 0.
                q_flux_sed = 0.
                dim_counter = 0
                omp = 0.
                fac = 0.

                !-----------------------------------------------------
                ! Combine all features into one vector

                ! Add temperature as input feature
                features(dim_counter+1:dim_counter + input_ver_dim) = real(t_i(i,j,1:input_ver_dim),4)
                dim_counter = dim_counter + input_ver_dim

                ! Add non-precipitating water mixing ratio as input feature using random forest (rf) approach from earlier paper
                ! TODO Currently we do not use rf_uses_rh option, but may add it back in later
                ! if (rf_uses_rh) then
                ! ! If using generalised relative humidity convert non-precip. water content to rel. hum
                !     do k=1,nzm
                !         omn = omegan(tabs(i,j,k))
                !         qsat(k) = omn * qsatw(tabs(i,j,k),pres(k)) + (1.-omn) * qsati(tabs(i,j,k),pres(k))
                !     end do
                !     features(dim_counter+1:dim_counter+input_ver_dim) = real(q_i(i,j,1:input_ver_dim)/qsat(1:input_ver_dim),4)
                !     dim_counter = dim_counter + input_ver_dim
                ! else
                ! ! if using non-precipitating water as water content
                    features(dim_counter+1:dim_counter+input_ver_dim) = real(q_i(i,j,1:input_ver_dim),4)
                    dim_counter =  dim_counter + input_ver_dim
                ! endif

                ! Add distance to the equator as input feature
                ! y is a proxy for insolation and surface albedo as both are only a function of |y| in SAM
                ! TODO Will this be true in CAM/SCAM?
                features(dim_counter+1) = y_in(j)
                dim_counter = dim_counter+1

                !-----------------------------------------------------
                ! Call the forward method of the NN on the input features
                ! Scaling and normalisation done as layers in NN

                call net_forward(features, outputs, nrf)

                !-----------------------------------------------------
                ! Separate physical outputs from NN output vector and apply physical constraints

                ! TODO: The following code applies various physical constraints to the output of the NN
                ! These are likely to be required for any application, so should they be in the NN rather than the interface?
                ! If in the NN can they be generalised?

                ! Temperature rest tendency
                t_rad_rest_tend(1:nrf) = outputs(1:nrf)
                out_dim_counter = nrf

                ! Temperature flux
                ! BC: advection surface flux is zero
                t_flux_adv(1) = 0.0
                t_flux_adv(2:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrfq)
                out_dim_counter = out_dim_counter + nrfq

                ! Non-precip. water flux
                ! BC: advection surface flux is zero
                q_flux_adv(1) = 0.0
                q_flux_adv(2:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrfq)
                out_dim_counter = out_dim_counter + nrfq

                ! Non-precip. water must be >= 0, so ensure advective flux will not reduce it below 0
                do k=2,nrf
                    if (q_flux_adv(k).lt.0) then
                        if ( q(i,j,k).lt.-q_flux_adv(k)* irhoadzdz(k)) then
                            q_flux_adv(k) = -q(i,j,k)/irhoadzdz(k)
                        end if
                    else
                        ! If gaining water content ensure that we are not gaining more than is in cell below
                        ! TODO: Not sure I fully understant this constraint
                        if (q(i,j,k-1).lt.q_flux_adv(k)* irhoadzdz(k)) then
                            q_flux_adv(k) = q(i,j,k-1)/irhoadzdz(k)
                        end if
                    end if
                end do

                ! Calculate tendencies
                q_tendency_auto(1:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrf)
                out_dim_counter = out_dim_counter + nrf

                ! Convert flux to advective tendency via finite difference
                do k=1,nrf-1
                    t_tendency_adv(k) = - (t_flux_adv(k+1) - t_flux_adv(k)) * irhoadzdz(k)
                    q_tendency_adv(k) = - (q_flux_adv(k+1) - q_flux_adv(k)) * irhoadzdz(k)
                end do
                ! Apply finite difference boundary condition to advective tendencies
                t_tendency_adv(nrf) = - (0.0 - t_flux_adv(nrf)) * irhoadzdz(nrf)
                q_tendency_adv(nrf) = - (0.0 - q_flux_adv(nrf)) * irhoadzdz(nrf)
                ! q must be >= 0 so ensure tendency won't reduce below zero
                do k=1,nrf
                    if (q(i,j,k) .lt. -q_tendency_adv(k)) then
                        q_tendency_adv(k) = -q(i,j,k)
                    end if
                end do
                ! Apply advective tendencies to variables
                t(i,j,1:nrf) = t(i,j,1:nrf) + t_tendency_adv(1:nrf)
                q(i,j,1:nrf) = q(i,j,1:nrf) + q_tendency_adv(1:nrf)

                ! TODO ???
                do k=1,nrf
                    omp(k) = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
                    fac(k) = (fac_cond + fac_fus * (1.0 - omp(k)))
                    if (q_tendency_auto(k).lt.0) then
                        q_tend_tot(k) = min(-q_tendency_auto(k) * dtn, q(i,j,k))!q_tendency_auto(1:nrf) * dtn + q_tendency_adv(1:nrf) + q_tendency_sed(1:nrf)
                        q_tend_tot(k) = -q_tend_tot(k)
                    else
                        q_tend_tot(k) = q_tendency_auto(k) * dtn
                    endif
                end do

                ! Update q and t
                q(i,j,1:nrf) = q(i,j,1:nrf) + q_tend_tot(1:nrf)
                t(i,j,1:nrf) = t(i,j,1:nrf) - q_tend_tot(1:nrf)*fac(1:nrf)

                ! TODO ???
                q_flux_sed(1:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrf)
                ! q_flux must be >= 0
                do k=2,nrf
                    if (q_flux_sed(k).lt.0) then
                        if ( q(i,j,k).lt.-q_flux_sed(k)* irhoadzdz(k)) then
                            q_flux_sed(k) = -q(i,j,k)/irhoadzdz(k)
                        end if
                    else
                        if (q(i,j,k-1).lt.q_flux_sed(k)* irhoadzdz(k)) then
                            q_flux_sed(k) = q(i,j,k-1)/irhoadzdz(k)
                        end if
                    end if
                end do

                ! Calculate sed tendency via finite difference
                do k=1,nrf-1 ! One level less than I actually use
                    q_tendency_sed(k) = - (q_flux_sed(k+1) - q_flux_sed(k)) * irhoadzdz(k)
                end do
                ! Set value at top of nrf layr
                q_tendency_sed(nrf) = - (0.0 - q_flux_sed(nrf)) * irhoadzdz(nrf)

                ! If q sed tendency < 0 ensure it will not reduce q below 0
                do k=1,nrf
                    if (q_tendency_sed(k).lt.0) then
                        q_tendency_sed(k) = min(-q_tendency_sed(k), q(i,j,k))
                        q_tendency_sed(k) = -q_tendency_sed(k)
                    end if
                end do

                ! Apply sed tendency to variables
                t(i,j,1:nrf) = t(i,j,1:nrf) - q_tendency_sed(1:nrf)*(fac_fus+fac_cond)
                q(i,j,1:nrf) = q(i,j,1:nrf) +q_tendency_sed(1:nrf)

                ! Apply rest tendency to variables
                t(i,j,1:nrf) = t(i,j,1:nrf) + t_rad_rest_tend(1:nrf)*dtn

                ! Calculate surface precipitation
                ! Apply sed flux at surface
                precsfc(i,j) = precsfc(i,j)  - q_flux_sed(1)*dtn*rev_dz ! For statistics
                prec_xy(i,j) = prec_xy(i,j)  - q_flux_sed(1)*dtn*rev_dz ! For 2D output
                !
                do k=1, nrf
                    ! TODO: Both these lines have dz*(1/dz) can we just remove this?
                    precsfc(i,j) = precsfc(i,j) - q_tend_tot(k)*adz(k)*dz*rho(k)*rev_dz! removed the time step mult because q_tend_tot is already mult
                    prec_xy(i,j) = prec_xy(i,j) - q_tend_tot(k)*adz(k)*dz*rho(k)*rev_dz
                end do

                ! As a final check q must be >= 0.0, if not then set to 0.0, otherwise add tendencies
                do k = 1,nrf
                    q(i,j,k) = max(0.,q(i,j,k))
                end do

                ! qn (precip.) must be >= 0.0, if not then set to 0.0, otherwise add tendencies
                where (qn(i,j,1:nrf) .gt. 0.0)
                    qn(i,j,1:nrf) = qn(i,j,1:nrf) + q_tend_tot(1:nrf) + q_tendency_adv(1:nrf) + q_tendency_sed(1:nrf)
                end where
                where (qn(i,j,:) .lt. 0.0)
                    qn(i,j,:) = 0.0
                end where

            end do
        end do
        ! End of loop over x, y, columns

    end subroutine nn_convection_flux

    subroutine nn_convection_flux_finalize()
        !! Finalize the NN module

        ! Finalize the Neural Net deallocating arrays
        call nn_cf_net_finalize()

    end subroutine nn_convection_flux_finalize

    !###################################################################

    subroutine error_mesg (message)
      character(len=*), intent(in) :: message
          !! message to be written to output   (character string)

      ! TODO Since masterproc is a SAM variable adjust to print from all proc for now
      ! if(masterproc) print*, 'Neural network  module: ', message
      print*, 'Neural network  module: ', message
      stop

    end subroutine error_mesg

!     !###################################################################
!     ! Need omegan function to run with rf_uses_rh option (currently unused)
!     ! Ripped from SAM model:
!     ! https://github.com/yaniyuval/Neural_nework_parameterization/blob/f81f5f695297888f0bd1e0e61524590b4566bf03/sam_code_NN/omega.f90
! 
!     != unit 1 :: omegan
!     real function omegan(tabs)
!         != unit K :: tabs
!         real, intent(in) :: tabs
!             !! Absolute temperature
! 
!         omegan = max(0., min(1., (tabs-tbgmin)*a_bg))
! 
!         return
! 
!     end function omegan
! 
!     !###################################################################
!     ! Need qsatw functions to run with rf_uses_rh option (currently unused)
!     ! Ripped from SAM model:
!     ! https://github.com/yaniyuval/Neural_nework_parameterization/blob/f81f5f695297888f0bd1e0e61524590b4566bf03/sam_code_NN/sat.f90
! 
!     ! Saturation vapor pressure and mixing ratio.
!     ! Based on Flatau et.al, (JAM, 1992:1507)
! 
!     != unit mb :: esatw
!     real function esatw(t)
!       implicit none
!       != unit K :: t
!       real :: t  ! temperature (K)
! 
!       != unit :: a0
!       != unit :: mb / k :: a1, a2, a3, a4, a5, a6, a7, a8
!       real :: a0,a1,a2,a3,a4,a5,a6,a7,a8
!       data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
!               6.105851, 0.4440316, 0.1430341e-1, &
!               0.2641412e-3, 0.2995057e-5, 0.2031998e-7, &
!               0.6936113e-10, 0.2564861e-13,-0.3704404e-15/
!       !       6.11239921, 0.443987641, 0.142986287e-1, &
!       !       0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
!       !       0.640689451e-10, -0.952447341e-13,-0.976195544e-15/
! 
!       != unit K :: dt
!       real :: dt
!       dt = max(-80.,t-273.16)
!       esatw = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
!     end function esatw
! 
!     != unit 1 :: qsatw
!     real function qsatw(t,p)
!       implicit none
!       != unit K :: t
!       real :: t  ! temperature
! 
!       != unit mb :: p, esat
!       real :: p  ! pressure
!       real :: esat
! 
!       esat = esatw(t)
!       qsatw = 0.622 * esat/max(esat, p-esat)
!     end function qsatw
! 
!     ! real function dtesatw(t)
!     !   implicit none
!     !   real :: t  ! temperature (K)
!     !   real :: a0,a1,a2,a3,a4,a5,a6,a7,a8
!     !   data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
!     !             0.443956472, 0.285976452e-1, 0.794747212e-3, &
!     !             0.121167162e-4, 0.103167413e-6, 0.385208005e-9, &
!     !            -0.604119582e-12, -0.792933209e-14, -0.599634321e-17/
!     !   real :: dt
!     !   dt = max(-80.,t-273.16)
!     !   dtesatw = a0 + dt* (a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
!     ! end function dtesatw
!     !
!     ! real function dtqsatw(t,p)
!     !   implicit none
!     !   real :: t  ! temperature (K)
!     !   real :: p  ! pressure    (mb)
!     !   real :: dtesatw
!     !   dtqsatw=0.622*dtesatw(t)/p
!     ! end function dtqsatw
! 
!     real function esati(t)
!       implicit none
!       != unit K :: t
!       real :: t  ! temperature
!       real :: a0,a1,a2,a3,a4,a5,a6,a7,a8
!       data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
!               6.11147274, 0.503160820, 0.188439774e-1, &
!               0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
!               0.385852041e-9, 0.146898966e-11, 0.252751365e-14/
!       real :: dt
!       dt = max(-80.0, t-273.16)
!       esati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
!     end function esati
! 
!     != unit 1 :: qsati
!     real function qsati(t,p)
!       implicit none
!       != unit t :: K
!       real :: t  ! temperature
! 
!       != unit mb :: p
!       real :: p  ! pressure
! 
!       != unit mb :: esat
!       real :: esat
!       esat = esati(t)
!       qsati = 0.622 * esat/max(esat,p-esat)
!     end function qsati
! 
!     ! function dtesati(t)
!     !   implicit none
!     !   real :: t  ! temperature (K)
!     !   real :: a0,a1,a2,a3,a4,a5,a6,a7,a8
!     !   data a0,a1,a2,a3,a4,a5,a6,a7,a8 / &
!     !           0.503223089, 0.377174432e-1,0.126710138e-2, &
!     !           0.249065913e-4, 0.312668753e-6, 0.255653718e-8, &
!     !           0.132073448e-10, 0.390204672e-13, 0.497275778e-16/
!     !   real :: dtesati, dt
!     !   dt = max(-800. ,t-273.16)
!     !   dtesati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
!     !   return
!     ! end function dtesati
!     !
!     !
!     ! real function dtqsati(t,p)
!     !   implicit none
!     !   real :: t  ! temperature (K)
!     !   real :: p  ! pressure    (mb)
!     !   real :: dtesati
!     !   dtqsati = 0.622 * dtesati(t) / p
!     ! end function dtqsati

end module nn_convection_flux_mod
