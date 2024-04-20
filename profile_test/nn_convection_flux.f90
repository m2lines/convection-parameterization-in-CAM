module nn_convection_flux_mod
    !! Code to perform the Convection Flux parameterisation
    !! and interface to the nn_cf_net Neural Net
    !! Reference: https://doi.org/10.1029/2020GL091363
    !! Also see YOG20: https://doi.org/10.1038/s41467-020-17142-3

!---------------------------------------------------------------------
! Libraries to use
use nn_cf_net_mod, only: nn_cf_net_init, net_forward, nn_cf_net_finalize
use SAM_consts_mod, only: fac_cond, fac_fus, tprmin, a_pr, input_ver_dim, &
                          nrf, nrfq, dt_sam
use nf_out, only: nf_write_sam

implicit none
private


!---------------------------------------------------------------------
! public interfaces
public  nn_convection_flux, nn_convection_flux_init, nn_convection_flux_finalize, &
        esati, qsati, qsatw, dtqsatw, dtqsati


!---------------------------------------------------------------------
! local/private data

! Neural Net parameters used throughout module

integer :: n_inputs, n_outputs
    !! Length of input/output vector to the NN

logical :: do_init=.true.
    !! model initialisation is yet to be performed


!---------------------------------------------------------------------
! Functions and Subroutines

contains

    !-----------------------------------------------------------------
    ! Public Subroutines

    subroutine nn_convection_flux_init(nn_filename)
        !! Initialise the NN module

        character(len=136), intent(in) :: nn_filename
            !! NetCDF filename from which to read model weights

        ! Initialise the Neural Net from file and get info
        call nn_cf_net_init(nn_filename, n_inputs, n_outputs, nrf)

        ! Set init flag as complete
        do_init = .false.

    end subroutine nn_convection_flux_init


    subroutine nn_convection_flux(tabs_i, q_i, y_in, &
                                  tabs, &
                                  t, q, &
                                  rho, adz, dz, dtn, &
                                  prec_sed)
        !! Interface to the neural net that applies physical constraints and reshaping
        !! of variables.
        !! Operates on subcycle of timestep dtn to update t, q, precsfc, and prec_xy

        ! -----------------------------------
        ! Input Variables
        ! -----------------------------------

        ! ---------------------
        ! Fields from beginning of time step used as NN inputs
        ! ---------------------
        != unit s :: tabs_i
        real(8), intent(in) :: tabs_i(:, :)
            !! Temperature

        != unit 1 :: q_i
        real(8), intent(in) :: q_i(:, :)
            !! Non-precipitating water mixing ratio

        ! ---------------------
        ! Other fields from SAM
        ! ---------------------
        != unit m :: y_in
        real(8), intent(in) :: y_in(:)
            !! Distance of column from equator (proxy for insolation and sfc albedo)

        != unit K :: tabs
        real(8), intent(in) :: tabs(:, :)
            !! absolute temperature
        
        != unit (J / kg) :: t
        real(8), intent(inout) :: t(:, :)
            !! Liquid Ice static energy (cp*T + g*z − L(qliq + qice) − Lf*qice)
        
        != unit 1 :: q
        real(8), intent(inout) :: q(:, :)
            !! total water

        ! ---------------------
        ! reference vertical profiles:
        ! ---------------------
        != unit (kg / m**3) :: rho
        real(8), intent(in) :: rho(:)
            !! air density at pressure levels
        
        ! != unit mb :: pres
        ! real(8), intent(in) pres(nzm)
        !     !! pressure,mb at scalar levels
        
        != unit 1 :: adz
        real(8), intent(in) :: adz(:)
            !! ratio of the pressure level grid height spacing [m] to dz (lowest dz spacing)
        
        ! ---------------------
        ! Single value parameters from model/grid
        ! ---------------------
        != unit m :: dz
        real(8), intent(in) :: dz
            !! grid spacing in z direction for the lowest grid layer

        != unit s :: dtn
        real(8), intent(in) :: dtn
            !! current dynamical timestep (can be smaller than dt due to subcycling)

        ! -----------------------------------
        ! Output Variables
        ! -----------------------------------
        != unit (J / kg) :: t_delta_adv, t_delta_auto, t_delta_sed
        != unit 1 :: q_delta_adv, q_delta_auto, q_delta_sed
        real(8), dimension(size(tabs_i, 1), size(tabs_i, 2)) :: t_delta_adv, q_delta_adv, &
                                                             t_delta_auto, q_delta_auto, &
                                                             t_delta_sed, q_delta_sed
            !! delta values of t and q generated by the NN

        != unit kg :: t_rad_rest_tend
        real(8), dimension(size(tabs_i, 1), size(tabs_i, 2)) :: t_rad_rest_tend
            !! tendency of t generated by the NN
        
        != unit kg :: prec_sed
        real(8), intent(out), dimension(:)   :: prec_sed
            !! Surface precipitation due to sedimentation

        ! -----------------------------------
        ! Local Variables
        ! -----------------------------------
        integer  i, k, dim_counter, out_dim_counter
        integer nx
            !! Number of x points in a subdomain
        integer nzm
            !! Number of z points in a subdomain - 1
        ! real(8) :: omn
        !     !! variable to store result of omegan function
        !     !! Note that if using you will need to import omegan() from SAM_consts_mod

        ! Other variables
        real(8),   dimension(nrf) :: omp, fac
        real(8),   dimension(size(tabs_i, 2)) :: qsat, irhoadz, irhoadzdz

        ! -----------------------------------
        ! variables for NN
        ! -----------------------------------
        real(4), dimension(n_inputs) :: features
            !! Vector of input features for the NN
        real(4), dimension(n_outputs) :: outputs
            !! vector of output features from the NN
        ! NN outputs
        real(8),   dimension(nrf) :: t_flux_adv, q_flux_adv, q_tend_auto, &
                                  q_sed_flux
        ! Output variable t_rad_rest_tend is also an output from the NN (defined above)

        nx = size(tabs_i, 1)
        nzm = size(tabs_i, 2)

        ! Check that we have initialised all of the variables.
        if (do_init) call error_mesg('NN has not yet been initialised using nn_convection_flux_init.')

        ! Define useful variables relating to grid spacing to convert fluxes to tendencies
        do k=1,nzm
            irhoadz(k) = dtn/(rho(k)*adz(k)) !  Temporary factor for below
            irhoadzdz(k) = irhoadz(k)/dz ! 2.0 * dtn / (rho(k)*(z(k+1) - z(k-1))) [(kg.m/s)^-1]
        end do

        ! The NN operates on atmospheric columns which have been flattened into 2D
        do i=1,nx
            ! Initialize variables
            features = 0.
            dim_counter = 0
            outputs = 0.
            t_rad_rest_tend(i,:) = 0.
            t_flux_adv = 0.
            q_flux_adv = 0.
            t_delta_adv(i,:) = 0.
            q_delta_adv(i,:) = 0.
            q_tend_auto = 0.
            t_delta_auto(i,:) = 0.
            q_delta_auto(i,:) = 0.
            q_sed_flux = 0.
            t_delta_sed(i,:) = 0.
            q_delta_sed(i,:) = 0.
            prec_sed(i) = 0.
            omp = 0.
            fac = 0.

            !-----------------------------------------------------
            ! Combine all features into one vector

            ! Add temperature as input feature
            features(dim_counter+1:dim_counter + input_ver_dim) = real(tabs_i(i,1:input_ver_dim),4)
            dim_counter = dim_counter + input_ver_dim

            ! Add non-precipitating water mixing ratio as input feature using random forest (rf) approach from earlier paper
            ! Currently we do not use rf_uses_rh option, but may add it back in later
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
                features(dim_counter+1:dim_counter+input_ver_dim) = real(q_i(i,1:input_ver_dim),4)
                dim_counter =  dim_counter + input_ver_dim
            ! endif

            ! Add distance to the equator as input feature
            ! y is a proxy for insolation and surface albedo as both are only a function of |y| in SAM
            features(dim_counter+1) = y_in(i)
            dim_counter = dim_counter+1

            !-----------------------------------------------------
            ! Call the forward method of the NN on the input features
            ! Scaling and normalisation done as layers in NN

            call net_forward(features, outputs)

            !-----------------------------------------------------
            ! Separate physical outputs from NN output vector

            ! Moist Static Energy radiative + rest(microphysical changes) tendency
            t_rad_rest_tend(i,1:nrf) = outputs(1:nrf)
            out_dim_counter = nrf

            ! Moist Static Energy advective subgrid flux
            ! BC: advection surface flux is zero
            t_flux_adv(1) = 0.0
            t_flux_adv(2:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrfq)
            out_dim_counter = out_dim_counter + nrfq

            ! Total non-precip. water mix. ratio advective flux
            ! BC: advection surface flux is zero
            q_flux_adv(1) = 0.0
            q_flux_adv(2:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrfq)
            out_dim_counter = out_dim_counter + nrfq

            ! Total non-precip. water autoconversion tendency
            q_tend_auto(1:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrf)
            out_dim_counter = out_dim_counter + nrf

            ! total non-precip. water mix. ratio ice-sedimenting flux
            q_sed_flux(1:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrf)
            
            call nf_write_sam(q_flux_adv(:), "Q_FLX_ADV")
            call nf_write_sam(t_flux_adv(:), "T_FLX_ADV")
            call nf_write_sam(q_tend_auto(:), "Q_TEND_AUTO")
            call nf_write_sam(-q_tend_auto(:)*fac(1:nrf), "T_TEND_AUTO")
            call nf_write_sam(q_sed_flux(:), "Q_FLX_SED")
            call nf_write_sam(-q_sed_flux(:)*(fac_fus+fac_cond), "T_FLX_SED")

            !-----------------------------------------------------
            ! Apply physical constraints and update q and t

            ! Non-precip. water content must be >= 0, so ensure advective fluxes
            ! will not reduce it below 0 anywhere
            do k=2,nrf
                if (q_flux_adv(k).lt.0) then
                    ! If flux is negative ensure we don't lose more than is already present
                    if ( q(i,k).lt.-q_flux_adv(k)* irhoadzdz(k)) then
                        q_flux_adv(k) = -q(i,k)/irhoadzdz(k)
                    end if
                else
                    ! If flux is positive ensure we don't gain more than is in the box below
                    if (q(i,k-1).lt.q_flux_adv(k)* irhoadzdz(k)) then
                        q_flux_adv(k) = q(i,k-1)/irhoadzdz(k)
                    end if
                end if
            end do

            ! Convert advective fluxes to deltas
            do k=1,nrf-1
                t_delta_adv(i,k) = - (t_flux_adv(k+1) - t_flux_adv(k)) * irhoadzdz(k)
                q_delta_adv(i,k) = - (q_flux_adv(k+1) - q_flux_adv(k)) * irhoadzdz(k)
            end do
            ! Enforce boundary condition at top of column
            t_delta_adv(i,nrf) = - (0.0 - t_flux_adv(nrf)) * irhoadzdz(nrf)
            q_delta_adv(i,nrf) = - (0.0 - q_flux_adv(nrf)) * irhoadzdz(nrf)
            ! q must be >= 0 so ensure delta won't reduce it below zero
            do k=1,nrf
                if (q(i,k) .lt. -q_delta_adv(i,k)) then
                    q_delta_adv(i,k) = -q(i,k)
                end if
            end do

            ! Update q and t with delta values
            q(i,1:nrf) = q(i,1:nrf) + q_delta_adv(i,1:nrf)
            t(i,1:nrf) = t(i,1:nrf) + t_delta_adv(i,1:nrf)

            call nf_write_sam(q_delta_adv(1,:), "DQ_ADV")
            call nf_write_sam(t_delta_adv(1,:), "DT_ADV")

            ! ensure autoconversion tendency won't reduce q below 0
            do k=1,nrf
                omp(k) = max(0.,min(1.,(tabs(i,k)-tprmin)*a_pr))
                fac(k) = (fac_cond + fac_fus * (1.0 - omp(k)))
                if (q_tend_auto(k).lt.0) then
                    q_delta_auto(i,k) = - min(-q_tend_auto(k) * dtn, q(i,k))
                else
                    q_delta_auto(i,k) = q_tend_auto(k) * dtn
                endif
            end do

            ! Update with autoconversion q and t deltas (dt = -dq*(latent_heat/cp))
            q(i,1:nrf) = q(i,1:nrf) + q_delta_auto(i,1:nrf)
            t_delta_auto(i,1:nrf) = - q_delta_auto(i,1:nrf)*fac(1:nrf)
            t(i,1:nrf) = t(i,1:nrf) + t_delta_auto(i,1:nrf)

            call nf_write_sam(q_delta_auto(1,:), "DQ_AUTO")
            call nf_write_sam(t_delta_auto(1,:), "DT_AUTO")

            ! Ensure sedimenting ice will not reduce q below zero anywhere
            do k=2,nrf
                if (q_sed_flux(k).lt.0) then
                    ! If flux is negative ensure we don't lose more than is already present
                    if ( q(i,k).lt.-q_sed_flux(k)* irhoadzdz(k)) then
                        q_sed_flux(k) = -q(i,k)/irhoadzdz(k)
                    end if
                else
                    ! If flux is positive ensure we don't gain more than is in the box below
                    if (q(i,k-1).lt.q_sed_flux(k)* irhoadzdz(k)) then
                        q_sed_flux(k) = q(i,k-1)/irhoadzdz(k)
                    end if
                end if
            end do

            ! Convert sedimenting fluxes to deltas
            do k=1,nrf-1 ! One level less than I actually use
                q_delta_sed(i,k) = - (q_sed_flux(k+1) - q_sed_flux(k)) * irhoadzdz(k)
            end do
            ! Enforce boundary condition at top of column
            q_delta_sed(i,nrf) = - (0.0 - q_sed_flux(nrf)) * irhoadzdz(nrf)
            ! q must be >= 0 so ensure delta won't reduce it below zero
            do k=1,nrf
                if (q_delta_sed(i,k).lt.0) then
                    q_delta_sed(i,k) = min(-q_delta_sed(i,k), q(i,k))
                    q_delta_sed(i,k) = -q_delta_sed(i,k)
                end if
            end do

            ! Update q and t with sed q and t deltas (dt = -dq*(latent_heat/cp))
            q(i,1:nrf) = q(i,1:nrf) + q_delta_sed(i,1:nrf)
            t(i,1:nrf) = t(i,1:nrf) - q_delta_sed(i,1:nrf)*(fac_fus+fac_cond)

            call nf_write_sam(q_delta_sed(1,:), "DQ_SED")
            call nf_write_sam(-q_delta_sed(1,:)*(fac_fus+fac_cond), "DT_SED")

            ! Apply radiation rest tendency to variables (multiply by dtn to get dt)
            t(i,1:nrf) = t(i,1:nrf) + t_rad_rest_tend(i,1:nrf)*dtn

            call nf_write_sam(t_rad_rest_tend(i,1:nrf)*dtn, "DT_RAD")

            ! Calculate surface precipitation
            ! Apply sedimenting flux at surface to get rho*dq term
            prec_sed(i) = - q_sed_flux(1)*dtn/dz
            
            ! As a final check enforce q must be >= 0.0
            do k = 1,nrf
                q(i,k) = max(0.,q(i,k))
            end do
        end do
        ! End of loop over columns

    end subroutine nn_convection_flux


    subroutine nn_convection_flux_finalize()
        !! Finalize the NN module

        ! Finalize the Neural Net deallocating arrays
        call nn_cf_net_finalize()

    end subroutine nn_convection_flux_finalize

    
    !-----------------------------------------------------------------
    
    subroutine error_mesg (message)
      character(len=*), intent(in) :: message
          !! message to be written to output   (character string)

      ! Since masterproc is a SAM variable adjust to print from all proc for now
      ! if(masterproc) print*, 'Neural network  module: ', message
      print*, 'Neural network  module: ', message
      stop

    end subroutine error_mesg


    !-----------------------------------------------------------------
    ! Need qsatw functions to:
    !     - run with rf_uses_rh option (currently unused)
    !     - convert variables in CAM interface
    ! Ripped from SAM model:
    ! https://github.com/yaniyuval/Neural_nework_parameterization/blob/f81f5f695297888f0bd1e0e61524590b4566bf03/sam_code_NN/sat.f90

    ! Saturation vapor pressure and mixing ratio.
    ! Based on Flatau et.al, (JAM, 1992:1507)

    != unit mb :: esatw
    real(8) function esatw(t)
      implicit none
      != unit K :: t
      real(8) :: t  ! temperature (K)

      != unit :: a0
      != unit :: mb / k :: a1, a2, a3, a4, a5, a6, a7, a8
      real(8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
              6.105851, 0.4440316, 0.1430341e-1, &
              0.2641412e-3, 0.2995057e-5, 0.2031998e-7, &
              0.6936113e-10, 0.2564861e-13,-0.3704404e-15/
      !       6.11239921, 0.443987641, 0.142986287e-1, &
      !       0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
      !       0.640689451e-10, -0.952447341e-13,-0.976195544e-15/

      != unit K :: dt
      real(8) :: dt
      dt = max(-80.,t-273.16)
      esatw = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    end function esatw


    != unit 1 :: qsatw
    real(8) function qsatw(t,p)
      implicit none
      != unit K :: t
      real(8) :: t  ! temperature

      != unit mb :: p, esat
      real(8) :: p  ! pressure
      real(8) :: esat

      esat = esatw(t)
      qsatw = 0.622 * esat/max(esat, p-esat)
    end function qsatw


    real(8) function dtesatw(t)
      implicit none
      real(8) :: t  ! temperature (K)
      real(8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
                0.443956472, 0.285976452e-1, 0.794747212e-3, &
                0.121167162e-4, 0.103167413e-6, 0.385208005e-9, &
               -0.604119582e-12, -0.792933209e-14, -0.599634321e-17/
      real(8) :: dt
      dt = max(-80.,t-273.16)
      dtesatw = a0 + dt* (a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    end function dtesatw


    real(8) function dtqsatw(t,p)
      implicit none
      real(8) :: t  ! temperature (K)
      real(8) :: p  ! pressure    (mb)
      dtqsatw=0.622*dtesatw(t)/p
    end function dtqsatw


    real(8) function esati(t)
      implicit none
      != unit K :: t
      real(8) :: t  ! temperature
      real(8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
              6.11147274, 0.503160820, 0.188439774e-1, &
              0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
              0.385852041e-9, 0.146898966e-11, 0.252751365e-14/
      real(8) :: dt
      dt = max(-80.0, t-273.16)
      esati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    end function esati


    != unit 1 :: qsati
    real(8) function qsati(t,p)
      implicit none
      != unit t :: K
      real(8) :: t  ! temperature

      != unit mb :: p
      real(8) :: p  ! pressure

      != unit mb :: esat
      real(8) :: esat
      esat = esati(t)
      qsati = 0.622 * esat/max(esat,p-esat)
    end function qsati


    real(8) function dtesati(t)
      implicit none
      real(8) :: t  ! temperature (K)
      real(8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 / &
              0.503223089, 0.377174432e-1,0.126710138e-2, &
              0.249065913e-4, 0.312668753e-6, 0.255653718e-8, &
              0.132073448e-10, 0.390204672e-13, 0.497275778e-16/
      real(8) :: dt
      dt = max(-800. ,t-273.16)
      dtesati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
      return
    end function dtesati


    real(8) function dtqsati(t,p)
      implicit none
      real(8) :: t  ! temperature (K)
      real(8) :: p  ! pressure    (mb)
      dtqsati = 0.622 * dtesati(t) / p
    end function dtqsati


end module nn_convection_flux_mod