module nn_interface_CAM_mod
    !! Interface to convection parameterisation for the CAM model
    !! Reference: https://doi.org/10.1029/2020GL091363
    !! Also see YOG20: https://doi.org/10.1038/s41467-020-17142-3

!---------------------------------------------------------------------
! Libraries to use
use nn_convection_flux_mod, only: nn_convection_flux, &
                                  nn_convection_flux_init, nn_convection_flux_finalize
implicit none
private


!---------------------------------------------------------------------
! public interfaces
public  nn_convection_flux_CAM, &
        nn_convection_flux_CAM_init, nn_convection_flux_CAM_finalize


!---------------------------------------------------------------------
! local/private data

! Copied from nn_convection_flux.f90
! Outputs from NN are supplied at lowest 30 half-model levels
integer, parameter :: nrf = 30
    !! number of vertical levels the NN parameterisation uses

!---------------------------------------------------------------------
! Functions and Subroutines

contains

    !-----------------------------------------------------------------
    ! Public Subroutines

    subroutine nn_convection_flux_CAM_init(nn_filename)
        !! Initialise the NN module

        character(len=1024), intent(in) :: nn_filename
            !! NetCDF filename from which to read model weights

        ! Initialise the Neural Net from file and get info
        call nn_convection_flux_init(nn_filename)

    end subroutine nn_convection_flux_CAM_init


    subroutine nn_convection_flux_CAM(tabs_i, q_i, &
                                      tabs, &
                                      rho, adz, &
                                      dz, dtn, dy, &
                                      nx, ny, ny_gl, &
                                      nstep, nstatis, icycle, &
                                      t, q, precsfc, prec_xy)
        !! Interface to the nn_convection parameterisation for the SAM model

        integer :: j, k

        real, dimension(:,:,:) :: tabs_i, q_i, tabs, t, q
        real, dimension(:, :) :: precsfc, prec_xy
        real, dimension(:) :: rho, adz
        != unit m :: dz
        real, intent(in) :: dz
            !! grid spacing in z direction for the lowest grid layer
        real, intent(in) :: dtn
        real, intent(in) :: dy
        integer, intent(in) :: nx, ny, ny_gl, nstep, nstatis, icycle

        real :: y_in(nx, ny)
            !! Distance of column from equator (proxy for insolation and sfc albedo)

        real, dimension(nx, ny, nrf) :: t_rad_rest_tend, t_delta_adv, q_delta_adv, &
                                        t_delta_auto, t_delta_sed, &
                                        q_delta_auto, q_delta_sed
            !! deltas/tendencies returned by the parameterisation:
            !! radiation rest tendency, advective, autoconversion, sedimentation

        real, dimension(nx, ny)      :: prec_sed
            !! Sedimenting precipitation at surface
       
        ! TODO: Does CAM require precipitation?
        ! Initialise precipitation to 0 if required and at start of cycle
        if(mod(nstep-1,nstatis).eq.0 .and. icycle.eq.1) then
            precsfc(:,:)=0.
        end if

        ! distance to the equator
        ! y is a proxy for insolation and surface albedo as both are only a function of |y| in SAM
        do j=1, ny
            ! TODO: Set y_in as appropriate for CAM
            y_in(:,j) = 0.0
        enddo

        !-----------------------------------------------------
        
        ! TODO: Formulate the input variables to the parameterisation as required.

        !-----------------------------------------------------
        
        ! Run the neural net parameterisation
        call nn_convection_flux(tabs_i(:,:,1:nrf), q_i(:,:,1:nrf), y_in, &
                                tabs(:,:,1:nrf), &
                                t(:,:,1:nrf), q(:,:,1:nrf), &
                                rho, adz, dz, dtn, &
                                t_rad_rest_tend, &
                                t_delta_adv, q_delta_adv, &
                                t_delta_auto, q_delta_auto, &
                                t_delta_sed, q_delta_sed, prec_sed)

        !-----------------------------------------------------
        
        ! TODO: Update CAM variables with NN outputs
        ! q_delta_adv(:,:,:)
        ! q_delta_auto(:,:,:)
        ! q_delta_sed(:,:,:)
        ! t_delta_adv(:,:,:)
        ! t_delta_auto(:,:,:)
        ! t_delta_sed(:,:,:)
        ! t_rad_rest_tend(:,:,:)*dtn

        ! TODO: Update precipitation if required


    end subroutine nn_convection_flux_CAM


    subroutine nn_convection_flux_CAM_finalize()
        !! Finalize the module deallocating arrays

        call nn_convection_flux_finalize()

    end subroutine nn_convection_flux_CAM_finalize
    

    !-----------------------------------------------------------------
    ! Private Subroutines
    
    ! TODO: Check if CAM needs to allocate MPI Subdomains in some way
    
    ! TODO: Check for redundant variables once refactored (icycle, nstatis, precip etc.)

end module nn_interface_SAM_mod
