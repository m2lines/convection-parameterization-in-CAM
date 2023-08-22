module nn_interface_CAM_mod
    !! Interface to convection parameterisation for the CAM model
    !! Reference: https://doi.org/10.1029/2020GL091363
    !! Also see YOG20: https://doi.org/10.1038/s41467-020-17142-3

    ! TODO: Check if CAM needs to allocate MPI Subdomains in some way
    
    ! TODO: Check for redundant variables once refactored (icycle, nstatis, precip etc.)

    ! TODO: Remove commented lines for variables not required from SAM Sounding
    ! TODO: Check precisions when defined in python, saved as netcdf, and read to Fortran

!---------------------------------------------------------------------
! Libraries to use
use netcdf
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

real, allocatable, dimension(:) :: pres, presi, rho, adz
!= unit m :: dz
real :: dz
    !! grid spacing in z direction for the lowest grid layer

!---------------------------------------------------------------------
! Functions and Subroutines

contains

    !-----------------------------------------------------------------
    ! Public Subroutines

    subroutine nn_convection_flux_CAM_init(nn_filename, sounding_filename)
        !! Initialise the NN module

        character(len=1024), intent(in) :: nn_filename
            !! NetCDF filename from which to read model weights
        character(len=1024), intent(in) :: sounding_filename
            !! NetCDF filename from which to read SAM sounding and grid data

        ! Initialise the Neural Net from file and get info
        call nn_convection_flux_init(nn_filename)

        ! Read in the SAM sounding and grid data required
        call sam_sounding_init(sounding_filename)

    end subroutine nn_convection_flux_CAM_init


    subroutine nn_convection_flux_CAM(tabs_i, q_i, &
                                      tabs, &
                                      dtn, dy, &
                                      nx, ny, ny_gl, &
                                      nstep, nstatis, icycle, &
                                      t, q, precsfc, prec_xy)
        !! Interface to the nn_convection parameterisation for the SAM model

        integer :: j, k

        real, dimension(:,:,:) :: tabs_i, q_i, tabs, t, q
        real, dimension(:, :) :: precsfc, prec_xy
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
        call sam_sounding_finalize()

    end subroutine nn_convection_flux_CAM_finalize
    

    !-----------------------------------------------------------------
    ! Private Subroutines
    
    subroutine sam_sounding_init(filename)
        !! Read various profiles in from SAM sounding file

        ! This will be the netCDF ID for the file and data variable.
        integer :: ncid
        integer :: z_dimid, dz_dimid
        integer :: nz
        integer :: z_varid, dz_varid, pres_varid, presi_varid, rho_varid, adz_varid

        character(len=1024), intent(in) :: filename
            !! NetCDF filename from which to read data


        !-------------allocate arrays and read data-------------------

        ! Open the file. NF90_NOWRITE tells netCDF we want read-only
        ! access
        ! Get the varid or dimid for each variable or dimension based
        ! on its name.

        call check( nf90_open(trim(filename),NF90_NOWRITE,ncid ))

        call check( nf90_inq_dimid(ncid, 'z', z_dimid))
        call check( nf90_inquire_dimension(ncid, z_dimid, len=nz))
        ! call check( nf90_inq_dimid(ncid, 'dz', dz_dimid))
        ! call check( nf90_inquire_dimension(ncid, dz_dimid, len=ndz))

        ! Note that nz of sounding may be longer than nrf
        allocate(pres(nz))
        allocate(presi(nz))
        allocate(rho(nz))
        allocate(adz(nz))

        ! Read data in from nc file
        ! call check( nf90_inq_varid(ncid, "z", z_varid))
        ! call check( nf90_get_var(ncid, z_varid, z))
        call check( nf90_inq_varid(ncid, "pressure", pres_varid))
        call check( nf90_get_var(ncid, pres_varid, pres))
        call check( nf90_inq_varid(ncid, "interface_pressure", presi_varid))
        call check( nf90_get_var(ncid, presi_varid, presi))
        call check( nf90_inq_varid(ncid, "rho", rho_varid))
        call check( nf90_get_var(ncid, rho_varid, rho))
        call check( nf90_inq_varid(ncid, "adz", adz_varid))
        call check( nf90_get_var(ncid, adz_varid, adz))
        call check( nf90_inq_varid(ncid, "dz", dz_varid))
        call check( nf90_get_var(ncid, dz_varid, dz))

        ! Close the nc file
        call check( nf90_close(ncid))

        write(*,*) 'Finished reading SAM sounding file.'

    end subroutine sam_sounding_init


    subroutine sam_sounding_finalize()
        !! Deallocate module variables read from sounding

        deallocate(pres, presi, rho, adz)

    end subroutine sam_sounding_finalize


    ! TODO Duplicated from nn_cf_net.f90 - consider placing 1 copy in shared location?
    subroutine check(err_status)
        !! Check error status after netcdf call and print message for
        !! error codes.

        integer, intent(in) :: err_status
            !! error status from nf90 function

        if(err_status /= nf90_noerr) then
             write(*, *) trim(nf90_strerror(err_status))
        end if

    end subroutine check


end module nn_interface_CAM_mod
