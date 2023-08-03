module nn_interface_SAM_mod
    !! Interface to convection parameterisation for the SAM model
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
public  nn_convection_flux_SAM, &
        nn_convection_flux_SAM_init, nn_convection_flux_SAM_finalize


!---------------------------------------------------------------------
! local/private data

integer :: it, jt
    !! indices corresponding to the start of the grid domain for
    !! current MPI rank

! Copied from nn_convection_flux.f90
! Outputs from NN are supplied at lowest 30 half-model levels
integer, parameter :: nrf = 30
    !! number of vertical levels the NN parameterisation uses

!---------------------------------------------------------------------
! Functions and Subroutines

contains

    !-----------------------------------------------------------------
    ! Public Subroutines

    subroutine nn_convection_flux_SAM_init(nn_filename, rank, nx, ny, nsubdomains_x)
        !! Initialise the NN module

        character(len=1024), intent(in) :: nn_filename
            !! NetCDF filename from which to read model weights

        integer, intent(in) :: rank
            !! rank of the current subdomain task (default 0)

        integer, intent(in) :: nx, ny, nsubdomains_x
            !! information about the grid and subdomains
        
        ! Get indices for grid section on this MPI rank
        call task_rank_to_index(rank,it,jt,nx,ny,nsubdomains_x)

        ! Initialise the Neural Net from file and get info
        call nn_convection_flux_init(nn_filename)

    end subroutine nn_convection_flux_SAM_init


    subroutine nn_convection_flux_SAM(tabs_i, q_i, &
                                      tabs, &
                                      rho, adz, &
                                      dz, dtn, dy, &
                                      nx, ny, ny_gl, &
                                      nstep, nstatis, icycle, YES3D, &
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
        integer, intent(in) :: nx, ny, ny_gl, nstep, nstatis, icycle, YES3D

        real :: y_in(nx, ny)
            !! Distance of column from equator (proxy for insolation and sfc albedo)

        real, dimension(nx, ny, nrf) :: t_rad_rest_tend, t_delta_adv, q_delta_adv, &
                                        t_delta_auto, t_delta_sed, &
                                        q_delta_auto, q_delta_sed
        real, dimension(nx, ny)      :: prec_sed
            !! Sedimenting precipitation at surface
        
        ! Initialise precipitation to 0 if required and at start of cycle
        if(mod(nstep-1,nstatis).eq.0 .and. icycle.eq.1) then
            precsfc(:,:)=0.
        end if

        ! temperature - Note that this variable is called t_i in the main SAM code!!!
        ! tabs_i(i,j,1:input_ver_dim)

        ! non-precipitating water mixing ratio
        ! if (rf_uses_rh) then
        ! ! If using generalised relative humidity convert non-precip. water content to rel. hum
        !     do k=1,nzm
        !         omn = omegan(tabs(i,j,k))
        !         qsat(k) = omn * qsatw(tabs(i,j,k),pres(k)) + (1.-omn) * qsati(tabs(i,j,k),pres(k))
        !     end do
        !     ! non-precipitating water from mixing ratio
        !     real(q_i(i,j,1:input_ver_dim)/qsat(1:input_ver_dim),4)
        ! else
        !
        !     ! non-precipitating water
        !     q_i(i,j,1:input_ver_dim)
        ! end if

        ! distance to the equator
        ! y is a proxy for insolation and surface albedo as both are only a function of |y| in SAM
        do j=1, ny
            y_in(:,j) = real(abs(dy*(j+jt-(ny_gl+YES3D-1)/2-0.5)))
        enddo

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
        ! Update q and t with delta values
        ! advective, autoconversion (dt = -dq*(latent_heat/cp)),
        ! sedimentation (dt = -dq*(latent_heat/cp)),
        ! radiation rest tendency (multiply by dtn to get dt)
        q(:,:,1:nrf) = q(:,:,1:nrf) + q_delta_adv(:,:,:) &
                                    + q_delta_auto(:,:,:) &
                                    + q_delta_sed(:,:,:)
        t(:,:,1:nrf) = t(:,:,1:nrf) + t_delta_adv(:,:,:) &
                                    + t_delta_auto(:,:,:) &
                                    + t_delta_sed(:,:,:) &
                                    + t_rad_rest_tend(:,:,:)*dtn

        !-----------------------------------------------------
        ! Calculate surface precipitation

        ! Apply sedimentation at surface
        precsfc(:,:) = precsfc(:,:) + prec_sed(:,:) ! For statistics
        prec_xy(:,:) = prec_xy(:,:) + prec_sed(:,:) ! For 2D output
        
        ! Apply all autoconversion in the column and multiply by rho*adz for
        ! precip (rho*dq*adz) i.e. all precip. falls out on large model timescale
        do k=1, nrf
            precsfc(:,:) = precsfc(:,:) - q_delta_auto(:,:,k)*adz(k)*rho(k)
            prec_xy(:,:) = prec_xy(:,:) - q_delta_auto(:,:,k)*adz(k)*rho(k)
        end do

        ! As a final check enforce q must be >= 0.0 (should be redundant)
        where (q(:,:,1:nrf) .lt. 0)
          q = 0.0
        end where

    end subroutine nn_convection_flux_SAM


    subroutine nn_convection_flux_SAM_finalize()
        !! Finalize the module deallocating arrays

        call nn_convection_flux_finalize()

    end subroutine nn_convection_flux_SAM_finalize
    

    !-----------------------------------------------------------------
    ! Private Subroutines
    
    subroutine task_rank_to_index (rank, i, j, nx, ny, nsubdomains_x)
        !! returns the pair of  beginning indices for the subdomain on the
        !! global grid given the subdomain's rank.
        
        integer, intent(in) :: rank
            !! rank of MPI process
        integer, intent(in) :: nx, ny, nsubdomains_x
            !! information about the grid and subdomains
        integer :: i, j
            !! indices at which subdomain starts

        j = rank/nsubdomains_x
        i = rank - j*nsubdomains_x

        i = i * (nx)
        j = j * (ny)

    end subroutine task_rank_to_index

end module nn_interface_SAM_mod
