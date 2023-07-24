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
                                      dz, dtn, dy, ny, ny_gl, &
                                      nstep, nstatis, icycle, YES3D, &
                                      t, q, precsfc, prec_xy)
        !! Interface to the nn_convection parameterisation for the SAM model

        integer :: j, k

        real, dimension(:,:,:) :: tabs_i, q_i, tabs, t, q
        real, dimension(:, :) :: precsfc, prec_xy
        real, dimension(:) :: rho, adz
        real, intent(in) :: dz, dtn
        real, intent(in) :: dy
        integer, intent(in) :: ny, ny_gl, nstep, nstatis, icycle, YES3D

        real :: y_in(ny)
            !! Distance of column from equator (proxy for insolation and sfc albedo)

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
            y_in(j) = real(abs(dy*(j+jt-(ny_gl+YES3D-1)/2-0.5)))
        enddo

        ! Run the parameterization
        call nn_convection_flux(tabs_i, q_i, y_in, &
                                tabs, &
                                rho, adz, &
                                dz, dtn, &
                                t, q, precsfc, prec_xy)

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
