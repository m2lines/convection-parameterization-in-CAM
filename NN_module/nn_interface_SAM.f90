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


! Parameters from SAM that are used here
! From domain.f90
integer, parameter :: YES3D = 1
    !! Domain dimensionality: 1 - 3D, 0 - 2D
integer, parameter :: nx_gl = 72
    !! Number of grid points in X - Yani changed to 36 from 32
integer, parameter :: ny_gl = 180
    !! Number of grid points in Y
integer, parameter :: nz_gl = 48
    !! Number of pressure (scalar) levels
integer, parameter :: nsubdomains_x  = 3
    !! No of subdomains in x
integer, parameter :: nsubdomains_y  = 12
    !! No of subdomains in y

! From grid.f90
integer, parameter :: nx = nx_gl/nsubdomains_x
    !! Number of x points in a subdomain
integer, parameter :: ny = ny_gl/nsubdomains_y
    !! Number of y points in a subdomain
integer, parameter :: nz = nz_gl+1
    !! Number of z points in a subdomain
! Store useful variations on these values
integer, parameter :: nzm = nz-1
integer, parameter :: nxp3 = nx + 3
integer, parameter :: nyp3 = ny + 3 * YES3D
integer, parameter :: dimx1_s = -2
integer, parameter :: dimx2_s = nxp3
integer, parameter :: dimy1_s = 1-3*YES3D
integer, parameter :: dimy2_s = nyp3

real dy
    !! grid spacing in y direction
real dz
    !! grid spacing in z direction for the lowest grid layer, input to parameterisation


! Subcycling variables
integer nstep
    !! current number of performed time steps
integer icycle
    !! current subcycle
integer nstatis
    !! the interval between substeps to compute statistics

! Multitasking stuff
integer rank
    !! rank of the current subdomain task (default 0)
! logical masterproc
!     !! logical variable .true. if rank .eq. 0





! From vars.f90
! prognostic variables:

!= unit J :: t
real t(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    !! moist static energy

real q(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    !! total water
! diagnostic variables:
real qn(nx, ny, nzm)
    !! cloud water+cloud ice
! fluxes at the top and bottom of the domain:
real precsfc(nx,ny)
    !! surface precip. rate
!  Horizontally varying stuff (as a function of xy)
real prec_xy(nx,ny)
    !! surface precipitation rate
! reference vertical profiles:

!= unit (kg / m**3) :: rho
real rho(nzm)
    !! air density at pressure levels

real :: tabs(nx, ny, nzm)
    !! absolute temperature

! Input Variables

! Fields from beginning of time step used as NN inputs
real t_i(nx,ny,nzm)
    !! Temperature
real q_i(nx,ny,nzm)
    !! Non-precipitating water mixing ratio
real qp_i (nx,ny,nzm)
    !! Precipitating water mixing ratio
real adz(nzm)
    !! ratio of the grid spacing to dz for pressure levels
!= unit s :: dtn
real dtn
    !! current dynamical timestep (can be smaller than dt)
! ONLY Needed if using rf_uses_rh
! != unit mb :: pres
! real pres(nzm)
!     !! pressure,mb at scalar levels

!---------------------------------------------------------------------
! Functions and Subroutines

contains

    subroutine nn_convection_flux_SAM_init(nn_filename)
        !! Initialise the NN module

        character(len=1024), intent(in) :: nn_filename
            !! NetCDF filename from which to read model weights

        ! Get indices for grid section on this MPI rank
        call task_rank_to_index(rank,it,jt)

        ! Initialise the Neural Net from file and get info
        call nn_convection_flux_init(nn_filename)

    end subroutine nn_convection_flux_SAM_init

    subroutine nn_convection_flux_SAM()
        !! Interface to the nn_convection parameterisation for the SAM model

        integer :: j, k

        real :: y_in(ny)
            !! Distance of column from equator (proxy for insolation and sfc albedo)

        ! Initialise precipitation to 0 if required and at start of cycle
        if(mod(nstep-1,nstatis).eq.0 .and. icycle.eq.1) then
            precsfc(:,:)=0.
        end if

        ! temperature
        ! t_i(i,j,1:input_ver_dim)

        ! non-precipitating water mixing ratio
        ! if (rf_uses_rh) then
        ! ! If using generalised relative humidity convert non-precip. water content to rel. hum
        !     do k=1,nzm
        !         omn = omegan(tabs(i,j,k))
        !         qsat(k) = omn * qsatw(tabs(i,j,k),pres(k)) + (1.-omn) * qsati(tabs(i,j,k),pres(k))
        !     end do
        !     features(dim_counter+1:dim_counter+input_ver_dim) = real(q_i(i,j,1:input_ver_dim)/qsat(1:input_ver_dim),4)
        !     dim_counter = dim_counter + input_ver_dim
        ! else

        ! non-precipitating water
        ! q_i(i,j,1:input_ver_dim)

        ! precipitating water (rain+snow) mixing ratio content
        ! qp_i(i,j,1:input_ver_dim)

        ! distance to the equator
        ! y is a proxy for insolation and surface albedo as both are only a function of |y| in SAM
        do j=1, ny
            y_in(j) = real(abs(dy*(j+jt-(ny_gl+YES3D-1)/2-0.5)))
        enddo


        call nn_convection_flux(t_i, q_i, qp_i, y_in, &
                                rho, adz, tabs, dz, dtn, &
                                t, q, qn, precsfc, prec_xy)


    end subroutine nn_convection_flux_SAM

    subroutine nn_convection_flux_SAM_finalize()
        !! Finalize the module deallocating arrays

        call nn_convection_flux_finalize()

    end subroutine nn_convection_flux_SAM_finalize
    
    !###################################################################

    subroutine task_rank_to_index (rank,i,j)

      ! returns the pair of  beginning indices for the subdomain on the
      ! global grid given the subdomain's rank.
      integer, intent(in) :: rank
          !! rank of MPI process
      integer :: i, j
          !! indices at which subdomain starts

      j = rank/nsubdomains_x
      i = rank - j*nsubdomains_x

      i = i * (nx_gl/nsubdomains_x)
      j = j * (ny_gl/nsubdomains_y)

    end subroutine task_rank_to_index


end module nn_interface_SAM_mod
