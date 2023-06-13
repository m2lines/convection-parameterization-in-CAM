program run_tests
    use nn_cf_net_mod, only: relu, nn_cf_net_init, nn_cf_net_finalize

    use nn_convection_flux_mod, only: nn_convection_flux, nn_convection_flux_init, nn_convection_flux_finalize

    implicit none

    real(4), dimension(4) :: test_array = (/ -1.0, 0.0, 0.5, 1.0 /)
    integer :: nin, nout




    ! Parameters from SAM that are used here
    ! From domain.f90
    integer, parameter :: YES3D = 1
        !! Domain dimensionality: 1 - 3D, 0 - 2D
    integer, parameter :: nx_gl = 1
        !! Number of grid points in X - Yani changed to 36 from 32
    integer, parameter :: ny_gl = 1
        !! Number of grid points in Y
    integer, parameter :: nz_gl = 48
        !! Number of pressure (scalar) levels
    integer, parameter :: nsubdomains_x  = 1
        !! No of subdomains in x
    integer, parameter :: nsubdomains_y  = 1
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
    real tabs_i(nx,ny,nzm)
        !! Temperature
    real q_i(nx,ny,nzm)
        !! Non-precipitating water mixing ratio
    real :: y_in(ny)
        !! Distance of column from equator (proxy for insolation and sfc albedo)
    real adz(nzm)
        !! ratio of the grid spacing to dz for pressure levels
    != unit s :: dtn
    real dtn
        !! current dynamical timestep (can be smaller than dt)
    real dz
        !! current dynamical timestep (can be smaller than dt)

    ! Test NN routines

    call relu(test_array)

    write (*,*) test_array

    call nn_cf_net_init("./NN_weights_YOG_convection.nc", nin, nout, 30)
    call nn_cf_net_finalize()


    ! Test convection flux routines


    t = 0.
    q = 0.
    qn = 0.
    precsfc = 0.
    prec_xy = 0.
    rho = 1.
    tabs = 1.
    tabs_i = 287.15
    q_i = 0.2
    adz = 1.
    y_in = 1.
    dz = 1.
    dtn = 1.


    call nn_convection_flux_init("./NN_weights_YOG_convection.nc")

    call nn_convection_flux(tabs_i, q_i, y_in, &
                            tabs, &
                            rho, adz, &
                            dz, dtn, &
                            t, q, qn, precsfc, prec_xy)

    call nn_convection_flux_finalize()

    
    write (*,*) t(-2:2, 0, 1)
    write (*,*) t(0, -2:2, 1)
    write (*,*) t(-2, -2, 1:48)
    ! write (*,*) t(-1, -1, 1:48)

end program run_tests
