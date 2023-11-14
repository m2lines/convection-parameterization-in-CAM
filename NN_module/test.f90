module tests
  !! Module containing individual tests for the CAM ML code

  !--------------------------------------------------------------------------
  ! Libraries to use
  use netcdf
  use nn_cf_net, only: relu, net_forward, nn_cf_net_init, nn_cf_net_finalize
    use nn_cf_net_torch, only: nn_cf_net_torch_init, nn_cf_net_torch_finalize  use nn_convection_flux, only:   nn_convection_flux_forward, nn_convection_flux_init, nn_convection_flux_finalize
    use nn_convection_flux_torch, only: nn_convection_flux_forward_torch, nn_convection_flux_init_torch, nn_convection_flux_finalize_torch
  use test_utils, only: assert_array_equal

  implicit none

    real(4), dimension(4) :: test_array = (/ -1.0, 0.0, 0.5, 1.0 /)
    integer :: nin, nout

    character(len=1024) :: nn_filename = "../NN_weights_YOG_convection.nc"

    ! Parameters from SAM that are used here
    ! From domain.f90
    integer, parameter :: YES3D = 1
        !! Domain dimensionality: 1 - 3D, 0 - 2D
    integer, parameter :: nx_gl = 1
        !! Number of grid points in X - Yani changed to 36 from 32
    integer, parameter :: ny_gl = 1
        !! Number of grid points in Y
    integer, parameter :: nz_gl = 30
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

    integer :: nrf

    != unit J :: t
    real t(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
        !! moist static energy
    real q(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
        !! total water
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
    real tabs_i(nx, ny, nzm)
        !! Temperature
    real q_i(nx, ny, nzm)
        !! Non-precipitating water mixing ratio
    real :: y_in(nx, ny)
        !! Distance of column from equator (proxy for insolation and sfc albedo)
    real adz(nzm)
        !! ratio of the grid spacing to dz for pressure levels
    != unit s :: dtn
    real dtn
        !! current dynamical timestep (can be smaller than dt)
    real dz
        !! current dynamical timestep (can be smaller than dt)


    real, dimension(nx,ny,30) :: t_rad_rest_tend, &
                                 t_delta_adv, q_delta_adv, &
                                 t_delta_auto, q_delta_auto, &
                                 t_delta_sed, q_delta_sed
    real, dimension(nx,ny) :: prec_sed

    write (*,*) "Test array", test_array

    write(*,*) "=== Testing routines from nn_cf_net ==="
    call test_nn_net_routines()
    write(*,*) "======================================="

    write(*,*) "=== Testing routines from nn_cf_net_torch ==="
    call test_nn_net_torch_routines()
    write(*,*) "============================================="

    write(*,*) "=== Testing routines from nn_convection_flux ==="
    call test_nn_flux_routines()
    write(*,*) "================================================"

    write(*,*) "=== Testing routines from nn_convection_flux_torch ==="
    call test_nn_flux_torch_routines()
    write(*,*) "======================================================"

    contains

        subroutine test_nn_net_routines()

            write (*,*) "Testing relu"
            call relu(test_array)
            write (*,*) "relu ok"

            write (*,*) "Testing neural network initialisation"
            call nn_cf_net_init(nn_filename, nin, nout, 30)
            write (*,*) "nn_cf_net_init ok"

            write (*,*) "Testing neural network finalisation"
            call nn_cf_net_finalize()
            write (*,*) "nn_cf_net_finalize ok"

        end subroutine test_nn_net_routines

        subroutine test_nn_net_torch_routines()

            write (*,*) "Testing FTorch neural network initialisation"
            call nn_cf_net_torch_init(nn_filename, nin, nout, 30)
            write (*,*) "nn_cf_net_torch_init ok"

            write (*,*) "Testing FTorch neural network finalisation"
            call nn_cf_net_torch_finalize()
            write (*,*) "nn_cf_net_torch_finalize ok"

        end subroutine test_nn_net_torch_routines

        subroutine test_nn_flux_routines()

            ! Test convection flux routines
            write (*,*) "Testing convection flux initialisation"
            call nn_convection_flux_init(nn_filename)
            write (*,*) "nn_convection_flux_init ok"

            t = 0.
            q = 0.4
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

            nrf = 30

            write (*,*) "Testing convection flux forward"
            call nn_convection_flux_forward(tabs_i, q_i, y_in, &
                                    tabs, &
                                    t, q, &
                                    rho, adz, dz, dtn, &
                                    t_rad_rest_tend, &
                                    t_delta_adv, q_delta_adv, &
                                    t_delta_auto, q_delta_auto, &
                                    t_delta_sed, q_delta_sed, prec_sed)
            write (*,*) "nn_convection_flux_forward ok"

            q(:, :, 1:nrf) = q(:, :, 1:nrf) + q_delta_adv(:,:,:) &
                                        + q_delta_auto(:,:,:) &
                                        + q_delta_sed(:,:,:)
            t(:, :, 1:nrf) = t(:, :, 1:nrf) + t_delta_adv(:,:,:) &
                                        + t_delta_auto(:,:,:) &
                                        + t_delta_sed(:,:,:) &
                                        + t_rad_rest_tend(:,:,:)*dtn

            write(*,*) "Example output: "
            write (*,*) t(-2:2, 0, 1)
            write (*,*) t(0, -2:2, 1)
            write (*,*) t(-2, -2, 1:nz_gl)
            write (*,*) t(-1, -1, 1:nz_gl)

            write (*,*) "Testing convection flux finalisation"
            call nn_convection_flux_finalize()
            write (*,*) "nn_convection_flux_finalise ok"

        end subroutine test_nn_flux_routines

        subroutine test_nn_flux_torch_routines()

            write (*,*) "Testing FTorch convection flux initialisation"
            call nn_convection_flux_init_torch(nn_filename)
            write (*,*) "nn_convection_flux_init ok"

            t = 0.
            q = 0.4
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

            nrf = 30

            write (*,*) "Testing convection flux forward"
            call nn_convection_flux_forward_torch(tabs_i, q_i, y_in, &
                                    tabs, &
                                    t, q, &
                                    rho, adz, dz, dtn, &
                                    t_rad_rest_tend, &
                                    t_delta_adv, q_delta_adv, &
                                    t_delta_auto, q_delta_auto, &
                                    t_delta_sed, q_delta_sed, prec_sed)
            write (*,*) "nn_convection_flux_forward_torch ok"

            q(:, :, 1:nrf) = q(:, :, 1:nrf) + q_delta_adv(:,:,:) &
                                        + q_delta_auto(:,:,:) &
                                        + q_delta_sed(:,:,:)
            t(:, :, 1:nrf) = t(:, :, 1:nrf) + t_delta_adv(:,:,:) &
                                        + t_delta_auto(:,:,:) &
                                        + t_delta_sed(:,:,:) &
                                        + t_rad_rest_tend(:,:,:)*dtn

            write(*,*) "Example output: "
            write (*,*) t(-2:2, 0, 1)
            write (*,*) t(0, -2:2, 1)
            write (*,*) t(-2, -2, 1:nz_gl)
            write (*,*) t(-1, -1, 1:nz_gl)

            write (*,*) "Testing FTorch convection flux finalization"
            call nn_convection_flux_finalize_torch()
            write (*,*) "nn_convection_flux_finalize ok"

        end subroutine test_nn_flux_torch_routines

end program run_tests
