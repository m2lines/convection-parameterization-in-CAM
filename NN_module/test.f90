module tests
    !! Module containing individual tests for the CAM ML code
  
    !--------------------------------------------------------------------------
    ! Libraries to use
    use netcdf
    use nn_cf_net_mod, only: relu, nn_cf_net_init, nn_cf_net_forward, nn_cf_net_finalize
    use nn_cf_net_torch, only: nn_cf_net_torch_init, nn_cf_net_torch_forward, &
        nn_cf_net_torch_finalize
    use nn_convection_flux_mod, only: nn_convection_flux_forward, nn_convection_flux_init, &
        nn_convection_flux_finalize
    use test_utils, only: assert_array_equal
  
    implicit none

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
    integer, parameter :: nrf = 30

    != unit J :: t
    real t(nx, nzm)
    real t_torch(nx, nzm)
        !! moist static energy
    real q(nx, nzm)
    real q_torch(nx, nzm)
        !! total water
    ! fluxes at the top and bottom of the domain:
    real :: precsfc(nx) = 0.
        !! surface precip. rate
    !  Horizontally varying stuff (as a function of xy)
    real :: prec_xy(nx) = 0.
        !! surface precipitation rate
    ! reference vertical profiles:

    != unit (kg / m**3) :: rho
    real :: rho(nzm) = 1.2
        !! air density at pressure levels

    real :: tabs(nx, nzm) = 293.15
        !! absolute temperature

    ! Input Variables

    ! Fields from beginning of time step used as NN inputs
    real :: tabs_i(nx, nzm) = 293.15
        !! Temperature
    real :: q_i(nx, nzm) = 0.5
        !! Non-precipitating water mixing ratio
    real :: y_in(nx) = 0.0
        !! Distance of column from equator (proxy for insolation and sfc albedo)
    real :: adz(nzm) = 1.0
        !! ratio of the grid spacing to dz for pressure levels
    != unit s :: dtn
    real :: dtn = 2.0
        !! current dynamical timestep (can be smaller than dt)
    real :: dz = 100.0
        !! current dynamical timestep (can be smaller than dt)


    real, dimension(nx, nzm) :: t_rad_rest_tend, &
                                t_delta_adv, q_delta_adv, &
                                t_delta_auto, q_delta_auto, &
                                t_delta_sed, q_delta_sed
    real :: prec_sed(1)

    real(4), allocatable, dimension(:) :: features, features_torch
    real(4), allocatable, dimension(:) :: logits, logits_torch

    contains

        subroutine test_nn_net_routines()

            real(4), dimension(4) :: test_array = (/ -1.0, 0.0, 0.5, 1.0 /)

            write (*,*) "Testing relu"
            call relu(test_array)
            write (*,*) "relu ok"

            write (*,*) "Testing neural network initialisation"
            call nn_cf_net_init(nn_filename, nin, nout, 30)
            write (*,*) "nn_cf_net_init ok"

            write (*,*) "Testing neural network forward pass"
            allocate(features(nin))
            allocate(logits(nout))
            features = 1.
            logits = 0.
            call nn_cf_net_forward(features, logits)
            write (*,*) "nn_cf_net_forward ok"

            write (*,*) "Testing neural network finalisation"
            call nn_cf_net_finalize()
            write (*,*) "nn_cf_net_finalize ok"

        end subroutine test_nn_net_routines

        subroutine test_nn_net_torch_routines()

            write (*,*) "Testing FTorch neural network initialisation"
            call nn_cf_net_torch_init(nn_filename, nin, nout, 30)
            write (*,*) "nn_cf_net_torch_init ok"

            write (*,*) "Testing FTorch neural network forward pass"
            allocate(features_torch(nin))
            allocate(logits_torch(nout))
            features_torch = 1.
            logits_torch = 0.
            call nn_cf_net_torch_forward(features_torch, logits_torch)
            call assert_real_1d(logits, logits_torch, "Test NN output",  rtol_opt=1.0e-4)
            write (*,*) "nn_cf_net_torch_forward ok"

            write (*,*) "Testing FTorch neural network finalisation"
            call nn_cf_net_torch_finalize()
            write (*,*) "nn_cf_net_torch_finalize ok"

        end subroutine test_nn_net_torch_routines

        subroutine test_nn_flux_routines()

            ! Test convection flux routines
            write (*,*) "Testing convection flux initialisation"
            call nn_convection_flux_init(nn_cf_net_init, nn_filename)
            write (*,*) "nn_convection_flux_init ok"

            t = 0.
            q = 0.4

            write (*,*) "Testing convection flux forward"
            call nn_convection_flux_forward(nn_cf_net_forward, &
                                    tabs_i, q_i, y_in, &
                                    tabs, &
                                    t, q, &
                                    rho, adz, dz, dtn, &
                                    prec_sed)
            write (*,*) "nn_convection_flux_forward ok"

            ! q(:, 1:nrf) = q(:, 1:nrf) + q_delta_adv &
            !                             + q_delta_auto &
            !                             + q_delta_sed
            ! t(:, 1:nrf) = t(:, 1:nrf) + t_delta_adv &
            !                             + t_delta_auto &
            !                             + t_delta_sed &
            !                             + t_rad_rest_tend*dtn

            write (*,*) "Testing convection flux finalisation"
            call nn_cf_net_finalize()
            write (*,*) "nn_convection_flux_finalise ok"

        end subroutine test_nn_flux_routines

        subroutine test_nn_flux_torch_routines()

            write (*,*) "Testing FTorch convection flux initialisation"
            call nn_convection_flux_init(nn_cf_net_torch_init, nn_filename)
            write (*,*) "nn_convection_flux_init ok"

            t_torch = 0.
            q_torch = 0.4

            write (*,*) "Testing convection flux forward"
            call nn_convection_flux_forward(nn_cf_net_torch_forward, &
                                    tabs_i, q_i, y_in, &
                                    tabs, &
                                    t_torch, q_torch, &
                                    rho, adz, dz, dtn, prec_sed)

            ! q_torch(:, 1:nrf) = q_torch(:, 1:nrf) + q_delta_adv &
            !                             + q_delta_auto &
            !                             + q_delta_sed
            ! t_torch(:, 1:nrf) = t_torch(:, 1:nrf) + t_delta_adv &
            !                             + t_delta_auto &
            !                             + t_delta_sed &
            !                             + t_rad_rest_tend*dtn

            call assert_real_2d(t, t_torch, "Test t", rtol_opt=1.0e-4)
            call assert_real_2d(q, q_torch, "Test q", rtol_opt=1.0e-4)
            write (*,*) "nn_convection_flux_forward_torch ok"

            write (*,*) "Testing FTorch convection flux finalization"
            call nn_cf_net_torch_finalize()
            write (*,*) "nn_convection_flux_finalize ok"

        end subroutine test_nn_flux_torch_routines

        subroutine assert_real_1d(a, b, test_name, rtol_opt)

            character(len=*), intent(in) :: test_name
            real(4), dimension(:), intent(in) :: a, b
            real(4), intent(in), optional :: rtol_opt
            real(4) :: relative_error, rtol

            if (.not. present(rtol_opt)) then
                rtol = 1.0e-5
            else
                rtol = rtol_opt
            end if

            relative_error = maxval(abs(a/b - 1.0))
            call print_assert_real(test_name, (rtol > relative_error), relative_error)

        end subroutine assert_real_1d

        subroutine assert_real_2d(a, b, test_name, rtol_opt)

            character(len=*), intent(in) :: test_name
            real(4), dimension(:, :), intent(in) :: a, b
            real(4), intent(in), optional :: rtol_opt
            real(4) :: relative_error, rtol

            if (.not. present(rtol_opt)) then
                rtol = 1.0e-5
            else
                rtol = rtol_opt
            end if

            relative_error = maxval(abs(a/b - 1.0))
            call print_assert_real(test_name, (rtol > relative_error), relative_error)

        end subroutine assert_real_2d

        subroutine print_assert_real(test_name, is_close, relative_error)

            character(len=*), intent(in) :: test_name
            logical, intent(in) :: is_close
            real(4), intent(in) :: relative_error
            character(len=15) :: pass, fail

            fail = char(27)//'[31m'//'FAILED'//char(27)//'[0m'
            pass = char(27)//'[32m'//'PASSED'//char(27)//'[0m'

            if (is_close) then
            write(*, '(A, " :: [", A, "] maximum relative error = ", E11.4)') pass, trim(test_name), relative_error
            else
            write(*, '(A, " :: [", A, "] maximum relative error = ", E11.4)') fail, trim(test_name), relative_error
            end if

        end subroutine print_assert_real

end module tests


program run_tests

    use tests
    implicit none

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

    deallocate(features)
    deallocate(logits)
    deallocate(features_torch)
    deallocate(logits_torch)

end program run_tests