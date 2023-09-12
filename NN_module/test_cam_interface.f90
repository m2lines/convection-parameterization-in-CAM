program run_cam_tests
    use nn_cf_net_mod, only: relu, nn_cf_net_init, nn_cf_net_finalize

    use nn_convection_flux_mod, only: nn_convection_flux, nn_convection_flux_init, nn_convection_flux_finalize

    use nn_interface_CAM_mod, only: nn_convection_flux_CAM, nn_convection_flux_CAM_init, nn_convection_flux_CAM_finalize, &
    interp_to_sam

    implicit none

    real, dimension(4, 3) :: p_cam
    real, dimension(4, 3) :: var_cam
    real, dimension(4) :: ps_cam
    real, dimension(4, 30) :: var_sam

    ! Test interpolation by interpolating pressure to pressure => should match pres from SAM
    p_cam = reshape((/ 1000.0, 1000.0, 1000.0, 1000.0, 500.0, 500.0, 500.0, 500.0, 10.0, 10.0, 10.0, 10.0 /), (/ 4, 3 /))
    ps_cam = (/ 1000.0, 1000.0, 1000.0, 1000.0/)
    var_cam = reshape((/ 1000.0, 1000.0, 1000.0, 1000.0, 500.0, 500.0, 500.0, 500.0, 10.0, 10.0, 10.0, 10.0 /), (/ 4, 3 /))
    var_sam = -1.0

    call nn_convection_flux_CAM_init(trim("NN_weights_YOG_convection.nc"), trim("resources/SAM_sounding.nc"))
    call interp_to_sam(p_cam, ps_cam, var_cam, var_sam)
    call nn_convection_flux_CAM_finalize()
    do k=1,30
        write(*,*) var_sam(2, k)
    end do




end program run_cam_tests
