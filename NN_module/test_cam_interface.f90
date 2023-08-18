program run_cam_tests
    use nn_cf_net_mod, only: relu, nn_cf_net_init, nn_cf_net_finalize

    use nn_convection_flux_mod, only: nn_convection_flux, nn_convection_flux_init, nn_convection_flux_finalize

    use nn_interface_CAM_mod, only: nn_convection_flux_CAM, nn_convection_flux_CAM_init, nn_convection_flux_CAM_finalize
    implicit none

    call nn_convection_flux_CAM_init(trim("NN_weights_YOG_convection.nc"), trim("resources/SAM_sounding.nc"))
    call nn_convection_flux_CAM_finalize()



end program run_cam_tests
