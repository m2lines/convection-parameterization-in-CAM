program run_tests
    use nn_convection_flux_mod, only: relu, nn_convection_flux_init, nn_convection_flux_finalize

    implicit none

    real(4), dimension(4) :: test_array = (/ -1.0, 0.0, 0.5, 1.0 /)
    integer :: nin, nout

    call relu(test_array)

    write (*,*) test_array

    call nn_convection_flux_init("../torch_nets/qobsTTFFFFFTF30FFTFTF30TTFTFTFFF80FFTFTTF2699FFFF_X01_no_qp_no_adv_surf_F_Tin_qin_disteq_O_Trad_rest_Tadv_qadv_qout_qsed_RESCALED_7epochs_no_drop_REAL_NN_layers5in61out148_BN_F_te70.nc", nin, nout)
    call nn_convection_flux_finalize()

end program run_tests
