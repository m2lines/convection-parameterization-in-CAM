program run_tests
    use nn_convection_flux_mod, only: relu

    implicit none

    real(4), dimension(4) :: test_array = (/ -1.0, 0.0, 0.5, 1.0 /)

    call relu(test_array)

    write (*,*) test_array


end program run_tests
