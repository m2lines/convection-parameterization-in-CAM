program run_cam_tests
    use nn_cf_net, only: relu, nn_cf_net_init, nn_cf_net_finalize

    use nn_convection_flux, only: nn_convection_flux_forward, nn_convection_flux_init, nn_convection_flux_finalize

    use nn_interface_CAM, only: nn_convection_flux_CAM, nn_convection_flux_CAM_init, nn_convection_flux_CAM_finalize, &
    interp_to_sam, interp_to_cam, fetch_sam_data

    implicit none

    integer :: i, k

    real, dimension(4, 3) :: p_cam, p_int_cam
    real, dimension(4, 3) :: var_cam
    real, dimension(4) :: var_cam_surface
    real, dimension(4) :: ps_cam
    real, dimension(4, 30) :: var_sam

    real, dimension(48) :: pres_sam, presi_sam

    real, dimension(4, 31) :: p_cam2, p_int_cam2, var_cam2
    real, dimension(4, 93) :: p_cam3, p_int_cam3, var_cam3
    real, dimension(4, 11) :: p_cam4, p_int_cam4, var_cam4

    ! Initialise the NN module
    call nn_convection_flux_CAM_init(trim("NN_weights_YOG_convection.nc"), trim("resources/SAM_sounding.nc"))
    
    ! Fetch SAM grid data
    call fetch_sam_data(pres_sam, presi_sam)


    ! Check interpolation to SAM grid by interpolating variable of 1.0
    p_cam = reshape((/ 1000.0, 1000.0, 1000.0, 1000.0, 500.0, 500.0, 500.0, 500.0, 10.0, 10.0, 10.0, 10.0 /), (/ 4, 3 /))
    ps_cam = (/ 1000.0, 1000.0, 1000.0, 1000.0/) / 0.9

    var_cam = 1.0
    var_cam_surface = 1.0
    write(*,*) "Test interpolation with array of ones"
    call interp_to_sam(p_cam, ps_cam, var_cam, var_sam, var_cam_surface)
    do k=1,30
        write(*,*) k, pres_sam(k), var_sam(2, k)
    end do


    ! Check interpolation to SAM grid by interpolating pressure to pressure
    ! Set top of CAM to 1.0d-4
    ! => should match pres from SAM
    p_cam(:, 1) = pres_sam(5)
    p_cam(:, 2) = pres_sam(10)
    p_cam(:, 3) = 1.0d-4
    ps_cam = presi_sam(1)
    var_cam(:, 1) = pres_sam(5)
    var_cam(:, 2) = pres_sam(10)
    var_cam(:, 3) = 1.0d-4
    var_cam_surface = presi_sam(1)

    call interp_to_sam(p_cam, ps_cam, var_cam, var_sam, var_cam_surface)
    write(*,*) "Test interpolation mapping pressure to pressure"
    do k=1,30
        write(*,*) k, pres_sam(k), var_sam(2, k)
    end do


    ! Test interpolation from SAM grid to CAM via conservative regridding
    ! Match SAM and CAM grids
    write(*,*) "Test interpolation to CAM with column of density 1.0 and equal grids."
    do i=1,4
        p_cam2(i, 1:30) = pres_sam(1:30)
        p_int_cam2(i, 1:30) = presi_sam(1:30)
        ! Set SAM variable equal to cell size
        var_sam(i, :) = (presi_sam(1:30) - presi_sam(2:31))
        var_sam(i, :) = 1.0
    enddo
    p_cam2(:,31) = 0.0
    p_int_cam2(:,31) = 1.0
    do k=1,29
        write(*,*) k, pres_sam(k), var_sam(2, k) / (presi_sam(k)-presi_sam(k+1))
    end do
    write(*,*) 30, pres_sam(30), var_sam(2, 30) / (2.0*(presi_sam(30)-pres_sam(30)))

    call interp_to_cam(p_cam2, p_int_cam2, var_sam, var_cam2)
    do k=1,30
        write(*,*) k, p_cam2(2, k), var_sam(2, k), var_cam2(2, k), var_cam2(2, k) / (p_int_cam2(2, k)-p_int_cam2(2, k+1))
    end do
    write(*,*) 31, p_cam2(2, 31), 1.0, var_cam2(2, 31), var_cam2(2, 31) / (2.0*(p_cam2(2, 31) - p_int_cam2(2, 31)))
    do i=1,4
        write(*,*) sum(var_sam(i,:)), sum(var_cam2(i,:))
    enddo

    ! CAM grid finer than SAM grid
    do k=0, 29
    do i=1,4
        p_cam3(i, 3*k+1) = pres_sam(k+1)
        p_cam3(i, 3*k+2) = pres_sam(k+1) + (pres_sam(k+2)-pres_sam(k+1))/3
        p_cam3(i, 3*k+3) = pres_sam(k+1) + (pres_sam(k+2)-pres_sam(k+1))/2
        p_int_cam3(i, 3*k+1) = presi_sam(k+1)
        p_int_cam3(i, 3*k+2) = presi_sam(k+1) + (presi_sam(k+2)-presi_sam(k+1))/3
        p_int_cam3(i, 3*k+3) = presi_sam(k+1) + (presi_sam(k+2)-presi_sam(k+1))/2
        ! Set SAM variable equal to cell size
        var_sam(i, :) = (presi_sam(1:30) - presi_sam(2:31))
    enddo
    enddo
    p_cam3(:,91) = 20.0
    p_cam3(:,92) = 10.0
    p_cam3(:,93) = 0.0
    p_int_cam3(:,91) = 25.0
    p_int_cam3(:,92) = 12.0
    p_int_cam3(:,93) = 1.0

    call interp_to_cam(p_cam3, p_int_cam3, var_sam, var_cam3)
    write(*,*) "Test interpolation to CAM with column of density 1.0 and finer grid."
    do k=1,29
        write(*,*) k, pres_sam(k), var_sam(2, k), var_sam(2, k) / (presi_sam(k)-presi_sam(k+1))
    end do
    write(*,*) 30, pres_sam(30), var_sam(2, 30), var_sam(2, 30) / (2.0*(presi_sam(30)-pres_sam(30)))
    do k=1,92
        write(*,*) k, p_cam3(2, k), var_cam3(2, k), var_cam3(2, k) / (p_int_cam3(2, k)-p_int_cam3(2, k+1))
    end do
    write(*,*) 93, p_cam3(2, 93), var_cam3(2, 93), var_cam3(2, 93) / (2.0*(p_cam3(2, 93) - p_int_cam3(2, 93)))
    do i=1,4
        write(*,*) sum(var_sam(i,:)), sum(var_cam3(i,:))
    enddo

    ! CAM grid coarser than SAM grid
    do k=0, 10
    do i=1,4
        p_cam4(i, k+1) = pres_sam(3*k+1)
        p_int_cam4(i, k+1) = presi_sam(3*k+1)
        ! Set SAM variable equal to cell size
        var_sam(i, :) = (presi_sam(1:30) - presi_sam(2:31))
    enddo
    enddo
    do i=1,4
        p_cam4(i, 11) = pres_sam(30)
        p_int_cam4(i, 11) = presi_sam(30)
    enddo
    do k=1,11
        write(*,*) k, p_cam4(2, k), p_int_cam4(2, k), var_cam4(2, k)
    end do

    call interp_to_cam(p_cam4, p_int_cam4, var_sam, var_cam4)
    write(*,*) "Test interpolation to CAM with column of density 1.0 and coarser grid."
    do k=1,10
        write(*,*) k, p_cam4(2, k), var_cam4(2, k), var_cam4(2, k) / (p_int_cam4(2, k)-p_int_cam4(2, k+1))
        write(*,*) k, p_cam4(1, k), var_cam4(1, k), var_cam4(1, k) / (p_int_cam4(1, k)-p_int_cam4(1, k+1))
    end do
    write(*,*) 11, p_cam4(2, 11), var_cam4(2, 11), var_cam4(2, 11) / (2.0*(p_int_cam4(2, 11) - p_cam4(2, 11)))
    do i=1,4
        write(*,*) sum(var_sam(i,:)), sum(var_cam4(i,:))
    enddo





    ! Test interpolation from SAM grid to CAM via conservative regridding
    p_cam = reshape((/ 1000.0, 1000.0, 1000.0, 1000.0, 500.0, 500.0, 500.0, 500.0, 0.0, 0.0, 0.0, 0.0 /), (/ 4, 3 /))
    p_cam(:, 3) = pres_sam(30)
    p_int_cam(:, 1) = 1005.0
    p_int_cam(:, 2) = 750.0
    p_int_cam(:, 3) = 100.0
    do i=1,4
        ! Set SAM variable equal to cell size
        var_sam(i, :) = (presi_sam(1:30) - presi_sam(2:31))
    enddo
    do k=1,30
!        write(*,*) k, pres_sam(k), var_sam(2, k) / (presi_sam(k)-presi_sam(k+1))
    end do

    call interp_to_cam(p_cam, p_int_cam, var_sam, var_cam)
!    write(*,*) "Test interpolation to CAM"
    do k=1,3
!        write(*,*) k, p_cam(2, k), var_cam(2, k), var_cam(2, k) / (p_int_cam(2, k)-p_int_cam(2, k+1))
    end do


    ! Clean up
    call nn_convection_flux_CAM_finalize()

end program run_cam_tests
