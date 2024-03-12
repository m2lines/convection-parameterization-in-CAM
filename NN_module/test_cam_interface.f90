module cam_tests
  !! Module containing individual tests for the CAM ML code

  !--------------------------------------------------------------------------
  ! Libraries to use
  use netcdf
  use nn_cf_net_mod, only: relu, nn_cf_net_init, nn_cf_net_finalize, net_forward
  use nn_convection_flux_mod, only:   nn_convection_flux, nn_convection_flux_init, nn_convection_flux_finalize
  use nn_interface_CAM, only: nn_convection_flux_CAM, nn_convection_flux_CAM_init, nn_convection_flux_CAM_finalize, &
  interp_to_sam, interp_to_cam, fetch_sam_data, SAM_var_conversion, CAM_var_conversion
  use test_utils, only: assert_array_equal

  implicit none

  character(len=15) :: pass = char(27)//'[32m'//'PASSED'//char(27)//'[0m'
  character(len=15) :: fail = char(27)//'[31m'//'FAILED'//char(27)//'[0m'
  integer, parameter :: nrf = 30
  integer, parameter :: n_nn_out = 148
  real(4), dimension(n_nn_out) :: nn_out_ones

  contains

    subroutine test_interp_to_sam_match(test_name)
      !! Check interpolation to SAM grid by interpolating variable of value 1.0
      !! Define a CAM grid from 1111.0 to 10.0

      character(len=*), intent(in) :: test_name
      
      integer :: i

      real, dimension(4, 48) :: p_cam, p_int_cam
      real, dimension(4, 48) :: var_cam
      real, dimension(4) :: var_cam_surface
      real, dimension(4) :: ps_cam
      real, dimension(4, 48) :: var_sam, var_sam_exp

      real, dimension(48) :: pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam
          !! Data from the SAM soundings used in tests

      ! Fetch SAM grid data
      call fetch_sam_data(pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam)

      do i=1,4
          p_cam(i, 1:48) = pres_sam(1:48)
          p_int_cam(i, 1:47) = presi_sam(2:48)
          ! Set SAM variable equal to cell size (density 1.0)
          var_cam(i, :) = pres_sam(1:48)
      enddo
      ! Set interface of top of CAM grid
      p_int_cam(:,48) = p_cam(:, 48) + (p_int_cam(:, 47)-p_cam(:, 48))

      ps_cam(:) = presi_sam(1)
      var_cam_surface(:) = presi_sam(1)

      call interp_to_sam(p_cam, ps_cam, var_cam, var_sam, var_cam_surface)

      ! Set anything above 30 elems to zero as per interpolation routine
      var_sam_exp = var_cam
      var_sam_exp(:, 31:48) = 0.0
      call assert_array_equal(var_sam, var_sam_exp, test_name)

    end subroutine test_interp_to_sam_match

    subroutine test_interp_to_sam_one(test_name)
      !! Check interpolation to SAM grid by interpolating variable of value 1.0
      !! Define a CAM grid from 1111.0 to 10.0

      character(len=*), intent(in) :: test_name

      real, dimension(4, 3) :: p_cam
      real, dimension(4, 3) :: var_cam
      real, dimension(4) :: var_cam_surface
      real, dimension(4) :: ps_cam
      real, dimension(4, 30) :: var_sam, var_sam_exp

      p_cam = reshape((/ 1000.0, 1000.0, 1000.0, 1000.0, 500.0, 500.0, 500.0, 500.0, 10.0, 10.0, 10.0, 10.0 /), (/ 4, 3 /))
      ps_cam = (/ 1111.0, 1111.0, 1111.0, 1111.0/)

      var_cam = 1.0
      var_cam_surface = 1.0
      var_sam_exp = 1.0

      call interp_to_sam(p_cam, ps_cam, var_cam, var_sam, var_cam_surface)
      
      call assert_array_equal(var_sam, var_sam_exp, test_name)

    end subroutine test_interp_to_sam_one

    subroutine test_interp_to_sam_pres(test_name)
      !! Check interpolation to SAM grid by interpolating pressure to pressure
      !! Set top of CAM to 1.0d-4
      !! => should match pres from SAM

      character(len=*), intent(in) :: test_name

      integer :: i

      real, dimension(4, 3) :: p_cam
      real, dimension(4, 3) :: var_cam
      real, dimension(4) :: var_cam_surface
      real, dimension(4) :: ps_cam
      real, dimension(4, 30) :: var_sam, var_sam_exp

      real, dimension(48) :: pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam
          !! Data from the SAM soundings used in tests

      ! Fetch SAM grid data
      call fetch_sam_data(pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam)

      p_cam(:, 1) = pres_sam(5)
      p_cam(:, 2) = pres_sam(10)
      p_cam(:, 3) = 1.0d-4
      ps_cam = presi_sam(1)
      var_cam(:, 1) = pres_sam(5)
      var_cam(:, 2) = pres_sam(10)
      var_cam(:, 3) = 1.0d-4
      var_cam_surface = presi_sam(1)

      do i = 1,4
        var_sam_exp(i, :) = pres_sam(1:30)
      end do

      call interp_to_sam(p_cam, ps_cam, var_cam, var_sam, var_cam_surface)

      call assert_array_equal(var_sam_exp, var_sam, test_name)

    end subroutine test_interp_to_sam_pres

    subroutine test_interp_to_cam_match(test_name)
      !! Check conservative regridding to CAM grid by interpolating constant density
      !! => should conserve density
      !! For case CAM and SAM grids match and CAM grid inside SAM grid

      character(len=*), intent(in) :: test_name

      integer :: i

      real, dimension(4, 30) :: p_cam, var_cam, rho_cam, rho_cam_exp
      real, dimension(4, 31) :: p_int_cam
      real, dimension(4) :: ps_cam, var_cam_surface
      real, dimension(4, 30) :: var_sam

      real, dimension(48) :: pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam
          !! Data from the SAM soundings used in tests

      ! Fetch SAM grid data
      call fetch_sam_data(pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam)

      do i=1,4
          p_cam(i, 1:30) = pres_sam(1:30)
          p_int_cam(i, 1:31) = presi_sam(1:31)
          ! Set SAM variable equal to cell size ("density" 1.0)
          var_sam(i, :) = (presi_sam(1:30) - presi_sam(2:31))
          ! var_sam(i, :) = 1.0
      enddo

      call interp_to_cam(p_cam, p_int_cam, var_sam, var_cam)

      do i=1,30
          rho_cam(:, i) = var_cam(:, i) / (p_int_cam(:, i)-p_int_cam(:, i+1))
      end do
      
      rho_cam_exp = 1.0

      call assert_array_equal(rho_cam_exp, rho_cam, test_name)

      ! Could do an additional check on the sums of the variables to check
      ! variable has been conserved
      ! do i=1,4
      !     write(*,*) sum(var_sam(i,:)), sum(var_cam(i,:))
      ! enddo

    end subroutine test_interp_to_cam_match

    subroutine test_interp_to_cam_fine(test_name)
      !! Check conservative regridding to CAM grid by interpolating constant density
      !! => should conserve density
      !! For case CAM and SAM grids match and CAM grid inside SAM grid

      character(len=*), intent(in) :: test_name

      integer :: i, j

      real, dimension(4, 90) :: p_cam, var_cam, rho_cam, rho_cam_exp
      real, dimension(4, 91) :: p_int_cam
      real, dimension(4) :: var_cam_surface
      real, dimension(4) :: ps_cam, sum_sam, sum_cam
      real, dimension(4, 30) :: var_sam

      real, dimension(48) :: pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam
          !! Data from the SAM soundings used in tests

      ! Fetch SAM grid data
      call fetch_sam_data(pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam)

      ! CAM grid finer than SAM grid
      do j=0, 29
      do i=1,4
          p_int_cam(i, 3*j+1) = presi_sam(j+1)
          p_int_cam(i, 3*j+2) = presi_sam(j+1) + (presi_sam(j+2)-presi_sam(j+1))/3
          p_int_cam(i, 3*j+3) = presi_sam(j+1) + 2*(presi_sam(j+2)-presi_sam(j+1))/3
      enddo
      enddo
      p_int_cam(:, 91) = presi_sam(31)

      do j=1, 90
          p_cam(:, j) = (p_int_cam(:, j+1)+p_int_cam(:, j)) / 2.0
      end do

      do j=1, 4
          ! Set SAM variable equal to cell size
          var_sam(j, :) = (presi_sam(1:30) - presi_sam(2:31))
      end do

      call interp_to_cam(p_cam, p_int_cam, var_sam, var_cam)

      do i=1,92
          rho_cam(:, i) = var_cam(:, i) / (p_int_cam(:, i)-p_int_cam(:, i+1))
      end do

      rho_cam_exp = 1.0

!       ! Compare individual cells of regridded data
!       do i = 1,30
!         write(*,*) i, var_sam(1,i), presi_sam(i)
!         ! write(*,*) i, p_int_cam(1,i)
!         ! write(*,*) i, p_cam(1,i)
!       enddo
!       write(*,*) presi_sam(31)
!       do i = 1,90
!         write(*,*) i, var_cam(1,i), rho_cam(1,i), rho_cam_exp(1,i), p_int_cam(1, i)
!         ! write(*,*) i, p_int_cam(1,i)
!         ! write(*,*) i, p_cam(1,i)
!       enddo
!       write(*,*) p_int_cam(1, 91)
!       do j=0, 30
!           write(*,*) var_sam(1, j+1), sum(var_cam(1, 3*j+1:3*j+3)), var_cam(1, 3*j+1:3*j+3)
!       enddo
      call assert_array_equal(rho_cam, rho_cam_exp, test_name, rtol_opt=2.0E-5)

      ! Compare sums of the variables to check conservation between grids
      do i=1,4
        sum_sam(i) = sum(var_sam(i,:))
        sum_cam(i) = sum(var_cam(i,:))
      enddo
      call assert_array_equal(sum_sam, sum_cam, test_name//": integrated sum")

    end subroutine test_interp_to_cam_fine

    subroutine test_var_conv_sam_zero(test_name)
      !! Check variable conversion SAM->CAM with 0.0 results in 0.0 on the other side

      character(len=*), intent(in) :: test_name

      real, dimension(4, 48) :: t, q, tabs, qv, qc, qi
      real, dimension(4, 48) :: tabs_exp, qv_exp, qc_exp, qi_exp
      integer :: i
      real, dimension(48) :: pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam
          !! Data from the SAM soundings used in tests

      ! Fetch SAM grid data
      call fetch_sam_data(pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam)

      t = 0.0
      q = 0.0
      tabs = -100.0
      qv = -100.0
      qc = -100.0
      qi = -100.0

      call SAM_var_conversion(t, q, tabs, qv, qc, qi)
        !! Convert SAM t and q to tabs, qv, qc, qi used by CAM

      ! Add 1.0 to q values (they should be 0.0) to avoid division by 0.0 errors in check
      qv = qv + 1.0
      qc = qc + 1.0
      qi = qi + 1.0

      ! Check that, without moisture, tabs is (t - gamaz)
      do i=1,48
        tabs_exp(:,i) = 0.0 - gamaz_sam(i)
      end do
      qv_exp = 1.0
      qc_exp = 1.0
      qi_exp = 1.0

      call assert_array_equal(tabs, tabs_exp, test_name)
      call assert_array_equal(qv, qv_exp, test_name//": qv")
      call assert_array_equal(qc, qc_exp, test_name//": qc")
      call assert_array_equal(qi, qi_exp, test_name//": qi")

    end subroutine test_var_conv_sam_zero

    subroutine test_rev_var_conv_sat(test_name)
      !! Check variable conversion SAM->CAM with 0.0 results in 0.0 on the other side

      character(len=*), intent(in) :: test_name

      real, dimension(4, 48) :: t, q, tabs, qv, qc, qi
      real, dimension(4, 48) :: tabs_exp, qv_exp, qc_exp, qi_exp
      integer :: i
      real, dimension(48) :: pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam

      real, parameter :: rgas = 287.0
          !! Gas constant for dry air used in SAM [j / kg / K]

      ! Fetch SAM grid data
      call fetch_sam_data(pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam)

      ! Set up dry saturated air
      do i = 1,48
        tabs(:, i) = 100.0 * pres_sam(i) / (rho_sam(i) * rgas)
      end do
      qv = 0.0
      qc = 0.0
      qi = 0.0

      tabs_exp = -100.0
      qv_exp = -100.0
      qc_exp = -100.0
      qi_exp = -100.0

      call CAM_var_conversion(qv, qc, qi, q, tabs, t)
      call SAM_var_conversion(t, q, tabs_exp, qv_exp, qc_exp, qi_exp)

      ! Add 1.0 to q values (they should be 0.0) to avoid division by 0.0 errors in check
      qv = qv + 1.0
      qc = qc + 1.0
      qi = qi + 1.0

      qv_exp = 1.0
      qc_exp = 1.0
      qi_exp = 1.0

      call assert_array_equal(tabs, tabs_exp, test_name//": tabs")
      call assert_array_equal(qv, qv_exp, test_name//": qv")
      call assert_array_equal(qc, qc_exp, test_name//": qc")
      call assert_array_equal(qi, qi_exp, test_name//": qi")

    end subroutine test_rev_var_conv_sat

    subroutine test_rev_var_conv_moist(test_name)
      !! Check variable conversion SAM->CAM with 0.0 results in 0.0 on the other side

      character(len=*), intent(in) :: test_name

      real, dimension(4, 48) :: t, q, tabs, qv, qc, qi
      real, dimension(4, 48) :: r, p_sat, q_sat
      real, dimension(4, 48) :: tabs_exp, qv_exp, qc_exp, qi_exp
      integer :: i
      real, dimension(48) :: pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam

      real, parameter :: rgas = 287.0
          !! Gas constant for dry air used in SAM [j / kg / K]
      real, parameter :: rvap = 461.0
          !! Gas constant for vapour used in SAM [j / kg / K]

      ! Fetch SAM grid data
      call fetch_sam_data(pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam)

      ! Set up moist air
      do i = 1,48
        qv(:, i) = 0.01 * exp(-5.0d-4 * z_sam(i))
        r(:, i) = qv(:, i) / (1.0 - qv(:, i))

        ! WRONG FROM REFERENCE!!!!!tabs(:, i) = (1.0 + (rgas/rvap)) * 100.0 * pres_sam(i) / (rho_sam(i) * rgas * (1.0 + r(:, i) * (rvap/rgas)))
        ! CRUDE APPROX (actually pretty close) tabs(:, i) =   100.0 * pres_sam(i) / (rho_sam(i) * rgas * (1.0 + 0.61 * r(:, i)))
        tabs(:, i) = 100.0 * pres_sam(i) * (1.0 + r(:, i)) * (rgas/rvap) / (rgas * rho_sam(i) * (r(:, i) + (rgas/rvap)))

        ! Check that we are below saturation using Bolton (1980)
        p_sat(:, i) = 6.112 * exp(17.67 * (tabs(:,i) - 273.15) / (tabs(:,i) - 29.65))
        q_sat(:, i) = 0.622 * (p_sat(:, i) / (pres_sam(i) - p_sat(:, i)))
        ! write(*,*) tabs(1, i), p_sat(1, i), pres_sam(i), q_sat(1, i), qv(1, i)
        if (q_sat(1, i) < qv(1, i)) then
          write(*,*) "Specific humidity above saturation, please review q-profile. STOPPING."
          stop
        endif
      end do

      qc = 0.0
      qi = 0.0

      tabs_exp = -100.0
      qv_exp = -100.0
      qc_exp = -100.0
      qi_exp = -100.0

      call CAM_var_conversion(qv, qc, qi, q, tabs, t)
      call SAM_var_conversion(t, q, tabs_exp, qv_exp, qc_exp, qi_exp)

      ! Add 1.0 to q values (they should be 0.0) to avoid division by 0.0 errors in check
      qc = qc + 1.0
      qi = qi + 1.0

      qc_exp = 1.0
      qi_exp = 1.0

      call assert_array_equal(tabs, tabs_exp, test_name//": tabs")
      call assert_array_equal(qv, qv_exp, test_name//": qv")
      call assert_array_equal(qc, qc_exp, test_name//": qc")
      call assert_array_equal(qi, qi_exp, test_name//": qi")

    end subroutine test_rev_var_conv_moist

end module cam_tests

program run_cam_tests

  use cam_tests

  implicit none

  real, dimension(4) :: a = [1.0, 2.0, 3.0, 4.0]
  real, dimension(5) :: b = [1.0, 2.0, 3.0, 4.0, 6.0]
  integer, dimension(4) :: aa = [1, 2, 3, 4]
  integer, dimension(4) :: bb = [1, 2, 3, 4]

  character(len=1024) :: nn_file = "NN_weights_YOG_convection.nc"
  character(len=1024) :: sounding_file = "resources/SAM_sounding.nc"

  ! Initialise the NN module
  call nn_convection_flux_CAM_init(nn_file, sounding_file)
  
  ! Run tests
  ! CAM to SAM grid interpolation
  call test_interp_to_sam_match("Test interpolation to SAM from matching grid")
  call test_interp_to_sam_one("Test interpolation to SAM with array of ones")
  call test_interp_to_sam_pres("Test interpolation to SAM mapping pressure to pressure")

  ! SAM to CAM grid interpolation
  call test_interp_to_cam_match("Test interpolation to CAM mapping density 1:1 match grid")
  call test_interp_to_cam_fine("Test interpolation to CAM mapping density 1:1 fine grid")

  ! Variable conversion
  call test_var_conv_sam_zero("Test variable conversion SAM->CAM for 0.0")
  call test_rev_var_conv_sat("Test reversible variable conversion for dry air")
  call test_rev_var_conv_moist("Test reversible variable conversion for moist unsaturated air")

  ! Clean up
  call nn_convection_flux_CAM_finalize()

end program run_cam_tests
