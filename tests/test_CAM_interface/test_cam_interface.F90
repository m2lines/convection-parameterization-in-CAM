module cam_tests
  !! Module containing individual tests for the CAM ML code

  !--------------------------------------------------------------------------
  ! Libraries to use
  use netcdf
  use precision, only: dp

  use SAM_consts_mod, only: nrf, num_cols, sam_sounding, num_cells, num_cam_cells_fine, num_sam_cells, & 
  num_cam_cells_coarse
  use nn_interface_CAM, only: nn_convection_flux_CAM, nn_convection_flux_CAM_init, nn_convection_flux_CAM_finalize, &
  interp_to_sam, interp_to_cam, fetch_sam_data, SAM_var_conversion, CAM_var_conversion
  use test_utils, only: assert_array_equal

  implicit none

  character(len=15) :: pass = char(27)//'[32m'//'PASSED'//char(27)//'[0m'
  character(len=15) :: fail = char(27)//'[31m'//'FAILED'//char(27)//'[0m'
  integer, parameter :: n_nn_out = 148
  real(dp), dimension(n_nn_out) :: nn_out_ones

  contains

    subroutine test_interp_to_sam_match(test_name)
      !! Check interpolation to SAM grid by defining an idential CAM grid and
      !! interpolating a variable equal to the pressure
      !! Define a CAM grid consiting of 4 atmospheric columns
      !! from 1111.0 to 10.0

      character(len=*), intent(in) :: test_name
      
      integer :: i

      real(dp), dimension(num_cols, nrf) :: p_cam, p_int_cam
      real(dp), dimension(num_cols, nrf) :: var_cam
      real(dp), dimension(num_cols) :: var_cam_surface
      real(dp), dimension(num_cols) :: ps_cam
      real(dp), dimension(num_cols, nrf) :: var_sam, var_sam_exp

      real(dp), dimension(sam_sounding) :: pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam
          !! Data from the SAM soundings used in tests

      ! Fetch SAM grid data
      call fetch_sam_data(pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam)

      do i = 1, num_cols
          p_cam(i, :) = pres_sam(1 : nrf)
          p_int_cam(i, :) = presi_sam(1 : nrf)
          ! Set SAM variable equal to cell size (density 1.0)
          var_cam(i, :) = pres_sam(1 : nrf)
      enddo

      ps_cam(:) = presi_sam(1)
      var_cam_surface(:) = presi_sam(1)

      call interp_to_sam(p_cam, ps_cam, var_cam, var_sam, var_cam_surface)

      ! Compare the results of the interpolation scheme to expected output
      ! Set anything above 30 elems to zero as the parameterization and interpolation 
      ! code only uses the bottom 30 cells on the SAM grid
      var_sam_exp = var_cam
      call assert_array_equal(var_sam, var_sam_exp, test_name)

    end subroutine test_interp_to_sam_match

    subroutine test_interp_to_sam_one(test_name)
      !! Check interpolation to SAM grid by setting up a coarse CAM grid with every variable being 1.0
      !! interpolating to a new grid should also have of value 1.0 everywhere.

      character(len=*), intent(in) :: test_name

      integer :: i

      real(dp), dimension(num_cols, num_cells) :: p_cam
      real(dp), dimension(num_cols, num_cells) :: var_cam
      real(dp), dimension(num_cols) :: var_cam_surface
      real(dp), dimension(num_cols) :: ps_cam
      real(dp), dimension(num_cols, num_sam_cells) :: var_sam, var_sam_exp

      ! Set up a coarse CAM grid of 4 columns of pressures [1000, 500, 10] hPa with surface pressure 1111 hPa
      do i = 1, num_cols
          p_cam(i, 1) = 1000.0
          p_cam(i, 2) = 500.0
          p_cam(i, 3) = 10.0
      end do
      ps_cam = (/ 1111.0, 1111.0, 1111.0, 1111.0/)
      
      var_cam = 1.0
      var_cam_surface = 1.0
      var_sam_exp = 1.0

      call interp_to_sam(p_cam, ps_cam, var_cam, var_sam, var_cam_surface)
      
      call assert_array_equal(var_sam, var_sam_exp, test_name)

    end subroutine test_interp_to_sam_one

    subroutine test_interp_to_sam_pres(test_name)
      !! Check interpolation to SAM grid by interpolating pressure to pressure
      !! Use a coarse CAM grid of 3 cells and 4 columns
      !! => expected variable on SAM grid should be pressure at that point

      character(len=*), intent(in) :: test_name

      integer :: i

      real(dp), dimension(num_cols, num_cells) :: p_cam
      real(dp), dimension(num_cols, num_cells) :: var_cam
      real(dp), dimension(num_cols) :: var_cam_surface
      real(dp), dimension(num_cols) :: ps_cam
      real(dp), dimension(num_cols, num_sam_cells) :: var_sam, var_sam_exp

      real(dp), dimension(sam_sounding) :: pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam
          !! Data from the SAM soundings used in tests

      ! Fetch SAM grid data
      call fetch_sam_data(pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam)

      ! Set up coarse CAM grid with surface pressure equal to SAM surface,
      ! a top pressure of 1.0d-4 above the SAM grid, and other at 10 and 20 from the SAM grid
      ps_cam = presi_sam(1)
      p_cam(:, 1) = pres_sam(10)
      p_cam(:, 2) = pres_sam(20)
      p_cam(:, 3) = 1.0d-4

      ! Set the variable on the CAM grid equal to the pressure
      var_cam_surface = ps_cam
      var_cam = p_cam

      ! Expected variable value on the SAM grid will be equal to pressure at that point
      do i = 1, num_cols
        var_sam_exp(i, :) = pres_sam(1 : num_sam_cells)
      end do

      call interp_to_sam(p_cam, ps_cam, var_cam, var_sam, var_cam_surface)

      call assert_array_equal(var_sam_exp, var_sam, test_name, rtol_opt = 1.0d-14)

    end subroutine test_interp_to_sam_pres

    subroutine test_interp_to_cam_match(test_name)
      !! Check conservative regridding to CAM grid by interpolating constant density
      !! => should conserve density
      !! For case CAM and SAM grids match and CAM grid inside SAM grid

      character(len=*), intent(in) :: test_name

      integer :: i

      real(dp), dimension(num_cols, nrf) :: p_cam, var_cam, var_cam_exp
      real(dp), dimension(num_cols, nrf + 1) :: p_int_cam
      real(dp), dimension(num_cols) :: ps_cam, var_cam_surface
      real(dp), dimension(num_cols, nrf) :: var_sam

      real(dp), dimension(sam_sounding) :: pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam
          !! Data from the SAM soundings used in tests

      ! Fetch SAM grid data
      call fetch_sam_data(pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam)

      do i = 1, num_cols
          ! Set up a CAM grid that matches the lower nrf cells of the SAM grid
          p_cam(i, :) = pres_sam(1 : num_sam_cells)
          p_int_cam(i, :) = presi_sam(1 : nrf + 1)
          ! Set SAM (variable) density equal to 1.0
          var_sam(i, :) = 1.0
      enddo

      call interp_to_cam(p_cam, p_int_cam, p_int_cam(:, 1), var_sam, var_cam)
      
      var_cam_exp = 1.0

      call assert_array_equal(var_cam_exp, var_cam, test_name)
    
    end subroutine test_interp_to_cam_match

    subroutine test_interp_to_cam_coarse(test_name)
      !! Check conservative regridding to CAM coarse grid
      !! => should conserve density

      character(len=*), intent(in) :: test_name
      integer :: i, j

      real(dp), dimension(num_cols, num_cam_cells_coarse) :: p_cam, var_cam, var_cam_exp
      real(dp), dimension(num_cols, num_cam_cells_coarse + 1) :: p_int_cam
      real(dp), dimension(num_cols) :: ps_cam, sum_sam, sum_cam, var_cam_surface
      real(dp), dimension(num_cols, num_sam_cells) :: var_sam 

      real(dp), dimension(sam_sounding) :: pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam
        !! Data from the SAM soundings used in tests

      ! Fetch SAM grid data
      call fetch_sam_data(pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam)

      ! Define CAM grid coarser than SAM grid
      do j = 1, num_cam_cells_coarse + 1
          p_int_cam(:,j) = presi_sam(3 * (j - 1) + 1)
      enddo

      ! Get the CAM pressures from average interface pressures
      do j = 1, num_cam_cells_coarse 
        p_cam(:, j) = (p_int_cam(:, j + 1) + p_int_cam(:, j)) / 2.0
      end do

      do j = 1, num_cols
        ! Set SAM variable (density) equal to 1.0
        var_sam(j, :) = 1.0
      end do

      call interp_to_cam(p_cam, p_int_cam, p_int_cam(:, 1), var_sam, var_cam)

      ! Expected density 
      var_cam_exp = 1.0

      call assert_array_equal(var_cam, var_cam_exp, test_name, rtol_opt = 2.0D-5)

      ! Compare sums of the variables to check conservation between grids
      do i = 1, num_cols
        sum_sam(i) = sum(var_sam(i, :) * (presi_sam(1 : num_sam_cells) - presi_sam(2 : num_sam_cells + 1)))
        sum_cam(i) = sum(var_cam(i, :) * (p_int_cam(i, 1 : num_cam_cells_coarse) - p_int_cam(i, 2 : num_cam_cells_coarse + 1)))
      enddo

      call assert_array_equal(sum_sam, sum_cam, test_name//": integrated sum")

    end subroutine test_interp_to_cam_coarse

    subroutine test_interp_to_cam_coarse_variable_density(test_name)
      !! Check conservative regridding to CAM coarse grid with variable density
      !! => With integrated sum, this test should conserve density

      character(len=*), intent(in) :: test_name
      integer :: i, j

      real(dp), dimension(num_cols, num_cam_cells_coarse) :: p_cam, var_cam, var_cam_exp
      real(dp), dimension(num_cols, num_cam_cells_coarse + 1) :: p_int_cam
      real(dp), dimension(num_cols) :: ps_cam, sum_sam, sum_cam, var_cam_surface
      real(dp), dimension(num_cols, num_sam_cells) :: var_sam 

      real(dp), dimension(sam_sounding) :: pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam
        !! Data from the SAM soundings used in tests

      ! Fetch SAM grid data
      call fetch_sam_data(pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam)

      ! Define CAM grid coarser than SAM grid
      do j = 1, num_cam_cells_coarse + 1
          p_int_cam(:,j) = presi_sam(3 * (j - 1) + 1)
      enddo

      ! Get the CAM pressures from average interface pressures
      do j = 1, num_cam_cells_coarse 
        p_cam(:, j) = (p_int_cam(:, j + 1) + p_int_cam(:, j)) / 2.0
      end do

      do j = 1, num_cols
          do i = 1, num_sam_cells
            ! Set SAM variable density not equal to 1.0
            var_sam(j, i) = tan(real(i * j))
          end do
      end do

      call interp_to_cam(p_cam, p_int_cam, p_int_cam(:, 1), var_sam, var_cam)

      ! Compare sums of the variables to check conservation between grids
      do i = 1, num_cols
        sum_sam(i) = sum(var_sam(i, :) * (presi_sam(1 : num_sam_cells) - presi_sam(2 : num_sam_cells + 1)))
        sum_cam(i) = sum(var_cam(i, :) * (p_int_cam(i, 1 : num_cam_cells_coarse) - p_int_cam(i, 2 : num_cam_cells_coarse + 1)))
      enddo

      call assert_array_equal(sum_sam, sum_cam, test_name//": integrated sum")

    end subroutine test_interp_to_cam_coarse_variable_density

    subroutine test_interp_to_cam_fine(test_name)
      !! Check conservative regridding to CAM grid by interpolating constant density
      !! => should conserve density
      !! For case CAM and SAM grids match and CAM grid inside SAM grid

      character(len=*), intent(in) :: test_name

      integer :: i, j

      real(dp), dimension(num_cols, num_cam_cells_fine) :: p_cam, var_cam, var_cam_exp
      real(dp), dimension(num_cols, num_cam_cells_fine + 1) :: p_int_cam
      real(dp), dimension(num_cols) :: var_cam_surface
      real(dp), dimension(num_cols) :: ps_cam, sum_sam, sum_cam
      real(dp), dimension(num_cols, num_sam_cells) :: var_sam

      real(dp), dimension(sam_sounding) :: pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam
          !! Data from the SAM soundings used in tests

      ! Fetch SAM grid data
      call fetch_sam_data(pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam)

      ! CAM grid finer than SAM grid
      do j = 0, num_sam_cells - 1
        do i = 1, num_cols
          p_int_cam(i, 3*j+1) = presi_sam(j+1)
          p_int_cam(i, 3*j+2) = presi_sam(j+1) + (presi_sam(j+2) - presi_sam(j+1)) / 3
          p_int_cam(i, 3*j+3) = presi_sam(j+1) + 2*(presi_sam(j+2) - presi_sam(j+1)) / 3
        enddo
      enddo
      p_int_cam(:, num_cam_cells_fine + 1) = presi_sam(num_sam_cells + 1)
     
      do j = 1, num_cam_cells_fine
        p_cam(:, j) = (p_int_cam(:, j+1) + p_int_cam(:, j)) / 2.0
      end do

      do j = 1, num_cols
          ! Set SAM density equal to 1.0
          var_sam(j, :) = 1.0
      end do

      call interp_to_cam(p_cam, p_int_cam, p_int_cam(:, 1), var_sam, var_cam)

      var_cam_exp = 1.0

      call assert_array_equal(var_cam, var_cam_exp, test_name, rtol_opt = 2.0D-5)

      ! Compare sums of the variables to check conservation between grids
      do i = 1, num_cols
        sum_sam(i) = sum(var_sam(i, :) * (presi_sam(1 : num_sam_cells) - presi_sam(2 : num_sam_cells + 1)))
        sum_cam(i) = sum(var_cam(i, :) * (p_int_cam(i, 1 : num_cam_cells_fine) - p_int_cam(i, 2 : num_cam_cells_fine + 1)))
      enddo
    
      call assert_array_equal(sum_sam, sum_cam, test_name//": integrated sum")

    end subroutine test_interp_to_cam_fine

    subroutine test_interp_to_cam_fine_variable_density(test_name)
      !! Check conservative regridding to CAM fine grid with variable density
      !! => With integrated sum, this test should conserve density

      character(len=*), intent(in) :: test_name

      integer :: i, j

      real(dp), dimension(num_cols, num_cam_cells_fine) :: p_cam, var_cam, var_cam_exp
      real(dp), dimension(num_cols, num_cam_cells_fine + 1) :: p_int_cam
      real(dp), dimension(num_cols) :: var_cam_surface
      real(dp), dimension(num_cols) :: ps_cam, sum_sam, sum_cam
      real(dp), dimension(num_cols, num_sam_cells) :: var_sam

      real(dp), dimension(sam_sounding) :: pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam
          !! Data from the SAM soundings used in tests

      ! Fetch SAM grid data
      call fetch_sam_data(pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam)

      ! CAM grid finer than SAM grid
      do j = 0, num_sam_cells - 1
        do i = 1, num_cols
          p_int_cam(i, 3*j+1) = presi_sam(j+1)
          p_int_cam(i, 3*j+2) = presi_sam(j+1) + (presi_sam(j+2) - presi_sam(j+1)) / 3
          p_int_cam(i, 3*j+3) = presi_sam(j+1) + 2*(presi_sam(j+2) - presi_sam(j+1)) / 3
        enddo
      enddo
      p_int_cam(:, num_cam_cells_fine + 1) = presi_sam(num_sam_cells + 1)
     
      do j = 1, num_cam_cells_fine
        p_cam(:, j) = (p_int_cam(:, j+1) + p_int_cam(:, j)) / 2.0
      end do

      do j = 1, num_cols
        do i = 1, num_sam_cells
          ! Set SAM to variable density not equal to 1.0
          var_sam(j, i) = exp(real(i * j))
        end do
      end do

      call interp_to_cam(p_cam, p_int_cam, p_int_cam(:, 1), var_sam, var_cam)

      ! Compare sums of the variables to check conservation between grids
      do i = 1, num_cols
        sum_sam(i) = sum(var_sam(i, :) * (presi_sam(1 : num_sam_cells) - presi_sam(2 : num_sam_cells + 1)))
        sum_cam(i) = sum(var_cam(i, :) * (p_int_cam(i, 1 : num_cam_cells_fine) - p_int_cam(i, 2 : num_cam_cells_fine + 1)))
      enddo
    
      call assert_array_equal(sum_sam, sum_cam, test_name//": integrated sum")

    end subroutine test_interp_to_cam_fine_variable_density      

    subroutine test_var_conv_sam_zero(test_name)
      !! Check variable conversion SAM->CAM with 0.0 results in 0.0 on the other side

      character(len=*), intent(in) :: test_name

      real(dp), dimension(num_cols, sam_sounding) :: t, q, tabs, qv, qc, qi
      real(dp), dimension(num_cols, sam_sounding) :: tabs_exp, qv_exp, qc_exp, qi_exp
      integer :: i
      real(dp), dimension(sam_sounding) :: pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam
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
      do i = 1, sam_sounding
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

      real(dp), dimension(num_cols, sam_sounding) :: t, q, tabs, qv, qc, qi
      real(dp), dimension(num_cols, sam_sounding) :: tabs_exp, qv_exp, qc_exp, qi_exp
      integer :: i
      real(dp), dimension(sam_sounding) :: pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam

      real(dp), parameter :: rgas = 287.0
          !! Gas constant for dry air used in SAM [j / kg / K]

      ! Fetch SAM grid data
      call fetch_sam_data(pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam)

      ! Set up dry saturated air
      do i = 1, sam_sounding
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

      real(dp), dimension(num_cols, sam_sounding) :: t, q, tabs, qv, qc, qi
      real(dp), dimension(num_cols, sam_sounding) :: r, p_sat, q_sat
      real(dp), dimension(num_cols, sam_sounding) :: tabs_exp, qv_exp, qc_exp, qi_exp
      integer :: i
      real(dp), dimension(sam_sounding) :: pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam

      real(dp), parameter :: rgas = 287.0
          !! Gas constant for dry air used in SAM [j / kg / K]
      real(dp), parameter :: rvap = 461.0
          !! Gas constant for vapour used in SAM [j / kg / K]

      ! Fetch SAM grid data
      call fetch_sam_data(pres_sam, presi_sam, gamaz_sam, rho_sam, z_sam)

      ! Set up moist air
      do i = 1, sam_sounding
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

  real(dp), dimension(num_cols) :: a = [1.0, 2.0, 3.0, 4.0]
  real(dp), dimension(5) :: b = [1.0, 2.0, 3.0, 4.0, 6.0]
  integer, dimension(num_cols) :: aa = [1, 2, 3, 4]
  integer, dimension(num_cols) :: bb = [1, 2, 3, 4]

  character(len=1024) :: nn_file = "../../YOG_convection/NN_weights_YOG_convection.nc"
  character(len=1024) :: sounding_file = "../../YOG_convection/resources/SAM_sounding.nc"

  ! Initialise the NN module
  call nn_convection_flux_CAM_init(nn_file, sounding_file)
  
  ! Run tests
  ! CAM to SAM grid interpolation
  call test_interp_to_sam_match("Test interpolation to SAM from matching grid")
  call test_interp_to_sam_one("Test interpolation to SAM with array of ones")
  call test_interp_to_sam_pres("Test interpolation to SAM mapping pressure to pressure")

  ! SAM to CAM grid interpolation
  call test_interp_to_cam_match("Test interpolation to CAM mapping density 1:1 match grid")
  call test_interp_to_cam_coarse("Test interpolation to CAM mapping density 1:1 coarse grid")
  call test_interp_to_cam_fine("Test interpolation to CAM mapping density 1:1 fine grid")

  ! SAM to CAM grid interpolation with variable densities
  call test_interp_to_cam_coarse_variable_density("Test interpolation to CAM mapping with variable density coarse grid")
  call test_interp_to_cam_fine_variable_density("Test interpolation to CAM mapping with variable density fine grid")

  ! Variable conversion
  call test_var_conv_sam_zero("Test variable conversion SAM->CAM for 0.0")
  call test_rev_var_conv_sat("Test reversible variable conversion for dry air")
  call test_rev_var_conv_moist("Test reversible variable conversion for moist unsaturated air")

  ! Clean up
  call nn_convection_flux_CAM_finalize()

end program run_cam_tests
