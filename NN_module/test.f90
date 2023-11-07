module tests
  !! Module containing individual tests for the CAM ML code

  !--------------------------------------------------------------------------
  ! Libraries to use
  use netcdf
  use nn_cf_net, only: relu, net_forward, nn_cf_net_init, nn_cf_net_finalize
    use nn_cf_net_torch, only: nn_cf_net_torch_init, nn_cf_net_torch_finalize  use nn_convection_flux, only:   nn_convection_flux, nn_convection_flux_init, nn_convection_flux_finalize
  use test_utils, only: assert_array_equal

  implicit none

  character(len=15) :: pass = char(27)//'[32m'//'PASSED'//char(27)//'[0m'
  character(len=15) :: fail = char(27)//'[31m'//'FAILED'//char(27)//'[0m'
  integer, parameter :: nrf = 30
  integer, parameter :: n_nn_out = 148
  real(4), dimension(n_nn_out) :: nn_out_ones

  contains

    subroutine test_relu(test_name)
      !! Test relu function is working as expected

      character(len=*), intent(in) :: test_name
      real(4), dimension(4) :: test_array = (/ -1.0, 0.0, 0.5, 1.0 /)
      real(4), dimension(4) :: res_array = (/ 0.0, 0.0, 0.5, 1.0 /)

      call relu(test_array)

      call assert_array_equal(res_array, test_array, test_name)

    end subroutine test_relu

    subroutine test_nn_cf_init(test_name)
      !! Test NN initialisation is working as expected
      !! Checks that nin and nout are read in as expected
      !! Checking deeper requires interaction with nn_cf_net_mod module vars

      character(len=*), intent(in) :: test_name
      character(len=1024) :: nn_filename
      integer :: nin, nout

      nn_filename = "./NN_weights_YOG_convection.nc"

      call nn_cf_net_init(nn_filename, nin, nout, nrf)

      call nn_cf_net_finalize()

      if (nin == 61) then
        write(*, '(A, " :: [", A, " - nin]")') pass, trim(test_name)
        if (nout == n_nn_out) then
          write(*, '(A, " :: [", A, " - nout]")') pass, trim(test_name)
        else
          write(*, '(A, " :: [", A, "] with nout = ", I3)') fail, trim(test_name), nout
        end if
      else
        write(*, '(A, " :: [", A, "] with nin = ", I3)') fail, trim(test_name), nin
      end if
    
    end subroutine test_nn_cf_init

    subroutine test_nn(test_name)
      !! Test NN is producing the expected results.

      integer :: i, nin, nout
      character(len=*), intent(in) :: test_name
      character(len=1024) :: nn_filename
      real(4) :: nn_in(61)
      real(4) :: nn_out(n_nn_out)

      nn_in = 1.0

      nn_filename = "./NN_weights_YOG_convection.nc"
      call nn_cf_net_init(nn_filename, nin, nout, nrf)

      call net_forward(nn_in, nn_out)

      call nn_cf_net_finalize()

      call assert_array_equal(nn_out, nn_out_ones, test_name, 1.0e-4)

    end subroutine test_nn

    subroutine test_param(test_name)
      !! Test Parameterisation is producing the same results as it initially did.
      !! Run for a single column with physically plausible parameters.

      integer :: i, io, stat
      character(len=*), intent(in) :: test_name
      character(len=1024) :: nn_filename
      character(len=512) :: msg

      real :: tabs_i(1, 1, 48) = 293.15
      real :: q_i(1, 1, 48) = 0.5
      real :: y_in(1, 48) = 0.0
      real :: tabs(1, 1, 48) = 293.15
      real :: t_0(1, 1, 48) =  1.0e4
      real :: q_0(1, 1, 48) = 0.5
      real :: rho(48) = 1.2
      real :: adz(48) = 1.0
      real :: dz = 100.0
      real :: dtn = 2.0
      real, dimension(1,1,nrf) :: t_delta_adv, q_delta_adv, &
                                t_delta_auto, q_delta_auto, &
                                t_delta_sed, q_delta_sed
      real :: t_rad_rest_tend(1,1,nrf)
      real :: prec_sed(1,1)

      real, dimension(1,1,nrf) :: t_delta_adv_dat, q_delta_adv_dat, &
                                t_delta_auto_dat, q_delta_auto_dat, &
                                t_delta_sed_dat, q_delta_sed_dat
      real :: t_rad_rest_tend_dat(1,1,nrf)
      real :: prec_sed_dat(1,1)
      
      nn_filename = "./NN_weights_YOG_convection.nc"
      call nn_convection_flux_init(nn_filename)

      call nn_convection_flux(tabs_i, q_i, y_in, &
                                  tabs, &
                                  t_0, q_0, &
                                  rho, adz, dz, dtn, &
                                  t_rad_rest_tend, &
                                  t_delta_adv, q_delta_adv, &
                                  t_delta_auto, q_delta_auto, &
                                  t_delta_sed, q_delta_sed, prec_sed)

      call nn_convection_flux_finalize()
      nn_filename = "param_test.txt"

      ! Writing data out to file from original code runniing
      ! open(newunit=io, file=trim(nn_filename), status="replace", action="write", &
      !      iostat=stat, iomsg=msg)
      ! if (stat /= 0) then
      !   print *, trim(msg)
      !   stop
      ! end if
      ! do i = 1,nrf
      !   write(io, '(7E18.8)') t_delta_adv(1,1,i), q_delta_adv(1,1,i), &
      !                t_delta_auto(1,1,i), q_delta_auto(1,1,i), &
      !                t_delta_sed(1,1,i), q_delta_sed(1,1,i), t_rad_rest_tend(1,1,i)
      ! enddo
      ! write(io, '(E18.8)') prec_sed(1,1)
      ! close(io)

      open(newunit=io, file=trim(nn_filename), status="old", action="read", &
           iostat=stat, iomsg=msg)
      if (stat /= 0) then
        print *, trim(msg)
        stop
      end if
      do i = 1,nrf
        read(io, '(7E18.8)') t_delta_adv_dat(1,1,i), q_delta_adv_dat(1,1,i), &
                     t_delta_auto_dat(1,1,i), q_delta_auto_dat(1,1,i), &
                     t_delta_sed_dat(1,1,i), q_delta_sed_dat(1,1,i), &
                     t_rad_rest_tend_dat(1,1,i)
      enddo
      read(io, '(E18.8)') prec_sed_dat(1,1)
      close(io)

      call assert_array_equal(t_delta_adv, t_delta_adv_dat, test_name//" t adv", 1.0e-6)
      call assert_array_equal(q_delta_adv, q_delta_adv_dat, test_name//" q adv", 1.0e-6)
      call assert_array_equal(t_delta_auto, t_delta_auto_dat, test_name//" t auto", 1.0e-6)
      call assert_array_equal(q_delta_auto, q_delta_auto_dat, test_name//" q auto", 1.0e-6)
      call assert_array_equal(t_delta_sed, t_delta_sed_dat, test_name//" t sed", 1.0e-6)
      call assert_array_equal(q_delta_sed, q_delta_sed_dat, test_name//" q sed", 1.0e-6)
      call assert_array_equal(t_rad_rest_tend, t_rad_rest_tend_dat, test_name//" t rad", 1.0e-6)
      call assert_array_equal(prec_sed, prec_sed_dat, test_name//" prec", 1.0e-6)

    end subroutine test_param

    subroutine load_nn_out_ones(nn_ones_file)
      !! Load the result of running the NN with ones

      integer :: io, stat, i
      character(len=512) :: msg
      character(len=*) :: nn_ones_file

      open(newunit=io, file=trim(nn_ones_file), status="old", action="read", &
           iostat=stat, iomsg=msg)
      if (stat /= 0) then
        print *, trim(msg)
        stop
      end if
      do i = 1,n_nn_out
        read(io, *) nn_out_ones(i)
      enddo
      close(io)

    end subroutine load_nn_out_ones

end module tests



program run_tests

  use tests

  implicit none

  character(len=1024) :: nn_filename

  ! Test NN routines

  call test_relu("test_relu")
  call test_nn_cf_init("test_nn_cf_init")

  call load_nn_out_ones("nn_ones.txt")
  call test_nn("Test NN ones")

  call test_param("Test param simple")

end program run_tests
