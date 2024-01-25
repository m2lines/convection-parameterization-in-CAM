module tests
  !! Module containing individual tests for the CAM ML code

  !--------------------------------------------------------------------------
  ! Libraries to use
  use netcdf
  use nn_cf_net_mod, only: relu, net_forward, nn_cf_net_init, nn_cf_net_finalize
  use test_utils, only: assert_array_equal

  implicit none

  character(len=15) :: pass = char(27)//'[32m'//'PASSED'//char(27)//'[0m'
  character(len=15) :: fail = char(27)//'[31m'//'FAILED'//char(27)//'[0m'

  real(4), dimension(148) :: nn_out_ones

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

      call nn_cf_net_init(nn_filename, nin, nout, 30)

      call nn_cf_net_finalize()

      if (nin == 61) then
        write(*, '(A, " :: [", A, " - nin]")') pass, trim(test_name)
        if (nout == 148) then
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
      real(4) :: nn_out(148)

      nn_in = 1.0

      nn_filename = "./NN_weights_YOG_convection.nc"
      call nn_cf_net_init(nn_filename, nin, nout, 30)

      call net_forward(nn_in, nn_out)

      call nn_cf_net_finalize()

      call assert_array_equal(nn_out, nn_out_ones, test_name, 1.0e-4)

    end subroutine test_nn

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
      do i = 1,148
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

end program run_tests
