module test_utils

  use :: precision, only: wp

  implicit none

  character(len=15) :: pass = char(27)//'[32m'//'PASSED'//char(27)//'[0m'
  character(len=15) :: fail = char(27)//'[31m'//'FAILED'//char(27)//'[0m'

  interface assert_array_equal
    module procedure assert_array_equal_1d, assert_array_equal_2d, assert_array_equal_3d
  end interface

  interface print_assert
    module procedure print_assert
  end interface

  contains

  subroutine print_assert(test_name, is_close, relative_error)

    character(len=*), intent(in) :: test_name
    logical, intent(in) :: is_close
    real(wp), intent(in) :: relative_error

    if (is_close) then
      write(*, '(A, " :: [", A, "] maximum relative error = ", E11.4)') pass, trim(test_name), relative_error
    else
      write(*, '(A, " :: [", A, "] maximum relative error = ", E11.4)') fail, trim(test_name), relative_error
    end if

  end subroutine print_assert

  subroutine assert_array_equal_1d(a, b, test_name, rtol_opt)

    character(len=*), intent(in) :: test_name
    real(wp), intent(in), dimension(:) :: a, b
    real(wp), intent(in), optional :: rtol_opt
    real(wp) :: relative_error, rtol

    if (.not. present(rtol_opt)) then
      rtol = 1.0e-5
    else
      rtol = rtol_opt
    end if

    relative_error = maxval(abs(a/b - 1.0))
    
    call print_assert(test_name, (rtol > relative_error), relative_error)

  end subroutine assert_array_equal_1d

  subroutine assert_array_equal_2d(a, b, test_name, rtol_opt)

    character(len=*), intent(in) :: test_name
    real(wp), intent(in), dimension(:,:) :: a, b
    real(wp), intent(in), optional :: rtol_opt
    real(wp) :: relative_error, rtol

    if (.not. present(rtol_opt)) then
      rtol = 1.0e-5
    else
      rtol = rtol_opt
    end if

    relative_error = maxval(abs(a/b - 1.0))
    call print_assert(test_name, (rtol > relative_error), relative_error)

  end subroutine assert_array_equal_2d

  subroutine assert_array_equal_3d(a, b, test_name, rtol_opt)

    character(len=*), intent(in) :: test_name
    real(wp), intent(in), dimension(:,:,:) :: a, b
    real(wp), intent(in), optional :: rtol_opt
    real(wp) :: relative_error, rtol

    if (.not. present(rtol_opt)) then
      rtol = 1.0e-5
    else
      rtol = rtol_opt
    end if

    relative_error = maxval(abs(a/b - 1.0))
    call print_assert(test_name, (rtol > relative_error), relative_error)

  end subroutine assert_array_equal_3d

end module test_utils
