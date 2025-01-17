module test_utils

  use :: precision, only: sp, dp

  implicit none

  character(len=15) :: pass = char(27)//'[32m'//'PASSED'//char(27)//'[0m'
  character(len=15) :: fail = char(27)//'[31m'//'FAILED'//char(27)//'[0m'

  interface assert_array_equal
    module procedure &
      assert_array_equal_1d_sp, assert_array_equal_2d_sp, assert_array_equal_3d_sp, &
      assert_array_equal_1d_dp, assert_array_equal_2d_dp, assert_array_equal_3d_dp, &
      assert_array_equal_1d_int, assert_array_equal_2d_int, assert_array_equal_3d_int
  end interface

  interface print_assert
    module procedure print_assert_sp, print_assert_dp
  end interface

  contains

  subroutine print_assert_sp(test_name, is_close, relative_error)

    character(len=*), intent(in) :: test_name
    logical, intent(in) :: is_close
    real(sp), intent(in) :: relative_error

    if (is_close) then
      write(*, '(A, " :: [", A, "] maximum relative error = ", E11.4)') pass, trim(test_name), relative_error
    else
      write(*, '(A, " :: [", A, "] maximum relative error = ", E11.4)') fail, trim(test_name), relative_error
      stop 999
    end if

  end subroutine print_assert_sp

  subroutine print_assert_dp(test_name, is_close, relative_error)

    character(len=*), intent(in) :: test_name
    logical, intent(in) :: is_close
    real(dp), intent(in) :: relative_error

    if (is_close) then
      write(*, '(A, " :: [", A, "] maximum relative error = ", E11.4)') pass, trim(test_name), relative_error
    else
      write(*, '(A, " :: [", A, "] maximum relative error = ", E11.4)') fail, trim(test_name), relative_error
      stop 999
    end if

  end subroutine print_assert_dp

  subroutine print_assert_int(test_name, match, abs_error)

    character(len=*), intent(in) :: test_name
    logical, intent(in) :: match
    integer, intent(in) :: abs_error

    if (match) then
      write(*, '(A, " :: [", A, "] maximum absolute error = ", I8)') pass, trim(test_name), abs_error
    else
      write(*, '(A, " :: [", A, "] maximum absolute error = ", I8)') fail, trim(test_name), abs_error
      stop 999
    end if

  end subroutine print_assert_int

  subroutine assert_array_equal_1d_sp(a, b, test_name, rtol_opt)

    character(len=*), intent(in) :: test_name
    real(sp), intent(in), dimension(:) :: a, b
    real(sp), intent(in), optional :: rtol_opt
    real(sp) :: relative_error, rtol
    integer :: shape_error

    shape_error = maxval(abs(shape(a) - shape(b)))

    if (shape_error > 0) then
      write(*, '(A, " :: [", A, "] Arrays have mismatching shapes.")') fail, trim(test_name)
    else
      if (.not. present(rtol_opt)) then
        rtol = 1.0e-5
      else
        rtol = rtol_opt
      end if

      relative_error = maxval(abs(a/b - 1.0))
      call print_assert(test_name, (rtol > relative_error), relative_error)
    endif

  end subroutine assert_array_equal_1d_sp

  subroutine assert_array_equal_2d_sp(a, b, test_name, rtol_opt)

    character(len=*), intent(in) :: test_name
    real(sp), intent(in), dimension(:,:) :: a, b
    real(sp), intent(in), optional :: rtol_opt
    real(sp) :: relative_error, rtol
    integer :: shape_error

    shape_error = maxval(abs(shape(a) - shape(b)))

    if (shape_error > 0) then
      write(*, '(A, " :: [", A, "] Arrays have mismatching shapes.")') fail, trim(test_name)
    else
      if (.not. present(rtol_opt)) then
        rtol = 1.0e-5
      else
        rtol = rtol_opt
      end if
      relative_error = maxval(abs(a/b - 1.0))
      call print_assert(test_name, (rtol > relative_error), relative_error)
    end if

  end subroutine assert_array_equal_2d_sp

  subroutine assert_array_equal_3d_sp(a, b, test_name, rtol_opt)

    character(len=*), intent(in) :: test_name
    real(sp), intent(in), dimension(:,:,:) :: a, b
    real(sp), intent(in), optional :: rtol_opt
    real(sp) :: relative_error, rtol
    integer :: shape_error

    shape_error = maxval(abs(shape(a) - shape(b)))

    if (shape_error > 0) then
      write(*, '(A, " :: [", A, "] Arrays have mismatching shapes.")') fail, trim(test_name)
    else
      if (.not. present(rtol_opt)) then
        rtol = 1.0e-5
      else
        rtol = rtol_opt
      end if
      relative_error = maxval(abs(a/b - 1.0))
      call print_assert(test_name, (rtol > relative_error), relative_error)
    end if

  end subroutine assert_array_equal_3d_sp

  subroutine assert_array_equal_1d_dp(a, b, test_name, rtol_opt)

    character(len=*), intent(in) :: test_name
    real(dp), intent(in), dimension(:) :: a, b
    real(dp), intent(in), optional :: rtol_opt
    real(dp) :: relative_error, rtol
    integer :: shape_error

    shape_error = maxval(abs(shape(a) - shape(b)))

    if (shape_error > 0) then
      write(*, '(A, " :: [", A, "] Arrays have mismatching shapes.")') fail, trim(test_name)
    else
      if (.not. present(rtol_opt)) then
        rtol = 1.0e-5
      else
        rtol = rtol_opt
      end if
      relative_error = maxval(abs(a/b - 1.0))
      call print_assert(test_name, (rtol > relative_error), relative_error)
    end if

  end subroutine assert_array_equal_1d_dp

  subroutine assert_array_equal_2d_dp(a, b, test_name, rtol_opt)

    character(len=*), intent(in) :: test_name
    real(dp), intent(in), dimension(:,:) :: a, b
    real(dp), intent(in), optional :: rtol_opt
    real(dp) :: relative_error, rtol
    integer :: shape_error

    shape_error = maxval(abs(shape(a) - shape(b)))

    if (shape_error > 0) then
      write(*, '(A, " :: [", A, "] Arrays have mismatching shapes.")') fail, trim(test_name)
    else
      if (.not. present(rtol_opt)) then
        rtol = 1.0e-5
      else
        rtol = rtol_opt
      end if
      relative_error = maxval(abs(a/b - 1.0))
      call print_assert(test_name, (rtol > relative_error), relative_error)
    end if

  end subroutine assert_array_equal_2d_dp

  subroutine assert_array_equal_3d_dp(a, b, test_name, rtol_opt)

    character(len=*), intent(in) :: test_name
    real(dp), intent(in), dimension(:,:,:) :: a, b
    real(dp), intent(in), optional :: rtol_opt
    real(dp) :: relative_error, rtol
    integer :: shape_error

    shape_error = maxval(abs(shape(a) - shape(b)))

    if (shape_error > 0) then
      write(*, '(A, " :: [", A, "] Arrays have mismatching shapes.")') fail, trim(test_name)
    else
      if (.not. present(rtol_opt)) then
        rtol = 1.0e-5
      else
        rtol = rtol_opt
      end if
      relative_error = maxval(abs(a/b - 1.0))
      call print_assert(test_name, (rtol > relative_error), relative_error)
    end if

  end subroutine assert_array_equal_3d_dp

  subroutine assert_array_equal_1d_int(a, b, test_name)

    character(len=*), intent(in) :: test_name
    integer, intent(in), dimension(:) :: a, b
    integer :: error, shape_error

    shape_error = maxval(abs(shape(a) - shape(b)))

    if (shape_error > 0) then
      write(*, '(A, " :: [", A, "] Arrays have mismatching shapes.")') fail, trim(test_name)
    else
      error = maxval(abs(a - b))
      if (error > 0) then
        write(*, '(A, " :: [", A, "] maximum absolute error = ", I8)') fail, trim(test_name), error
      else
        write(*, '(A, " :: [", A, "] maximum absolute error = ", I8)') pass, trim(test_name), error
      end if
    end if

  end subroutine assert_array_equal_1d_int

  subroutine assert_array_equal_2d_int(a, b, test_name)

    character(len=*), intent(in) :: test_name
    integer, intent(in), dimension(:,:) :: a, b
    integer :: error, shape_error

    shape_error = maxval(abs(shape(a) - shape(b)))

    if (shape_error > 0) then
      write(*, '(A, " :: [", A, "] Arrays have mismatching shapes.")') fail, trim(test_name)
    else
      error = maxval(abs(a - b))
      if (error > 0) then
        write(*, '(A, " :: [", A, "] maximum absolute error = ", I8)') fail, trim(test_name), error
      else
        write(*, '(A, " :: [", A, "] maximum absolute error = ", I8)') pass, trim(test_name), error
      end if
    end if

  end subroutine assert_array_equal_2d_int

  subroutine assert_array_equal_3d_int(a, b, test_name)

    character(len=*), intent(in) :: test_name
    integer, intent(in), dimension(:,:,:) :: a, b
    integer :: error, shape_error

    shape_error = maxval(abs(shape(a) - shape(b)))

    if (shape_error > 0) then
      write(*, '(A, " :: [", A, "] Arrays have mismatching shapes.")') fail, trim(test_name)
    else
      error = maxval(abs(a - b))
      if (error > 0) then
        write(*, '(A, " :: [", A, "] maximum absolute error = ", I8)') fail, trim(test_name), error
      else
        write(*, '(A, " :: [", A, "] maximum absolute error = ", I8)') pass, trim(test_name), error
      end if
    end if

  end subroutine assert_array_equal_3d_int

end module test_utils
