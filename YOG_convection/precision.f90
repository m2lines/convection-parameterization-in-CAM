module precision

  ! Sets single and double precision to be used throughout the code.
  ! We use iso_fortran_env which is current Fortran standard.
  ! These may need changing by users to interface with a larger, older codes.
  ! For example, setting sp = 4 and dp = 8 etc.

  use, intrinsic :: iso_fortran_env, only : real32, real64

  implicit none

  public
  integer, parameter :: sp = real32
  integer, parameter :: dp = real64

end module precision
