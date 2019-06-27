program test

  ! All tests need the Fortran kinds corresponding to the C floating types.
  use, intrinsic :: iso_c_binding, only : C_FLOAT, C_DOUBLE, C_FLOAT_COMPLEX, &
       & C_DOUBLE_COMPLEX
  use bml
  use TEST_MODULE

  implicit none

  integer, parameter :: N = 7, M = 7
  type(TEST_TYPE) :: tester

  write(*, "(A)") "Testing "//MATRIX_TYPE//":"//MATRIX_PRECISION
  if(.not. tester%test_function(MATRIX_TYPE, REAL_NAME, REAL_KIND, N, M)) then
     write(*, "(A)") "Test failed"
     error stop
  end if

end program test
