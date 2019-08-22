module test_m
  ! All tests need the Fortran kinds corresponding to the C floating types.
  use, intrinsic :: iso_c_binding, only : C_FLOAT, C_DOUBLE, C_FLOAT_COMPLEX, &
       & C_DOUBLE_COMPLEX

  implicit none

  type, abstract :: test_t
  contains
    procedure(test_function_interface), nopass, deferred :: test_function
  end type test_t

  abstract interface
    function test_function_interface(matrix_type, element_type, &
         & element_precision, n, m) result(test_result)
      character(len=*), intent(in) :: matrix_type, element_type
      integer, intent(in) :: element_precision
      integer, intent(in) :: n, m
      logical :: test_result
    end function test_function_interface
  end interface

end module test_m
