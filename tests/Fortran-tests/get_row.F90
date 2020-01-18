module get_row

  use bml
  use prec
  use get_row_single_real
  use get_row_double_real
#ifdef BML_COMPLEX
  use get_row_single_complex
  use get_row_double_complex
#endif

  implicit none

  public :: test_get_row

contains

  function test_get_row(matrix_type, element_type, n, m) &
       & result(test_result)

    character(len=*), intent(in) :: matrix_type, element_type
    integer, intent(in) :: n, m
    character(20) :: element_kind
    logical :: test_result
    integer :: element_precision

    write(*,*)"Im in test_get_row"
    select case(element_type)
    case("single_real")
      element_kind = bml_real
      element_precision = sp
      test_result = test_get_row_single_real(matrix_type, element_kind,&
           &element_precision, n, m)
    case("double_real")
      element_kind = bml_real
      element_precision = dp
      test_result = test_get_row_double_real(matrix_type, element_kind,&
           &element_precision, n, m)
#ifdef BML_COMPLEX
    case("single_complex")
      element_kind = bml_complex
      element_precision = sp
      test_result = test_get_row_single_complex(matrix_type, element_kind,&
           &element_precision, n, m)
    case("double_complex")
      element_kind = bml_complex
      element_precision = dp
      test_result = test_get_row_double_complex(matrix_type, element_kind,&
           &element_precision, n, m)
#endif
    end select

  end function test_get_row

end module get_row
