module example_template_typed

  use bml
  use prec
  use bml_example_template_m

  implicit none

  public :: test_example_template_typed

contains

  function test_example_template_typed(matrix_type, element_kind, element_precision,&
       & n, m) result(test_result)

    character(len=*), intent(in) :: matrix_type, element_kind
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result
    DUMMY_KIND(DUMMY_PREC) :: variable

    variable = 1.0_DM

    ! ADD YOUR TEST HERE

  end function test_example_template_typed

end module example_template_typed
