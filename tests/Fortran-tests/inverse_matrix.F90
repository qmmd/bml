module inverse_matrix_m

  use bml
  use test_m

  implicit none

  private

  type, public, extends(test_t) :: inverse_matrix_t
   contains
     procedure, nopass :: test_function
  end type inverse_matrix_t

contains

  function test_function(matrix_type, element_type, element_precision, n, m) &
       & result(test_result)

    character(len=*), intent(in) :: matrix_type, element_type
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result
    double precision :: threshold

    type(bml_matrix_t) :: a
    type(bml_matrix_t) :: a_inverse
    type(bml_matrix_t) :: aux

    threshold = 0.0
    test_result = .false.

    call bml_random_matrix(matrix_type, element_type, element_precision, n, m, &
         & a)
    call bml_zero_matrix(matrix_type, element_type, element_precision, n, m, &
         & a_inverse)
    call bml_zero_matrix(matrix_type, element_type, element_precision, n, m, &
         & aux)
    call bml_print_matrix("A", a, 1, n, 1, n)
    call bml_inverse(a, a_inverse)
    call bml_print_matrix("A_inverse", a_inverse, 1, n, 1, n)

    call bml_multiply(a, a_inverse, aux, 1.0, 0.0, threshold)
    call bml_print_matrix("A*A_inverse", aux, 1, n, 1, n)
    test_result = .true.

    call bml_deallocate(a)
    call bml_deallocate(a_inverse)
    call bml_deallocate(aux)

  end function test_function

end module inverse_matrix_m
