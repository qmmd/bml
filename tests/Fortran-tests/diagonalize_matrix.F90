module diagonalize_matrix_m

  use bml
  use test_m

  implicit none

  private

  type, public, extends(test_t) :: diagonalize_matrix_t
  contains
    procedure, nopass :: test_function
  end type diagonalize_matrix_t

contains

  function test_function(matrix_type, element_type, element_precision, n, m) &
       & result(test_result)

    character(len=*), intent(in) :: matrix_type, element_type
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result
    double precision :: threshold

    type(bml_matrix_t) :: a
    type(bml_matrix_t) :: a_t
    type(bml_matrix_t) :: b
    type(bml_matrix_t) :: c
    type(bml_matrix_t) :: eigenvectors
    type(bml_matrix_t) :: eigenvectors_t
    double precision, allocatable :: eigenvalues(:)

    threshold = 0.0
    test_result = .false.

    call bml_random_matrix(matrix_type, element_type, element_precision, n, m, &
         & a)
    call bml_print_matrix("A", a, 1, n, 1, n)
    call bml_transpose(a, a_t)
    call bml_print_matrix("A_t", a_t, 1, n, 1, n)
    call bml_add(0.5d0, a, 0.5d0, a_t, threshold)
    call bml_print_matrix("A", a, 1, n, 1, n)
    allocate(eigenvalues(n))
    call bml_zero_matrix(matrix_type, element_type, element_precision, n, m, &
         & eigenvectors)
    call bml_diagonalize(a, eigenvalues, eigenvectors)
    call bml_transpose(eigenvectors, eigenvectors_t)
    call bml_zero_matrix(matrix_type, element_type, element_precision, n, m, b)
    call bml_zero_matrix(matrix_type, element_type, element_precision, n, m, c)
    call bml_multiply(eigenvectors_t, eigenvectors, b)
    !call bml_multiply(b, eigenvectors, c)
    write(*, *) eigenvalues
    call bml_print_matrix("eigenvectors", eigenvectors, 1, n, 1, n)
    call bml_print_matrix("U^t U", b, 1, n, 1, n)
    test_result = .true.

    call bml_deallocate(a)
    call bml_deallocate(a_t)
    call bml_deallocate(b)
    call bml_deallocate(c)
    call bml_deallocate(eigenvectors)
    call bml_deallocate(eigenvectors_t)

    deallocate(eigenvalues)

  end function test_function

end module diagonalize_matrix_m
