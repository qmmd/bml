module get_element_typed

  use bml
  use bml_elemental_m
  use prec

  implicit none

  public :: test_get_element_typed

contains

  function test_get_element_typed(matrix_type, element_kind, element_precision, n, m) &
       & result(test_result)

    character(len=*), intent(in) :: matrix_type, element_kind
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result

    type(bml_matrix_t) :: a

    integer :: i, j

    DUMMY_KIND(DUMMY_PREC), allocatable :: a_dense(:, :)
    DUMMY_KIND(DUMMY_PREC) :: a_ij
    real(DUMMY_PREC) :: tol

    if(element_precision == sp)then
      tol = 1e-6
    elseif(element_precision == dp)then
      tol = 1d-12
    endif

    call bml_random_matrix(matrix_type, element_kind, element_precision, n, m, &
         & a)

    test_result = .true.

    call bml_export_to_dense(a, a_dense)
    call bml_print_matrix("A", a_dense, 1, n, 1, n)
    do i = 1, n
      do j = 1, n
        call bml_get_element(a_ij, a, i, j)
        if(abs(a_ij-a_dense(i, j)) > tol) then
          test_result = .false.
          print *, "matrix element mismatch"
          print *, "got ", a_ij
          print *, "expected ", a_dense(i, j)
          return
        end if
      end do
    end do
    call bml_deallocate(a)

  end function test_get_element_typed

end module get_element_typed
