module transpose_matrix_typed

  use bml
  use prec
  use bml_transpose_m

  implicit none

  public :: test_transpose_matrix_typed

contains

  function test_transpose_matrix_typed(matrix_type, element_kind, element_precision, n, m) &
       & result(test_result)

    character(len=*), intent(in) :: matrix_type, element_kind
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result

    type(bml_matrix_t) :: a
    type(bml_matrix_t) :: b
    type(bml_matrix_t) :: c

    DUMMY_KIND(DUMMY_PREC), allocatable :: a_dense(:, :)
    DUMMY_KIND(DUMMY_PREC), allocatable :: b_dense(:, :)

    real(DUMMY_PREC) :: tol

    if(element_precision == sp)then
      tol = 1e-6
    elseif(element_precision == dp)then
      tol = 1d-12
    endif

    call bml_random_matrix(matrix_type, element_kind, element_precision, n, m, &
         & a)
    call bml_transpose(a, b)
    call bml_copy_new(a, c)

    call bml_export_to_dense(a, a_dense)
    call bml_export_to_dense(b, b_dense)

    if(maxval(abs(a_dense-transpose(b_dense))) > tol) then
      test_result = .false.
      print *, "matrices are not transposes"
    else
      test_result = .true.
    end if

    call bml_deallocate(a)
    call bml_deallocate(b)
    call bml_deallocate(c)

    deallocate(a_dense)
    deallocate(b_dense)

  end function test_transpose_matrix_typed

end module transpose_matrix_typed
