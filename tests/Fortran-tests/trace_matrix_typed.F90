module trace_matrix_typed

  use bml
  use prec
  use bml_trace_m

  implicit none

  public :: test_trace_matrix_typed

contains

  function test_trace_matrix_typed(matrix_type, element_kind, element_precision, n, m) &
       & result(test_result)

    character(len=*), intent(in) :: matrix_type, element_kind
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result

    type(bml_matrix_t) :: a

    DUMMY_KIND(DUMMY_PREC), allocatable :: a_dense(:, :)

    real(dp) :: tr_a
    real(dp) :: tr_reference
    real(DUMMY_PREC) :: rel_tol

    integer :: i

    if(element_precision == sp)then
      rel_tol = 1e-6
    elseif(element_precision == dp)then
      rel_tol = 1d-12
    endif

    call bml_random_matrix(matrix_type, element_kind, element_precision, n, m, &
         & a)
    tr_a = bml_trace(a)

    tr_reference = 0_MP
    call bml_export_to_dense(a, a_dense)
    do i = 1, n
      tr_reference = tr_reference+a_dense(i, i)
    end do

    if(abs((tr_a - tr_reference)/tr_reference) > rel_tol) then
      test_result = .false.
      print *, "trace incorrect"
    else
      test_result = .true.
    end if

    call bml_deallocate(a)
    deallocate(a_dense)

  end function test_trace_matrix_typed

end module trace_matrix_typed
