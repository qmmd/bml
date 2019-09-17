module trace_matrix_m

  use bml
  use test_m

  implicit none

  private

  type, public, extends(test_t) :: trace_matrix_t
  contains
    procedure, nopass :: test_function
  end type trace_matrix_t

contains

  function test_function(matrix_type, element_type, element_precision, n, m) &
       & result(test_result)

    character(len=*), intent(in) :: matrix_type, element_type
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result

#if defined(SINGLE_REAL) || defined(SINGLE_COMPLEX)
    double precision :: rel_tol = 1e-6
#else
    double precision :: rel_tol = 1d-12
#endif

    type(bml_matrix_t) :: a

    REAL_TYPE, allocatable :: a_dense(:, :)

    double precision :: tr_a
    double precision :: tr_reference

    integer :: i

    call bml_random_matrix(matrix_type, element_type, element_precision, n, m, &
         & a)
    tr_a = bml_trace(a)

    tr_reference = 0
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

  end function test_function

end module trace_matrix_m
