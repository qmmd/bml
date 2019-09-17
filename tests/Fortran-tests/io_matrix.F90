module io_matrix_m

  use bml
  use test_m

  implicit none

  private

  type, public, extends(test_t) :: io_matrix_t
  contains
    procedure, nopass :: test_function
  end type io_matrix_t

contains

  function test_function(matrix_type, element_type, element_precision, n, m) &
       & result(test_result)

    character(len=*), intent(in) :: matrix_type, element_type
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result

    character(len=*), parameter :: fname="ftest_matrix.mtx"

    type(bml_matrix_t) :: a
    type(bml_matrix_t) :: b

    REAL_TYPE, allocatable :: a_dense(:, :)
    REAL_TYPE, allocatable :: b_dense(:, :)

    allocate(a_dense(n,n))
    allocate(b_dense(n,n))

    call bml_random_matrix(matrix_type, element_type, element_precision, n, m, &
         & a)
    call bml_write_matrix(a, fname)
    call bml_zero_matrix(matrix_type, element_type, element_precision, n, m, b)
    call bml_read_matrix(b, fname)

    call bml_export_to_dense(a, a_dense)
    call bml_export_to_dense(b, b_dense)
    call bml_print_matrix("A", a_dense, 1, n, 1, n)
    call bml_print_matrix("B", b_dense, 1, n, 1, n)
    if(maxval(abs(a_dense - b_dense)) > 1e-12) then
      test_result = .false.
      print *, "matrices are not identical"
    else
      test_result = .true.
    end if

  end function test_function

end module io_matrix_m
