module bandwidth_m

  use bml
  use test_m

  implicit none

  private

  type, public, extends(test_t) :: bandwidth_t
  contains
    procedure, nopass :: test_function
  end type bandwidth_t

contains

  function test_function(matrix_type, element_type, element_precision, n, m) &
       & result(test_result)

    character(len=*), intent(in) :: matrix_type, element_type
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result

    type(bml_matrix_t) :: a
    REAL_TYPE :: a_dense(n, n)
    integer :: i, j
    integer :: bandwidth

    test_result = .true.

    a_dense = 0
    do i = 1, n
      do j = i, min(i + 4, n)
        a_dense(i, j) = 1
      end do
    end do
    call bml_print_matrix("A", a_dense, 1, n, 1, n)
    call bml_import_from_dense(matrix_type, a_dense, a, 0d0, n)
    call bml_print_matrix("A", a, 1, n, 1, n)

    do i = 1, n
      bandwidth = bml_get_row_bandwidth(a, i)
      if(bandwidth /= min(i + 4, n) - i + 1) then
        print *, "row ", i, " incorrect bandwidth, ", bandwidth, &
             ", expected ", min(i + 4, n) - i + 1
        test_result = .false.
      else
        print *, "row ", i, " correct bandwidth, ", bandwidth
      end if
    end do
    bandwidth = bml_get_bandwidth(a)
    print *, "matrix banddwidth = ", bandwidth
    if(bandwidth /= 5) then
      print *, "total bandwidth incorrect"
      test_result = .false.
    end if

    call bml_deallocate(a)

    call bml_banded_matrix(matrix_type, element_type, element_precision, n, &
         & m/2, a)
    call bml_print_matrix("A", a, 1, n, 1, n)

    call bml_deallocate(a)

  end function test_function

end module bandwidth_m
