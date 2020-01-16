module threshold_matrix_typed

  use bml
  use prec
  use bml_threshold_m

  implicit none

  public :: test_threshold_matrix_typed

contains

  function test_threshold_matrix_typed(matrix_type, element_kind, element_precision, n, m) &
       & result(test_result)

    character(len=*), intent(in) :: matrix_type, element_kind
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result

    type(bml_matrix_t) :: a
    DUMMY_KIND(DUMMY_PREC), allocatable :: a_dense(:, :)
    integer :: i, j

    call bml_random_matrix(matrix_type, element_kind, element_precision, n, m, &
         & a)
    call bml_print_matrix("A", a, 1, n, 1, n)
    call bml_threshold(a, 0.5_dp)
    call bml_print_matrix("A", a, 1, n, 1, n)
    call bml_export_to_dense(a, a_dense)

    test_result = .true.
    do i = 1, n
      do j = 1, n
        if(abs(a_dense(i, j)) > 0.0_MP .and. abs(a_dense(i, j)) < 0.5_MP) then
          test_result = .false.
          call bml_print_matrix("A", a_dense, 1, n, 1, n)
          print *, "matrix not thresholded"
          return
        end if
      end do
    end do

    call bml_deallocate(a)

  end function test_threshold_matrix_typed

end module threshold_matrix_typed
