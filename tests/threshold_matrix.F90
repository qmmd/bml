module threshold_matrix_m

  use bml
  use test_m

  implicit none

  private

  type, public, extends(test_t) :: threshold_matrix_t
   contains
     procedure, nopass :: test_function
  end type threshold_matrix_t

contains

  function test_function(matrix_type, element_type, element_precision, n, m) &
       & result(test_result)

    character(len=*), intent(in) :: matrix_type, element_type
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result

    type(bml_matrix_t) :: a
    REAL_TYPE, allocatable :: a_dense(:, :)
    integer :: i, j

    call bml_random_matrix(matrix_type, element_type, element_precision, n, m, &
         & a)
    call bml_print_matrix("A", a, 1, n, 1, n)
    call bml_threshold(a, 0.5d0)
    call bml_print_matrix("A", a, 1, n, 1, n)
    call bml_export_to_dense(a, a_dense)

    test_result = .true.
    do i = 1, n
       do j = 1, n
          if(abs(a_dense(i, j)) > 0.0 .and. abs(a_dense(i, j)) < 0.5) then
             test_result = .false.
             call bml_print_matrix("A", a_dense, 1, n, 1, n)
             print *, "matrix not thresholded"
             return
          end if
       end do
    end do

    call bml_deallocate(a)

  end function test_function

end module threshold_matrix_m
