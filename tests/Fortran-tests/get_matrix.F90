module get_matrix_m

  use bml
  use test_m

  implicit none

  private

  type, public, extends(test_t) :: get_matrix_t
   contains
     procedure, nopass :: test_function
  end type get_matrix_t

contains

  function test_function(matrix_type, element_type, element_precision, n, m) &
       & result(test_result)

    character(len=*), intent(in) :: matrix_type, element_type
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result

    type(bml_matrix_t) :: a

    integer :: i, j

    REAL_TYPE, allocatable :: a_dense(:, :)
    REAL_TYPE :: a_ij

    call bml_random_matrix(matrix_type, element_type, element_precision, n, m, &
         & a)

    test_result = .true.

    call bml_export_to_dense(a, a_dense)
    call bml_print_matrix("A", a_dense, 1, n, 1, n)
    do i = 1, n
       do j = 1, n
          call bml_get(a_ij, a, i, j)
          if(abs(a_ij-a_dense(i, j)) > 1e-12) then
             test_result = .false.
             print *, "matrix element mismatch"
             print *, "got ", a_ij
             print *, "expected ", a_dense(i, j)
             return
          end if
       end do
    end do
    call bml_deallocate(a)

  end function test_function

end module get_matrix_m
