module convert_matrix_m

  use bml
  use test_m

  implicit none

  private

  type, public, extends(test_t) :: convert_matrix_t
   contains
     procedure, nopass :: test_function
  end type convert_matrix_t

contains

  function test_function(matrix_type, element_type, element_precision, n, m) &
       & result(test_result)

    character(len=*), intent(in) :: matrix_type, element_type
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result

    type(bml_matrix_t) :: a
    type(bml_matrix_t) :: b

    double precision, allocatable :: a_random(:, :)

    REAL_TYPE, allocatable :: a_dense(:, :)
    REAL_TYPE, allocatable :: b_dense(:, :)

    integer :: i, j

    test_result = .true.

    allocate(a_random(n, n))
    call random_number(a_random)
    allocate(a_dense(n, n))
    do i = 1, n
       do j = 1, n
          a_dense(i, j) = a_random(i, j)
       end do
    end do
    call bml_import_from_dense(matrix_type, a_dense, a, 0.0d0, m)
    call bml_export_to_dense(a, b_dense)
    call bml_print_matrix("A", a_dense, 1, n, 1, n)
    call bml_print_matrix("B", b_dense, 1, n, 1, n)
    if (maxval(abs(a_dense - b_dense)) > 1e-12) then
       print *, "Matrix element mismatch"
       test_result = .false.
    end if
    if(test_result) then
       print *, "Test passed"
    end if

    call bml_deallocate(a)
    call bml_deallocate(b)

    deallocate(a_dense)
    deallocate(b_dense)

  end function test_function

end module convert_matrix_m
