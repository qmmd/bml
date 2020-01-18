module set_element_typed

  use bml
  use bml_elemental_m
  use prec

  implicit none

  public :: test_set_element_typed

contains

  function test_set_element_typed(matrix_type, element_kind, element_precision, n, m) &
       & result(test_result)

    character(len=*), intent(in) :: matrix_type, element_kind
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result

    type(bml_matrix_t) :: a

    integer :: i, j

    DUMMY_KIND(DUMMY_PREC), allocatable :: a_dense(:, :)
    DUMMY_KIND(DUMMY_PREC) :: a_ij

    call bml_zero_matrix(matrix_type, element_kind, element_precision, n, m, &
         & a)

    test_result = .true.

    do i = 1, n
      do j = 1, n
        a_ij = real(i,DUMMY_PREC)*real(j,DUMMY_PREC)
        call bml_set_element(a, i, j, a_ij)
      end do
    end do

    call bml_print_matrix("A", a_dense, 1, n, 1, n)
    call bml_export_to_dense(a, a_dense)

    do i = 1, n
      do j = 1, n
        a_ij = real(i,DUMMY_PREC)*real(j,DUMMY_PREC)
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
    deallocate(a_dense)

  end function test_set_element_typed

end module set_element_typed
