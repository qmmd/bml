module allocate_matrix_typed

  use bml
  use prec
  use bml_allocate_m


  implicit none

  public :: test_allocate_matrix_typed

contains

  function test_allocate_matrix_typed(matrix_type, element_kind, element_precision,&
       & n, m) result(test_result)

    character(len=*), intent(in) :: matrix_type, element_kind
    integer, intent(in) :: element_precision
    integer, intent(in) :: n, m
    logical :: test_result
    real(DUMMY_PREC) :: tol

    type(bml_matrix_t) :: a

    DUMMY_KIND(DUMMY_PREC), allocatable :: a_dense(:, :)
    integer :: i, j

    test_result = .true.

    if(element_precision == sp)then
      tol = 1e-6
    elseif(element_precision == dp)then
      tol = 1d-12
    endif

    call bml_random_matrix(matrix_type, element_kind, element_precision, n, m, &
         & a)
    call bml_export_to_dense(a, a_dense)
    if(lbound(a_dense, 1) /= 1 .or. lbound(a_dense, 2) /= 1) then
      print *, "incorrect lbound"
      print *, "lbound(a_dense, 1) = ", lbound(a_dense, 1)
      print *, "lbound(a_dense, 2) = ", lbound(a_dense, 2)
      test_result = .false.
      return
    end if
    if(ubound(a_dense, 1) /= n .or. ubound(a_dense, 2) /= n) then
      print *, "incorrect ubound"
      print *, "ubound(a_dense, 1) = ", ubound(a_dense, 1)
      print *, "ubound(a_dense, 2) = ", ubound(a_dense, 2)
      test_result = .false.
      return
    end if
    call bml_print_matrix("A", a, 1, n, 1, n)

    deallocate(a_dense)

    call bml_deallocate(a)
    call bml_identity_matrix(matrix_type, element_kind, element_precision, &
         & n, m, a)
    call bml_export_to_dense(a, a_dense)
    call bml_print_matrix("A", a, lbound(a_dense, 1), ubound(a_dense, 1), &
         lbound(a_dense, 2), ubound(a_dense, 2))
    write(*,*)a_dense
    do i = 1, n
      do j = 1, n
        if(i == j) then
          if(abs(a_dense(i, j)-1.0_MP) > tol) then
            print *, "Incorrect value on diagonal", a_dense(i, j)
            test_result = .false.
            return
          end if
        else
          if(abs(a_dense(i, j)) > tol) then
            print *, "Incorrect value off diagonal", a_dense(i, j)
            test_result = .false.
            return
          end if
        end if
      end do
    end do
    print *, "Identity matrix test passed"

    call bml_zero_matrix(matrix_type, element_kind, element_precision, n, m, a)
    call bml_deallocate(a)

  end function test_allocate_matrix_typed

end module allocate_matrix_typed
